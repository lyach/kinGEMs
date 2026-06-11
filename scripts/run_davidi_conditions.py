#!/usr/bin/env python3
"""
Davidi Multi-Condition Pipeline
================================

Run enzyme-constrained FBA for each experimental condition defined in a CSV
file (e.g. ``davidi_media.csv``).  Data preparation (steps 1-3) is performed
once and the prepared model is deep-copied for each condition.

Usage:
    python scripts/run_davidi_conditions.py configs/davidi_iML1515_GEM.json
    python scripts/run_davidi_conditions.py configs/davidi_iML1515_GEM.json --force

The CSV must have columns:
  - ``condition_id``  – label for the condition
  - ``avg_growth``    – growth rate to fix on the biomass reaction
  - ``EX_*``          – exchange reactions; 1000 = open, 0 = closed, NaN = closed
"""

from __future__ import annotations

from copy import deepcopy
from datetime import datetime
import json
import os
import random
import sys
import time
import warnings

import cobra
import pandas as pd

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from kinGEMs.dataset import (
    annotate_model_with_kcat_and_gpr,
    load_model,
    merge_substrate_sequences,
    prepare_model_data,
    process_kcat_predictions,
    convert_to_irreversible,
)
from kinGEMs.dataset_modelseed import prepare_modelseed_model_data
from kinGEMs.modeling.optimize import (
    apply_medium_irreversible,
    run_optimization_with_dataframe,
)

warnings.filterwarnings('ignore')
try:
    import gurobipy
    gurobipy.setParam('OutputFlag', 0)
except ImportError:
    pass


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def load_config(config_path: str) -> dict:
    with open(config_path, 'r') as f:
        return json.load(f)


def is_modelseed_model(model_name: str) -> bool:
    return '_genome_' in model_name.lower()


def determine_biomass_reaction(model):
    obj_rxns = {
        rxn.id: rxn.objective_coefficient
        for rxn in model.reactions
        if rxn.objective_coefficient != 0
    }
    if not obj_rxns:
        raise ValueError("No objective reaction found in model")
    return next(iter(obj_rxns.keys()))


def find_predictions_file(model_name: str, cpipred_dir: str) -> str:
    import glob
    patterns = [
        f"X06A_kinGEMs_{model_name}_predictions.csv",
        f"*{model_name}*predictions.csv",
        f"*{model_name.replace('_GEM', '')}*predictions.csv",
    ]
    if '_GEM' in model_name:
        base = model_name.replace('_GEM', '')
        patterns.append(f"*ecoli_{base}*predictions.csv")
        patterns.append(f"*{base}*predictions.csv")

    for pat in patterns:
        matches = glob.glob(os.path.join(cpipred_dir, pat))
        if matches:
            print(f"  Found predictions file: {os.path.basename(matches[0])}")
            return matches[0]

    raise FileNotFoundError(
        f"No CPI-Pred predictions file found for '{model_name}' in {cpipred_dir}"
    )


def get_enzyme_constraints(config: dict) -> dict:
    """Read enzyme constraint flags from config with pipeline defaults."""
    ec = config.get('enzyme_constraints', {})
    return {
        'multi_enzyme_off': ec.get('multi_enzyme_off', False),
        'isoenzymes_off': ec.get('isoenzymes_off', False),
        'promiscuous_off': ec.get('promiscuous_off', False),
        'complexes_off': ec.get('complexes_off', False),
    }


def parse_condition_medium(row: pd.Series):
    """Extract exchange-reaction medium dict and growth value from a CSV row.

    Returns
    -------
    condition_id : str
    avg_growth : float
    medium : dict
        ``{reaction_id: flux_value}`` for every ``EX_*`` column with a
        non-NaN value.
    """
    condition_id = row['condition_id']
    avg_growth = float(row['avg_growth'])

    medium = {}
    for col, val in row.items():
        if col.startswith('EX_') and pd.notna(val):
            medium[col] = float(val)

    return condition_id, avg_growth, medium


# ---------------------------------------------------------------------------
# Data preparation (steps 1-3, run once)
# ---------------------------------------------------------------------------

def prepare_data(config: dict, force_regenerate: bool = False):
    """Run pipeline steps 1-3 and return the prepared model + processed data.

    Returns
    -------
    model : cobra.Model
        Irreversible model annotated with kcat / GPR info.
    processed_data : pd.DataFrame
    biomass_reaction : str
    """
    model_name = config['model_name']
    organism = config.get('organism', 'Unknown')
    solver_name = config.get('solver', 'glpk')
    results_subdir = config.get('results_subdir', None)

    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    data_dir = os.path.join(project_root, "data")
    raw_data_dir = os.path.join(data_dir, "raw")
    interim_data_dir = os.path.join(data_dir, "interim", model_name)
    processed_data_dir = os.path.join(data_dir, "processed", model_name)
    cpipred_dir = os.path.join(data_dir, "interim", "CPI-Pred predictions")
    os.makedirs(interim_data_dir, exist_ok=True)
    os.makedirs(processed_data_dir, exist_ok=True)

    if results_subdir == "BiGG_models":
        model_path = os.path.join(raw_data_dir, "BiGG_models", f"{model_name}.xml")
    else:
        model_path = os.path.join(raw_data_dir, f"{model_name}.xml")

    substrates_output = os.path.join(interim_data_dir, f"{model_name}_substrates.csv")
    sequences_output = os.path.join(interim_data_dir, f"{model_name}_sequences.csv")
    merged_data_output = os.path.join(interim_data_dir, f"{model_name}_merged_data.csv")
    processed_data_output = os.path.join(processed_data_dir, f"{model_name}_processed_data.csv")

    # Determine biomass reaction from the original (reversible) model
    temp_model = cobra.io.read_sbml_model(model_path)
    temp_model.solver = solver_name
    biomass_reaction = config.get('biomass_reaction') or determine_biomass_reaction(temp_model)
    print(f"  Biomass reaction: {biomass_reaction}")
    print(f"  Baseline growth (unconstrained): {temp_model.slim_optimize():.4f}")

    # --- Step 1: Prepare model data ---
    print("=== Step 1: Preparing model data ===")
    if (not force_regenerate
            and os.path.exists(substrates_output)
            and os.path.exists(sequences_output)):
        print("  Found existing files, loading cached data")
        substrate_df = pd.read_csv(substrates_output)
        sequences_df = pd.read_csv(sequences_output)
        model = load_model(model_path)
        model.solver = solver_name
        model = convert_to_irreversible(model)
    else:
        if is_modelseed_model(model_name):
            metadata_dir = config.get('metadata_dir',
                                      os.path.join(data_dir, "Biolog experiments"))
            model, substrate_df, sequences_df = prepare_modelseed_model_data(
                model_path=model_path,
                substrates_output=substrates_output,
                sequences_output=sequences_output,
                organism=organism,
                metadata_dir=metadata_dir,
            )
        else:
            model, substrate_df, sequences_df = prepare_model_data(
                model_path=model_path,
                substrates_output=substrates_output,
                sequences_output=sequences_output,
                organism=organism,
                convert_irreversible=True,
            )

    model.solver = solver_name
    print(f"  Model: {len(model.genes)} genes, {len(model.reactions)} reactions")

    # --- Step 2: Merge substrate and sequence data ---
    print("=== Step 2: Merging substrate and sequence data ===")
    if not force_regenerate and os.path.exists(merged_data_output):
        merged_data = pd.read_csv(merged_data_output)
    else:
        merged_data = merge_substrate_sequences(
            substrate_df=substrate_df,
            sequences_df=sequences_df,
            model=model,
            output_path=merged_data_output,
        )
    print(f"  Merged data: {len(merged_data)} rows")

    # --- Step 3: Process kcat predictions ---
    print("=== Step 3: Processing CPI-Pred kcat values ===")
    if not force_regenerate and os.path.exists(processed_data_output):
        processed_data = pd.read_csv(processed_data_output)
    else:
        predictions_path = find_predictions_file(model_name, cpipred_dir)
        processed_data = process_kcat_predictions(
            merged_df=merged_data,
            predictions_csv_path=predictions_path,
            output_path=processed_data_output,
        )
    print(f"  Processed data: {len(processed_data)} rows")

    if 'kcat_mean' in processed_data.columns and 'kcat' not in processed_data.columns:
        processed_data['kcat'] = processed_data['kcat_mean']
    elif 'kcat_y' in processed_data.columns and 'kcat' not in processed_data.columns:
        processed_data['kcat'] = processed_data['kcat_y']

    model = annotate_model_with_kcat_and_gpr(model=model, df=processed_data)

    return model, processed_data, biomass_reaction


# ---------------------------------------------------------------------------
# Per-condition optimization (step 4 only)
# ---------------------------------------------------------------------------

def run_condition(
    model,
    processed_data,
    biomass_reaction: str,
    condition_id: str,
    avg_growth: float,
    medium: dict,
    enzyme_upper_bound: float,
    solver_name: str,
    output_dir: str,
    enzyme_constraints: dict,
):
    """Run enzyme-constrained FBA for a single experimental condition.

    Returns a dict with summary metrics.
    """
    t0 = time.time()

    model_copy = deepcopy(model)
    model_copy.solver = solver_name

    # Apply medium + fix growth
    apply_medium_irreversible(
        model_copy,
        medium,
        growth_reaction=biomass_reaction,
        growth_value=None, # to constrain growth, use avg_growth
        verbose=True,
    )

    # COBRApy baseline (media-constrained, no enzyme constraints)
    cobra_sol = model_copy.optimize()
    cobra_biomass = cobra_sol.objective_value if cobra_sol.status == 'optimal' else None
    if cobra_biomass is not None:
        print(f"  Baseline growth (media-constrained): {cobra_biomass:.4f}")
    else:
        print(f"  Baseline growth (media-constrained): infeasible ({cobra_sol.status})")

    # Enzyme-constrained optimization
    sol_val, df_FBA, gene_seq_dict, _ = run_optimization_with_dataframe(
        model=model_copy,
        processed_df=processed_data,
        objective_reaction=biomass_reaction,
        enzyme_upper_bound=enzyme_upper_bound,
        enzyme_ratio=True,
        maximization=True,
        multi_enzyme_off=enzyme_constraints['multi_enzyme_off'],
        isoenzymes_off=enzyme_constraints['isoenzymes_off'],
        promiscuous_off=enzyme_constraints['promiscuous_off'],
        complexes_off=enzyme_constraints['complexes_off'],
        output_dir=None,
        save_results=False,
        print_reaction_conditions=False,
        verbose=False,
        solver_name=solver_name,
        medium=None,
        medium_upper_bound=False,
    )

    feasible = sol_val is not None
    elapsed = time.time() - t0

    # Save per-condition results
    cond_dir = os.path.join(output_dir, condition_id)
    os.makedirs(cond_dir, exist_ok=True)
    if df_FBA is not None and not df_FBA.empty:
        df_FBA.to_csv(os.path.join(cond_dir, "df_FBA.csv"), index=False)

    return {
        'condition_id': condition_id,
        'avg_growth': avg_growth,
        'feasible': feasible,
        'ec_biomass': round(float(sol_val), 6) if feasible else None,
        'cobra_biomass': round(float(cobra_biomass), 6) if cobra_biomass is not None else None,
        'n_medium_open': sum(1 for v in medium.values() if pd.notna(v) and v > 0),
        'n_medium_closed': sum(1 for v in medium.values() if pd.notna(v) and v == 0),
        'elapsed_s': round(elapsed, 1),
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    config_path = sys.argv[1]
    force_regenerate = '--force' in sys.argv or '-f' in sys.argv

    config = load_config(config_path)
    model_name = config['model_name']
    enzyme_upper_bound = config.get('enzyme_upper_bound', 0.15)
    solver_name = config.get('solver', 'glpk')
    enzyme_constraints = get_enzyme_constraints(config)

    # Resolve conditions CSV path relative to project root
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    conditions_csv = config.get('conditions_csv')
    if conditions_csv is None:
        print("Error: config must include 'conditions_csv' pointing to the media CSV")
        sys.exit(1)
    if not os.path.isabs(conditions_csv):
        conditions_csv = os.path.join(project_root, conditions_csv)

    run_id = f"{model_name}_davidi_{datetime.today().strftime('%Y%m%d')}_{random.randint(1000, 9999)}"
    results_subdir = config.get('results_subdir', None)
    if results_subdir:
        output_dir = os.path.join(project_root, "results", "davidi_conditions", results_subdir, run_id)
    else:
        output_dir = os.path.join(project_root, "results", "davidi_conditions", run_id)
    os.makedirs(output_dir, exist_ok=True)

    print("\n" + "=" * 70)
    print(f"=== Davidi Multi-Condition Pipeline for {model_name} ===")
    print("=" * 70)
    print(f"Run ID:         {run_id}")
    print(f"Conditions CSV: {conditions_csv}")
    print(f"Results dir:    {output_dir}")
    print(f"Enzyme UB:      {enzyme_upper_bound}")
    print(f"Solver:         {solver_name}")
    print(
        "Enzyme mode:    "
        f"multi_enzyme_off={enzyme_constraints['multi_enzyme_off']}, "
        f"isoenzymes_off={enzyme_constraints['isoenzymes_off']}, "
        f"promiscuous_off={enzyme_constraints['promiscuous_off']}, "
        f"complexes_off={enzyme_constraints['complexes_off']}"
    )
    print("=" * 70)

    # ------------------------------------------------------------------
    # 1) Data preparation (once)
    # ------------------------------------------------------------------
    t_prep = time.time()
    print("\n--- Data Preparation (steps 1-3, run once) ---")
    model, processed_data, biomass_reaction = prepare_data(config, force_regenerate)
    print(f"  Data prep finished in {time.time() - t_prep:.1f}s\n")

    # ------------------------------------------------------------------
    # 2) Read conditions CSV
    # ------------------------------------------------------------------
    conditions_df = pd.read_csv(conditions_csv)
    n_conditions = len(conditions_df)
    print(f"Loaded {n_conditions} conditions from {os.path.basename(conditions_csv)}")

    # ------------------------------------------------------------------
    # 3) Loop over conditions
    # ------------------------------------------------------------------
    summary_rows = []
    for idx, row in conditions_df.iterrows():
        condition_id, avg_growth, medium = parse_condition_medium(row)
        print(f"\n{'='*60}")
        print(f"[{idx + 1}/{n_conditions}] Condition: {condition_id}  |  growth={avg_growth}")
        print(f"{'='*60}")

        result = run_condition(
            model=model,
            processed_data=processed_data,
            biomass_reaction=biomass_reaction,
            condition_id=condition_id,
            avg_growth=avg_growth,
            medium=medium,
            enzyme_upper_bound=enzyme_upper_bound,
            solver_name=solver_name,
            output_dir=output_dir,
            enzyme_constraints=enzyme_constraints,
        )
        summary_rows.append(result)

        status = "FEASIBLE" if result['feasible'] else "INFEASIBLE"
        ec = result['ec_biomass'] if result['ec_biomass'] is not None else "N/A"
        cobra = result['cobra_biomass'] if result['cobra_biomass'] is not None else "N/A"
        print(f"  Result: {status}  |  COBRApy: {cobra}  |  EC biomass: {ec}  |  {result['elapsed_s']}s")

    # ------------------------------------------------------------------
    # 4) Save summary
    # ------------------------------------------------------------------
    summary_df = pd.DataFrame(summary_rows)
    summary_path = os.path.join(output_dir, "summary.csv")
    summary_df.to_csv(summary_path, index=False)

    # Save a copy of the config used
    with open(os.path.join(output_dir, "config_used.json"), 'w') as f:
        json.dump(config, f, indent=2)

    # ------------------------------------------------------------------
    # Final report
    # ------------------------------------------------------------------
    n_feasible = summary_df['feasible'].sum()
    print("\n" + "=" * 70)
    print("=== Pipeline Complete ===")
    print("=" * 70)
    print(f"Conditions tested:  {n_conditions}")
    print(f"Feasible:           {n_feasible}/{n_conditions}")
    print(f"Summary saved to:   {summary_path}")
    print(f"Results directory:  {output_dir}")
    print("=" * 70)


if __name__ == '__main__':
    main()
