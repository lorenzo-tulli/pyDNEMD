#!/usr/bin/env python3

"""
Analyse equilibrium simulations: compute and plot Cα RMSD and RMSF
for all runs, using the energy-minimised structure as reference.

Outputs (written to output_dir/results/equilibrium/):
    rmsd_run<N>.xvg      — raw GROMACS RMSD data
    rmsf_run<N>.xvg      — raw GROMACS RMSF data
    rmsd_all_runs.png    — overlaid RMSD plot
    rmsf_all_runs.png    — overlaid RMSF plot
    summary.csv          — per-run mean/max RMSD and RMSF

Usage:
    python scripts/analyse_equilibrium.py --config examples/config.yaml
    python scripts/analyse_equilibrium.py --config examples/config.yaml --run 3
"""

import argparse
from pathlib import Path
from modules.config import Config
from modules.equilibrium_analysis import EquilibriumAnalyser
from modules.analysis import plot_rmsd, plot_rmsf, write_summary_csv
from modules.utils import ensure_dir, get_logger

logger = get_logger("analyse_equilibrium")


def main():
    parser = argparse.ArgumentParser(
        description="Analyse equilibrium simulations (Cα RMSD + RMSF)."
    )
    parser.add_argument("--config", required=True)
    parser.add_argument("--run",    type=int, default=None)
    args = parser.parse_args()

    cfg         = Config.from_yaml(args.config)
    results_dir = ensure_dir(Path(cfg.output_dir) / "results" / "equilibrium")
    run_ids     = [args.run] if args.run else list(range(1, cfg.n_runs + 1))

    summary_rows, rmsd_xvgs, rmsf_xvgs, labels = [], [], [], []

    for run_id in run_ids:
        analyser = EquilibriumAnalyser(cfg, run_id, results_dir)
        row = analyser.analyse()
        if row is None:
            continue
        summary_rows.append(row)
        rmsd_xvgs.append(analyser.rmsd_xvg)
        rmsf_xvgs.append(analyser.rmsf_xvg)
        labels.append(f"Run {run_id}")

    if not summary_rows:
        logger.error("No runs could be analysed.")
        return

    plot_rmsd(rmsd_xvgs, labels,
              out_png=results_dir / "rmsd_all_runs.png",
              title=f"{cfg.system_name} — Cα RMSD vs EM reference")
    plot_rmsf(rmsf_xvgs, labels,
              out_png=results_dir / "rmsf_all_runs.png",
              title=f"{cfg.system_name} — Cα RMSF per residue")
    write_summary_csv(summary_rows, results_dir / "summary.csv")

    print(f"\n{'Run':>5}  {'Mean RMSD (Å)':>14}  {'Max RMSD (Å)':>13}  "
          f"{'Mean RMSF (Å)':>14}  {'Max RMSF (Å)':>13}")
    print("-" * 66)
    for row in summary_rows:
        print(f"{row['run']:>5}  {row['mean_rmsd_A']:>14.3f}  {row['max_rmsd_A']:>13.3f}  "
              f"{row['mean_rmsf_A']:>14.3f}  {row['max_rmsf_A']:>13.3f}")
    print(f"\nResults written to: {results_dir}")


if __name__ == "__main__":
    main()
