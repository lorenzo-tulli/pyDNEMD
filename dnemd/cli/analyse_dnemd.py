#!/usr/bin/env python3
"""
Analyse DNEMD results: aggregate Cα displacement vectors (NE - NP)
across runs and equilibrium time points, apply SE significance filter,
and write statistics + B-factor PDBs.

Outputs (written to output_dir/results/dnemd/):
    stat_<system>_<tp>ps_<N>SE.txt        — per-CA statistics table
    vec_norm_<system>_<tp>ps_<N>SE.pdb    — CA-only PDB with B-factors = |displacement|
    summary.csv                            — aggregated stats across all time points

Usage:
    # All time points defined in config
    python scripts/analyse_dnemd.py --config examples/config.yaml

    # Single time point
    python scripts/analyse_dnemd.py --config examples/config.yaml --time-point 1000

    # SLURM array (one task per time point)
    python scripts/analyse_dnemd.py --config examples/config.yaml --task-id $SLURM_ARRAY_TASK_ID
"""
import argparse
import csv
import sys
from pathlib import Path

import numpy as np

from dnemd.config import Config
from dnemd.dnemd_analysis import (
    collect_vectors,
    compute_statistics,
    write_stats_txt,
    write_bfactor_pdb,
)
from dnemd.utils import ensure_dir, get_logger

logger = get_logger("analyse_dnemd")


# ---------------------------------------------------------------------------
# Per-timepoint analysis
# ---------------------------------------------------------------------------

def analyse_timepoint(
    cfg: Config,
    time_point_ps: int,
    dnemd_base: Path,
    results_dir: Path,
    ref_pdb: Path,
    runs: range,
    time_range_ns: range,
) -> dict | None:
    """
    Full analysis for one time point.
    Returns a summary dict, or None if no valid frames were found.
    """
    logger.info(f"--- Time point: {time_point_ps} ps ---")

    try:
        vectors_arr = collect_vectors(
            base_dir=dnemd_base,
            time_point_ps=time_point_ps,
            runs=runs,
            time_range_ns=time_range_ns,
        )
    except RuntimeError as e:
        logger.warning(str(e))
        return None

    (
        avg_vectors,
        avg_disp,
        se,
        se_mag,
        adjusted_avg_disp,
        _adjusted_avg_vectors,
    ) = compute_statistics(vectors_arr, se_threshold=cfg.se_threshold)

    n_samples = len(vectors_arr)
    label     = f"{cfg.system_name}_{time_point_ps}ps_{cfg.se_threshold}SE"

    # -------------------------------------------------------- Stats file --
    txt_path = results_dir / f"stat_{label}.txt"
    write_stats_txt(
        path=txt_path,
        avg_vectors=avg_vectors,
        avg_disp=avg_disp,
        se=se,
        se_mag=se_mag,
        n_samples=n_samples,
    )
    logger.info(f"Stats written: {txt_path}")

    # ------------------------------------------------------- B-factor PDB --
    if not ref_pdb.exists():
        logger.warning(
            f"Reference PDB not found: {ref_pdb}\n"
            "Skipping B-factor PDB. Provide a 'closest_to_average.pdb' "
            "in the DNEMD base directory."
        )
        pdb_path = None
    else:
        pdb_path = results_dir / f"vec_norm_{label}.pdb"
        write_bfactor_pdb(
            ref_pdb=ref_pdb,
            adjusted_avg_disp=adjusted_avg_disp,
            out_pdb=pdb_path,
        )

    # ---------------------------------------------------- Summary stats --
    n_sig     = int((adjusted_avg_disp > 0).sum())
    n_total   = len(avg_disp)
    nonzero   = avg_disp[avg_disp > 0]
    mean_disp = float(nonzero.mean()) if nonzero.size else 0.0
    max_disp  = float(avg_disp.max())

    logger.info(
        f"Significant CA atoms : {n_sig}/{n_total} "
        f"| mean |disp| = {mean_disp:.3f} Å "
        f"| max |disp|  = {max_disp:.3f} Å"
    )

    return {
        "time_point_ps":    time_point_ps,
        "n_samples":        n_samples,
        "n_ca_total":       n_total,
        "n_ca_significant": n_sig,
        "pct_significant":  round(100 * n_sig / n_total, 1) if n_total else 0.0,
        "mean_disp_A":      round(mean_disp, 3),
        "max_disp_A":       round(max_disp,  3),
        "txt_file":         str(txt_path),
        "pdb_file":         str(pdb_path) if pdb_path else "N/A",
    }


# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------

def write_dnemd_summary_csv(rows: list[dict], out_csv: Path):
    if not rows:
        return
    with open(out_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    logger.info(f"Summary CSV written: {out_csv}")


def print_summary_table(rows: list[dict], system_name: str):
    print(f"\n{'='*74}")
    print(f"  DNEMD Analysis Summary — {system_name}")
    print(f"{'='*74}")
    print(
        f"{'Time (ps)':>10}  {'Samples':>8}  {'Sig. CA':>8}  "
        f"{'Total CA':>9}  {'Sig. (%)':>9}  {'Mean |d| (Å)':>13}  {'Max |d| (Å)':>12}"
    )
    print("-" * 74)
    for r in rows:
        print(
            f"{r['time_point_ps']:>10}  "
            f"{r['n_samples']:>8}  "
            f"{r['n_ca_significant']:>8}  "
            f"{r['n_ca_total']:>9}  "
            f"{r['pct_significant']:>8.1f}%  "
            f"{r['mean_disp_A']:>13.3f}  "
            f"{r['max_disp_A']:>12.3f}"
        )
    print(f"{'='*74}\n")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Analyse DNEMD Cα displacement vectors (NE - NP)."
    )
    parser.add_argument("--config",      required=True, help="Path to config.yaml")
    parser.add_argument(
        "--time-point", type=int, default=None,
        help="Single time point in ps to analyse (default: all from config)",
    )
    parser.add_argument(
        "--task-id", type=int, default=None,
        help="SLURM array task ID — index into config time_points_ps list",
    )
    parser.add_argument(
        "--dnemd-dir", default=None,
        help=(
            "Root directory containing TRJDUMP_EQ/NE/NP subdirectories. "
            "Defaults to output_dir from config."
        ),
    )
    parser.add_argument(
        "--ref-pdb", default=None,
        help=(
            "Path to reference PDB for B-factor mapping "
            "(closest-to-average structure). "
            "Defaults to <dnemd-dir>/closest_to_average.pdb"
        ),
    )
    parser.add_argument(
        "--runs", type=int, default=None,
        help="Number of independent runs (default: n_runs from config)",
    )
    parser.add_argument(
        "--start-ns", type=int, default=None,
        help="First equilibrium time point to include in ns (default: from config extract_ns_start)",
    )
    parser.add_argument(
        "--end-ns", type=int, default=None,
        help="Last equilibrium time point to include in ns (default: from config extract_ns_end)",
    )
    parser.add_argument(
        "--interval-ns", type=int, default=None,
        help="Interval between equilibrium time points in ns (default: from config extract_ns_interval)",
    )
    args = parser.parse_args()

    # ---------------------------------------------------------------- Config --
    cfg = Config.from_yaml(args.config)

    start_ns    = args.start_ns    if args.start_ns    is not None else cfg.extract_ns_start
    end_ns      = args.end_ns      if args.end_ns      is not None else cfg.extract_ns_end
    interval_ns = args.interval_ns if args.interval_ns is not None else cfg.extract_ns_interval

    # Resolve which time points to analyse
    if args.time_point is not None:
        time_points = [args.time_point]
    elif args.task_id is not None:
        if args.task_id >= len(cfg.time_points_ps):
            logger.error(
                f"--task-id {args.task_id} is out of range "
                f"(config has {len(cfg.time_points_ps)} time points)."
            )
            sys.exit(1)
        time_points = [cfg.time_points_ps[args.task_id]]
    else:
        time_points = cfg.time_points_ps

    # Resolve directories
    dnemd_base = Path(args.dnemd_dir) if args.dnemd_dir else Path(cfg.output_dir)
    results_dir = ensure_dir(Path(cfg.output_dir) / "results" / "dnemd")

    ref_pdb = (
        Path(args.ref_pdb)
        if args.ref_pdb
        else dnemd_base / "closest_to_average.pdb"
    )

    n_runs      = args.runs if args.runs else cfg.n_runs
    runs        = range(1, n_runs + 1)
    time_range_ns = range(start_ns, end_ns + 1, interval_ns)

    logger.info(f"System      : {cfg.system_name}")
    logger.info(f"DNEMD base  : {dnemd_base}")
    logger.info(f"Results dir : {results_dir}")
    logger.info(f"Ref PDB     : {ref_pdb}")
    logger.info(f"Runs        : {list(runs)}")
    logger.info(f"EQ range    : {start_ns}–{end_ns} ns every {interval_ns} ns")
    logger.info(f"Time points : {time_points} ps")
    logger.info(f"SE threshold: {cfg.se_threshold}")

    # ------------------------------------------------------------ Analysis --
    summary_rows = []

    for tp in time_points:
        row = analyse_timepoint(
            cfg=cfg,
            time_point_ps=tp,
            dnemd_base=dnemd_base,
            results_dir=results_dir,
            ref_pdb=ref_pdb,
            runs=runs,
            time_range_ns=time_range_ns,
        )
        if row is not None:
            summary_rows.append(row)

    if not summary_rows:
        logger.error(
            "No time points could be analysed. "
            "Check that NE/NP/EQ trajectories are present under:\n"
            f"  {dnemd_base}"
        )
        sys.exit(1)

    # ------------------------------------------------- Write CSV + print --
    write_dnemd_summary_csv(
        rows=summary_rows,
        out_csv=results_dir / "summary.csv",
    )
    print_summary_table(summary_rows, cfg.system_name)

    logger.info("analyse_dnemd.py finished.")


if __name__ == "__main__":
    main()
