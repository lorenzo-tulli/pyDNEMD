#!/usr/bin/env python3

"""
Set up and optionally run the full GROMACS equilibration pipeline:
  em -> step1 (NVT heavy) -> step2 (NVT CA) -> step3 (NPT CA)
     -> step4 (NPT backbone) -> production

Usage:
    python scripts/run_equilibrium.py --config examples/config.yaml
    python scripts/run_equilibrium.py --config examples/config.yaml --run 2
    python scripts/run_equilibrium.py --config examples/config.yaml --setup-only
"""

import argparse
from dnemd.config import Config
from dnemd.equilibration import EquilibrationPipeline

def main():
    parser = argparse.ArgumentParser(description="Run GROMACS equilibration pipeline.")
    parser.add_argument("--config",     required=True)
    parser.add_argument("--run",        type=int, default=None)
    parser.add_argument("--setup-only", action="store_true")
    args = parser.parse_args()

    cfg     = Config.from_yaml(args.config)
    run_ids = [args.run] if args.run else list(range(1, cfg.n_runs + 1))

    for run_id in run_ids:
        EquilibrationPipeline(cfg, run_id, setup_only=args.setup_only).run_all()

if __name__ == "__main__":
    main()
