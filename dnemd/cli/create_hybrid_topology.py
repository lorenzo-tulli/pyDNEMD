#!/usr/bin/env python3

"""
Build a hybrid (dual-topology) system for mutation-perturbation D-NEMD.

Reads wt_gro/wt_topology and mutant_gro/mutant_topology from the config,
merges them via BioSimSpace, and writes a GROMACS-format hybrid topology
to <output_dir>/hybrid_topology/. Requires BioSimSpace in the current
environment — see the README's "Optional: protein mutation perturbations"
section for install instructions.

After this command finishes, point input_gro/topology in your config at
the .gro/.top files it writes, and mdp_dir at templates/protein_mutation/,
then run dnemd-equilibrium as usual.

Usage:
    dnemd-create-hybrid-topology --config config.yaml
"""

import argparse
from dnemd.config import Config
from dnemd.hybrid_topology import build_hybrid_topology
from dnemd.utils import get_logger

logger = get_logger("create_hybrid_topology")

REQUIRED_FIELDS = ("wt_gro", "wt_topology", "mutant_gro", "mutant_topology")


def main():
    parser = argparse.ArgumentParser(
        description="Build a hybrid topology for mutation perturbations."
    )
    parser.add_argument("--config", required=True, help="Path to config.yaml")
    args = parser.parse_args()

    cfg = Config.from_yaml(args.config)

    missing = [f for f in REQUIRED_FIELDS if not getattr(cfg, f)]
    if missing:
        logger.error(
            "Missing required config fields for hybrid topology generation: "
            + ", ".join(missing)
        )
        raise SystemExit(1)

    build_hybrid_topology(cfg)


if __name__ == "__main__":
    main()
