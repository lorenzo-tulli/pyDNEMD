"""
Build a hybrid (dual-topology) system for mutation-perturbation D-NEMD.

Merges a wild-type and mutant structure via BioSimSpace into a single
GROMACS free-energy topology (A-state/B-state atoms in one file), the
input EquilibrationPipeline needs for the mutation perturbation mode.

BioSimSpace is imported lazily, inside build_hybrid_topology(), so that
importing this module — or any dnemd CLI command — never requires it.
Only actually running dnemd-create-hybrid-topology does.
"""
import re
import shutil
from pathlib import Path
from dnemd.utils import ensure_dir, get_logger

logger = get_logger("hybrid_topology")

# BioSimSpace's FreeEnergy.Relative API requires num_lam >= 3 even though
# we only want a fixed lambda=0 topology (this is a non-equilibrium switch,
# not lambda-dynamics — see docs/design/mutation-perturbation.md). The
# lambda=0.5/1 windows it's forced to create get discarded afterwards.
N_LAMBDAS = 3


def build_hybrid_topology(cfg) -> Path:
    """
    Merge cfg.wt_gro/wt_topology and cfg.mutant_gro/mutant_topology into a
    hybrid topology, writing GROMACS input files to
    <output_dir>/hybrid_topology/.

    Runs a single BSS FreeEnergyMinimisation protocol with setup_only=True
    purely to get BSS's GROMACS free-energy topology writer to run once —
    nothing is actually simulated here; EquilibrationPipeline handles the
    real nvt1-npt2-production staging afterwards, using this output as its
    starting structure/topology.
    """
    try:
        import BioSimSpace as bss
    except ImportError as e:
        raise ImportError(
            "dnemd-create-hybrid-topology requires BioSimSpace, which isn't "
            "installed in this environment. See the README's 'Optional: "
            "protein mutation perturbations' section."
        ) from e

    out_dir = ensure_dir(Path(cfg.output_dir) / "hybrid_topology")

    logger.info("Reading wild-type system...")
    wt_system = bss.IO.readMolecules([cfg.wt_gro, cfg.wt_topology])

    logger.info("Reading mutant system...")
    mut_system = bss.IO.readMolecules([cfg.mutant_gro, cfg.mutant_topology])

    wt_protein = wt_system[0]
    mut_protein = mut_system[0]

    wt_residues = wt_protein.getResidues()
    mut_residues = mut_protein.getResidues()
    if len(wt_residues) != len(mut_residues):
        raise ValueError(
            f"WT and mutant structures have different residue counts "
            f"({len(wt_residues)} vs {len(mut_residues)}) — residue-by-residue "
            "mutation detection assumes matching numbering."
        )

    logger.info("Detecting mutated residues...")
    roi = []
    for i, res in enumerate(wt_residues):
        mut_res = mut_residues[i]
        if res.name() != mut_res.name():
            logger.info(f"  residue {i}: {res.name()} -> {mut_res.name()}")
            roi.append(i)
    if not roi:
        raise ValueError("No mutated residues detected between WT and mutant structures.")
    logger.info(f"Region of interest: {roi}")

    _warn_charge_and_volume_delta(wt_protein, mut_protein)

    logger.info("Mapping and aligning...")
    mapping = bss.Align.matchAtoms(
        molecule0=wt_protein, molecule1=mut_protein, roi=roi, complete_rings_only=True,
    )
    aligned_wt = bss.Align.rmsdAlign(
        molecule0=wt_protein, molecule1=mut_protein, roi=roi, mapping=mapping,
    )

    logger.info("Merging...")
    merged_protein = bss.Align.merge(
        molecule0=aligned_wt, molecule1=mut_protein, mapping=mapping, roi=roi,
    )

    merged_system = wt_system.copy()
    merged_system.removeMolecules(wt_protein)
    merged_system.addMolecules(merged_protein)

    logger.info("Writing GROMACS hybrid topology...")
    protocol = bss.Protocol.FreeEnergyMinimisation(num_lam=N_LAMBDAS)
    bss.FreeEnergy.Relative(
        system=merged_system,
        protocol=protocol,
        engine="GROMACS",
        work_dir=str(out_dir),
        setup_only=True,
        ignore_warnings=True,
    )

    kept = _keep_lambda_zero_only(out_dir)
    logger.info(f"Hybrid topology written to {kept or out_dir}")
    return kept or out_dir


def _warn_charge_and_volume_delta(wt_protein, mut_protein):
    """
    Informational check for the no-re-solvation constraint (see
    docs/design/mutation-perturbation.md): the NE run starts from an
    equilibrium frame solvated/ionised for the WT system. A net-charge
    change is handled gracefully by PME's neutralising background and
    isn't a stability concern; a heavy-atom-count change is a real,
    separate approximation (solvent cavity size around the mutated
    residue won't match the new equilibrium at t=0). Both are logged as
    warnings, not errors — this doesn't block anything.
    """
    wt_charge = wt_protein.charge().value()
    mut_charge = mut_protein.charge().value()
    if abs(wt_charge - mut_charge) > 1e-3:
        logger.warning(
            f"Net charge changes by {mut_charge - wt_charge:+.2f} e (WT -> mutant) "
            "— handled by PME's neutralising background, not expected to affect dynamics."
        )

    wt_heavy = sum(1 for a in wt_protein.getAtoms() if a.element().symbol() != "H")
    mut_heavy = sum(1 for a in mut_protein.getAtoms() if a.element().symbol() != "H")
    if wt_heavy != mut_heavy:
        logger.warning(
            f"Heavy-atom count changes by {mut_heavy - wt_heavy:+d} (WT -> mutant) "
            "— the NE run starts from a WT-solvated frame, so this mutation isn't "
            "fully volume-conservative. Known approximation of the no-re-solvation "
            "design, not an error."
        )


def _keep_lambda_zero_only(out_dir: Path) -> Path | None:
    """
    Discard the lambda=0.5 and lambda=1 output BSS is forced to create.

    Identifies the lambda=0 folder by reading each generated .mdp for
    init-lambda = 0, rather than assuming BSS's folder-naming convention
    (undocumented and not worth hardcoding a guess for). If no folder can
    be identified this way, everything is left in place for manual
    inspection rather than risking deleting the wrong one.
    """
    candidates = [d for d in out_dir.iterdir() if d.is_dir()]
    kept = None
    for d in candidates:
        for mdp in d.glob("*.mdp"):
            if re.search(r"init-lambda\s*=\s*0(\.0+)?\b", mdp.read_text()):
                kept = d
                break
        if kept:
            break

    if kept is None:
        logger.warning(
            "Could not identify the lambda=0 output folder by inspecting "
            "generated .mdp files — leaving all BSS output in place under "
            f"{out_dir} for manual inspection."
        )
        return None

    for d in candidates:
        if d != kept:
            shutil.rmtree(d)
    logger.info(f"Kept lambda=0 output: {kept.name}")
    return kept
