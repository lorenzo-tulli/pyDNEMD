"""
DNEMD vector analysis — aggregates Cα displacement vectors (NE - NP)
across runs and equilibrium time points.
"""
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import align
from pathlib import Path
from modules.utils import get_logger, ensure_dir

logger = get_logger("dnemd_analysis")


def load_and_align_ca(eq_path: str, ne_path: str, np_path: str) -> np.ndarray | None:
    try:
        EQ = mda.Universe(eq_path)
        NE = mda.Universe(ne_path)
        NP = mda.Universe(np_path)
        align.alignto(NE, EQ, select="name CA")
        align.alignto(NP, EQ, select="name CA")
        return NE.select_atoms("name CA").positions - NP.select_atoms("name CA").positions
    except Exception as e:
        logger.warning(f"Skipped EQ={eq_path}: {e}")
        return None


def collect_vectors(
    base_dir: str | Path,
    time_point_ps: int,
    runs: range = range(1, 6),
    time_range_ns: range = range(50, 500, 5),
) -> np.ndarray:
    """
    Collect all (NE - NP) Cα displacement vectors for one time point.
    Directory layout mirrors your existing TRJDUMP structure:
        base_dir/TRJDUMP_EQ/EQ_<run>/<time>ns/Run<run>EQ<time>ns<tp>ps.gro
        base_dir/TRJDUMP_NE/NE_<run>/<time>ns/Run<run>NE<time>ns<tp>ps.gro
        base_dir/TRJDUMP_NP/NP_<run>/<time>ns/Run<run>NP<time>ns<tp>ps.gro
    """
    base_dir = Path(base_dir)
    vectors = []

    for run in runs:
        for t in time_range_ns:
            eq = base_dir / f"TRJDUMP_EQ/EQ_{run}/{t}ns/Run{run}EQ{t}ns{time_point_ps}ps.gro"
            ne = base_dir / f"TRJDUMP_NE/NE_{run}/{t}ns/Run{run}NE{t}ns{time_point_ps}ps.gro"
            np_ = base_dir / f"TRJDUMP_NP/NP_{run}/{t}ns/Run{run}NP{t}ns{time_point_ps}ps.gro"

            diff = load_and_align_ca(str(eq), str(ne), str(np_))
            if diff is not None:
                vectors.append(diff)

    if not vectors:
        raise RuntimeError(f"No valid frames found at {time_point_ps} ps under {base_dir}")

    return np.array(vectors)   # (N_samples, N_CA, 3)


def compute_statistics(vectors_arr: np.ndarray):
    avg_vectors  = np.mean(vectors_arr, axis=0)          # (N_CA, 3)
    avg_disp     = np.linalg.norm(avg_vectors, axis=1)   # (N_CA,)
    sd           = np.std(vectors_arr, axis=0)
    se           = sd / np.sqrt(len(vectors_arr))

    ax, ay, az = avg_vectors[:, 0], avg_vectors[:, 1], avg_vectors[:, 2]
    denom    = np.where(avg_disp > 0, avg_disp, 1e-12)
    se_mag   = (1.0 / denom) * np.sqrt(
        ax**2 * se[:, 0]**2 + ay**2 * se[:, 1]**2 + az**2 * se[:, 2]**2
    )

    adjusted_avg_disp    = np.where(avg_disp - se_mag > 0, avg_disp, 0.0)
    adjusted_avg_vectors = np.where(np.abs(avg_vectors) - se > 0, avg_vectors, 0.0)

    return avg_vectors, avg_disp, se, se_mag, adjusted_avg_disp, adjusted_avg_vectors


def write_stats_txt(path: str | Path, avg_vectors, avg_disp, se, se_mag, n_samples):
    with open(path, "w") as f:
        f.write("#" * 128 + "\n")
        f.write("#                 Average Displacement                Sample Size"
                "                                   SE (1SD)                SE (1SD) of"
                "              SE (2SD)                SE (2SD) of       #\n")
        f.write("#              --------------------------       ---------------------"
                "      Average        --------------------------          Average"
                "        --------------------------          Average       #\n")
        f.write("#       CA     x-axis    y-axis    z-axis       x         y         z"
                "    Displacement     x-axis    y-axis    z-axis    Displacement vector"
                "  x-axis    y-axis    z-axis    Displacement vector #\n")
        f.write("#" * 128 + "\n")

        for i in range(len(avg_vectors)):
            f.write(
                f"{i+1:>9} {avg_vectors[i,0]:>8.3f} {avg_vectors[i,1]:>8.3f} "
                f"{avg_vectors[i,2]:>8.3f} "
                f"{n_samples:>8} {n_samples:>8} {n_samples:>8} {avg_disp[i]:>12.3f} "
                f"{se[i,0]:>10.3f} {se[i,1]:>10.3f} {se[i,2]:>10.3f} {se_mag[i]:>12.3f} "
                f"{2*se[i,0]:>10.3f} {2*se[i,1]:>10.3f} {2*se[i,2]:>10.3f} "
                f"{2*se_mag[i]:>12.3f}\n"
            )


def write_bfactor_pdb(ref_pdb: str | Path, adjusted_avg_disp: np.ndarray,
                      out_pdb: str | Path):
    u = mda.Universe(str(ref_pdb))
    ca = u.select_atoms("name CA")
    ca_u = mda.Merge(ca)
    ca_u.add_TopologyAttr("tempfactors")
    for atom, val in zip(ca_u.atoms, adjusted_avg_disp):
        atom.tempfactor = float(val)
    ca_u.atoms.write(str(out_pdb))
    logger.info(f"B-factor PDB written: {out_pdb}")
