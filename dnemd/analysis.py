"""
Parse GROMACS XVG files and plot RMSD / RMSF.
"""
import numpy as np
from pathlib import Path
from dnemd.utils import get_logger

logger = get_logger("analysis")


def parse_xvg(path: str | Path) -> tuple[np.ndarray, np.ndarray]:
    """Return (x, y) arrays from a GROMACS XVG file, skipping header lines."""
    x, y = [], []
    with open(path) as f:
        for line in f:
            if line.startswith(("#", "@")):
                continue
            cols = line.split()
            if len(cols) >= 2:
                x.append(float(cols[0]))
                y.append(float(cols[1]))
    return np.array(x), np.array(y)


def plot_rmsd(xvg_paths: list[str | Path], labels: list[str],
              out_png: str | Path, title: str = "Cα RMSD"):
    """Plot RMSD curves for multiple runs on one figure."""
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        logger.warning("matplotlib not installed — skipping RMSD plot.")
        return

    fig, ax = plt.subplots(figsize=(8, 4))
    for path, label in zip(xvg_paths, labels):
        t, r = parse_xvg(path)
        ax.plot(t, r * 10, label=label)   # nm -> Å

    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("RMSD (Å)")
    ax.set_title(title)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(str(out_png), dpi=150)
    plt.close(fig)
    logger.info(f"RMSD plot saved: {out_png}")


def plot_rmsf(xvg_paths: list[str | Path], labels: list[str],
              out_png: str | Path, title: str = "Cα RMSF per residue"):
    """Plot RMSF per residue for multiple runs."""
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        logger.warning("matplotlib not installed — skipping RMSF plot.")
        return

    fig, ax = plt.subplots(figsize=(10, 4))
    for path, label in zip(xvg_paths, labels):
        res, r = parse_xvg(path)
        ax.plot(res, r * 10, label=label)   # nm -> Å

    ax.set_xlabel("Residue")
    ax.set_ylabel("RMSF (Å)")
    ax.set_title(title)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(str(out_png), dpi=150)
    plt.close(fig)
    logger.info(f"RMSF plot saved: {out_png}")


def write_summary_csv(runs_data: list[dict], out_csv: str | Path):
    """Write a CSV table of per-run mean RMSD and mean RMSF."""
    import csv
    with open(out_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["run", "mean_rmsd_A", "max_rmsd_A",
                                                "mean_rmsf_A", "max_rmsf_A"])
        writer.writeheader()
        writer.writerows(runs_data)
    logger.info(f"Summary CSV written: {out_csv}")
