# pyDNEMD

A Python package for running and analysing **Dynamical Non-Equilibrium MD (D-NEMD)** simulations with GROMACS.

## Installation

create and activate a new environment.
```bash
mamba create -n pyDNEMD python=3.11
mamba activate pyDNEMD
```
Install the package.
```bash
pip install git+https://github.com/lorenzo-tulli/pyDNEMD.git
```

## Requirements

- Python ≥ 3.10
- GROMACS (`gmx` or `gmx_mpi` on your `PATH`)
- Python dependencies installed automatically: `numpy`, `pyyaml`, `MDAnalysis`, `matplotlib`

### Optional: protein mutation perturbations

The ligand-removal pipeline above needs nothing extra. If you also want to generate hybrid
topologies for **mutation** perturbations (`dnemd-create-hybrid-topology`), you need
[BioSimSpace](https://biosimspace.openbiosim.org/) in the same environment. BioSimSpace is
conda-only (no PyPI package) and only supports Linux and macOS — install it *before* pyDNEMD,
in the same environment:

```bash
conda install -n pyDNEMD -c openbiosim biosimspace
pip install git+https://github.com/lorenzo-tulli/pyDNEMD.git
```

Installing pyDNEMD after BioSimSpace (not before) avoids pip's resolver trying to touch
packages conda already manages. pyDNEMD's own dependencies (`numpy`, `MDAnalysis`,
`matplotlib`, `pyyaml`) are near-certainly already satisfied by any recent BioSimSpace
environment — worth a quick check afterwards:

```bash
python -c "import MDAnalysis; print(MDAnalysis.__version__)"  # should be >= 2.0
```

If `pip install` ever complains about those four packages specifically in a BioSimSpace
environment, `pip install --no-deps git+https://github.com/lorenzo-tulli/pyDNEMD.git` is a
safe escape hatch.

If you only need the ligand-removal pipeline, skip this section entirely — BioSimSpace is
never imported unless you actually run `dnemd-create-hybrid-topology`.

## Workflow overview

1. **Equilibration** — run independent equilibrium replicates
2. **NE/NP setup** — create non-equilibrium and null-perturbation input files
3. **Frame extraction** — dump trajectory frames for analysis
4. **Analysis** — compute Cα displacement vectors and generate statistics/B-factor PDBs

## Command-line tools

After installation the following commands are available:

| Command | Description |
|---|---|
| `dnemd-equilibrium` | Run GROMACS equilibration pipeline |
| `dnemd-create-ne-np` | Create NE/NP input files |
| `dnemd-run-ne` | Run non-equilibrium simulations |
| `dnemd-run-np` | Run null perturbation simulations |
| `dnemd-extract` | Extract EQ/NE/NP frames from trajectories |
| `dnemd-analyse-equilibrium` | Compute Cα RMSD and RMSF |
| `dnemd-analyse-dnemd` | Compute Cα displacement vectors (NE − NP) |

All commands require a `--config` argument pointing to a YAML configuration file.
See `examples/config.yaml` for a template.

## Quick start

```bash
# Copy and edit the config template
cp examples/config.yaml my_config.yaml

# 1. Run equilibration (all replicates)
dnemd-equilibrium --config my_config.yaml

# 2. Analyse equilibrium
dnemd-analyse-equilibrium --config my_config.yaml

# 3. Create NE/NP input files
dnemd-create-ne-np --config my_config.yaml --start 50000 --frequency 5000

# 4. Run NE/NP simulations
dnemd-run-ne --config my_config.yaml
dnemd-run-np --config my_config.yaml

# 5. Extract frames (run after NE/NP simulations finish)
dnemd-extract --config my_config.yaml --all

# 6. Analyse D-NEMD results
dnemd-analyse-dnemd --config my_config.yaml
```

## Python API

If you want to embed this pipeline inside a larger Python workflow, add custom logic between steps, or automate things programmatically, you could follow this usage example.

```python
from dnemd.config import Config
from dnemd.equilibration import EquilibrationPipeline

cfg = Config.from_yaml("my_config.yaml")
EquilibrationPipeline(cfg, run_id=1).run_all()
```
