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
