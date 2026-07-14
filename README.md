# pyDNEMD

A Python package for running and analysing **Dynamical Non-Equilibrium MD (D-NEMD)** simulations with GROMACS. A step-by-step guide using a Docker container is also available following the tutorial in [pyDNEMD-workshop](https://github.com/lorenzo-tulli/pyDNEMD-workshop.git).

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
| `dnemd-create-hybrid-topology` | Build a hybrid topology for mutation perturbations (requires BioSimSpace) |
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
Solvate your favourite Protein-ligand complex using gromacs, then:

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

## Mutation perturbation quick start

The steps above default to `perturbation: ligand_removal`. For a mutation perturbation,
set `perturbation: mutation` in your config, along with `wt_gro`/`wt_topology`/
`mutant_gro`/`mutant_topology` and `mdp_dir: templates/protein_mutation` (requires
BioSimSpace — see "Optional: protein mutation perturbations" above). There's one extra
step first, and `dnemd-create-ne-np`/`dnemd-run-ne` behave differently under the hood
(the NE leg runs as two chained simulations — a fast switch, then the response phase —
instead of one), but the commands themselves are the same:

```bash
cp examples/config.yaml my_mutation_config.yaml
# edit: perturbation: mutation, wt_gro/wt_topology, mutant_gro/mutant_topology,
#       mdp_dir: templates/protein_mutation

# 0. Build the hybrid topology (mutation-only step)
dnemd-create-hybrid-topology --config my_mutation_config.yaml
# then point input_gro/topology in my_mutation_config.yaml at the .gro/.top
# files it wrote to <output_dir>/hybrid_topology/

# 1-6. Same as the ligand-removal workflow above
dnemd-equilibrium --config my_mutation_config.yaml
dnemd-analyse-equilibrium --config my_mutation_config.yaml
dnemd-create-ne-np --config my_mutation_config.yaml --start 50000 --frequency 5000
dnemd-run-ne --config my_mutation_config.yaml
dnemd-run-np --config my_mutation_config.yaml
dnemd-extract --config my_mutation_config.yaml --all
dnemd-analyse-dnemd --config my_mutation_config.yaml
```

`dnemd-run-ne` runs the switch phase, then automatically builds and runs the response
phase from the switch phase's checkpoint — nothing extra to invoke. `dnemd-extract`/
`dnemd-analyse-dnemd` don't need to know which mode was used; they only ever look at the
final NE/NP trajectories, which come out under the same file names in both modes.

## Python API

If you want to embed this pipeline inside a larger Python workflow, add custom logic between steps, or automate things programmatically, you could follow this usage example.

```python
from dnemd.config import Config
from dnemd.equilibration import EquilibrationPipeline

cfg = Config.from_yaml("my_config.yaml")
EquilibrationPipeline(cfg, run_id=1).run_all()
```
