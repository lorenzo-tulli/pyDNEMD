# pyDNEMD Workshop Container
# Provides: GROMACS, PyMOL (open-source), pyDNEMD, JupyterLab
#
# Build:   docker build -f docker/Dockerfile -t pydnemd-workshop .
# Pull:    docker pull ghcr.io/lorenzo-tulli/pydnemd-workshop:latest
# Run:     docker run -p 8888:8888 ghcr.io/lorenzo-tulli/pydnemd-workshop:latest

FROM continuumio/miniconda3:24.1.2-0

# ── System dependencies ────────────────────────────────────────────────────────
RUN apt-get update && apt-get install -y --no-install-recommends \
        git \
        build-essential \
        libgl1 \
        libglu1-mesa \
    && rm -rf /var/lib/apt/lists/*

# ── Conda environment ──────────────────────────────────────────────────────────
RUN conda install -y -c conda-forge \
        python=3.11 \
        gromacs \
        pymol-open-source \
        "numpy>=1.20" \
        "pyyaml>=5.1" \
        "mdanalysis>=2.0" \
        "matplotlib>=3.5" \
        pandas \
        jupyterlab \
        ipykernel \
    && conda clean --all --yes

# ── gmx_mpi symlink ───────────────────────────────────────────────────────────
# conda-forge installs 'gmx'; symlink so the default config (gmx: gmx_mpi) works.
RUN ln -s "$(which gmx)" /usr/local/bin/gmx_mpi

# ── Install pyDNEMD from GitHub ───────────────────────────────────────────────
RUN pip install --no-cache-dir git+https://github.com/lorenzo-tulli/pyDNEMD.git

# ── Workshop directory layout ──────────────────────────────────────────────────
# /workshop/notebook          — workshop.ipynb + config.yaml
# /workshop/templates         — MDP templates
# /workshop/examples/inputs   — bundled input files (topology, coordinates, index)
# /workshop/examples/output   — pre-computed simulation output (ready to analyse)
# No volume mounts needed — everything is self-contained.

WORKDIR /workshop

COPY workshop.ipynb          notebook/workshop.ipynb
COPY templates/              templates/
COPY examples/inputs/        examples/inputs/
COPY examples/output/        examples/output/

# Pre-fill config.yaml with container paths so participants can run immediately
COPY examples/config.yaml    notebook/config.yaml
RUN sed -i \
        -e 's|input_gro:.*|input_gro: /workshop/examples/inputs/solv_ions.gro|' \
        -e 's|topology:.*|topology:  /workshop/examples/inputs/topol.top|' \
        -e 's|index_ndx:.*|index_ndx: /workshop/examples/inputs/index.ndx|' \
        -e 's|mdp_dir:.*|mdp_dir:   /workshop/templates|' \
        -e 's|output_dir:.*|output_dir: /workshop/examples/output|' \
        -e 's|eq_sim_dir:.*|eq_sim_dir: /workshop/examples/output|' \
        -e 's|ne_sim_dir:.*|ne_sim_dir: /workshop/examples/output|' \
        -e 's|np_sim_dir:.*|np_sim_dir: /workshop/examples/output|' \
        -e 's|ne_index:.*|ne_index: /workshop/examples/output/perturbed_topology/perturbed_index.ndx|' \
        notebook/config.yaml

# ── JupyterLab port ───────────────────────────────────────────────────────────
EXPOSE 8888

CMD ["jupyter", "lab", \
     "--ip=0.0.0.0", \
     "--port=8888", \
     "--no-browser", \
     "--allow-root", \
     "--notebook-dir=/", \
     "--ServerApp.token="]
