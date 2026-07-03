from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional
import yaml


@dataclass
class Config:
    # Input files
    input_gro:   str = "solv_ions.gro"
    topology:    str = "topol.top"
    index_ndx:   str = "index.ndx"

    # Labels and executables
    system_name: str = "SYSTEM"
    gmx:         str = "gmx_mpi"
    mdp_dir:     str = "templates"
    output_dir:  str = "output"

    # Mutation perturbation: hybrid topology generation (dnemd-create-hybrid-topology)
    wt_gro:          str = ""
    wt_topology:     str = ""
    mutant_gro:      str = ""
    mutant_topology: str = ""

    # Equilibrium
    n_runs: int = 5

    # Extraction defaults (overridden by CLI)
    extract_start_ps:     int = 50000
    extract_end_ps:       int = 500000
    extract_frequency_ps: int = 5000

    # Paths for frame extraction (used by dnemd-extract)
    eq_sim_dir:          str = ""
    ne_sim_dir:          str = ""
    np_sim_dir:          str = ""
    ne_index:            str = ""
    extract_ns_start:    int = 50
    extract_ns_end:      int = 250
    extract_ns_interval: int = 5

    # DNEMD analysis
    time_points_ps: list = field(default_factory=lambda: [0, 10, 100, 1000, 5000])
    se_threshold:   int  = 1

    @classmethod
    def from_yaml(cls, path: str) -> "Config":
        with open(path) as f:
            data = yaml.safe_load(f)
        fields = {k: v for k, v in data.items() if k in cls.__dataclass_fields__}
        return cls(**fields)

    def validate(self):
        errors = []
        for attr, label in [
            ("input_gro", "input_gro"),
            ("topology",  "topology"),
            ("index_ndx", "index_ndx"),
        ]:
            p = Path(getattr(self, attr))
            if not p.exists():
                errors.append(f"  '{label}' not found: {p}")
        if not Path(self.mdp_dir).exists():
            errors.append(f"  'mdp_dir' not found: {self.mdp_dir}")
        if errors:
            raise FileNotFoundError("Config validation failed:\n" + "\n".join(errors))
