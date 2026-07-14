import BioSimSpace as bss
import os

input_directory = "/group/chem/oliveira2t/vi24769/dnemd-nes-pipeline/inputs"
outputs = "/group/chem/oliveira2t/vi24769/dnemd-nes-pipeline"
os.makedirs(outputs, exist_ok=True)

n_lambdas = 3 

min_protocol = bss.Protocol.FreeEnergyMinimisation(num_lam=n_lambdas)

nvt1_protocol = bss.Protocol.FreeEnergyEquilibration(
    num_lam=n_lambdas,
    pressure=None,
    runtime=1*bss.Units.Time.nanosecond,
    temperature=250*bss.Units.Temperature.kelvin,
    restraint="heavy",
    force_constant=10
)

nvt2_protocol = bss.Protocol.FreeEnergyEquilibration(
    num_lam=n_lambdas,
    pressure=None,
    runtime=2*bss.Units.Time.nanosecond,
    temperature=310*bss.Units.Temperature.kelvin,
    restraint="heavy",
    force_constant=5
)

npt1_protocol = bss.Protocol.FreeEnergyEquilibration(
    num_lam=n_lambdas,
    runtime=2.5*bss.Units.Time.nanosecond,
    pressure=1*bss.Units.Pressure.atm,
    restraint="backbone",
    force_constant=2
)

npt2_protocol = bss.Protocol.FreeEnergyEquilibration(
    num_lam=n_lambdas,
    runtime=5*bss.Units.Time.nanosecond,
    pressure=1*bss.Units.Pressure.atm,
)

eq_protocol = bss.Protocol.FreeEnergyProduction(
    num_lam=n_lambdas,
    runtime=500*bss.Units.Time.nanosecond,
    temperature=310*bss.Units.Temperature.kelvin,
    pressure=1*bss.Units.Pressure.atm
)

print(f"reading in wt system")

wt_system = bss.IO.readMolecules([f"{input_directory}/bound/state_A/WT-DIMER_Mg-ATP_solvated.gro",
                                  f"{input_directory}/bound/state_A/topol.top"])

print(f"reading in mut system")

mut_system = bss.IO.readMolecules([f"{input_directory}/bound/state_B/R858-DIMER_Mg-ATP_solvated.gro",
                                   f"{input_directory}/bound/state_B/topol.top"])

print("extracting the proteins")

wt_protein = wt_system[0]
mut_protein = mut_system[0]

print("detecting the different residue between the two proteins")

roi = []
for i, res in enumerate(wt_protein.getResidues()):
    if res.name() != mut_protein.getResidues()[i].name():
        print(f"Found mutation at residue {i}: {res.name()} -> {mut_protein.getResidues()[i].name()}")
        roi.append(i)
print(f"Region of interest (mutated residues): {roi}")

print("creating mapping")
mapping = bss.Align.matchAtoms(
        molecule0=wt_protein,
        molecule1=mut_protein,
        roi=roi,
        complete_rings_only=True
)

print("aligning")
aligned_wt_protein = bss.Align.rmsdAlign(
        molecule0=wt_protein,
        molecule1=mut_protein,
        roi=roi,
        mapping=mapping
)

print("merging")
merged_protein = bss.Align.merge(
        molecule0=aligned_wt_protein,
        molecule1=mut_protein,
        mapping=mapping,
        roi=roi
)

print("updating system with merged protein")

merged_system = wt_system.copy()
merged_system.removeMolecules(wt_protein)
merged_system.addMolecules(merged_protein)

print("setting up bound states")

bss.FreeEnergy.Relative(
    system=merged_system,
    protocol=min_protocol,
    engine="GROMACS",
    work_dir=outputs + "/bound/min/",
    setup_only=True,
    ignore_warnings=True
)

bss.FreeEnergy.Relative(
    system=merged_system,
    protocol=nvt1_protocol,
    engine="GROMACS",
    work_dir=outputs + "/bound/nvt1/",
    setup_only=True,
    ignore_warnings=True
)

bss.FreeEnergy.Relative(
    system=merged_system,
    protocol=nvt2_protocol,
    engine="GROMACS",
    work_dir=outputs + "/bound/nvt2/",
    setup_only=True,
    ignore_warnings=True
)

bss.FreeEnergy.Relative(
    system=merged_system,
    protocol=npt1_protocol,
    engine="GROMACS",
    work_dir=outputs + "/bound/npt1/",
    setup_only=True,
    ignore_warnings=True
)

bss.FreeEnergy.Relative(
    system=merged_system,
    protocol=npt2_protocol,
    engine="GROMACS",
    work_dir=outputs + "/bound/npt2/",
    setup_only=True,
    ignore_warnings=True
)

bss.FreeEnergy.Relative(
    system=merged_system,
    protocol=eq_protocol,
    engine="GROMACS",
    work_dir=outputs + "/bound/",
    setup_only=True,
    ignore_warnings=True
)
