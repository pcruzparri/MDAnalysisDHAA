from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout

from openff.toolkit.topology import Molecule
from openmmforcefields.generators import GAFFTemplateGenerator, SystemGenerator


ligand = Molecule.from_file("SEP Fixed.mol") # Some how use the .Mol2 file for our "ligand" in here, nop
                                          # Rename this to change file
ligand.name = "LIG"
gaff = GAFFTemplateGenerator(molecules=ligand, forcefield='gaff-2.11')

psf = CharmmPsfFile('12Zn-8Ser pH.psf') # Rename these to change file
crd = CharmmCrdFile('12Zn-8Ser pH.crd') # Rename these to change file

#'''
forcefield = ForceField('amber19-all.xml', 'amber19/tip3pfb.xml')
forcefield.registerTemplateGenerator(gaff.generator)

psf.setBox(5*nanometer,5*nanometer,5*nanometer)


'''
for atom in list(psf.topology.residues())[0].atoms():
    print(f"PSF Atom: {atom.name}, Element: {atom.element.symbol}")

# Print the atoms in your OpenFF Molecule
for atom in ligand.atoms:
    print(f"MOL Atom: {atom.name}, Element: {atom.symbol}")
'''
    
system = forcefield.createSystem(psf.topology, nonbondedMethod=PME,
        nonbondedCutoff=1*nanometer, constraints=HBonds)


'''
system_generator = SystemGenerator(
    forcefields=['amber19-all.xml', 'amber19/tip3pfb.xml'],
    small_molecule_forcefield='gaff-2.11',
    molecules=[ligand],
    periodic_forcefield_kwargs={'nonbondedMethod': PME, 'nonbondedCutoff': 1*nanometer, 'constraints': HBonds}
)
psf.setBox(5*nanometer, 5*nanometer, 5*nanometer)
system = system_generator.create_system(psf.topology)
'''

system.addForce(MonteCarloBarostat(1*bar, 300*kelvin)) # IDK about the frequency that this should be checked -> so I'm leaving this as default frequency=25
integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.001*picoseconds) # 2fs is standard? I think so..., the paper mentioned 1fs, but I don't really have the time
simulation = Simulation(psf.topology, system, integrator)
simulation.context.setPositions(crd.positions)
simulation.minimizeEnergy()

# simulation.reporters.append(DCDReporter('output_long.dcd', 100))
# simulation.reporters.append(StateDataReporter(stdout, 100, step=True,
#        potentialEnergy=True, temperature=True)) # Do you want a print statement?... YES I DO
# simulation.step(630000)

simulation.reporters.append(StateDataReporter(stdout, 1000, step=True, # This 1000 is just an arbitrary number, I just want to see that the sim is working
                                              potentialEnergy=True, temperature=True, density=True, volume=True))

print("Starting Eq Step 1: 100ps const temp") # 100ps is equal to 100,000fs so a total of 50k steps
                                              # Around 20sec per 1000 steps
                                              # Actually wrong, I wasn't using GPU because CUDA version mis-match, now it is ~1sec per 1000 steps, give or take the same for the other ones
simulation.step(100_000)

print("Starting Eq Step 2: 100ps System compression 100MPa") # Around 24sec per 1000 steps; Full run 20:33, avg checks out give or take
#simulation.addforce(MonteCarloBarostat(1*bar,300*kelvin)) #  This code is wrong, check above I added barostat earlier so I don't have to reparamatrize the simulation
simulation.context.setParameter(MonteCarloBarostat.Pressure(), 1000*bar)
simulation.step(100_000)

print("Starting Eq Step 3: 100ps System relaxation 0.1MPa (~1atm)") # Around 22sec per 1000 steps; Full run 20:54, avg is 25sec
simulation.context.setParameter(MonteCarloBarostat.Pressure(), 1*bar)
simulation.step(100_000)


# Saving equilibrated system to a state file, see other python file to run actual simulation
#simulation.saveState(f"{input("Name this equilibrated state: ")}.xml")
simulation.saveState("12Zn-8Ser pH 1fs State.xml")

