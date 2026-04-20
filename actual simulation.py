from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout

from openff.toolkit.topology import Molecule
from openmmforcefields.generators import GAFFTemplateGenerator, SystemGenerator


ligand = Molecule.from_file("SEP Fixed.mol")
ligand.name = "LIG"
gaff = GAFFTemplateGenerator(molecules=ligand, forcefield='gaff-2.11')

psf = CharmmPsfFile('4Zn-8Ser pH.psf')
crd = CharmmCrdFile('4Zn-8Ser pH.crd')


forcefield = ForceField('amber19-all.xml', 'amber19/tip3pfb.xml')
forcefield.registerTemplateGenerator(gaff.generator)

psf.setBox(5*nanometer,5*nanometer,5*nanometer)

system = forcefield.createSystem(psf.topology, nonbondedMethod=PME,
        nonbondedCutoff=1*nanometer, constraints=HBonds)
integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.001*picoseconds)
simulation = Simulation(psf.topology, system, integrator)
# simulation.context.setPositions(crd.positions)
# simulation.minimizeEnergy()

simulation.loadState("4Zn-8Ser pH 1fs State.xml") # Load the state that was gotten from "equilibrating sys.py"
                                            # Have to manually "clean"/remove the motecarlo line from the xml


print("Starting Sim")

simulation.reporters.append(DCDReporter('Traj 4Zn-8Ser pH 1fs.dcd', 1000)) # Appreantly I'm supposed to write every 1000 steps to avoid the dcd file from being too big, we'll see
simulation.reporters.append(StateDataReporter(stdout, 1000, step=True,
        potentialEnergy=True, temperature=True, density=True, volume=True))

#simulation.step(5000) # Test that volume is constant and try to get a time estimate
simulation.step(10_000_000)
print("Finished :D")
