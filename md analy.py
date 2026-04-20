import MDAnalysis as mda
import MDAnalysis.analysis.rdf as rdf
import os
import matplotlib.pyplot as plt
import math
import numpy as np

script_dir = os.path.dirname(__file__)

pdb_path = os.path.join(script_dir, '4Zn-8Ser pH', 'step4_pbcsetup.pdb')
dcd_path = os.path.join(script_dir, "Traj 4Zn-8Ser pH 1fs.dcd")

u = mda.Universe(pdb_path, dcd_path)
print(u)
print(len(u.trajectory))

print(u.select_atoms('name ZN'))
zinc = u.select_atoms('name ZN') # All of our Zn2+ ions

print(u.select_atoms('segid HAA*'))

i=1
ligands = [] # All of our SEP
while len(u.select_atoms(f'segid HAA{i}')) != 0:
    ligands.append(u.select_atoms(f'segid HAA{i}'))
    i += 1

print(ligands)

def moving_average(data, window_size):
    cumsum = np.cumsum([0] + data)
    moving_avg = (cumsum[window_size:] - cumsum[:-window_size]) / window_size
    return moving_avg

window = 100
colors = ['blue','green','pink','red']
i=1 #for labeling
for ligand in ligands:
    #print(ligand.center_of_mass())
    j=1
    if i!=6: # HAA6 is the best one; shows most interaction with Zn2+ ions
        i+=1
        continue
    for zn_mol in zinc:
        #print(zn_mol.position)
        y = []
        for frame in u.trajectory[0:10001]:
            y.append(np.linalg.norm(ligand.center_of_mass()-zn_mol.position)) 
        y = moving_average(y, window)
        x = np.linspace(0, 10, len(y))    
        plt.plot(x, y, label=f"HAA{i} with Zn #{j}", color=colors[j-1], linewidth=2)
        j+=1
    i+=1    
    plt.xlabel("time(ns)")
    plt.ylabel("distance(A)") 
    plt.legend()
    plt.show()

plt.show()


