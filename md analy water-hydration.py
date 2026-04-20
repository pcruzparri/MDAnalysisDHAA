import MDAnalysis as mda
import matplotlib.pyplot as plt
import os
import numpy as np
import math
import matplotlib.pyplot as plt


script_dir = os.path.dirname(__file__)

pdb_path = os.path.join(script_dir, '4Zn-8Ser pH', 'step4_pbcsetup.pdb')
dcd_path = os.path.join(script_dir, "Traj 4Zn-8Ser pH 1fs.dcd")

u = mda.Universe(pdb_path, dcd_path)

zincs = u.select_atoms('name ZN')
print(len(zincs))

center_atom = zincs[0]

radius = 2.08
oxy_hydrated = u.select_atoms(f"name OH2 and around {radius} index {center_atom.index}")



i=1
for zinc in zincs:
    y=[]
    for frame in u.trajectory[0:10001]:
        y.append(len(u.select_atoms(f"name OH2 and around {radius} index {zinc.index}")))
    x = list(range(1,10001)) #time
    print(f'y:{len(y)}, x:{len(x)}')
    plt.scatter(x,y,label=f"Zn #{i}")
    plt.xlabel("Time (ps)")
    plt.ylabel("Num O")
    plt.show()
    i+=1

plt.xlabel("Time (ps)")
plt.ylabel("Num O")
plt.show()
