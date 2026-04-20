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

# This proves that the auto select for angles won't work
# The phi_selection() is reuturning None :(
for res in u.residues[0:8]:
    phi = res.phi_selection()
    if phi is None:
        names = None
    else:
        names = phi.names
    print('{}: {} '.format(res.resname, names))

print(u.residues[0])

# Really helpful visulizer https://proteopedia.org/w/Tutorial:Ramachandran_principle_and_phi_psi_angles

# C3 is alpha carbon, C2 is carboxyl carbon, only for our case
c_alpha = u.select_atoms('name C3')
c_carboxyl = u.select_atoms('name C2')
c_nitrogen = u.select_atoms('name N2')
carboxyl_meas = u.select_atoms('name C4') #For Phi
nitrogen_meas = u.select_atoms('name N1') #For Psi

print(type(c_nitrogen.positions[0]))

print(type(c_alpha))
print(len(c_carboxyl))



'''
print()
#example for the first lig
#Phi - Nitrogen C-alpha
alpha_nitro_vect = c_nitrogen.positions[0] - c_alpha.positions[0] #normal for phi, ref end for psi
#alpha_nitro_vect /= np.linalg.norm(alpha_nitro_vect)
nitro_carbo_vect = carboxyl_meas.positions[0] - c_nitrogen.positions[0] #moving end phi
alpha_carbo_vect = c_carboxyl.positions[0] - c_alpha.positions[0] #ref end phi, normal for psi
carbo_nitro_vect = nitrogen_meas.positions[0] - c_carboxyl.positions[0] #moving end psi

norm_moving_phi = nitro_carbo_vect - (np.dot(nitro_carbo_vect, alpha_nitro_vect))/((np.linalg.norm(alpha_nitro_vect)**2))*alpha_nitro_vect
print(norm_moving_phi)
print(np.linalg.norm(norm_moving_phi))
norm_ref_phi = alpha_carbo_vect - (np.dot(alpha_carbo_vect, alpha_nitro_vect))/((np.linalg.norm(alpha_nitro_vect)**2))*alpha_nitro_vect

phi_angle = math.acos(np.dot(norm_moving_phi, norm_ref_phi)/(np.linalg.norm(norm_moving_phi)*np.linalg.norm(norm_ref_phi)))*180/math.pi
print(phi_angle)
#pos or neg
print(np.cross(norm_moving_phi, norm_ref_phi))
print(alpha_nitro_vect)

a = np.cross(norm_moving_phi, norm_ref_phi)
b = alpha_nitro_vect
print(f'{a[0]/b[0]} {a[1]/b[1]} {a[2]/b[2]}')

if np.cross(norm_moving_phi, norm_ref_phi)[0]/alpha_nitro_vect[0] < 0:
    phi_angle *= -1

print(phi_angle)

#Psi angles
norm_moving_psi = carbo_nitro_vect - (np.dot(carbo_nitro_vect, alpha_carbo_vect))/((np.linalg.norm(alpha_carbo_vect)**2))*alpha_carbo_vect
norm_ref_psi = alpha_nitro_vect - (np.dot(alpha_nitro_vect, alpha_carbo_vect))/((np.linalg.norm(alpha_carbo_vect)**2))*alpha_carbo_vect

psi_angle = math.acos(np.dot(norm_moving_psi, norm_ref_psi)/(np.linalg.norm(norm_moving_psi)*np.linalg.norm(norm_ref_psi)))*180/math.pi

if np.cross(norm_moving_psi, norm_ref_psi)[0]/alpha_carbo_vect[0] > 0:
    psi_angle *= -1
'''

phi_x = []
psi_y = []
# TODO: Make into function, check for optimization.
for i in range(0,8):
    for frame in u.trajectory[0:10001]:

        
        # Calculate the vectors needed for phi and psi dihedral angles for this frame
        alpha_nitro_vect = c_nitrogen.positions[i] - c_alpha.positions[i] #normal for phi, ref end for psi
        nitro_carbo_vect = carboxyl_meas.positions[i] - c_nitrogen.positions[i] #moving end phi
        alpha_carbo_vect = c_carboxyl.positions[i] - c_alpha.positions[i] #ref end phi, normal for psi
        carbo_nitro_vect = nitrogen_meas.positions[i] - c_carboxyl.positions[i] #moving end psi
            
        #Phi angles
        norm_moving_phi = nitro_carbo_vect - (np.dot(nitro_carbo_vect, alpha_nitro_vect))/((np.linalg.norm(alpha_nitro_vect)**2))*alpha_nitro_vect #Projection of nitro_carbo_vect onto the plane of alpha_nitro
        norm_ref_phi = alpha_carbo_vect - (np.dot(alpha_carbo_vect, alpha_nitro_vect))/((np.linalg.norm(alpha_nitro_vect)**2))*alpha_nitro_vect

        phi_value = np.dot(norm_moving_phi, norm_ref_phi)/(np.linalg.norm(norm_moving_phi)*np.linalg.norm(norm_ref_phi))
        phi_angle = math.acos(np.clip(phi_value,-1.0,1.0))*180/math.pi

        if np.cross(norm_moving_phi, norm_ref_phi)[0]/alpha_carbo_vect[0] < 0:
            phi_angle *= -1
        #Psi angles
        norm_moving_psi = carbo_nitro_vect - (np.dot(carbo_nitro_vect, alpha_carbo_vect))/((np.linalg.norm(alpha_carbo_vect)**2))*alpha_carbo_vect
        norm_ref_psi = alpha_nitro_vect - (np.dot(alpha_nitro_vect, alpha_carbo_vect))/((np.linalg.norm(alpha_carbo_vect)**2))*alpha_carbo_vect

        psi_value = np.dot(norm_moving_psi, norm_ref_psi)/(np.linalg.norm(norm_moving_psi)*np.linalg.norm(norm_ref_psi))
        psi_angle = math.acos(np.clip(psi_value,-1.0,1.0))*180/math.pi

        if np.cross(norm_moving_psi, norm_ref_psi)[0]/alpha_carbo_vect[0] > 0:
            psi_angle *= -1

        phi_x.append(phi_angle)
        psi_y.append(psi_angle)

plt.scatter(phi_x, psi_y, alpha=0.2, marker='.', edgecolors='none')
plt.xlabel("$\\phi$")
plt.ylabel("$\\psi$")
plt.xlim(-180, 180)
plt.ylim(-180, 180)
plt.show()