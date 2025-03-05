# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 18:48:32 2024

@author: Oxide
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Sep  2 15:46:32 2024

@author: Oxide
"""

import vtk
from vtk.util.numpy_support import vtk_to_numpy

# Charger le fichier VTK
file_path = 'D:\\W10 Documents\\Documents\\Le programme PydyCo\\Le programme PydyCo\\IN_vc.vtk'
reader = vtk.vtkStructuredGridReader()
reader.SetFileName(file_path)
reader.Update()

# Extraire les données du fichier VTK
data = reader.GetOutput()

# Extraire les vecteurs de courant électrique (Ix, Iy, Iz)
vectors_vtk = data.GetPointData().GetVectors()
current_vectors = vtk_to_numpy(vectors_vtk)

"""# Extraire les coordonnées des points (x, y, z)
points_vtk = data.GetPoints()
points = vtk_to_numpy(points_vtk)
"""
# Extraire les coordonnées des points (x, y, z)
points_vtk = data.GetPoints()
points = vtk_to_numpy(points_vtk.GetData())  # Utilisation de GetData() pour obtenir les coordonnées


print("Exemple de données de courant :")
print(current_vectors[:5])
print("Exemple de coordonnées :")
print(points[:5])


import numpy as np

# Décomposer les vecteurs de courant
Ix = current_vectors[:, 0]
Iy = current_vectors[:, 1]
Iz = current_vectors[:, 2]

ex=points[:, 0]
ey= points[:, 1]
ez=points[:, 2]
ex=ex.reshape(106,23)
ey=ey.reshape(106,23)

# Supposons que les points sont sur une grille régulière
Nx, Ny = 106, 23 # Ajuste en fonction des dimensions réelles
Lx, Ly = 106.0, -23.0  # Longueur des axes (à ajuster selon tes données)

# Créer les grilles de nombres d'onde
kx = np.fft.fftfreq(Nx, d=Lx/Nx) * 2 * np.pi
ky = np.fft.fftfreq(Ny, d=Ly/Ny) * 2 * np.pi

kx, ky = np.meshgrid(kx, ky)


# Appliquer la transformée de Fourier sur Ix et Iy
Ix_fft = np.fft.fft2(Ix.reshape(Ny, Nx))
Iy_fft = np.fft.fft2(Iy.reshape(Ny, Nx))



# Calculer le champ magnétique Bz dans l'espace fréquentiel
Bz_fft = (1j * (kx * Iy_fft - ky * Ix_fft)) / (kx**2 + ky**2 + 1e-10)  # Éviter la division par zéro

# Revenir dans le domaine spatial
Bz = np.fft.ifft2(Bz_fft).real

print("Champ magnétique Bz calculé :")
print(Bz)



import matplotlib.pyplot as plt

plt.imshow(Bz, extent=(0, Lx, 0, Ly), origin='lower', cmap='jet')
plt.colorbar(label='Champ Magnétique Bz')
plt.title('Champ Magnétique Calculé à partir du Courant Électrique')
plt.xlabel('x')
plt.ylabel('y')
plt.show()