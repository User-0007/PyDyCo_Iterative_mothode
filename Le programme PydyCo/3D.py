# -*- coding: utf-8 -*-
"""
Created on Mon Sep  2 15:46:32 2024

@author: Oxide
"""

import vtk
from vtk.util.numpy_support import vtk_to_numpy
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # For 3D plotting

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

# Extraire les coordonnées des points (x, y, z)
points_vtk = data.GetPoints()
points = vtk_to_numpy(points_vtk.GetData())  # Utilisation de GetData() pour obtenir les coordonnées

print("Exemple de données de courant :")
print(current_vectors[:5])
print("Exemple de coordonnées :")
print(points[:5])

# Décomposer les vecteurs de courant
Ix = current_vectors[:, 0]
Iy = current_vectors[:, 1]
Iz = current_vectors[:, 2]

# Extraire les coordonnées des points (x, y, z)
ex = points[:, 0]
ey = points[:, 1]
ez = points[:, 2]  # For 3D plotting, we'll use ez as the z-coordinate

# Reshape the data for grid-like plotting
Nx, Ny = 106, 23  # Ajuste en fonction des dimensions réelles
ex = ex.reshape(Ny, Nx)
ey = ey.reshape(Ny, Nx)

# Supposons que les points sont sur une grille régulière
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

# Plotting the 3D surface
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Create a 3D surface plot for Bz
surf = ax.plot_surface(ex, ey, Bz, cmap='jet', edgecolor='none')

ax.set_xlabel('X axis')
ax.set_ylabel('Y axis')
ax.set_zlabel('Bz (Magnetic Field)')
plt.title('3D Magnetic Field (Bz) Calculated from Electric Current')

# Add a color bar to reflect Bz values
fig.colorbar(surf, ax=ax, label='Champ Magnétique Bz')

plt.show()
