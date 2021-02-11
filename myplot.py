from NanoTCAD_ViDES import *
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

def plot3d(grid, Z, n=1):
	X,Y = meshgrid(grid.gridx, grid.gridy)
	Phi = reshape(Z, X.shape)
	fig = plt.figure()
	ax = fig.gca(projection='3d')
	# Plot the surface.
	surf = ax.plot_surface(X[::n], Y[::n], Phi[::n] ,rstride=1, cstride=1, cmap='inferno', edgecolor='none')

	fig.colorbar(surf, shrink=0.5, aspect=5)
	plt.show()
	# return;

    