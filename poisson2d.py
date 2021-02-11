'''
Sat 29 Jun 2019 09:53:17 GMT
Working perfectly
'''

from NanoTCAD_ViDES import *
from myplot import plot3d

def dope_reservoir(grid,interface,delta,molar_fraction,bbox):
  # xmin=bbox[0];
  # xmax=bbox[1];
  # ymin=bbox[2];
  # ymax=bbox[3];
  # # index=nonzero((xmin<=grid.x2D[grid.swap])&(xmax>=grid.x2D[grid.swap])&(ymin<=grid.y2D[grid.swap])&(ymax>=grid.y2D[grid.swap]))
  # index=nonzero((xmin<=grid.x2D)&(xmax>=grid.x2D)&(ymin<=grid.y2D)&(ymax>=grid.y2D))
  # interface.fixed_charge[index] = molar_fraction/delta*1e9;

  xmin=bbox[0];
  xmax=bbox[1];
  ymin=bbox[2];
  ymax=bbox[3];

  index=nonzero((xmin<=grid.x2D[grid.swap])&(xmax>=grid.x2D[grid.swap])&(ymin<=grid.y2D[grid.swap])&(ymax>=grid.y2D[grid.swap]))
  interface.fixed_charge[grid.swap[index]]=molar_fraction/delta*1e9;


# Dimensions
tox = 3.0; tch = 5.0
Lx = 2*tox+tch
Lsd = 10.0; Lch = 10.0
Ly = Lch + 2.0*Lsd

Nx = 23; Ny = 61;
# x = linspace(0, Lx, Nx)
x = nonuniformgrid(array([0, .5, tox,0.25, tox+tch, 0.25, 2*tox+tch, 0.5]))
xc = x
y = linspace(0, Ly, Ny)
yc = y
delta = .15
# print(y.shape)
# exit()

x = linspace(0, 10, 11)
xc = linspace(0, 10, 6)
y = linspace(10, 20, 11)
yc = y
grid = grid2D(x,y,xc,yc)

print(grid.y2D)
exit()
# Gate definition
top_gate=gate("hex", grid.xmax, grid.xmax, grid.ymin+Lsd, grid.ymax-Lsd)
top_gate.Ef=-1;

bottom_gate=gate("hex", grid.xmin, grid.xmin, grid.ymin+Lsd, grid.ymax-Lsd)
bottom_gate.Ef=-1;

# right_gate=gate("hex", grid.xmin, grid.xmax, grid.ymax, grid.ymax)
# right_gate.Ef=-2;

# left_gate=gate("hex", grid.xmin, grid.xmax, grid.ymin, grid.ymin)
# left_gate.Ef=0;

# Region definition
oxide = region("hex", grid.xmin, grid.xmax, grid.ymin, grid.ymax)
oxide.eps = 3.9

channel = region("hex", grid.xmin+tox, grid.xmax-tox, grid.ymin, grid.ymax)
channel.eps = 11.7

# Interface 
# p=interface2D(grid, top_gate, bottom_gate, right_gate, left_gate, channel);
# p=interface2D(grid, top_gate, bottom_gate, left_gate, oxide, channel);
p=interface2D(grid, top_gate, bottom_gate, oxide, channel);
p.normpoisson=1e-5;

fraction=5e-3
dope_reservoir(grid,p,delta,fraction,array([tox, tox+tch, 0, Lsd]))
dope_reservoir(grid,p,delta,fraction,array([tox, tox+tch, Lsd+Lch, 2*Lsd+Lch]))
# dope_reservoir(grid,p,delta, fraction,array([grid.xmin, grid.xmax, grid.ymin, grid.ymin+10]));
# dope_reservoir(grid,p,delta,-fraction,array([grid.xmin, grid.xmax, grid.ymax-10, grid.ymax]));

solve_Poisson(grid,p);

# X,Y = meshgrid(grid.gridx, grid.gridy)
# Phi = reshape(p.Phi, X.shape)
# fig = plt.figure()
# ax = fig.gca(projection='3d')
# # Plot the surface.
# surf = ax.plot_surface(X, Y, Phi ,rstride=1, cstride=1, cmap='inferno', edgecolor='none')

# fig.colorbar(surf, shrink=0.5, aspect=5)
# plt.show()
# exit()

plot3d(grid, p.fixed_charge, n=5)

