'''
Working perfectly
x-axis is the confinement direction
y-axis is the transport direction
'''

from NanoTCAD_ViDES import *
from myplot import plot3d

def dope_reservoir(grid,interface,delta,molar_fraction,bbox):
  xmin=bbox[0];
  xmax=bbox[1];
  ymin=bbox[2];
  ymax=bbox[3];

  ## Originally it had the following line, I had to change it
  # index=nonzero((xmin<=grid.x2D[grid.swap])&(xmax>=grid.x2D[grid.swap])&(ymin<=grid.y2D[grid.swap])&(ymax>=grid.y2D[grid.swap]))
  index=nonzero((xmin<=grid.x2D)&(xmax>=grid.x2D)&(ymin<=grid.y2D)&(ymax>=grid.y2D))
  # print(grid.x2D.reshape(len(grid.gridy), len(grid.gridx)))
  # print(index)
  ## Originally it had the following line, I had to change it
  # interface.fixed_charge[grid.swap[index]]=molar_fraction/delta*1e9;
  interface.fixed_charge[index]=molar_fraction/delta*1e9;


Lsd = 10.; Lch = 10.; Ny = 61
tox = 2.; tch = 3.; Nx = 29

x = linspace(-0.5*(2*tox+tch), 0.5*(2*tox+tch), Nx)
y = linspace(0, 2.*Lsd+Lch, Ny)
xc = x
yc = y
delta = 0.15

grid = grid2D(x,y,xc,yc)

print(grid.xmin, grid.xmax, grid.ymin, grid.ymax)
print(grid.gridx[1]-grid.gridx[0], grid.gridy[1]-grid.gridy[0])

top_gate=gate("hex", grid.xmax, grid.xmax, grid.ymin+Lsd, grid.ymax-Lsd); top_gate.Ef = -1
bottom_gate=gate("hex", grid.xmin, grid.xmin, grid.ymin+Lsd, grid.ymax-Lsd); bottom_gate.Ef = -1

# Region definition
top_oxide = region("hex", grid.xmax-tox, grid.xmax, grid.ymin, grid.ymax)
top_oxide.eps = 3.9

bottom_oxide = region("hex", grid.xmin, grid.xmin+tox, grid.ymin, grid.ymax)
bottom_oxide.eps = 3.9

channel = region("hex", grid.xmin+tox, grid.xmax-tox, grid.ymin, grid.ymax)
channel.eps = 11.9

p = interface2D(grid, top_gate, bottom_gate, bottom_oxide, top_oxide, channel);
p.normpoisson = 1e-5;

fraction = -1e-3
dope_reservoir(grid, p, delta, fraction, array([grid.xmin+tox, grid.xmax-tox, grid.ymin, grid.ymin+Lsd]))
dope_reservoir(grid, p, delta, fraction, array([grid.xmin+tox, grid.xmax-tox, grid.ymax-Lsd, grid.ymax]))

solve_Poisson(grid,p);

# plot3d(grid, p.Phi, n=1)
X,Y = meshgrid(grid.gridx, grid.gridy)
dy = grid.gridy[1]-grid.gridy[0]
Phi2d = reshape(p.Phi, X.shape)

Phi_x = sum(Phi2d, axis=0) * dy/(grid.ymax-grid.ymin)
plot(grid.gridx, Phi_x)
show()

exit()