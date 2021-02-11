from NanoTCAD_ViDES import * 
import sys, os
from module_TMD import *
from myplot import plot3d

rank = 0
def solve_self_consistent(grid,interface,channel):
    normad=1e30;
    interface.Phiold=interface.Phi.copy();
    counter=1;
    rank=0;

    while (normad>interface.normd):
        channel.Phi=interface.Phi[grid.swap] # I pass the potential in correspondence of the atoms of the material for which I compute the NEGF
        channel.charge_T();

        writeout("--------------------------------------------")
        string="            CURRENT = %s A/m" %(channel.current());
        writeout(string);
        writeout("--------------------------------------------")
        
        interface.free_charge[grid.swap]=channel.charge # I pass back the free_charge term to the 3D domain
        
        # savetxt("ncar.ini",interface.free_charge);
        # savetxt("Phi.ini",interface.Phi);
        
        solve_Poisson(grid,interface); # I solve Poisson
        normad=max(abs(interface.Phiold-interface.Phi))
        interface.Phi=interface.Phi+(interface.underel)*(interface.Phiold-interface.Phi)
        del interface.Phiold;
        interface.Phiold=interface.Phi.copy();
        
        print()
        string="Iteration # %s; ||Phi-Phiold||2 = %s > %s" %(counter,normad, interface.normd)
        writeout(string) 
        print() 
        counter=counter+1;
        if (counter>600):
            return;

# I create the grid
# xg=nonuniformgrid(array([-2.0,1,0,0.05,2.0,1]))
xg = linspace(-2, 2, 21)

FLAKE = TMD(30.0,"n");
acc = FLAKE.acc;
kF = 2*pi/(3*sqrt(3)*acc);
kymax = pi/FLAKE.delta;
Nky = 32.0;
dk = kymax/Nky;
FLAKE.kmax = pi/FLAKE.delta;
FLAKE.kmin = 0;
FLAKE.dk = dk;

grid = grid2D(xg, FLAKE.y, FLAKE.x, FLAKE.y);

# I take care of the solid
Oxide_Bottom=region("hex",grid.xmin,-1,grid.ymin,grid.ymax)
Oxide_Bottom.eps=3.9;

Oxide_Top=region("hex",1,grid.xmax,grid.ymin,grid.ymax)
Oxide_Top.eps=3.9;

Channel=region("hex",-1,1,grid.ymin,grid.ymax)
Channel.eps=10;

top_gate_2=gate("hex",grid.xmax,grid.xmax,grid.ymin,grid.ymax);
top_gate=gate("hex",grid.xmax,grid.xmax,10.0,20.0);
bottom_gate_2=gate("hex",grid.xmin,grid.xmin,grid.ymin,grid.ymax);
bottom_gate=gate("hex",grid.xmin,grid.xmin,10.0,20.0);

p=interface2D(grid,Oxide_Bottom,Oxide_Top, Channel,top_gate_2, top_gate,bottom_gate, bottom_gate_2);

fraction_source=0.001
fraction_drain=0.001
dope_reservoir(grid,p,FLAKE,fraction_source,array([-1,1,grid.ymin,10.0]));
dope_reservoir(grid,p,FLAKE,fraction_drain,array([-1,1,20.0,grid.ymax]));
# dope_reservoir(grid,p,FLAKE,-0.01,array([-1,1,10.0,20.0]));

Vgs = 0.5
bottom_gate_2.Ef=-0; set_gate(p,bottom_gate_2)
top_gate_2.Ef=-0; set_gate(p,top_gate_2)
bottom_gate.Ef=-Vgs; set_gate(p,bottom_gate)
top_gate.Ef=-Vgs; set_gate(p,top_gate)

solve_Poisson(grid,p);
# plot3d(grid, p.fixed_charge, 5)
Phi = p.Phi.reshape((len(grid.gridy), len(grid.gridx)))
Phi_mean = (1/grid.ymax)*sum(Phi, axis=0)*(grid.gridy[1]-grid.gridy[0])
# Phi_mean = sum(Phi, axis=0)
plot3d(grid, p.Phi)
# plot(grid.gridx, Phi[len(grid.gridy)/2, :])
# plot(Phi_mean)
# plot(grid.gridx, Phi_mean)
show()
exit()




# plot(p.Phi[grid.swap]); show()
# solve_init(grid,p,FLAKE);

Vgmin=0.0;
Vgmax=1.0;
Vgstep=0.05;

Np=int(abs(Vgmin-Vgmax)/Vgstep)+1;
vg=zeros(Np);
current=zeros(Np);
p.underel=0.1;

counter = 0;
Vgs = 0.1;
FLAKE.mu1 = -0.0
FLAKE.mu2 = -0.1

try:
    os.mkdir('datiout')
except:
    pass

while (Vgs<=Vgmax):
    bottom_gate.Ef=-Vgs; set_gate(p,bottom_gate)
    top_gate.Ef=-Vgs; set_gate(p,top_gate)
    p.normpoisson=1e-1;
    p.normd=5e-2;
    solve_self_consistent(grid,p,FLAKE);
    vg[counter]=Vgs;
    current[counter]=FLAKE.current();
    # I save the output files
    string="./datiout/Phi%s.out" %Vgs;
    savetxt(string,p.Phi);
    string="./datiout/ncar%s.out" %Vgs;
    savetxt(string,p.free_charge);
    a=[FLAKE.E,FLAKE.T];
    string="./datiout/T%s.out" %Vgs;
    savetxt(string,transpose(a));
    string="./datiout/jayn%s.out" %Vgs;
    fp=open(string,"w");
    string2="%s" %current[counter];
    fp.write(string2);
    fp.close();
    counter=counter+1;
    Vgs=Vgs+Vgstep;

tempo=[vg,current]
savetxt("./datiout/idvgs.out",transpose(tempo));
