# Figure 3(b): equilibria of system (10) continued in mu for eta = -1.5
# Run: auto run.py (with AUTO-07p set up); writes eq_*.dat (branches) and eq_*_special.dat (Hopf and saddle-node points)
from auto import *

# Equilibria, continued in mu in both directions
r1 = rl(r(e='varRun', c='genericRun.t', NMX=50000, NTST=0, NCOL=4, JAC=1, IPS=1, ISP=2, NPR=1000, ISW=1, ICP=[1,2,11,12,13,14,15], DS=1e-3, DSMAX=2e-2))
r2 = rl(r(r1(1), DS='-'))

# Save first branch: branch (mu, Psi, stable) and special points (mu, Psi, type)
stab = r1[0].stability()
with open('eq_1.dat', 'w') as f:
    for i, p in enumerate(r1[0]):
        s = [k for k in stab if abs(k) >= i + 1][0]
        f.write('%.10e %.10e %d\n' % (p['PAR(1)'], p['PAR(14)'], s < 0))
with open('eq_1_special.dat', 'w') as f:
    for p in r1[0]:
        if p['TY name'] in ('LP', 'HB'):
            f.write('%.10e %.10e %s\n' % (p['PAR(1)'], p['PAR(14)'], p['TY name']))

# Save second branch
stab = r2[0].stability()
with open('eq_2.dat', 'w') as f:
    for i, p in enumerate(r2[0]):
        s = [k for k in stab if abs(k) >= i + 1][0]
        f.write('%.10e %.10e %d\n' % (p['PAR(1)'], p['PAR(14)'], s < 0))
with open('eq_2_special.dat', 'w') as f:
    for p in r2[0]:
        if p['TY name'] in ('LP', 'HB'):
            f.write('%.10e %.10e %s\n' % (p['PAR(1)'], p['PAR(14)'], p['TY name']))
