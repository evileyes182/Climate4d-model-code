# Figure 3(b): equilibria of system (10) continued in mu for eta = -1.5
# Run: python3 run.py (with AUTO-07p set up); writes eq_*.dat (branches) and eq_*_special.dat (Hopf and saddle-node points)
from auto import *

# Equilibria, continued in mu in both directions
r1 = rl(r(e='varRun', c='genericRun.t', NMX=50000, NTST=0, NCOL=4, JAC=1, IPS=1, ISP=2, NPR=1000, ISW=1, ICP=[1,2,11,12,13,14,15], DS=1e-3, DSMAX=2e-2))
r2 = rl(r(r1(1), DS='-'))

# Save: branch (mu, Psi, stable) and special points (mu, Psi, type)
for run, name in [(r1, 'eq_1'), (r2, 'eq_2')]:
    stab = run[0].stability()
    with open(name + '.dat', 'w') as f:
        for i, p in enumerate(run[0]):
            s = [k for k in stab if abs(k) >= i + 1][0]
            f.write('%.10e %.10e %d\n' % (p['PAR(1)'], p['PAR(14)'], s < 0))
    with open(name + '_special.dat', 'w') as f:
        for p in run[0]:
            if p['TY name'] in ('LP', 'HB'):
                f.write('%.10e %.10e %s\n' % (p['PAR(1)'], p['PAR(14)'], p['TY name']))
