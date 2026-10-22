# Figure 4: one Hopf curve H and one saddle-node curve SN of system (10) in the (mu, eta)-plane,
# continued from the equilibria at eta = -1.5
# Run: auto run.py (with AUTO-07p set up); writes hb_*.dat (mu, eta, l1) and lp_*.dat (mu, eta)
from auto import *

# Equilibria, continued in mu
r1 = rl(r(e='varRun', c='genericRun.t', NMX=50000, NTST=0, NCOL=4, JAC=1, IPS=1, ISP=2, NPR=1000, ISW=1, ICP=[1,2,11,12,13,14,15], DS=1e-3, DSMAX=2e-2))

# Hopf curve from HB2 in both directions; save (mu, eta, first Lyapunov coefficient from the d-file)
for ds, name in [(1e-3, 'hb_1'), (-1e-3, 'hb_2')]:
    hb = r(r1('HB2'), ISW=2, ICP=[1,2,11], JAC=1, NMX=20000, NPR=20000, DS=ds, UZSTOP={1:[-5e-3,0.0], 2:[-6.2,0.0]})
    l1 = {}
    for line in open('fort.9'):
        if 'GH   Function' in line:
            l1[int(line.split()[1])] = float(line.split()[4])
    l1[1] = l1[2]
    with open(name + '.dat', 'w') as f:
        for i, p in enumerate(hb[0]):
            f.write('%.10e %.10e %.10e\n' % (p['PAR(1)'], p['PAR(2)'], l1[i + 1]))

# Saddle-node curve from LP1 in both directions; save (mu, eta)
for ds, name in [(1e-3, 'lp_1'), (-1e-3, 'lp_2')]:
    lp = r(r1('LP1'), ISW=2, ICP=[1,2], JAC=1, NMX=20000, NPR=20000, DS=ds, UZSTOP={1:[-5e-3,0.0], 2:[-6.2,0.0]})
    with open(name + '.dat', 'w') as f:
        for p in lp[0]:
            f.write('%.10e %.10e\n' % (p['PAR(1)'], p['PAR(2)']))
