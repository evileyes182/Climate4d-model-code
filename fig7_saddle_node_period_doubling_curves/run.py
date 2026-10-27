# Figure 7: saddle-node curve S and period-doubling curve PD of periodic orbits of system (10) in the (mu, eta)-plane,
# continued from the S and PD points bounding the second stable interval of the branch from H1 at eta = -4.08
# Run: auto run.py (with AUTO-07p set up; about 10 min); writes s_*.dat and pd_*.dat (mu, eta)
from auto import *

# Equilibria, continued in mu
r1 = rl(r(e='varRun', c='genericRun.t', NMX=50000, NTST=0, NCOL=4, JAC=1, IPS=1, ISP=2, NPR=1000, ISW=1, ICP=[1,2,11,12,13,14,15], DS=1e-3, DSMAX=2e-2))

# PO from H1 (HB5 on the branch), up to the second PD point
po1 = rl(r(r1('HB5'), NMX=750, NTST=500, NCOL=4, STOP={'BP5'}, UZSTOP={11:1e6}, JAC=1, IPS=2, ISP=2, ISW=1, ICP=[1,11,2,12,13,14,15], DS=1e-4, DSMAX=2e-3))

# S curve from LP5 (the second stable interval): generate starting data,
# then continue in both directions; save (mu, eta)
st_s = r(po1('LP5'), ISW=2, ICP=[1,2,11], JAC=1, NPR=5000)
s1 = r(st_s, NMX=5000, DS=1e-4, UZSTOP={1:[-2.17e-3,-2.04e-3], 2:[-4.3,-3.65]})
with open('s_1.dat', 'w') as f:
    for p in s1[0]:
        f.write('%.10e %.10e\n' % (p['PAR(1)'], p['PAR(2)']))
s2 = r(st_s, NMX=5000, DS=-1e-4, UZSTOP={1:[-2.17e-3,-2.04e-3], 2:[-4.3,-3.65]})
with open('s_2.dat', 'w') as f:
    for p in s2[0]:
        f.write('%.10e %.10e\n' % (p['PAR(1)'], p['PAR(2)']))

# PD curve from PD2 (the second stable interval): generate starting data,
# then continue in both directions; save (mu, eta)
st_pd = r(po1('PD2'), ISW=2, ICP=[1,2,11], JAC=1, NPR=5000)
pd1 = r(st_pd, NMX=5000, DS=1e-4, UZSTOP={1:[-2.17e-3,-2.04e-3], 2:[-4.3,-3.65]})
with open('pd_1.dat', 'w') as f:
    for p in pd1[0]:
        f.write('%.10e %.10e\n' % (p['PAR(1)'], p['PAR(2)']))
pd2 = r(st_pd, NMX=5000, DS=-1e-4, UZSTOP={1:[-2.17e-3,-2.04e-3], 2:[-4.3,-3.65]})
with open('pd_2.dat', 'w') as f:
    for p in pd2[0]:
        f.write('%.10e %.10e\n' % (p['PAR(1)'], p['PAR(2)']))
