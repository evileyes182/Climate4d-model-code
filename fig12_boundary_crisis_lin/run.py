# Figure 12: boundary-crisis curve BC of system (10), a fold of homoclinic orbits of a saddle periodic orbit,
# computed with Lin's method from the periodic orbit of H1 at (mu, eta) = (-2.0672e-3, -4.1)
from auto import *

## First run to get periodic solution 
r1 = rl(r(e='hom_po_lin', c='genericRun.t', NMX=50000, NTST=25, NCOL=4, IPS=1, ISP=2, NPR=1000, ISW=1, ICP=[1,2,11,12,13,14,15], DS=1e-3, DSMAX=2e-2))
r2 = rl(r(r1(1), DS='-'))
r1 = rl(r1 + r2)
# Comptute PO 1
po1 = rl(r(r1('HB4'), NMX=2000, NTST=500, NCOL=4, STOP={'BP5'},UZSTOP={1:-2.0672e-3, 11:1e6}, IPS=2, ISP=2, ISW=1, ICP=[1,11,2,12,13,14,15], DS=1e-4, DSMAX=2e-3))
## Second run to get unstable eignvalues 
floq1 = rl(r(po1('UZ1'), ICP=['lambda1','h1',11], IPS=4, NDIM=8, NBC=9, NINT=1,STOP=['BP3'],UZSTOP={}, NMX=9000, NTST=50, DS=1e-3, DSMAX=2e-3))
floq2 = r(floq1(1), DS='-', STOP=['BP1'])
floq = rl(floq1 + floq2)

## Fourth run to get stable eignfunctions  
stEF1 = rl(r(floq('BP1'), ICP=['lambda1','h1',11], ISW=-1, ILP=0, UZSTOP={'h1':1.0}, NMX=1000))
stEF2 = rl(r(stEF1('UZ1'), DS=1e-4, DSMAX=1e-3, ICP=['lambda2','h2',11,'lambda1'], NDIM=12, NBC=14, ISW=1, ILP=0, UZSTOP={'h2':1.0}, NMX=1000))
stEF2 = rl(r(stEF2('BP1'), ISW=-1, NMX=3000)) 
## Third run to get eignfunction unstabl
unST1 = rl(r(stEF2('UZ1'), ICP=['lambda3','h3',11,'lambda2','lambda1'], PAR={'lambda3': 0.95}, DS=-1e-3, NDIM=16, NBC=19, ISW=1, ILP=0, UZSTOP={'h3':1.0}))
unST1 = rl(r(unST1('BP1'), ISW=-1))

## Unstable manifold
runWu = rl(r(unST1('UZ1'), ICP=['sigma+', 'T+', 'T','lambda2','lambda1', 'lambda3','LIN'], 
            ISW=1, NMX=25000, NDIM=20, NBC=24, UZR={'sigma+':0}, UZSTOP={},STOP=['UZ44'],
            NCOL=4, NTST=400, DSMAX=1e3, DS=1e1))

## Integrate backwards from periodic solution 
runWs = rl(r(runWu('UZ44'), 
            ICP=['T-','sigma-','T','T+','lambda2','lambda1', 'lambda3'],NDIM=24,NBC=29,
            STOP=['UZ5'],UZSTOP={'sigma-':0}, DS=1e1, NPR=50, NMX=500, NTST=800))  
## Close lins gap
gap = run(runWs('UZ1'),ICP=['LIN','T+','T-','eps2','eps1','T','sigma+','lambda3','lambda2','lambda1'],
        NBC=33, NDIM=24, NINT=0, THL={'T':0}, STOP=['UZ1'],UZR={'LIN':0}, NMX=2000, ISW=1, NTST=500,  DSMAX=2e2, DS=1e2)

" ~~Continue homoclinic curves~~ "
hom1 = rl(r(gap('UZ1'), NINT=1, ILP=1, ICP=[1,2,'T+','T-','epu1', 'eps2','eps1','T','sigma+','lambda3','lambda2','lambda1'], NMX=1000, STOP=['LP1'], THL={'T':0,'T+':0,'T-':0}))

#~~~ closed loop lp
# LP1
lp1 = rl(r(hom1('LP1'), NINT=1, NPR=500, NMX=2750, NTST=500, ISW=2, SP={'BP0'}, ICP=[1,2,'T+','T-','epu1', 'eps2','eps1','T','sigma+','lambda3','lambda2','lambda1'], THL={'T':0,'T+':0,'T-':0}))
lp1 = rl(r(lp1('EP1'), NMX=50))

# Save: BC curve (mu, eta)
with open('bc.dat', 'w') as f:
    for p in lp1[0]:
        f.write('%.10e %.10e\n' % (p['PAR(1)'], p['PAR(2)']))
