# Numerical methods and examples

This is Appendix A of J. Bailie, P. Subramanian and B. Krauskopf, *Regular and chaotic Welander oscillations in a four-dimensional conceptual model for the Atlantic Meridional Overturning Circulation*, [arXiv:2603.11577](https://arxiv.org/abs/2603.11577). Equation and figure numbers refer to the paper.

The computations in the paper were carried out with the continuation package AUTO-07p [1] and with Julia [2]. This repository contains six self-contained folders that each reproduce one of the computations described below; each contains the right-hand side of system (2.10) as a Fortran file for AUTO-07p (or as a Julia function), a driver `run.py` (or `run.jl`) that performs the computation and writes the results as plain data files, and a plotting script `plot.jl` based on Makie [3]. A folder is run by executing the driver and then the plotting script. Note that AUTO-07p, Python 3, and Julia must be installed.

For example, the branches of equilibria in Figure 3(b) are computed and plotted by running, from the top folder of this repository,

```bash
cd fig3b_equilibria
auto run.py
julia plot.jl
```

The branches of equilibria in Figure 3(b), with their saddle-node and Hopf bifurcations $\mathrm{SN}$ and $\mathrm{H}$, are computed in the folder [`fig3b_equilibria`](fig3b_equilibria) by continuation in $\mu$ for $\eta=-1.5$. The Hopf curve and the saddle-node curve in Figure 4 are computed in the folder [`fig4_hopf_saddle_node_curves`](fig4_hopf_saddle_node_curves) by continuing a point $\mathrm{H}$ and a point $\mathrm{SN}$ in the $(\mu,\eta)$-plane. Along the Hopf curve, the first Lyapunov coefficient $\ell_1$ [4] is the test function for generalized Hopf points $\mathrm{GH}$ and distinguishes supercritical ($\ell_1<0$, solid) from subcritical ($\ell_1>0$, dashed) segments.

The branch of periodic orbits in Figure 6(c), with its saddle-node, period-doubling and torus bifurcations $\mathrm{S}$, $\mathrm{PD}$ and $\mathrm{T}$, is computed in the folder [`fig6c_periodic_orbits`](fig6c_periodic_orbits) by continuation in $\mu$ from $\mathrm{H}_1$ for $\eta=-4.08$. The branch is plotted as the minimum and maximum of $\Psi$ along each orbit. We follow the bifurcations $\mathrm{S}$ and $\mathrm{PD}$ as curves in the $(\mu,\eta)$-plane in the folder [`fig7_saddle_node_period_doubling_curves`](fig7_saddle_node_period_doubling_curves); they bound the second region of stable periodic orbits in Figure 7.

The chaotic attractor in Figure 11 is computed in the folder [`fig11_chaotic_attractor_lyapunov`](fig11_chaotic_attractor_lyapunov) at $(\mu,\eta)=(-2.128\times10^{-3},-3.99)$ by integrating system (2.10) with `DifferentialEquations.jl` [5] from a point in the two-dimensional unstable eigenspace of the northward equilibrium. The largest Lyapunov exponent $\lambda_{\max}$ of the attractor is computed with the QR method as implemented in `ChaosTools.jl` [6]. We present this computation for a single parameter point; the same computation on a grid in the $(\mu,\eta)$-plane gives the values of $\lambda_{\max}$ in Figure 12.

The boundary-crisis curve $\mathrm{BC}_2$ in Figure 12 is computed with Lin's method [7] in the folder [`fig12_boundary_crisis_lin`](fig12_boundary_crisis_lin) as a homoclinic tangency of a saddle periodic orbit. The periodic orbit from $\mathrm{H}_1$ for $\eta=-4.1$ is continued to $\mu=-2.0672\times10^{-3}$, where orbit segments in its unstable and stable manifolds are grown up to a common section. The two segments form a homoclinic orbit once Lin's gap between them is closed. This orbit is continued to a fold, which is then followed in the $(\mu,\eta)$-plane.

## References

1. E. J. Doedel and B. E. Oldeman, *AUTO-07P: Continuation and bifurcation software for ordinary differential equations*, Concordia University, Montreal, Canada, 2012.
2. J. Bezanson, A. Edelman, S. Karpinski and V. B. Shah, Julia: A fresh approach to numerical computing, *SIAM Review*, 59 (2017), pp. 65–98.
3. S. Danisch and J. Krumbiegel, Makie.jl: Flexible high-performance data visualization for Julia, *Journal of Open Source Software*, 6 (2021), p. 3349.
4. Y. A. Kuznetsov, *Elements of Applied Bifurcation Theory*, 2nd ed., Springer, 1998.
5. C. Rackauckas and Q. Nie, DifferentialEquations.jl – A performant and feature-rich ecosystem for solving differential equations in Julia, *Journal of Open Research Software*, 5 (2017), p. 15.
6. G. Datseris and U. Parlitz, *Nonlinear Dynamics: A Concise Introduction Interlaced with Code*, Springer, 2022.
7. B. Krauskopf and T. Rieß, A Lin's method approach to finding and continuing heteroclinic connections involving periodic orbits, *Nonlinearity*, 21 (2008), pp. 1655–1690.
