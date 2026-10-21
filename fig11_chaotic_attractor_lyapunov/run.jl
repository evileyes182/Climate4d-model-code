# Figure 11: chaotic attractor of system (10) at (mu, eta) = (-2.128e-3, -3.99), reached by shooting off the
# 2-D unstable manifold of the northward equilibrium, and its largest Lyapunov exponent lambda_max (Figure 12)
# Run: julia run.jl (about 1 min); writes equilibrium.dat, eigenvalues.dat, attractor.dat, events.dat and lyapunov.dat
using StaticArrays, DifferentialEquations, OrdinaryDiffEq, ChaosTools, DynamicalSystemsBase, LinearAlgebra, DelimitedFiles

include("4D_climate_define.jl")
AMOCfunc = ODEFunction(AutoSysODE2!, jac=AutoSysJac!)

SV = 1e6
σ  = 2.1e4*1e6
p  = @SVector[-2.128e-3, -3.99, 8SV/σ, SV, σ, 0.02]

# Northward equilibrium: Newton from a hard-coded guess
xeq = [1.21, 2.41, -3.76, 1.60]
F = zeros(4)
J = zeros(4, 4)
for i in 1:20
    AutoSysODE2!(F, xeq, p, 0.0)
    AutoSysJac!(J, xeq, p, 0.0)
    xeq .-= J \ F
end
AutoSysJac!(J, xeq, p, 0.0)
λ, V = eigen(J)       # λ[3:4] = a ∓ ib with a > 0, λ[1:2] < 0
v = V[:, 4]

# Shooting point on the unstable eigenplane
δ = 1e-4
θ = 6.2
ic0 = xeq .+ δ .* (cos(θ) .* real.(v) .+ sin(θ) .* imag.(v))
writedlm("equilibrium.dat", [xeq ic0])
writedlm("eigenvalues.dat", [real.(λ) imag.(λ)])

# Two-stage solve: transient [0, TMIN] not stored, window [TMIN, TMAX] stored
thresh = 5.0
TMIN = 2e6
TMAX = 4e6
tau_slow = 3.0e4
tau_fast = 1.0e3

# 1) Transient stage [0, TMIN]: run with divergence guard, do not store steps
div_cb = DiscreteCallback((u,t,int)->maximum(abs,u) > thresh, terminate!)
prob1 = ODEProblem(AMOCfunc, ic0, (0.0, TMIN), p)
sol1  = solve(prob1, AutoVern7(Rodas5P());
              callback=div_cb, abstol=1e-9, reltol=1e-9,
              maxiters=10_000_000, save_everystep=false,
              save_start=false, save_end=true, dtmax=0.07*tau_fast)

stateT = sol1.u[end]
prob2  = ODEProblem(AMOCfunc, stateT, (TMIN, TMAX), p)
sol    = solve(prob2, AutoVern7(Rodas5P());
               callback=div_cb, abstol=1e-10, reltol=1e-9,
               maxiters=10_000_000, save_everystep=true, dtmax=0.01*tau_fast)

state = sol.u[end]

# Time series over the Figure 11 window t in [0, 4e5] after the transient: t, y_N, y_E, x_B, x_N, K_N, Psi
ϕ   = 1.7e-4*(7.0 - 2.65)
x_E = (25.0 - 2.65)/(7.0 - 2.65)
y_B = 0.8e-3/1.7e-4*(34.538 - 35.0)/(7.0 - 2.65)
k1  = 1e-2*SV/σ
t   = TMIN:10.0:TMIN + 4e5
U   = reduce(hcat, sol(t).u)'
KN  = k1 .+ 0.5 .* (p[3] - k1) .* (1 .+ tanh.(((U[:, 1] .- y_B) .- (U[:, 4] .- U[:, 3]) .- p[2]) ./ p[6]))
Psi = ϕ .* ((U[:, 1] .- U[:, 4]) .- (U[:, 2] .- x_E))
writedlm("attractor.dat", [t .- TMIN U KN Psi])

# Convective shutdown events: minima of K_N below κ_star, less a small tolerance
κ_star = k1 + 0.5*(p[3] - k1)*(1 - inv(sqrt(2.0)))
ev = [i for i in 2:length(KN)-1 if KN[i-1] > KN[i] <= KN[i+1] && KN[i] < κ_star - max(1e-7, 1e-4*κ_star)]
writedlm("events.dat", [t[ev] .- TMIN U[ev, :] KN[ev] Psi[ev]])

# Largest Lyapunov exponent from the state at TMAX (QR method): transient Ttr, then N steps of Δt
Ttr   = tau_slow
Δt    = 0.2*tau_fast
Tavg  = 250.0*tau_slow
N     = round(Int, Tavg/Δt)

alg = AutoVern7(Rodas5P(autodiff = AutoFiniteDiff()))

ds  = CoupledODEs(AutoSysODE2!, state, p;
    diffeq=(alg=alg, abstol=1e-10, reltol=1e-9,
            save_everystep=false, maxiters=10_000_000))

# k = 1: largest exponent only; the tangent dynamics use a ForwardDiff Jacobian of AutoSysODE2!
λ1 = lyapunovspectrum(ds, N, 1; Δt=Δt, Ttr=Ttr, show_progress=false)[1]
writedlm("lyapunov.dat", λ1)
