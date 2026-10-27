# System (10): right-hand side AutoSysODE2! and Jacobian AutoSysJac!, p = (mu, eta, kappa_N, SV, sigma, epsilon)
# Included by run.jl
using StaticArrays
# ---------- RHS -------------------------------------------------
function AutoSysODE2!(du, u, p, t)
       # use ordinary scalar unpacking; u and du can be any AbstractVector
       y_N, y_E, x_B, x_N = u
       m, η, k2, SV, σ, ep = p
   
       αT=1.7e-4; αS=0.8e-3; TaN=7.0; TaL=25.0; T0=2.65
       S0=35.0; SB=34.538; g=3.17e-8
       VN=7.2106e15; VL=6.3515e16; VB=18VN
       ϕ=αT*(TaN-T0); x_E=(TaL-T0)/(TaN-T0)
       y_B=αS/αT*(SB-S0)/(TaN-T0)
       W=5.456*SV/σ; δN=g*VN/σ; k1=1e-2*SV/σ
       k2n=k2; k2l=k2/2; Ψ=ϕ*((y_N-x_N)-(y_E-x_E))
       κN=k1+0.5*(k2n-k1)*(1+tanh(((y_N-y_B)-(x_N-x_B)-η)/ep))
       κE=k1+0.5*(k2l-k1)*(1+tanh(((y_E-y_B)-(x_E-x_B)-η)/ep))
       Ψh=Ψ*tanh(Ψ/ep)
   
       TVE=VL/VN; TVB=VB/VN
       du[1]=m+Ψ/2*(y_E-y_B)+Ψh/2*(y_E+y_B-2y_N)+W*(y_E-y_N)-κN*(y_N-y_B)
       du[2]=(-m+Ψ/2*(y_B-y_N)+Ψh/2*(y_B+y_N-2y_E)-W*(y_E-y_N)-κE*(y_E-y_B))/TVE
       du[3]=(Ψ/2*(x_N-x_E)+Ψh/2*(x_N+x_E-2x_B)+κN*(x_N-x_B)+κE*(x_E-x_B))/TVB
       du[4]=-δN*(x_N-1)+Ψ/2*(x_E-x_B)+Ψh/2*(x_E+x_B-2x_N)+W*(x_E-x_N)-κN*(x_N-x_B)
end

# ---------- Jacobian (in-place) --------------------------------

function AutoSysJac!(J::AbstractMatrix, u::AbstractVector, p::SVector{6,Float64}, t)
        m, η, k2, SV, σ, ep = p
        y_N, y_E, x_B, x_N = u

        αT=1.7e-4; αS=0.8e-3; TaN=7; TaL=25; T0=2.65; S0=35; SB=34.538
        y_B=αS/αT*(SB-S0)/(TaN-T0)
        g=3.17e-8; VN=7.2106e15; VL=6.3515e16; VB=18VN
        ϕ=αT*(TaN-T0); x_E=(TaL-T0)/(TaN-T0)
        W=5.456*SV/σ; δN=g*VN/σ; k1=1e-2*SV/σ; k2n=k2; k2l=k2/2
        Ψ=ϕ*((y_N-x_N)-(y_E-x_E)); TVE=VL/VN; TVB=VB/VN

        c1=cosh((η-x_B+x_N+y_B-y_N)/ep)
        c2=cosh(ϕ*(x_E-x_N-y_E+y_N)/ep)
        t1=tanh((η-x_B+x_N+y_B-y_N)/ep)
        t2=tanh(ϕ*(x_E-x_N-y_E+y_N)/ep)

        J[1,1]=(-2k1-2W+ϕ*(-y_B+y_E)-((k1-k2n)*(y_B-y_N)/c1^2)/ep+
               (ϕ^2*(y_B+y_E-2y_N)*(x_E-x_N-y_E+y_N)/c2^2)/ep-
               (k1-k2n)*(-1+t1)+ϕ*(y_B+y_E-2y_N)*t2-
               2ϕ*(x_E-x_N-y_E+y_N)*t2)/2

        J[1,2]=(2W+ϕ*(y_B-y_E)+ϕ*(x_E-x_N-y_E+y_N)+
               (ϕ^2*(y_B+y_E-2y_N)*(-x_E+x_N+y_E-y_N)/c2^2)/ep-
               ϕ*(y_B+y_E-2y_N)*t2+ϕ*(x_E-x_N-y_E+y_N)*t2)/2

        J[1,3]=-((k1-k2n)*(y_B-y_N)/c1^2)/(2ep)

        J[1,4]=(ep*ϕ*(y_B-y_E)+(k1-k2n)*(y_B-y_N)/c1^2+
               ϕ^2*(y_B+y_E-2y_N)*(-x_E+x_N+y_E-y_N)/c2^2-
               ep*ϕ*(y_B+y_E-2y_N)*t2)/(2ep)

        J[2,1]=(2W+ϕ*(y_B-y_N)-ϕ*(x_E-x_N-y_E+y_N)+
               (ϕ^2*(y_B-2y_E+y_N)*(x_E-x_N-y_E+y_N)/c2^2)/ep+
               ϕ*(y_B-2y_E+y_N)*t2+ϕ*(x_E-x_N-y_E+y_N)*t2)/(2TVE)

        J[2,2]=-(ep*k1+ep*W+0.5ep*ϕ*(y_B-y_N)+
                0.5(k1-k2l)*(y_B-y_E)/cosh((η-x_B+x_E+y_B-y_E)/ep)^2-
                0.5ϕ^2*(-x_E+x_N+y_E-y_N)*(y_B-2y_E+y_N)/c2^2+
                0.5ep*(k1-k2l)*(-1+tanh((η-x_B+x_E+y_B-y_E)/ep))+
                0.5ep*ϕ*(y_B-2y_E+y_N)*t2+
                ep*ϕ*(x_E-x_N-y_E+y_N)*t2)/(ep*TVE)

        J[2,3]=-((k1-k2l)*(y_B-y_E)/cosh((η-x_B+x_E+y_B-y_E)/ep)^2)/(2ep*TVE)

        J[2,4]=ϕ*(-y_B+y_N+(ϕ*(-x_E+x_N+y_E-y_N)*(y_B-2y_E+y_N)/c2^2)/ep-
                 (y_B-2y_E+y_N)*t2)/(2TVE)

        J[3,1]=(ϕ*(-x_E+x_N)+((k1-k2n)*(x_B-x_N)/c1^2)/ep+
               (ϕ^2*(-2x_B+x_E+x_N)*(x_E-x_N-y_E+y_N)/c2^2)/ep+
               ϕ*(-2x_B+x_E+x_N)*t2)/(2TVB)

        J[3,2]=-0.5(ϕ*(-x_E+x_N)-((k1-k2l)*(x_B-x_E)/cosh((η-x_B+x_E+y_B-y_E)/ep)^2)/ep+
                    (ϕ^2*(-2x_B+x_E+x_N)*(x_E-x_N-y_E+y_N)/c2^2)/ep+
                    ϕ*(-2x_B+x_E+x_N)*t2)/TVB

        J[3,3]=-(2k1-((0.5(k1-k2l)*(x_B-x_E)/cosh((η-x_B+x_E+y_B-y_E)/ep)^2)/ep)+
                 ((0.5(k1-k2n)*(x_B-x_N)/c1^2)/ep)+
                 0.5(k1-k2l)*(-1+tanh((η-x_B+x_E+y_B-y_E)/ep))+
                 0.5(k1-k2n)*(-1+t1)+ϕ*(x_E-x_N-y_E+y_N)*t2)/TVB

        J[3,4]=(2k1+ϕ*(x_E-x_N)+ϕ*(x_E-x_N-y_E+y_N)-
               ((k1-k2n)*(x_B-x_N)/c1^2)/ep-
               (ϕ^2*(-2x_B+x_E+x_N)*(x_E-x_N-y_E+y_N)/c2^2)/ep+
               (k1-k2n)*(-1+t1)-ϕ*(-2x_B+x_E+x_N)*t2+
               ϕ*(x_E-x_N-y_E+y_N)*t2)/(2TVB)

        J[4,1]=(ϕ*(-x_B+x_E)-((k1-k2n)*(x_B-x_N)/c1^2)/ep+
               (ϕ^2*(x_B+x_E-2x_N)*(x_E-x_N-y_E+y_N)/c2^2)/ep+
               ϕ*(x_B+x_E-2x_N)*t2)/2

        J[4,2]=ϕ*(x_B-x_E-(ϕ*(x_B+x_E-2x_N)*(x_E-x_N-y_E+y_N)/c2^2)/ep-
                 (x_B+x_E-2x_N)*t2)/2

        J[4,3]=k1-ϕ*(x_E-x_N-y_E+y_N)/2-((0.5(k1-k2n)*(x_B-x_N)/c1^2)/ep)+
               0.5(k1-k2n)*(-1+t1)+ϕ*(x_E-x_N-y_E+y_N)*t2/2

        J[4,4]=-(δN*ep+ep*k1+ep*W-0.5ep*ϕ*(x_B-x_E)-
                 ((0.5(k1-k2n)*(x_B-x_N)/c1^2))+
                 0.5ϕ^2*(x_B+x_E-2x_N)*(x_E-x_N-y_E+y_N)/c2^2+
                 0.5ep*(k1-k2n)*(-1+t1)+0.5ep*ϕ*(x_B+x_E-2x_N)*t2+
                 ep*ϕ*(x_E-x_N-y_E+y_N)*t2)/ep
    end