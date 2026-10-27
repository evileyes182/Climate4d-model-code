!----------------------------------------------------------------------
!----------------------------------------------------------------------
!   varRun :  system (10) for AUTO-07p, U = (y_N, y_E, x_B, x_N),
!             PAR(1:4) = (mu, eta, kappa_N, epsilon), PAR(14:15) = min, max Psi.
!             Compiled and run by run.py (python3 run.py).
!----------------------------------------------------------------------
!----------------------------------------------------------------------
SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP)
!     ---------- ----
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)  
      DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM, NDIM), DFDP(NDIM, *)

      ! Declare variables with double precision
      DOUBLE PRECISION :: y_N, y_E, x_B, x_N, m, eta, ep, k2
      DOUBLE PRECISION :: TaN, TaL, T0, S0, SB, sig_noSV
      DOUBLE PRECISION :: alphaT, alphaS, y_B, g, x_E, TildeVE, TildeVB
      DOUBLE PRECISION :: W, ps, delta_N, Psi, k2n, k2l, k1, kappaN, kappaE 
      DOUBLE PRECISION :: SV, sig, phi, VN, VL, VB

      ! Constants and parameters
      SV     = 1.0D6
      sig    = 2.1D4 * SV
      sig_noSV = 2.1D4
      alphaT = 1.7D-4
      alphaS = 0.8D-3
      TaN    = 7.0D0
      TaL    = 25.0D0
      T0     = 2.65D0
      S0     = 35.0D0
      SB     = 34.538D0
      VN     = 7.2106D15
      VL     = 6.3515D16
      VB     = 18.0D0 * VN
      g      = 3.17D-8
          
      ! Computed parameters
      y_B    = alphaS / alphaT * (SB - S0) / (TaN - T0)
      phi    = alphaT * (TaN - T0)
      x_E    = (TaL - T0) / (TaN - T0)
      TildeVE = 6.3515D1 / 7.2106D0  ! VL / VN
      TildeVB = 18.0D0
      W      = 5.456D0 / sig_noSV
      delta_N = g * VN / sig
    
      ! State variables 
      y_N = U(1)
      y_E = U(2) 
      x_B = U(3)
      x_N = U(4)
    
      ! Model parameters 
      m   = PAR(1)
      eta = PAR(2)
      k2  = PAR(3)
      ep  = PAR(4)
    
      ! Advective functions
      ps = phi * ((y_N - x_N) - (y_E - x_E))
      Psi = ps * TANH(ps / ep)
    
      k2n = k2
      k2l = k2 / 2.0D0
      k1  = 1.0D-2 / sig_noSV
      kappaN = k1 + 0.5D0 * (k2n - k1) * &
               (1.0D0 + TANH(((y_N - y_B) - (x_N - x_B) - eta) / ep))
      kappaE = k1 + 0.5D0 * (k2l - k1) * &
               (1.0D0 + TANH(((y_E - y_B) - (x_E - x_B) - eta) / ep))
    
      F(1) = m + ps / 2.0D0 * (y_E - y_B) + &
             Psi / 2.0D0 * (y_E + y_B - 2.0D0 * y_N) + &
             W * (y_E - y_N) - kappaN * (y_N - y_B)
      F(2) = 1.0D0 / TildeVE * (-m + ps / 2.0D0 * (y_B - y_N) + &
             Psi / 2.0D0 * (y_B + y_N - 2.0D0 * y_E) - &
             W * (y_E - y_N) - kappaE * (y_E - y_B))
      F(3) = 1.0D0 / TildeVB * (ps / 2.0D0 * (x_N - x_E) + &
             Psi / 2.0D0 * (x_N + x_E - 2.0D0 * x_B) + &
             kappaN * (x_N - x_B) + kappaE * (x_E - x_B))
      F(4) = -delta_N * (x_N - 1.0D0) + ps / 2.0D0 * (x_E - x_B) + &
             Psi / 2.0D0 * (x_E + x_B - 2.0D0 * x_N) + &
             W * (x_E - x_N) - kappaN * (x_N - x_B)

!------------------------------------------------------------
!------- STATE JACOBIAN
!------------------------------------------------------------
 IF(IJAC.EQ.0)RETURN

DFDU(1,1) = (-2*k1 - 2*W + phi*(-y_B + y_E) - (1.*(k1 - 1.*k2n)*(y_B - 1.*y_N)* &
           (1 - tanh((eta - x_B + x_N + y_B - y_N)/ep)**2))/ep + &
           (phi**2*(y_B + y_E - 2*y_N)*(x_E - x_N - y_E + y_N)* &
           (1 - tanh((phi*(x_E - x_N - y_E + y_N))/ep)**2))/ep - &
           1.*(k1 - 1.*k2n)*(-1. + tanh((eta - x_B + x_N + y_B - y_N)/ep)) + &
           phi*(y_B + y_E - 2*y_N)*tanh((phi*(x_E - x_N - y_E + y_N))/ep) - &
           2*phi*(x_E - x_N - y_E + y_N)*tanh((phi*(x_E - x_N - y_E + y_N))/ep))/2.
!------------------------------------------------------------
DFDU(1,2) = (2*W + phi*(y_B - y_E) + phi*(x_E - x_N - y_E + y_N) + &
             (phi**2*(y_B + y_E - 2*y_N)*(-x_E + x_N + y_E - y_N)* &
             (1 - tanh((phi*(x_E - x_N - y_E + y_N))/ep)**2))/ep - &
             phi*(y_B + y_E - 2*y_N)*tanh((phi*(x_E - x_N - y_E + y_N))/ep) + &
             phi*(x_E - x_N - y_E + y_N)*tanh((phi*(x_E - x_N - y_E + y_N))/ep))/2.
!------------------------------------------------------------
DFDU(1,3) = (-0.5*(k1 - 1.*k2n)*(y_B - 1.*y_N)* &
            (1 - tanh((eta - x_B + x_N + y_B - y_N)/ep)**2))/ep
!------------------------------------------------------------
DFDU(1,4) = (ep*phi*(y_B - y_E) + 1.*(k1 - 1.*k2n)*(y_B - 1.*y_N)* &
             (1 - tanh((eta - x_B + x_N + y_B - y_N)/ep)**2) + &
             1.*phi**2*(y_B + y_E - 2.*y_N)*(-1.*x_E + x_N + y_E - 1.*y_N)* &
             (1 - tanh((phi*(x_E - x_N - y_E + y_N))/ep)**2) - &
             1.*ep*phi*(y_B + y_E - 2.*y_N)*tanh((phi*(x_E - x_N - y_E + y_N))/ep))/(2.*ep)
!------------------------------------------------------------
!------------------------------------------------------------
!------------------------------------------------------------
!------------------------------------------------------------
DFDU(2,1) = (2*W + phi*(y_B - y_N) - phi*(x_E - x_N - y_E + y_N) + &
             (phi**2*(y_B - 2*y_E + y_N)*(x_E - x_N - y_E + y_N)* &
             (1 - tanh((phi*(x_E - x_N - y_E + y_N))/ep)**2))/ep + &
             phi*(y_B - 2*y_E + y_N)*tanh((phi*(x_E - x_N - y_E + y_N))/ep) + &
             phi*(x_E - x_N - y_E + y_N)*tanh((phi*(x_E - x_N - y_E + y_N))/ep))/(2.*TildeVE)
!------------------------------------------------------------
DFDU(2,2) = (-1.*(ep*k1 + ep*W + 0.5*ep*phi*(y_B - y_N) + &
             0.5*(k1 - k2l)*(y_B - y_E)*(1 - tanh((eta - x_B + x_E + y_B - y_E)/ep)**2) - &
             0.5*phi**2*(-x_E + x_N + y_E - y_N)*(y_B - 2*y_E + y_N)* &
             (1 - tanh((phi*(x_E - x_N - y_E + y_N))/ep)**2) + 0.5*ep*(k1 - k2l)* &
             (-1. + tanh((eta - x_B + x_E + y_B - y_E)/ep)) + 0.5*ep*phi*(y_B - 2*y_E + y_N)* &
             tanh((phi*(x_E - x_N - y_E + y_N))/ep) + ep*phi*(x_E - x_N - y_E + y_N)* &
             tanh((phi*(x_E - x_N - y_E + y_N))/ep)))/(ep*TildeVE)
!------------------------------------------------------------
DFDU(2,3) = (-0.5*(k1 - k2l)*(y_B - y_E)* &
            (1 - tanh((eta - x_B + x_E + y_B - y_E)/ep)**2))/(ep*TildeVE)
!------------------------------------------------------------
DFDU(2,4) = (phi*(-y_B + y_N + (phi*(-x_E + x_N + y_E - y_N)*(y_B - 2*y_E + y_N)* &
             (1 - tanh((phi*(x_E - x_N - y_E + y_N))/ep)**2))/ep - &
             (y_B - 2*y_E + y_N)*tanh((phi*(x_E - x_N - y_E + y_N))/ep)))/(2.*TildeVE)
!------------------------------------------------------------
!------------------------------------------------------------
!------------------------------------------------------------
!------------------------------------------------------------
DFDU(3,1) = (phi*(-x_E + x_N) + ((k1 - k2n)*(x_B - x_N)* &
            (1 - tanh((eta - x_B + x_N + y_B - y_N)/ep)**2))/ep + &
            (phi**2*(-2*x_B + x_E + x_N)*(x_E - x_N - y_E + y_N)* &
            (1 - tanh((phi*(x_E - x_N - y_E + y_N))/ep)**2))/ep + &
            phi*(-2*x_B + x_E + x_N)*tanh((phi*(x_E - x_N - y_E + y_N))/ep))/(2.*TildeVB)
!------------------------------------------------------------
DFDU(3,2) = -0.5*(phi*(-x_E + x_N) - ((k1 - k2l)*(x_B - x_E)* &
             (1 - tanh((eta - x_B + x_E + y_B - y_E)/ep)**2))/ep + &
             (phi**2*(-2*x_B + x_E + x_N)*(x_E - x_N - y_E + y_N)* &
             (1 - tanh((phi*(x_E - x_N - y_E + y_N))/ep)**2))/ep + &
             phi*(-2*x_B + x_E + x_N)*tanh((phi*(x_E - x_N - y_E + y_N))/ep))/TildeVB
!------------------------------------------------------------
DFDU(3,3) = -((2*k1 - (0.5*(k1 - k2l)*(x_B - x_E)* &
             (1 - tanh((eta - x_B + x_E + y_B - y_E)/ep)**2))/ep - &
             (0.5*(k1 - k2n)*(x_B - x_N)* &
             (1 - tanh((eta - x_B + x_N + y_B - y_N)/ep)**2))/ep + &
             0.5*(k1 - k2l)*(-1. + tanh((eta - x_B + x_E + y_B - y_E)/ep)) + &
             0.5*(k1 - k2n)*(-1. + tanh((eta - x_B + x_N + y_B - y_N)/ep)) + &
             phi*(x_E - x_N - y_E + y_N)*tanh((phi*(x_E - x_N - y_E + y_N))/ep))/TildeVB)
!------------------------------------------------------------
DFDU(3,4) = (2*k1 + phi*(x_E - x_N) + phi*(x_E - x_N - y_E + y_N) - &
             ((k1 - k2n)*(x_B - x_N)* &
             (1 - tanh((eta - x_B + x_N + y_B - y_N)/ep)**2))/ep - &
             (phi**2*(-2*x_B + x_E + x_N)*(x_E - x_N - y_E + y_N)* &
             (1 - tanh((phi*(x_E - x_N - y_E + y_N))/ep)**2))/ep + &
             (k1 - k2n)*(-1. + tanh((eta - x_B + x_N + y_B - y_N)/ep)) - &
             phi*(-2*x_B + x_E + x_N)*tanh((phi*(x_E - x_N - y_E + y_N))/ep) + &
             phi*(x_E - x_N - y_E + y_N)*tanh((phi*(x_E - x_N - y_E + y_N))/ep))/(2.*TildeVB)
!------------------------------------------------------------
!------------------------------------------------------------
!------------------------------------------------------------
!------------------------------------------------------------
DFDU(4,1) = (phi*(-x_B + x_E) - ((k1 - k2n)*(x_B - x_N)* &
            (1 - tanh((eta - x_B + x_N + y_B - y_N)/ep)**2))/ep + &
            (phi**2*(x_B + x_E - 2*x_N)*(x_E - x_N - y_E + y_N)* &
            (1 - tanh((phi*(x_E - x_N - y_E + y_N))/ep)**2))/ep + &
            phi*(x_B + x_E - 2*x_N)*tanh((phi*(x_E - x_N - y_E + y_N))/ep))/2.
!------------------------------------------------------------
DFDU(4,2) = (phi*(x_B - x_E - (phi*(x_B + x_E - 2*x_N)*(x_E - x_N - y_E + y_N)* &
             (1 - tanh((phi*(x_E - x_N - y_E + y_N))/ep)**2))/ep - &
             (x_B + x_E - 2*x_N)*tanh((phi*(x_E - x_N - y_E + y_N))/ep)))/2.
!------------------------------------------------------------
DFDU(4,3) = k1 - (phi*(x_E - x_N - y_E + y_N))/2. - (0.5*(k1 - k2n)*(x_B - x_N)* &
            (1 - tanh((eta - x_B + x_N + y_B - y_N)/ep)**2))/ep + &
            0.5*(k1 - k2n)*(-1. + tanh((eta - x_B + x_N + y_B - y_N)/ep)) + &
            (phi*(x_E - x_N - y_E + y_N)*tanh((phi*(x_E - x_N - y_E + y_N))/ep))/2.
!------------------------------------------------------------
DFDU(4,4) = (-1.*(delta_N*ep + ep*k1 + ep*W - 0.5*ep*phi*(x_B - x_E) - &
             0.5*(k1 - k2n)*(x_B - x_N)* &
             (1 - tanh((eta - x_B + x_N + y_B - y_N)/ep)**2) + &
             0.5*phi**2*(x_B + x_E - 2.*x_N)*(x_E - x_N - y_E + y_N)* &
             (1 - tanh((phi*(x_E - x_N - y_E + y_N))/ep)**2) + &
             0.5*ep*(k1 - k2n)*(-1. + tanh((eta - x_B + x_N + y_B - y_N)/ep)) + &
             0.5*ep*phi*(x_B + x_E - 2.*x_N)*tanh((phi*(x_E - x_N - y_E + y_N))/ep) + &
             ep*phi*(x_E - x_N - y_E + y_N)*tanh((phi*(x_E - x_N - y_E + y_N))/ep)))/ep
!------------------------------------------------------------

 IF(IJAC.EQ.1)RETURN

!------------------------------------------------------------
!------- PARAMETER JACOBIAN
!------------------------------------------------------------
DFDP(1,1) = 1
!------------------------------------------------------------
DFDP(1,2) = (0.5*(k1 - 1.*k2)*(y_B - 1.*y_N)* &
            (1 - tanh((eta - x_B + x_N + y_B - y_N)/ep)**2))/ep
!------------------------------------------------------------
DFDP(1,3) = -0.5*(-y_B + y_N)*(1 - tanh((eta - x_B + x_N + y_B - y_N)/ep))
!------------------------------------------------------------
DFDP(1,4) = -0.5*(1.*(k1 - 1.*k2)*(y_B - 1.*y_N)*(eta - 1.*x_B + x_N + &
            y_B - 1.*y_N)*(1 - tanh((eta - x_B + x_N + y_B - y_N)/ep)**2) + &
            phi**2*(y_B + y_E - 2*y_N)*(x_E - x_N - y_E + y_N)**2*(1 - &
            tanh((phi*(x_E - x_N - y_E + y_N))/ep)**2))/ep**2
!------------------------------------------------------------
!------------------------------------------------------------
!------------------------------------------------------------
!------------------------------------------------------------
DFDP(2,1) = -(1/TildeVE)
!------------------------------------------------------------
DFDP(2,2) = (0.5*(k1 - 0.5*k2)*(y_B - 1.*y_E)* &
            (1 - tanh((eta - x_B + x_E + y_B - y_E)/ep)**2))/(ep*TildeVE)
!------------------------------------------------------------
DFDP(2,3) = (-0.25*(-y_B + y_E)*(1 - tanh((eta - x_B + x_E + y_B - y_E)/ep)))/TildeVE
!------------------------------------------------------------
DFDP(2,4) = -0.5*(1.*(k1 - 0.5*k2)*(y_B - 1.*y_E)*(eta - 1.*x_B + x_E + y_B - 1.*y_E)* &
            (1 - tanh((eta - x_B + x_E + y_B - y_E)/ep)**2) + &
            phi**2*(y_B - 2*y_E + y_N)*(x_E - x_N - y_E + y_N)**2*(1 - &
            tanh((phi*(x_E - x_N - y_E + y_N))/ep)**2))/(ep**2*TildeVE)
!------------------------------------------------------------
!------------------------------------------------------------
!------------------------------------------------------------
!------------------------------------------------------------
DFDP(3,1) = 0
!------------------------------------------------------------
DFDP(3,2) = (-0.5*(k1 - 0.5*k2)*(x_B - 1.*x_E)* &
            (1 - tanh((eta - x_B + x_E + y_B - y_E)/ep)**2) - &
            0.5*(k1 - 1.*k2)*(x_B - 1.*x_N)* &
            (1 - tanh((eta - x_B + x_N + y_B - y_N)/ep)**2))/(ep*TildeVB)
!------------------------------------------------------------
DFDP(3,3) = (0.25*(-x_B + x_E)*(1 - tanh((eta - x_B + x_E + y_B - y_E)/ep)) + &
             0.5*(-x_B + x_N)*(1 - tanh((eta - x_B + x_N + y_B - y_N)/ep)))/TildeVB
!------------------------------------------------------------
DFDP(3,4) = -0.5*(1.*(k1 - 0.5*k2)*(x_B - 1.*x_E)*(-1.*eta + x_B - 1.*x_E - 1.*y_B + y_E)* &
            (1 - tanh((eta - x_B + x_E + y_B - y_E)/ep)**2) + &
            1.*(k1 - 1.*k2)*(x_B - 1.*x_N)*(-1.*eta + x_B - 1.*x_N - 1.*y_B + y_N)* &
            (1 - tanh((eta - x_B + x_N + y_B - y_N)/ep)**2) + &
            phi**2*(-2*x_B + x_E + x_N)*(x_E - x_N - y_E + y_N)**2*(1 - &
            tanh((phi*(x_E - x_N - y_E + y_N))/ep)**2))/(ep**2*TildeVB)
!------------------------------------------------------------
!------------------------------------------------------------
!------------------------------------------------------------
!------------------------------------------------------------
DFDP(4,1) = 0
!------------------------------------------------------------
DFDP(4,2) = (0.5*(k1 - 1.*k2)*(x_B - 1.*x_N)* &
            (1 - tanh((eta - x_B + x_N + y_B - y_N)/ep)**2))/ep
!------------------------------------------------------------
DFDP(4,3) = -0.5*(-x_B + x_N)*(1 - tanh((eta - x_B + x_N + y_B - y_N)/ep))
!------------------------------------------------------------
DFDP(4,4) = -0.5*(-1.*(k1 - 1.*k2)*(x_B - 1.*x_N)*(-1.*eta + x_B - 1.*x_N - 1.*y_B + y_N)* &
            (1 - tanh((eta - x_B + x_N + y_B - y_N)/ep)**2) + &
            phi**2*(x_B + x_E - 2*x_N)*(x_E - x_N - y_E + y_N)**2*(1 - &
            tanh((phi*(x_E - x_N - y_E + y_N))/ep)**2))/ep**2

END SUBROUTINE FUNC
!----------------------------------------------------------------------

      SUBROUTINE STPNT(NDIM,U,PAR,T)
!     ---------- -----

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: T

      PAR(1) = -0.002006D0
      PAR(2) = -1.5D0
      PAR(3) = 8.0 / (2.1D4)
      PAR(4) = 0.02D0

      U(1) = -3.5435D0
      U(2) = 0.9597D0
      U(3) = 3.8512D0
      U(4) = 1.1678D0


      END SUBROUTINE STPNT
!----------------------------------------------------------------------

      SUBROUTINE BCND
      END SUBROUTINE BCND
!----------------------------------------------------------------------

      SUBROUTINE ICND
      END SUBROUTINE ICND
!----------------------------------------------------------------------
  

      SUBROUTINE FOPT
      END SUBROUTINE FOPT
!----------------------------------------------------------------------

SUBROUTINE PVLS(NDIM,U,PAR)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NDIM
    DOUBLE PRECISION, INTENT(IN) :: U(NDIM)
    DOUBLE PRECISION, INTENT(INOUT) :: PAR(*)
    DOUBLE PRECISION, EXTERNAL :: GETP
    INTEGER :: NDX, NCOL, NTST
    
    ! Get mesh information
    NDX = NINT(GETP('NDX',0,U))
    NTST = NINT(GETP('NTST',0,U))
    NCOL = NINT(GETP('NCOL',0,U))
    
    ! Set standard solution measures
    PAR(12) = GETP('NDIM',0,U)
    PAR(13) = GETP('STA',0,U)
    
    ! Minimum and maximum of PS over all NTST*NCOL+1 mesh points (one point when NTST = 0)
    CALL PSMINMAX(NDX,NTST*NCOL,U,PAR(14),PAR(15))
END SUBROUTINE PVLS


! Solution dimensioned as U(NDX,0:NTST*NCOL), as in the AUTO demo pvl
SUBROUTINE PSMINMAX(NDX,M,U,PS_MIN,PS_MAX)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NDX, M
    DOUBLE PRECISION, INTENT(IN) :: U(NDX,0:M)
    DOUBLE PRECISION, INTENT(OUT) :: PS_MIN, PS_MAX
    
    ! Declare model variables
    DOUBLE PRECISION :: alphaT, TaN, TaL, T0, phi, x_E
    DOUBLE PRECISION :: PS(0:M)
    
    ! Initialize key parameters needed for PS calculation
    alphaT = 1.7D-4
    TaN = 7.0D0
    TaL = 25.0D0
    T0 = 2.65D0
    phi = alphaT * (TaN - T0)
    x_E = (TaL - T0) / (TaN - T0)
    
    PS = phi * ((U(1,:) - U(4,:)) - (U(2,:) - x_E))
    PS_MIN = MINVAL(PS)
    PS_MAX = MAXVAL(PS)
END SUBROUTINE PSMINMAX
