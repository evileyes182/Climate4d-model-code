SUBROUTINE RHS(U,PAR,F,JAC,A)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: U(4), PAR(*)
    LOGICAL, INTENT(IN) :: JAC
    DOUBLE PRECISION, INTENT(OUT) :: F(4), A(4,4)
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

    IF(JAC)THEN
        A(1,1) = (-2*k1 - 2*W + phi*(-y_B + y_E) - (1.*(k1 - 1.*k2n)*(y_B - 1.*y_N)* &
                (1 - tanh((eta - x_B + x_N + y_B - y_N)/ep)**2))/ep + &
                (phi**2*(y_B + y_E - 2*y_N)*(x_E - x_N - y_E + y_N)* &
                (1 - tanh((phi*(x_E - x_N - y_E + y_N))/ep)**2))/ep - &
                1.*(k1 - 1.*k2n)*(-1. + tanh((eta - x_B + x_N + y_B - y_N)/ep)) + &
                phi*(y_B + y_E - 2*y_N)*tanh((phi*(x_E - x_N - y_E + y_N))/ep) - &
                2*phi*(x_E - x_N - y_E + y_N)*tanh((phi*(x_E - x_N - y_E + y_N))/ep))/2.
        !------------------------------------------------------------
        A(1,2) = (2*W + phi*(y_B - y_E) + phi*(x_E - x_N - y_E + y_N) + &
                    (phi**2*(y_B + y_E - 2*y_N)*(-x_E + x_N + y_E - y_N)* &
                    (1 - tanh((phi*(x_E - x_N - y_E + y_N))/ep)**2))/ep - &
                    phi*(y_B + y_E - 2*y_N)*tanh((phi*(x_E - x_N - y_E + y_N))/ep) + &
                    phi*(x_E - x_N - y_E + y_N)*tanh((phi*(x_E - x_N - y_E + y_N))/ep))/2.
        !------------------------------------------------------------
        A(1,3) = (-0.5*(k1 - 1.*k2n)*(y_B - 1.*y_N)* &
                    (1 - tanh((eta - x_B + x_N + y_B - y_N)/ep)**2))/ep
        !------------------------------------------------------------
        A(1,4) = (ep*phi*(y_B - y_E) + 1.*(k1 - 1.*k2n)*(y_B - 1.*y_N)* &
                    (1 - tanh((eta - x_B + x_N + y_B - y_N)/ep)**2) + &
                    1.*phi**2*(y_B + y_E - 2.*y_N)*(-1.*x_E + x_N + y_E - 1.*y_N)* &
                    (1 - tanh((phi*(x_E - x_N - y_E + y_N))/ep)**2) - &
                    1.*ep*phi*(y_B + y_E - 2.*y_N)*tanh((phi*(x_E - x_N - y_E + y_N))/ep))/(2.*ep)
        !------------------------------------------------------------
        !------------------------------------------------------------
        !------------------------------------------------------------
        !------------------------------------------------------------
        A(2,1) = (2*W + phi*(y_B - y_N) - phi*(x_E - x_N - y_E + y_N) + &
                    (phi**2*(y_B - 2*y_E + y_N)*(x_E - x_N - y_E + y_N)* &
                    (1 - tanh((phi*(x_E - x_N - y_E + y_N))/ep)**2))/ep + &
                    phi*(y_B - 2*y_E + y_N)*tanh((phi*(x_E - x_N - y_E + y_N))/ep) + &
                    phi*(x_E - x_N - y_E + y_N)*tanh((phi*(x_E - x_N - y_E + y_N))/ep))/(2.*TildeVE)
        !------------------------------------------------------------
        A(2,2) = (-1.*(ep*k1 + ep*W + 0.5*ep*phi*(y_B - y_N) + &
                    0.5*(k1 - k2l)*(y_B - y_E)*(1 - tanh((eta - x_B + x_E + y_B - y_E)/ep)**2) - &
                    0.5*phi**2*(-x_E + x_N + y_E - y_N)*(y_B - 2*y_E + y_N)* &
                    (1 - tanh((phi*(x_E - x_N - y_E + y_N))/ep)**2) + 0.5*ep*(k1 - k2l)* &
                    (-1. + tanh((eta - x_B + x_E + y_B - y_E)/ep)) + 0.5*ep*phi*(y_B - 2*y_E + y_N)* &
                    tanh((phi*(x_E - x_N - y_E + y_N))/ep) + ep*phi*(x_E - x_N - y_E + y_N)* &
                    tanh((phi*(x_E - x_N - y_E + y_N))/ep)))/(ep*TildeVE)
        !------------------------------------------------------------
        A(2,3) = (-0.5*(k1 - k2l)*(y_B - y_E)* &
                    (1 - tanh((eta - x_B + x_E + y_B - y_E)/ep)**2))/(ep*TildeVE)
        !------------------------------------------------------------
        A(2,4) = (phi*(-y_B + y_N + (phi*(-x_E + x_N + y_E - y_N)*(y_B - 2*y_E + y_N)* &
                    (1 - tanh((phi*(x_E - x_N - y_E + y_N))/ep)**2))/ep - &
                    (y_B - 2*y_E + y_N)*tanh((phi*(x_E - x_N - y_E + y_N))/ep)))/(2.*TildeVE)
        !------------------------------------------------------------
        !------------------------------------------------------------
        !------------------------------------------------------------
        !------------------------------------------------------------
        A(3,1) = (phi*(-x_E + x_N) + ((k1 - k2n)*(x_B - x_N)* &
                    (1 - tanh((eta - x_B + x_N + y_B - y_N)/ep)**2))/ep + &
                    (phi**2*(-2*x_B + x_E + x_N)*(x_E - x_N - y_E + y_N)* &
                    (1 - tanh((phi*(x_E - x_N - y_E + y_N))/ep)**2))/ep + &
                    phi*(-2*x_B + x_E + x_N)*tanh((phi*(x_E - x_N - y_E + y_N))/ep))/(2.*TildeVB)
        !------------------------------------------------------------
        A(3,2) = -0.5*(phi*(-x_E + x_N) - ((k1 - k2l)*(x_B - x_E)* &
                    (1 - tanh((eta - x_B + x_E + y_B - y_E)/ep)**2))/ep + &
                    (phi**2*(-2*x_B + x_E + x_N)*(x_E - x_N - y_E + y_N)* &
                    (1 - tanh((phi*(x_E - x_N - y_E + y_N))/ep)**2))/ep + &
                    phi*(-2*x_B + x_E + x_N)*tanh((phi*(x_E - x_N - y_E + y_N))/ep))/TildeVB
        !------------------------------------------------------------
        A(3,3) = -((2*k1 - (0.5*(k1 - k2l)*(x_B - x_E)* &
                    (1 - tanh((eta - x_B + x_E + y_B - y_E)/ep)**2))/ep - &
                    (0.5*(k1 - k2n)*(x_B - x_N)* &
                    (1 - tanh((eta - x_B + x_N + y_B - y_N)/ep)**2))/ep + &
                    0.5*(k1 - k2l)*(-1. + tanh((eta - x_B + x_E + y_B - y_E)/ep)) + &
                    0.5*(k1 - k2n)*(-1. + tanh((eta - x_B + x_N + y_B - y_N)/ep)) + &
                    phi*(x_E - x_N - y_E + y_N)*tanh((phi*(x_E - x_N - y_E + y_N))/ep))/TildeVB)
        !------------------------------------------------------------
        A(3,4) = (2*k1 + phi*(x_E - x_N) + phi*(x_E - x_N - y_E + y_N) - &
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
        A(4,1) = (phi*(-x_B + x_E) - ((k1 - k2n)*(x_B - x_N)* &
                    (1 - tanh((eta - x_B + x_N + y_B - y_N)/ep)**2))/ep + &
                    (phi**2*(x_B + x_E - 2*x_N)*(x_E - x_N - y_E + y_N)* &
                    (1 - tanh((phi*(x_E - x_N - y_E + y_N))/ep)**2))/ep + &
                    phi*(x_B + x_E - 2*x_N)*tanh((phi*(x_E - x_N - y_E + y_N))/ep))/2.
        !------------------------------------------------------------
        A(4,2) = (phi*(x_B - x_E - (phi*(x_B + x_E - 2*x_N)*(x_E - x_N - y_E + y_N)* &
                    (1 - tanh((phi*(x_E - x_N - y_E + y_N))/ep)**2))/ep - &
                    (x_B + x_E - 2*x_N)*tanh((phi*(x_E - x_N - y_E + y_N))/ep)))/2.
        !------------------------------------------------------------
        A(4,3) = k1 - (phi*(x_E - x_N - y_E + y_N))/2. - (0.5*(k1 - k2n)*(x_B - x_N)* &
                    (1 - tanh((eta - x_B + x_N + y_B - y_N)/ep)**2))/ep + &
                    0.5*(k1 - k2n)*(-1. + tanh((eta - x_B + x_N + y_B - y_N)/ep)) + &
                    (phi*(x_E - x_N - y_E + y_N)*tanh((phi*(x_E - x_N - y_E + y_N))/ep))/2.
        !------------------------------------------------------------
        A(4,4) = (-1.*(delta_N*ep + ep*k1 + ep*W - 0.5*ep*phi*(x_B - x_E) - &
                    0.5*(k1 - k2n)*(x_B - x_N)* &
                    (1 - tanh((eta - x_B + x_N + y_B - y_N)/ep)**2) + &
                    0.5*phi**2*(x_B + x_E - 2.*x_N)*(x_E - x_N - y_E + y_N)* &
                    (1 - tanh((phi*(x_E - x_N - y_E + y_N))/ep)**2) + &
                    0.5*ep*(k1 - k2n)*(-1. + tanh((eta - x_B + x_N + y_B - y_N)/ep)) + &
                    0.5*ep*phi*(x_B + x_E - 2.*x_N)*tanh((phi*(x_E - x_N - y_E + y_N))/ep) + &
                    ep*phi*(x_E - x_N - y_E + y_N)*tanh((phi*(x_E - x_N - y_E + y_N))/ep)))/ep
        !------------------------------------------------------------
    ENDIF
END SUBROUTINE RHS


SUBROUTINE D2RHS(U,PAR,H,FP,AP)
    ! Second derivatives of the model: H(i,j,k) = d2 f_i / du_j du_k,
    ! FP(i,1:2) = d f_i / d(mu, eta), AP(i,j) = d2 f_i / du_j deta
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: U(4), PAR(*)
    DOUBLE PRECISION, INTENT(OUT) :: H(4,4,4), FP(4,2), AP(4,4)
    DOUBLE PRECISION :: y_N, y_E, x_B, x_N, m, eta, ep, k2
    DOUBLE PRECISION :: TaN, TaL, T0, S0, SB, sig_noSV
    DOUBLE PRECISION :: alphaT, alphaS, y_B, g, x_E, TildeVE, TildeVB
    DOUBLE PRECISION :: W, delta_N, k1
    DOUBLE PRECISION :: SV, sig, phi, VN, VL, VB
    DOUBLE PRECISION :: q0, q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11
    DOUBLE PRECISION :: q12, q13, q14, q15, q16, q17, q18, q19, q20, q21, q22, q23
    DOUBLE PRECISION :: q24, q25, q26, q27, q28, q29, q30, q31, q32, q33, q34, q35
    DOUBLE PRECISION :: q36, q37, q38, q39, q40, q41, q42, q43, q44, q45, q46, q47
    DOUBLE PRECISION :: q48, q49, q50, q51, q52, q53, q54, q55, q56, q57, q58, q59
    DOUBLE PRECISION :: q60, q61, q62, q63, q64, q65, q66, q67, q68, q69, q70, q71
    DOUBLE PRECISION :: q72, q73, q74, q75, q76, q77, q78, q79, q80, q81, q82, q83
    DOUBLE PRECISION :: q84, q85, q86, q87, q88, q89, q90, q91, q92, q93, q94, q95
    DOUBLE PRECISION :: q96, q97, q98, q99, q100, q101, q102, q103, q104, q105, q106, q107
    DOUBLE PRECISION :: q108, q109, q110, q111, q112, q113, q114, q115, q116, q117, q118, q119
    DOUBLE PRECISION :: q120, q121, q122, q123, q124, q125, q126, q127, q128, q129, q130, q131

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

      k1  = 1.0D-2 / sig_noSV

      q0 = 1d0/ep
      q1 = -k2
      q2 = k1 + q1
      q3 = y_B - y_N
      q4 = eta - x_B
      q5 = q3 + q4 + x_N
      q6 = tanh(q0*q5)**2 - 1
      q7 = q2*q6
      q8 = q0*q7
      q9 = ep**(-2)
      q10 = tanh(q0*q5)
      q11 = y_B + y_E - 2*y_N
      q12 = -x_N
      q13 = -y_E
      q14 = q12 + q13 + x_E + y_N
      q15 = phi*q0
      q16 = q14*q15
      q17 = tanh(q16)
      q18 = q17**2 - 1
      q19 = phi**2*q18
      q20 = q11*q19
      q21 = phi**3*q17
      q22 = q11*q21
      q23 = q14*q9
      q24 = q18*q23
      q25 = q0*q20 - q22*q24
      q26 = phi*q17
      q27 = q0*q19
      q28 = q14*q27
      q29 = 2*q26 - 2*q28
      q30 = q16*q18
      q31 = (3.0d0/2.0d0)*q30
      q32 = q15*q18
      q33 = q17*q23
      q34 = q11*q32 - q20*q33
      q35 = q34 + 1.0d0/2.0d0
      q36 = phi*((3.0d0/2.0d0)*q17 - q31 + q35)
      q37 = -q0*q10*q3 + 1.0d0/2.0d0
      q38 = -q37*q8
      q39 = q2*(q10**2 - 1)
      q40 = q0*q39
      q41 = (1.0d0/2.0d0)*q40
      q42 = q10*q9
      q43 = q3*q39*q42
      q44 = -q43
      q45 = q26 - q28
      q46 = q25 + q41 + q44 + q45
      q47 = q17 - q30
      q48 = (1.0d0/2.0d0)*q17
      q49 = (1.0d0/2.0d0)*phi
      q50 = q0*q14
      q51 = q18*q50
      q52 = q49*q51
      q53 = q48 - q52
      q54 = phi*(-q35 - q53)
      q55 = q0*q10
      q56 = -q17 + q30 + 1
      q57 = y_B - 2*y_E + y_N
      q58 = q32*q57
      q59 = q19*q33*q57
      q60 = q58 - q59
      q61 = 1d0/TildeVE
      q62 = phi*q61
      q63 = q60 + 1.0d0/2.0d0
      q64 = -q62*((3.0d0/2.0d0)*q17 - q31 - q63)
      q65 = -q62*(q48 - q52 - q63)
      q66 = 2*k1 + q1
      q67 = q13 + y_B
      q68 = q4 + q67 + x_E
      q69 = tanh(q0*q68)
      q70 = q66*(q69**2 - 1)
      q71 = (1.0d0/2.0d0)*q70
      q72 = q0*q71
      q73 = q27*q57
      q74 = q67*q69
      q75 = q71*q9
      q76 = q74*q75
      q77 = q21*q24
      q78 = q66*(tanh(q0*q68)**2 - 1)
      q79 = (1.0d0/4.0d0)*q0*q78
      q80 = -2*q0*q74 + 1
      q81 = q61*q80
      q82 = q62*(q47 - q58 + q59)
      q83 = (1.0d0/4.0d0)*q70
      q84 = q0*q83
      q85 = q61*q76
      q86 = -q26*q50 + 1
      q87 = q12 + x_B
      q88 = q10*q8*q87
      q89 = -2*x_B + x_E + x_N
      q90 = q19*q89
      q91 = q21*q51
      q92 = -q89*q91 + q90
      q93 = 1d0/TildeVB
      q94 = q0*q93
      q95 = -q86
      q96 = -q90*q94*q95
      q97 = q7*q87
      q98 = q42*q97
      q99 = q45 + q98
      q100 = (1.0d0/2.0d0)*q8 + q99
      q101 = -q100*q93
      q102 = (1.0d0/2.0d0)*q26
      q103 = -1.0d0/2.0d0*q28
      q104 = q39*q42*q87
      q105 = q104 + q41
      q106 = q0*q90 - q77*q89
      q107 = q93*(q102 + q103 + q105 + q106 + q49)
      q108 = (1.0d0/2.0d0)*q78
      q109 = x_B - x_E
      q110 = q109*q69
      q111 = q108*q110
      q112 = q0*q111
      q113 = -q93*(q111*q9 - q26 + q28 + q79)
      q114 = -phi*q93*(q32*q89 - q33*q90 + q53 + 1.0d0/2.0d0)
      q115 = q104 + q40
      q116 = q93*(q115 + q45)
      q117 = q55*q87
      q118 = q117*q39
      q119 = q117 + 1.0d0/2.0d0
      q120 = q119*q8
      q121 = x_B + x_E - 2*x_N
      q122 = q121*q19
      q123 = q0*q122
      q124 = q123*q86
      q125 = q102 + q103 - q49
      q126 = q105 + q125
      q127 = q121*q77 - q123
      q128 = -q100 - q127
      q129 = q49*q56
      q130 = phi*(-q121*q32 + q122*q33 + q47)
      q131 = -q125 - q8 - q98
      H(1,1,1) = q10*q2*q3*q6*q9 - q25 - q29 - q8
      H(1,1,2) = q36
      H(1,1,3) = q38
      H(1,1,4) = q46
      AP(1,1) = q37*q40
      H(1,2,1) = q36
      H(1,2,2) = phi*(-q34 - q47 - 1)
      H(1,2,3) = 0
      H(1,2,4) = q54
      AP(1,2) = 0
      H(1,3,1) = q38
      H(1,3,2) = 0
      H(1,3,3) = q43
      H(1,3,4) = q44
      AP(1,3) = q44
      H(1,4,1) = q46
      H(1,4,2) = q54
      H(1,4,3) = q44
      H(1,4,4) = q0*(-q20 + q22*q51 + q3*q39*q55)
      AP(1,4) = q43
      FP(1,1) = 1
      FP(1,2) = -q3*q41
      H(2,1,1) = -q62*(q56 + q60)
      H(2,1,2) = q64
      H(2,1,3) = 0
      H(2,1,4) = q65
      AP(2,1) = 0
      H(2,2,1) = q64
      H(2,2,2) = q61*(q29 + q57*q77 - q72 - q73 + q76)
      H(2,2,3) = -q79*q81
      H(2,2,4) = q82
      AP(2,2) = q81*q84
      H(2,3,1) = 0
      H(2,3,2) = -q61*q79*q80
      H(2,3,3) = q85
      H(2,3,4) = 0
      AP(2,3) = -q85
      H(2,4,1) = q65
      H(2,4,2) = q82
      H(2,4,3) = 0
      H(2,4,4) = -q61*q73*q86
      AP(2,4) = 0
      FP(2,1) = -q61
      FP(2,2) = -q61*q67*q84
      H(3,1,1) = -q94*(q88 + q92)
      H(3,1,2) = q96
      H(3,1,3) = q101
      H(3,1,4) = q107
      AP(3,1) = q104*q93
      H(3,2,1) = q96
      H(3,2,2) = -q94*(q112 + q92)
      H(3,2,3) = q113
      H(3,2,4) = q114
      AP(3,2) = q110*q75*q93
      H(3,3,1) = q101
      H(3,3,2) = q113
      H(3,3,3) = -q94*(q108 + q112 + q7 + q88)
      H(3,3,4) = q116
      AP(3,3) = q94*(q110*q72 + q118 + (1.0d0/2.0d0)*q39 + q83)
      H(3,4,1) = q107
      H(3,4,2) = q114
      H(3,4,3) = q116
      H(3,4,4) = -q93*(phi + q106 + q8 + q99)
      AP(3,4) = -q120*q93
      FP(3,1) = 0
      FP(3,2) = (1.0d0/4.0d0)*q94*(q109*q78 + 2*q97)
      H(4,1,1) = q0*(q118 + q121*q91 - q122)
      H(4,1,2) = q124
      H(4,1,3) = q126
      H(4,1,4) = q128
      AP(4,1) = -q104
      H(4,2,1) = q124
      H(4,2,2) = q123*q95
      H(4,2,3) = q129
      H(4,2,4) = q130
      AP(4,2) = 0
      H(4,3,1) = q126
      H(4,3,2) = q129
      H(4,3,3) = q40*(q117 + 1)
      H(4,3,4) = q131
      AP(4,3) = -q120
      H(4,4,1) = q128
      H(4,4,2) = q130
      H(4,4,3) = q131
      H(4,4,4) = q115 + q127 + q29
      AP(4,4) = q119*q40
      FP(4,1) = 0
      FP(4,2) = -q41*q87

END SUBROUTINE D2RHS


SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP) 
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NDIM, IJAC, ICP(*)
    DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
    DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
    DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM),DFDP(NDIM,*)
    ! Declare variables with double precision
    DOUBLE PRECISION :: y_N, y_E, x_B, x_N, m, eta, ep, k2
    DOUBLE PRECISION :: TaN, TaL, T0, S0, SB, sig_noSV
    DOUBLE PRECISION :: alphaT, alphaS, y_B, g, x_E, TildeVE, TildeVB
    DOUBLE PRECISION :: W, ps, delta_N, Psi, k2n, k2l, k1, kappaN, kappaE 
    DOUBLE PRECISION :: SV, sig, phi, VN, VL, VB
	DOUBLE PRECISION T
    DOUBLE PRECISION A(4,4)
    DOUBLE PRECISION :: h1, h2, h3
    DOUBLE PRECISION :: lambda1, lambda2, lambda3
    ! Second derivatives: H(i,j,k) = d2 f_i/du_j du_k, FP = df/d(mu,eta), AP = dA/deta
    DOUBLE PRECISION H(4,4,4), FP(4,2), AP(4,4)
    INTEGER :: i, k

    ! Eig values and eigfunctions norms + PO period
    h1      = PAR(20) 
    lambda1 = PAR(21)
    h2      = PAR(22)
    lambda2 = PAR(23)
    h3      = PAR(24)
    lambda3 = PAR(25)
    T = PAR(11)

    IF (IJAC > 0) DFDU(1:NDIM,1:NDIM) = 0
    IF (IJAC == 2) DFDP(1:NDIM,1:48) = 0

    ! PO 
    CALL RHS(U,PAR,F,.TRUE.,A)
    IF (IJAC > 0) CALL D2RHS(U,PAR,H,FP,AP)
	IF(NDIM==4)THEN
        ! IPS=2: AUTO scales by the period itself
        IF (IJAC > 0) DFDU(1:4,1:4) = A
        IF (IJAC == 2) DFDP(1:4,1:2) = FP
        RETURN
    ENDIF
    IF (IJAC > 0) DFDU(1:4,1:4) = T*A
    IF (IJAC == 2) THEN
        DFDP(1:4,11) = F(1:4)
        DFDP(1:4,1:2) = T*FP
    ENDIF

    ! Store (stable) floq 
    F(5:8) = MATMUL((A(:,:)),U(5:8))
    F(1:8) = T*F(1:8)
    IF (IJAC > 0) THEN
        DO i = 1, 4
            DO k = 1, 4
                DFDU(4+i,k) = T*SUM(H(i,:,k)*U(5:8))
            ENDDO
        ENDDO
        DFDU(5:8,5:8) = T*A
    ENDIF
    IF (IJAC == 2) THEN
        DFDP(5:8,11) = MATMUL(A,U(5:8))
        DFDP(5:8,2) = T*MATMUL(AP,U(5:8))
    ENDIF
    IF (NDIM==8) RETURN

    ! Store (stable) floq 
    F(9:12) = MATMUL(A(:,:),U(9:12))
    F(9:12) = T*F(9:12)
    IF (IJAC > 0) THEN
        DO i = 1, 4
            DO k = 1, 4
                DFDU(8+i,k) = T*SUM(H(i,:,k)*U(9:12))
            ENDDO
        ENDDO
        DFDU(9:12,9:12) = T*A
    ENDIF
    IF (IJAC == 2) THEN
        DFDP(9:12,11) = MATMUL(A,U(9:12))
        DFDP(9:12,2) = T*MATMUL(AP,U(9:12))
    ENDIF
    IF (NDIM==12) RETURN

    ! Store (unstable) floq 
    ! F(13:16) = MATMUL(A(:,:),U(13:16))
    F(13:16) = T*MATMUL(A,U(13:16)) + LOG(ABS(lambda3))*U(13:16)
    ! F(13:16) = T*F(13:16)
    IF (IJAC > 0) THEN
        DO i = 1, 4
            DO k = 1, 4
                DFDU(12+i,k) = T*SUM(H(i,:,k)*U(13:16))
            ENDDO
            DFDU(12+i,12+i) = LOG(ABS(lambda3))
        ENDDO
        DFDU(13:16,13:16) = DFDU(13:16,13:16) + T*A
    ENDIF
    IF (IJAC == 2) THEN
        DFDP(13:16,11) = MATMUL(A,U(13:16))
        DFDP(13:16,2) = T*MATMUL(AP,U(13:16))
        DFDP(13:16,25) = U(13:16)/lambda3
    ENDIF
    IF (NDIM==16) RETURN

    ! Segment Wu
    CALL RHS(U(17:20),PAR,F(17:20),.TRUE.,A)
    T = PAR(12)
    IF (IJAC > 0) DFDU(17:20,17:20) = T*A
    IF (IJAC == 2) THEN
        CALL D2RHS(U(17:20),PAR,H,FP,AP)
        DFDP(17:20,12) = F(17:20)
        DFDP(17:20,1:2) = T*FP
    ENDIF
    F(17:20) = F(17:20) * T  
    IF (NDIM==20) RETURN

    ! Segment Ws 
    CALL RHS(U(21:24),PAR,F(21:24),.TRUE.,A)
    T = PAR(13)
    IF (IJAC > 0) DFDU(21:24,21:24) = T*A
    IF (IJAC == 2) THEN
        CALL D2RHS(U(21:24),PAR,H,FP,AP)
        DFDP(21:24,13) = F(21:24)
        DFDP(21:24,1:2) = T*FP
    ENDIF
    F(21:24) = F(21:24) * T 
    IF (NDIM==24) RETURN
END SUBROUTINE FUNC


SUBROUTINE STPNT(NDIM,U,PAR,T)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NDIM
    DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
    DOUBLE PRECISION, INTENT(IN) :: T
    DOUBLE PRECISION, SAVE :: s(16)
    DOUBLE PRECISION :: h1, h2, h3
    DOUBLE PRECISION :: lambda1, lambda2, lambda3
    DOUBLE PRECISION :: eps1, eps2, epu1, gap, SEG
    h1      = PAR(20) 
    lambda1 = PAR(21)
    h2      = PAR(22)
    lambda2 = PAR(23)
    h3      = PAR(24)
    lambda3 = PAR(25)
    eps1    = PAR(34)
    eps2    = PAR(35)
    epu1    = PAR(36)
    gap     = PAR(29)

    ! Epsilon values
    PAR(34) = 1e-6 ! eps1
    PAR(35) = 1e-6 ! eps2
    PAR(36) = 1e-2 ! epu1

    ! Inital gueses for stable floq 
    PAR(23) = -0.1
    PAR(24) = -0.1
    
    ! System parameters
    PAR(1) = -0.002006D0
    PAR(2) = -4.1D0
    PAR(3) = 8.0 / (2.1D4)
    PAR(4) = 0.02D0

    PAR(12) = 0
    PAR(13) = 0

    ! Section definition. Here x-component = PAR(40)
    PAR(40) = 1.35

    ! Base equilibrium for intial continuation 
    U(1) = -3.5435D0
    U(2) = 0.9597D0
    U(3) = 3.8512D0
    U(4) = 1.1678D0

    ! Set up eignspaces for segment continuation
    IF(NDIM==20)THEN ! unstable eignspace (1D)
        IF(PAR(12)==0) THEN
            s(1:4) = U(1:4)
            s(9:12) = U(9:12)
            s(1:16) = U(1:16)
       ENDIF
       U(17:20) = s(1:4) + epu1*s(13:16)
       RETURN
    ELSEIF(NDIM==24)THEN ! stable eignsapce (2D)
       IF(PAR(13)==0)THEN
            s(5:8) = U(5:8)
            s(9:12) = U(9:12)
       ENDIF
       U(21:24) = s(1:4) + eps1*s(5:8) + eps2*s(9:12)
       RETURN
    ENDIF
END SUBROUTINE STPNT


    ! IF (FIRST) THEN ! initialization for BCND
    !    FIRST = .FALSE.
    !    IF (NDIM==20) THEN ! After WU continuation 
    !       PAR(27) = GETP("BV0",20,U) - SIG ! sigma+
    !    ELSEIF (NDIM == 24) THEN  ! After WS continuation 
    !       NBC = AINT(GETP("NBC",0,U)) 
    !       IF (NBC == 29) THEN ! Right as we do the Ws cont 
    !         PAR(28) = GETP("BV1",24,U) - SIG ! sigma-
    !       ELSEIF (NBC == 33) THEN ! When we start LIN 
    !         IF (DOT_PRODUCT(PAR(30:33),PAR(30:33)) > 0) RETURN
    !         DO i=1,4
    !             a(i) = GETP("BV0",17+i,U)
    !             b(i) = GETP("BV1",20+i,U)
    !             d(i) = GETP("BV0",16+i,U) - GETP("BV1",20+i,U)
    !         ENDDO
    !         normlv = sqrt(DOT_PRODUCT(d,d))
    !         PAR(29) = normlv  ! gap size in PAR(29)
    !         PAR(30:33) = d(1:4)/normlv ! Lin vector in PAR(30)-PAR(33)



SUBROUTINE BCND(NDIM,PAR,ICP,NBC,U0,U1,FB,IJAC,DBC)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NDIM, ICP(*), NBC, IJAC
    DOUBLE PRECISION, INTENT(IN) :: PAR(*), U0(NDIM), U1(NDIM)
    DOUBLE PRECISION, INTENT(OUT) :: FB(NBC)
    DOUBLE PRECISION, INTENT(INOUT) :: DBC(NBC,*)

    DOUBLE PRECISION :: h1, h2, h3
    DOUBLE PRECISION :: lambda1, lambda2, lambda3
    DOUBLE PRECISION :: eps1, eps2, epu1, gap, SEG
    INTEGER :: i, N, P
    h1      = PAR(20) 
    lambda1 = PAR(21)
    h2      = PAR(22)
    lambda2 = PAR(23)
    h3      = PAR(24)
    lambda3 = PAR(25)
    eps1    = PAR(34)
    eps2    = PAR(35)
    epu1    = PAR(36)
    gap     = PAR(29)

    !!! Section on 2nd coordinate 
    SEG = PAR(40)

    ! Jacobian columns: U0(i) -> i, U1(i) -> N+i, PAR(j) -> P+j
    N = NDIM
    P = 2*NDIM
    IF (IJAC > 0) DBC(1:NBC,1:2*NDIM) = 0
    IF (IJAC == 2) DBC(1:NBC,P+1:P+48) = 0

    !!! BOUNDARY CONDITIONS FOR PO
    FB(1:4) = U0(1:4) - U1(1:4)  ! State
    IF (IJAC > 0) THEN
        DO i = 1, 4
            DBC(i,i) = 1
            DBC(i,N+i) = -1
        ENDDO
    ENDIF

    !!! BOUNDARY CONDITIONS FOR FLOQ 
    !! FLOQ 1 
    FB(5:8) = U1(5:8) - lambda1*U0(5:8)  
    FB(9) = DOT_PRODUCT(U0(5:8),U0(5:8)) - h1
    IF (IJAC > 0) THEN
        DO i = 1, 4
            DBC(4+i,N+4+i) = 1
            DBC(4+i,4+i) = -lambda1
        ENDDO
        DBC(9,5:8) = 2*U0(5:8)
    ENDIF
    IF (IJAC == 2) THEN
        DBC(5:8,P+21) = -U0(5:8)
        DBC(9,P+20) = -1
    ENDIF
    IF (NBC==9) RETURN 

    !! FLOQ 2 
    FB(10:13) = U1(9:12) - lambda2*U0(9:12)  
    FB(14) = h2 - DOT_PRODUCT(U0(9:12),U0(9:12))  
    IF (IJAC > 0) THEN
        DO i = 1, 4
            DBC(9+i,N+8+i) = 1
            DBC(9+i,8+i) = -lambda2
        ENDDO
        DBC(14,9:12) = -2*U0(9:12)
    ENDIF
    IF (IJAC == 2) THEN
        DBC(10:13,P+23) = -U0(9:12)
        DBC(14,P+22) = 1
    ENDIF
    IF (NBC==14) RETURN 

    !! FLOQ 3 
    FB(15:18) =  U1(13:16) - SIGN(1.0d0,lambda1)*U0(13:16) !U1(13:16) - lambda3*U0(13:16)  
    FB(19) = h3 - DOT_PRODUCT(U0(13:16), U0(13:16))  
    IF (IJAC > 0) THEN
        DO i = 1, 4
            DBC(14+i,N+12+i) = 1
            DBC(14+i,12+i) = -SIGN(1.0d0,lambda1)
        ENDDO
        DBC(19,13:16) = -2*U0(13:16)
    ENDIF
    IF (IJAC == 2) DBC(19,P+24) = 1
    IF (NBC==19) RETURN 

    !! WU 
    FB(20:23) = U0(17:20) - (U0(1:4) + epu1*U0(13:16))  ! espace
    FB(24)    = U1(17) - SEG - PAR(27) ! on section 
    IF (IJAC > 0) THEN
        DO i = 1, 4
            DBC(19+i,16+i) = 1
            DBC(19+i,i) = -1
            DBC(19+i,12+i) = -epu1
        ENDDO
        DBC(24,N+17) = 1
    ENDIF
    IF (IJAC == 2) THEN
        DBC(20:23,P+36) = -U0(13:16)
        DBC(24,P+40) = -1
        DBC(24,P+27) = -1
    ENDIF
    IF (NBC==24) RETURN

    !! WS 
    FB(25:28) = U1(21:24) - (U0(1:4) + eps1*U0(5:8) + eps2*U0(9:12))  ! espace
    FB(29) = U0(21) - SEG - PAR(28) ! on section 
    IF (IJAC > 0) THEN
        DO i = 1, 4
            DBC(24+i,N+20+i) = 1
            DBC(24+i,i) = -1
            DBC(24+i,4+i) = -eps1
            DBC(24+i,8+i) = -eps2
        ENDDO
        DBC(29,21) = 1
    ENDIF
    IF (IJAC == 2) THEN
        DBC(25:28,P+34) = -U0(5:8)
        DBC(25:28,P+35) = -U0(9:12)
        DBC(29,P+40) = -1
        DBC(29,P+28) = -1
    ENDIF
    IF (NBC==29) RETURN

    !! LIN
    FB(30:33) = U1(17:20) - U0(21:24) - gap*PAR(30:33)
    IF (IJAC > 0) THEN
        DO i = 1, 4
            DBC(29+i,N+16+i) = 1
            DBC(29+i,20+i) = -1
        ENDDO
    ENDIF
    IF (IJAC == 2) THEN
        DBC(30:33,P+29) = -PAR(30:33)
        DO i = 1, 4
            DBC(29+i,P+29+i) = -gap
        ENDDO
    ENDIF
END SUBROUTINE BCND


SUBROUTINE ICND(NDIM,PAR,ICP,NINT,U,UOLD,UDOT,UPOLD,FI,IJAC,DINT)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NDIM, ICP(*), NINT, IJAC
    DOUBLE PRECISION, INTENT(IN) :: PAR(*)
    DOUBLE PRECISION, INTENT(IN) :: U(NDIM), UOLD(NDIM), UDOT(NDIM), UPOLD(NDIM)
    DOUBLE PRECISION, INTENT(OUT) :: FI(NINT)
    DOUBLE PRECISION, INTENT(INOUT) :: DINT(NINT,*)

    ! Integral phase condition
    FI(1) = DOT_PRODUCT(U(1:4),UPOLD(1:4))
    IF (IJAC > 0) THEN
        DINT(1,1:NDIM) = 0
        DINT(1,1:4) = UPOLD(1:4)
    ENDIF
    IF (IJAC == 2) DINT(1,NDIM+1:NDIM+48) = 0
    IF (NINT==1) RETURN

!    FI(2) = DOT_PRODUCT(UPOLD(21:24),U(21:24)-UOLD(21:24))
!    FI(3) = DOT_PRODUCT(UPOLD(17:20),U(17:20)-UOLD(17:20))
!    IF (NINT==3) RETURN
END SUBROUTINE ICND
  

SUBROUTINE FOPT
END SUBROUTINE FOPT


SUBROUTINE PVLS(NDIM,U,PAR)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NDIM
    DOUBLE PRECISION, INTENT(IN) :: U(NDIM)
    DOUBLE PRECISION, INTENT(INOUT) :: PAR(*)
    
    DOUBLE PRECISION, EXTERNAL :: GETP
    DOUBLE PRECISION d(4), normlv, a(4), b(4)
    DOUBLE PRECISION SIG
    INTEGER i, NBC
    LOGICAL, SAVE :: FIRST = .TRUE.

!   PAR(27) : sigma+ : U0(22)-SIG (x-distance W^u(P) from section x=SIG)
!   PAR(28) : sigma- : U1(18)-SIG (x-distance W^s(E) from section x=SIG)
!   PAR(29) : eta    : gap size for Lin vector
!   PAR(30:33) : lin vector 
!   PAR(34) = eps1 
!   PAR(35) = eps2
!   PAR(36) = epu1

    SIG = PAR(40)

    IF (FIRST) THEN ! initialization for BCND
       FIRST = .FALSE.
       IF (NDIM==20) THEN ! After WU continuation 
          PAR(27) = GETP("BV0",17,U) - SIG ! sigma+
       ELSEIF (NDIM == 24) THEN  ! After WS continuation 
          NBC = AINT(GETP("NBC",0,U)) 
          IF (NBC == 29) THEN ! Right as we do the Ws cont 
            PAR(28) = GETP("BV1",21,U) - SIG ! sigma-
          ELSEIF (NBC == 33) THEN ! When we start LIN 
            IF (DOT_PRODUCT(PAR(30:33),PAR(30:33)) > 0) RETURN
            DO i=1,4
                a(i) = GETP("BV0",16+i,U)
                b(i) = GETP("BV1",20+i,U)
                d(i) = GETP("BV1",16+i,U) - GETP("BV0",20+i,U)
            ENDDO
            normlv = sqrt(DOT_PRODUCT(d,d))
            PAR(29) = normlv  ! gap size in PAR(29)
            PAR(30:33) = d(1:4)/normlv ! Lin vector in PAR(30)-PAR(33)
          ENDIF
       ENDIF
       RETURN
    ENDIF

    PAR(45:48) =  a(1:4) - b(1:4) - PAR(29)*PAR(30:33)
END SUBROUTINE PVLS

    ! FB(30:33) = U0(17:20) - U1(21:24) - gap*PAR(30:33)

!  U0(17:20) - U1(21:24)