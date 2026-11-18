PROGRAM RK5
    IMPLICIT NONE
    
!!======Ingresamos la cuadruple precision====================
    INTEGER, PARAMETER :: ip = selected_int_kind(18)
    INTEGER, PARAMETER :: rp = selected_real_kind(33, 4931)

    INTEGER(ip) :: nvar, nok, nbad
    REAL(rp) :: eps, h1, hmin
    INTEGER(ip) :: i, m
    REAL(rp) :: t, tnext, dt
    REAL(rp) :: ystart(3), y2start(3)
    REAL(rp) :: eps1, eps2, eps3, sep_z
    REAL(rp) :: sigma, rho, beta 
    
    nvar = 3_ip
    eps = 1.e-30_rp       ! Tolerancia acorde a cuadruple precision (perturbacion ~1e-31)
    h1 = 1.e-7_rp         ! Paso inicial de prueba
    hmin = 1.e-32_rp      ! Minimo valor permitido para el paso
    
    t = 0.0_rp
    dt = 0.001_rp
    m = 1000000_ip        ! 1,000,000 pasos -> t en [0, 1000]
    
    ! Condiciones iniciales Trayectoria 1 (Inciso a)
    ystart(1) = 200.0_rp
    ystart(2) = -10.0_rp
    ystart(3) = 3.0_rp
    
    ! Perturbaciones iniciales (Inciso b)
    eps1 = 1.0e-30_rp
    eps2 = 1.0e-31_rp
    eps3 = -1.0e-29_rp
    
    ! Condiciones iniciales Trayectoria 2: (xo + eps1, yo + eps2, zo - eps3)
    y2start(1) = 200.0_rp + eps1
    y2start(2) = -10.0_rp + eps2
    y2start(3) = 3.0_rp - eps3
    
    ! Constantes en el Atractor de Lorenz
    sigma = 10.0_rp
    rho = 29.5_rp
    beta = 8.0_rp / 3.0_rp
   
    CALL condicones_utilizadas()
    
    ! Archivo de salida con las columnas:
    ! Col 1: t | Col 2: y(t) original | Col 3: y'(t) perturbada | Col 4: |z'(t) - z(t)|
    OPEN(10, FILE='resultados.dat')
    
    sep_z = ABS(y2start(3) - ystart(3))
    WRITE(10, *) t, ystart(2), y2start(2), sep_z
    
    DO i = 1, m
        tnext = REAL(i, rp) * dt   ! Evita acumulacion de error de redondeo en t
        
        CALL odeint(ystart, nvar, t, tnext, eps, h1, hmin, nok, nbad, derivs, rkqs)
        CALL odeint(y2start, nvar, t, tnext, eps, h1, hmin, nok, nbad, derivs, rkqs)
        
        t = tnext
        sep_z = ABS(y2start(3) - ystart(3))
        
        WRITE(10, *) t, ystart(2), y2start(2), sep_z
    END DO
    
    CALL puntos_fijos(sigma, rho, beta)
   
    CLOSE(10)
    PRINT *, 'Terminado.....'

CONTAINS !!=============================================================================

    SUBROUTINE rkqs(y, dydx, n, x, htry, eps, yscal, hdid, hnext, derivs)
        IMPLICIT NONE
        INTEGER(ip), INTENT(IN) :: n
        INTEGER(ip), PARAMETER :: NMAX = 50_ip
        REAL(rp), INTENT(IN) :: eps, htry, yscal(n), dydx(n)
        REAL(rp), INTENT(INOUT) :: x, y(n)
        REAL(rp), INTENT(OUT) :: hdid, hnext
        EXTERNAL derivs

        INTEGER(ip) :: i
        REAL(rp) :: errmax, h, htemp, xnew
        REAL(rp) :: yerr(NMAX), ytemp(NMAX)
        REAL(rp), PARAMETER :: SAFETY = 0.9_rp, PGROW = -0.2_rp, PSHRNK = -0.25_rp, ERRCON = 1.89e-4_rp

        h = htry

        DO
            CALL rkck(y, dydx, n, x, h, ytemp, yerr, derivs)
            
            errmax = 0.0_rp
            DO i = 1, n
                errmax = MAX(errmax, ABS(yerr(i) / yscal(i)))
            END DO
            
            errmax = errmax / eps
            
            IF (errmax > 1.0_rp) THEN
                htemp = SAFETY * h * (errmax ** PSHRNK)
                h = SIGN(MAX(ABS(htemp), 0.1_rp * ABS(h)), h)
                xnew = x + h
                IF (xnew == x) STOP 'underflow en h'
                CYCLE
            ELSE
                IF (errmax > ERRCON) THEN
                    hnext = SAFETY * h * (errmax ** PGROW)
                ELSE
                    hnext = 5.0_rp * h
                END IF
                hdid = h
                x = x + h
                DO i = 1, n
                    y(i) = ytemp(i)
                END DO
                EXIT
            END IF
        END DO

    END SUBROUTINE rkqs

    SUBROUTINE rkck(y, dydx, n, x, h, yout, yerr, derivs)
        IMPLICIT NONE
        INTEGER(ip), INTENT(IN) :: n
        INTEGER(ip), PARAMETER :: NMAX = 50_ip
        REAL(rp), INTENT(IN) :: h, x, y(n), dydx(n)
        REAL(rp), INTENT(OUT) :: yout(n), yerr(n)
        EXTERNAL derivs

        INTEGER(ip) :: i
        REAL(rp) :: ak2(NMAX), ak3(NMAX), ak4(NMAX), ak5(NMAX), ak6(NMAX), ytemp(NMAX)
        REAL(rp), PARAMETER :: A2 = 0.2_rp, A3 = 0.3_rp, A4 = 0.6_rp, A5 = 1.0_rp, A6 = 0.875_rp
        REAL(rp), PARAMETER :: B21 = 0.2_rp, B31 = 3.0_rp/40.0_rp, B32 = 9.0_rp/40.0_rp
        REAL(rp), PARAMETER :: B41 = 0.3_rp, B42 = -0.9_rp, B43 = 1.2_rp
        REAL(rp), PARAMETER :: B51 = -11.0_rp/54.0_rp, B52 = 2.5_rp, B53 = -70.0_rp/27.0_rp, B54 = 35.0_rp/27.0_rp
        REAL(rp), PARAMETER :: B61 = 1631.0_rp/55296.0_rp, B62 = 175.0_rp/512.0_rp, B63 = 575.0_rp/13824.0_rp
        REAL(rp), PARAMETER :: B64 = 44275.0_rp/110592.0_rp, B65 = 253.0_rp/4096.0_rp
        REAL(rp), PARAMETER :: C1 = 37.0_rp/378.0_rp, C3 = 250.0_rp/621.0_rp, C4 = 125.0_rp/594.0_rp, C6 = 512.0_rp/1771.0_rp
        REAL(rp), PARAMETER :: DC1 = C1 - 2825.0_rp/27648.0_rp, DC3 = C3 - 18575.0_rp/48384.0_rp
        REAL(rp), PARAMETER :: DC4 = C4 - 13525.0_rp/55296.0_rp, DC5 = -277.0_rp/14336.0_rp, DC6 = C6 - 0.25_rp

        DO i = 1, n
            ytemp(i) = y(i) + B21 * h * dydx(i)
        END DO
        CALL derivs(x + A2 * h, ytemp, ak2)
        
        DO i = 1, n
            ytemp(i) = y(i) + h * (B31 * dydx(i) + B32 * ak2(i))
        END DO
        CALL derivs(x + A3 * h, ytemp, ak3)
        
        DO i = 1, n
            ytemp(i) = y(i) + h * (B41 * dydx(i) + B42 * ak2(i) + B43 * ak3(i))
        END DO
        CALL derivs(x + A4 * h, ytemp, ak4)
        
        DO i = 1, n
            ytemp(i) = y(i) + h * (B51 * dydx(i) + B52 * ak2(i) + B53 * ak3(i) + B54 * ak4(i))
        END DO
        CALL derivs(x + A5 * h, ytemp, ak5)
        
        DO i = 1, n
            ytemp(i) = y(i) + h * (B61 * dydx(i) + B62 * ak2(i) + B63 * ak3(i) + B64 * ak4(i) + B65 * ak5(i))
        END DO
        CALL derivs(x + A6 * h, ytemp, ak6)
        
        DO i = 1, n
            yout(i) = y(i) + h * (C1 * dydx(i) + C3 * ak3(i) + C4 * ak4(i) + C6 * ak6(i))
        END DO
        
        DO i = 1, n
            yerr(i) = h * (DC1 * dydx(i) + DC3 * ak3(i) + DC4 * ak4(i) + DC5 * ak5(i) + DC6 * ak6(i))
        END DO

    END SUBROUTINE rkck

    SUBROUTINE odeint(ystart, nvar, x1, x2, eps, h1, hmin, nok, nbad, derivs, rkqs)
        IMPLICIT NONE
        INTEGER(ip), INTENT(IN) :: nvar
        INTEGER(ip), INTENT(OUT) :: nok, nbad
        REAL(rp), INTENT(IN) :: eps, h1, hmin, x1, x2
        REAL(rp), INTENT(INOUT) :: ystart(nvar)
        EXTERNAL derivs, rkqs

        INTEGER(ip), PARAMETER :: MAXSTP = 1000000_ip, NMAX = 50_ip
        REAL(rp), PARAMETER :: TINY = 1.e-35_rp
        
        INTEGER(ip) :: i, nstp
        REAL(rp) :: h, hdid, hnext, x, dydx(NMAX), y(NMAX), yscal(NMAX)

        x = x1
        h = SIGN(h1, x2 - x1)
        nok = 0_ip
        nbad = 0_ip
        
        DO i = 1, nvar
            y(i) = ystart(i)
        END DO

        DO nstp = 1, MAXSTP
            CALL derivs(x, y, dydx)
            
            DO i = 1, nvar
                yscal(i) = ABS(y(i)) + ABS(h * dydx(i)) + TINY
            END DO

            IF ((x + h - x2) * (x + h - x1) > 0.0_rp) h = x2 - x
            
            CALL rkqs(y, dydx, nvar, x, h, eps, yscal, hdid, hnext, derivs)

            IF ((x - x2) * (x2 - x1) >= 0.0_rp) THEN
                DO i = 1, nvar
                    ystart(i) = y(i)
                END DO
                RETURN
            END IF

            IF (ABS(hnext) < hmin) STOP 'paso menor que minimo'
            h = hnext
        END DO

        STOP 'demasiados pasos en odeint'
    END SUBROUTINE odeint

    SUBROUTINE derivs(x, y, dydt)
        IMPLICIT NONE
        REAL(rp), INTENT(IN) :: x, y(3)
        REAL(rp), INTENT(OUT) :: dydt(3)
        
        dydt(1) = sigma * (y(2) - y(1))
        dydt(2) = y(1) * (rho - y(3)) - y(2)
        dydt(3) = y(1) * y(2) - beta * y(3)
    END SUBROUTINE derivs
    
    SUBROUTINE puntos_fijos(sigma, rho, beta)
        IMPLICIT NONE
        REAL(rp), INTENT(IN) :: sigma, rho, beta
        REAL(rp) :: xx, yy, zz
         
        IF (rho > 1.0_rp) THEN
            xx = SQRT(beta * (rho - 1.0_rp))
            yy = SQRT(beta * (rho - 1.0_rp))
            zz = rho - 1.0_rp
            WRITE(*, *) 'Punto fijo P+: ', xx, yy, zz
            WRITE(*, *) 'Punto fijo P-: ', -xx, -yy, zz
        ELSE
            WRITE(*, *) 'Solo existe el punto fijo en el origen (0,0,0)'
        END IF
    END SUBROUTINE puntos_fijos

    SUBROUTINE condicones_utilizadas()
        IMPLICIT NONE
        REAL(rp) :: TTtime
        
        WRITE(*,*) "Sistema de Lorenz - Examen 1"
        WRITE(*,*) "Trayectoria 1 (x0, y0, z0):", ystart(1), ystart(2), ystart(3)
        WRITE(*,*) "Trayectoria 2 (x'0,y'0,z'0):", y2start(1), y2start(2), y2start(3)
        WRITE(*,*) "Parametros: sigma =", sigma, "rho =", rho, "beta =", beta
        TTtime = REAL(m, rp) * dt
        WRITE(*,*) "Tiempo simulado:", TTtime 
    END SUBROUTINE condicones_utilizadas

END PROGRAM RK5

