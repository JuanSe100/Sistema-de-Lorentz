PROGRAM RK5
    IMPLICIT NONE
    
!!======Ingresamos la cuadruple precision====================
    INTEGER, PARAMETER:: ip=selected_int_kind(18)
    INTEGER, PARAMETER:: rp=selected_real_kind(33,4931)
   

    INTEGER(ip) :: nvar, nok, nbad
    REAL(rp) :: eps, h1, hmin
    INTEGER(ip):: i, m
    REAL(rp) :: t, tnext, dt
    REAL(rp):: ystart(1)

    nvar = 1_ip
    eps = 1.e-6_rp
    h1 = 1.e-6_rp
    hmin = 1.e-15_rp

    
    t = 0.0_rp
    dt = 0.01_rp
    m = 1000_ip
    ystart(1) = 1.0_rp

    OPEN(10, FILE='datos.dat')
    WRITE(10, *) t, ystart(1)
    
    DO i = 1, m
        tnext = t + dt
        CALL odeint(ystart, nvar, t, tnext, eps, h1, hmin, nok, nbad, derivs, rkqs)
        t = tnext
  !     WRITE(*, *) tnext, ystart(1)
        WRITE(10, *) tnext, ystart(1)
    END DO

    CLOSE(10)
    PRINT *, 'Terminado.....'

CONTAINS !!=============================================================================

    SUBROUTINE rkqs(y, dydx, n, x, htry, eps, yscal, hdid, hnext, derivs)
        IMPLICIT NONE
        INTEGER, PARAMETER:: ip=selected_int_kind(18)
        INTEGER, PARAMETER:: rp=selected_real_kind(33,4931)
        
        INTEGER(ip), INTENT(IN) :: n
        INTEGER(ip), PARAMETER :: NMAX = 50_ip
        REAL(rp), INTENT(IN) :: eps, htry, yscal(n), dydx(n)
        REAL(rp), INTENT(INOUT) :: x, y(n)
        REAL(rp), INTENT(OUT) :: hdid, hnext
        EXTERNAL derivs

        ! Paso Runge-Kutta 5to orden con monitoreo de error de truncamiento
        ! para asegurar la precision y ajustar el tamano de paso.
        INTEGER(ip) :: i
        REAL(rp) :: errmax, h, htemp, xnew
        REAL(rp) :: yerr(NMAX), ytemp(NMAX)
        REAL(rp), PARAMETER :: SAFETY = 0.9_rp, PGROW = -0.2_rp, PSHRNK = -0.25_rp, ERRCON = 1.89e-4_rp
        ! El valor ERRCON es igual a (5/SAFETY)**(1/PGROW)

        h = htry ! Establece el paso al valor de prueba inicial

        DO
            CALL rkck(y, dydx, n, x, h, ytemp, yerr, derivs) ! Da un paso
            
            errmax = 0.0_rp ! Evalua la precision
            DO i = 1, n
                errmax = MAX(errmax, ABS(yerr(i) / yscal(i)))
            END DO
            
            errmax = errmax / eps ! Escala relativo a la tolerancia
            
            IF (errmax > 1.0_rp) THEN ! Error grande, reduce el paso
                htemp = SAFETY * h * (errmax ** PSHRNK)
                h = SIGN(MAX(ABS(htemp), 0.1_rp * ABS(h)), h) ! No mas de factor de 10
                xnew = x + h
                IF (xnew == x) STOP 'underflow en h'
                CYCLE ! Intenta de nuevo
            ELSE ! Paso exitoso. Calcula el proximo paso
                IF (errmax > ERRCON) THEN
                    hnext = SAFETY * h * (errmax ** PGROW)
                ELSE ! Incremento no mayor a un factor de 5
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
        INTEGER, PARAMETER:: ip=selected_int_kind(18)
        INTEGER, PARAMETER:: rp=selected_real_kind(33,4931)
        
        INTEGER(ip), INTENT(IN) :: n
        INTEGER(ip), PARAMETER :: NMAX = 50_ip
        REAL(rp), INTENT(IN) :: h, x, y(n), dydx(n)
        REAL(rp), INTENT(OUT) :: yout(n), yerr(n)
        EXTERNAL derivs

        ! Metodo Runge-Kutta de Cash-Karp de quinto orden.
        ! Avanza la solucion un intervalo h y devuelve variables en yout.
        ! Devuelve el error en yerr usando el metodo de 4to orden.
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

        DO i = 1, n ! Primer paso
            ytemp(i) = y(i) + B21 * h * dydx(i)
        END DO
        
        CALL derivs(x + A2 * h, ytemp, ak2) ! Segundo paso
        
        DO i = 1, n
            ytemp(i) = y(i) + h * (B31 * dydx(i) + B32 * ak2(i))
        END DO
        
        CALL derivs(x + A3 * h, ytemp, ak3) ! Tercer paso
        
        DO i = 1, n
            ytemp(i) = y(i) + h * (B41 * dydx(i) + B42 * ak2(i) + B43 * ak3(i))
        END DO
        
        CALL derivs(x + A4 * h, ytemp, ak4) ! Cuarto paso
        DO i = 1, n
            ytemp(i) = y(i) + h * (B51 * dydx(i) + B52 * ak2(i) + B53 * ak3(i) + &
                                   B54 * ak4(i))
        END DO
        
        CALL derivs(x + A5 * h, ytemp, ak5) ! Quinto paso
        
        DO i = 1, n
            ytemp(i) = y(i) + h * (B61 * dydx(i) + B62 * ak2(i) + B63 * ak3(i) + &
                                   B64 * ak4(i) + B65 * ak5(i))
        END DO
        
        CALL derivs(x + A6 * h, ytemp, ak6) ! Sexto paso
        DO i = 1, n ! Acumula incrementos con los pesos adecuados
            yout(i) = y(i) + h * (C1 * dydx(i) + C3 * ak3(i) + C4 * ak4(i) + &
                                  C6 * ak6(i))
        END DO
        
        DO i = 1, n
            ! Estima error como diferencia entre 4to y 5to orden
            yerr(i) = h * (DC1 * dydx(i) + DC3 * ak3(i) + DC4 * ak4(i) + DC5 * ak5(i) + &
                           DC6 * ak6(i))
        END DO

    END SUBROUTINE rkck

    SUBROUTINE odeint(ystart, nvar, x1, x2, eps, h1, hmin, nok, nbad, derivs, rkqs)
        IMPLICIT NONE
        INTEGER, PARAMETER:: ip=selected_int_kind(18)
        INTEGER, PARAMETER:: rp=selected_real_kind(33,4931)
        
        INTEGER(ip), INTENT(IN) :: nvar
        INTEGER(ip), INTENT(OUT) :: nok, nbad
        REAL(rp), INTENT(IN) :: eps, h1, hmin, x1, x2
        REAL(rp), INTENT(INOUT) :: ystart(nvar)
        EXTERNAL derivs, rkqs

        INTEGER(ip), PARAMETER :: MAXSTP = 10000_ip, NMAX = 50_ip
        REAL(rp), PARAMETER :: TINY = 1.e-30_rp

        ! Driver Runge-Kutta adaptativo. Integra ystart desde x1 a x2 con precision eps.
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
            
            DO i = 1, nvar ! Escalado para monitorear la precision
                yscal(i) = ABS(y(i)) + ABS(h * dydx(i)) + TINY
            END DO

            ! Si el paso sobrepasa el limite x2, lo ajusta
            IF ((x + h - x2) * (x + h - x1) > 0.0) h = x2 - x
            
            CALL rkqs(y, dydx, nvar, x, h, eps, yscal, hdid, hnext, derivs)

         !! IF (hdid == h) THEN
         !!    nok = nok + 1_ip
         !!  ELSE
         !!    nbad = nbad + 1_ip
         !!   END IF

            IF ((x - x2) * (x2 - x1) >= 0.0_rp) THEN ! Terminamos?
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
        INTEGER, PARAMETER:: ip=selected_int_kind(18)
        INTEGER, PARAMETER:: rp=selected_real_kind(33,4931)
        
        REAL(rp), INTENT(IN) :: x, y(1)
        REAL(rp), INTENT(OUT) :: dydt(1)
        REAL(rp) :: a
        a = 3.3_rp
        dydt(1) = a * y(1)
    END SUBROUTINE derivs

END PROGRAM RK5  !! Fin del programa
