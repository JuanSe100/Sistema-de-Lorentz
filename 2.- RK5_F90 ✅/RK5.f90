PROGRAM RK5
    IMPLICIT NONE

    INTEGER :: nvar, nok, nbad
    REAL :: eps, h1, hmin
    INTEGER :: i, m
    REAL :: t, tnext, dt
    REAL :: ystart(1)

    nvar = 1
    eps = 1.e-6
    h1 = 1.e-6
    hmin = 1.e-15

    t = 0.0
    dt = 0.01
    m = 10

    ystart(1) = 1.0

    OPEN(10, FILE='datos.dat')
    WRITE(10, *) t, ystart(1)

    DO i = 1, m
        tnext = t + dt
        CALL odeint(ystart, nvar, t, tnext, eps, h1, hmin, nok, nbad, derivs, rkqs)
        t = tnext
        WRITE(10, *) t, ystart(1)
    END DO

    CLOSE(10)
    PRINT *, 'Terminado. Revisa el archivo datos.dat'

CONTAINS !!=============================================================================

    SUBROUTINE rkqs(y, dydx, n, x, htry, eps, yscal, hdid, hnext, derivs)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: n
        INTEGER, PARAMETER :: NMAX = 50
        REAL, INTENT(IN) :: eps, htry, yscal(n), dydx(n)
        REAL, INTENT(INOUT) :: x, y(n)
        REAL, INTENT(OUT) :: hdid, hnext
        EXTERNAL derivs

        ! Paso Runge-Kutta 5to orden con monitoreo de error de truncamiento
        ! para asegurar la precision y ajustar el tamano de paso.
        INTEGER :: i
        REAL :: errmax, h, htemp, xnew
        REAL :: yerr(NMAX), ytemp(NMAX)
        REAL, PARAMETER :: SAFETY = 0.9, PGROW = -0.2, PSHRNK = -0.25, ERRCON = 1.89e-4
        ! El valor ERRCON es igual a (5/SAFETY)**(1/PGROW)

        h = htry ! Establece el paso al valor de prueba inicial

        DO
            CALL rkck(y, dydx, n, x, h, ytemp, yerr, derivs) ! Da un paso
            
            errmax = 0.0 ! Evalua la precision
            DO i = 1, n
                errmax = MAX(errmax, ABS(yerr(i) / yscal(i)))
            END DO
            
            errmax = errmax / eps ! Escala relativo a la tolerancia
            
            IF (errmax > 1.0) THEN ! Error grande, reduce el paso
                htemp = SAFETY * h * (errmax ** PSHRNK)
                h = SIGN(MAX(ABS(htemp), 0.1 * ABS(h)), h) ! No mas de factor de 10
                xnew = x + h
                IF (xnew == x) STOP 'underflow en h'
                CYCLE ! Intenta de nuevo
            ELSE ! Paso exitoso. Calcula el proximo paso
                IF (errmax > ERRCON) THEN
                    hnext = SAFETY * h * (errmax ** PGROW)
                ELSE ! Incremento no mayor a un factor de 5
                    hnext = 5.0 * h
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
        INTEGER, INTENT(IN) :: n
        INTEGER, PARAMETER :: NMAX = 50
        REAL, INTENT(IN) :: h, x, y(n), dydx(n)
        REAL, INTENT(OUT) :: yout(n), yerr(n)
        EXTERNAL derivs

        ! Metodo Runge-Kutta de Cash-Karp de quinto orden.
        ! Avanza la solucion un intervalo h y devuelve variables en yout.
        ! Devuelve el error en yerr usando el metodo de 4to orden.
        INTEGER :: i
        REAL :: ak2(NMAX), ak3(NMAX), ak4(NMAX), ak5(NMAX), ak6(NMAX), ytemp(NMAX)
        REAL, PARAMETER :: A2 = 0.2, A3 = 0.3, A4 = 0.6, A5 = 1.0, A6 = 0.875
        REAL, PARAMETER :: B21 = 0.2, B31 = 3.0/40.0, B32 = 9.0/40.0
        REAL, PARAMETER :: B41 = 0.3, B42 = -0.9, B43 = 1.2
        REAL, PARAMETER :: B51 = -11.0/54.0, B52 = 2.5, B53 = -70.0/27.0, B54 = 35.0/27.0
        REAL, PARAMETER :: B61 = 1631.0/55296.0, B62 = 175.0/512.0, B63 = 575.0/13824.0
        REAL, PARAMETER :: B64 = 44275.0/110592.0, B65 = 253.0/4096.0
        REAL, PARAMETER :: C1 = 37.0/378.0, C3 = 250.0/621.0, C4 = 125.0/594.0, C6 = 512.0/1771.0
        REAL, PARAMETER :: DC1 = C1 - 2825.0/27648.0, DC3 = C3 - 18575.0/48384.0
        REAL, PARAMETER :: DC4 = C4 - 13525.0/55296.0, DC5 = -277.0/14336.0, DC6 = C6 - 0.25

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
        INTEGER, INTENT(IN) :: nvar
        INTEGER, INTENT(OUT) :: nok, nbad
        REAL, INTENT(IN) :: eps, h1, hmin, x1, x2
        REAL, INTENT(INOUT) :: ystart(nvar)
        EXTERNAL derivs, rkqs

        INTEGER, PARAMETER :: MAXSTP = 10000, NMAX = 50
        REAL, PARAMETER :: TINY = 1.e-30

        ! Driver Runge-Kutta adaptativo. Integra ystart desde x1 a x2 con precision eps.
        INTEGER :: i, nstp
        REAL :: h, hdid, hnext, x, dydx(NMAX), y(NMAX), yscal(NMAX)

        x = x1
        h = SIGN(h1, x2 - x1)
        nok = 0
        nbad = 0
        
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

            IF (hdid == h) THEN
                nok = nok + 1
            ELSE
                nbad = nbad + 1
            END IF

            IF ((x - x2) * (x2 - x1) >= 0.0) THEN ! Terminamos?
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
        REAL, INTENT(IN) :: x, y(1)
        REAL, INTENT(OUT) :: dydt(1)
        REAL :: a
        a = 3.3
        dydt(1) = a * y(1)
    END SUBROUTINE derivs

END PROGRAM RK5
