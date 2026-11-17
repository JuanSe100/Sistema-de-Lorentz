PROGRAM RK5
      IMPLICIT NONE
!! Estamos trabajando con simple precision
!! Declaración de variables del programa
      INTEGER nvar, nok, nbad, i, m  
      REAL eps, h1, hmin, t, tnext, dt
      REAL ystart(1)
      

      nvar = 1          !! Numero de variables
      eps = 1.e-6       !!Tolerancia del programa
      h1 = 1.e-6        !! Tamaño del paso 
      hmin = 1.e-15     !! Paso Minimo permitido
      
      t = 0.0           !! Tiempo inicial de simulacion
      dt = 0.01         !! Salto de tiempo 
      m = 1000          !! Numero de saltos
      ystart(1) = 1.0   !! Condicion inicial
     
      OPEN(10, FILE='datos.dat')   !! Apertura del archivo de datos (10)
      WRITE(10, *) t, ystart(1)    !! Escritura del primer elementos CI
                                       
      DO i = 1, m                  !! m= 1000, mil pasos de 0.01
        tnext = t + dt             !! Actualización del paso
        CALL odeint(ystart,nvar,t,tnext,eps,h1,hmin,nok,nbad,derivs,rkqs) !! Llamada a la función
        t = tnext                   !!Actualizacion de t
        WRITE(10, *) t, ystart(1)   !! Impresión de resultados
      END DO
      
      CLOSE(10)                      !! Cerramos el archivo de datos
      PRINT *, 'Terminado. Revisa el archivo datos.dat' !! Avisa de finalizacion

      CONTAINS  !!Subrrutinas del programa(--FUNCIONES--)
      !!====================PRIMERA SUBRUTINA=================
      SUBROUTINE rkqs(y,dydx,n,x,htry,eps,yscal,hdid,hnext,derivs)
      INTEGER n,NMAX
      REAL eps,hdid,hnext,htry,x,dydx(n),y(n),yscal(n)
      EXTERNAL derivs
      PARAMETER (NMAX=50) ! Numero maximo de ecuaciones
!     Paso Runge-Kutta 5to orden con monitoreo de error de truncamiento
!     para asegurar la precision y ajustar el tamano de paso.
      INTEGER i
      REAL errmax,h,htemp,xnew,yerr(NMAX),ytemp(NMAX),SAFETY,PGROW,
     * PSHRNK,ERRCON
      PARAMETER (SAFETY=0.9,PGROW=-.2,PSHRNK=-.25,ERRCON=1.89e-4)
!     El valor ERRCON es igual a (5/SAFETY)**(1/PGROW)
      h=htry ! Establece el paso al valor de prueba inicial
    1 call rkck(y,dydx,n,x,h,ytemp,yerr,derivs) ! Da un paso
      
      errmax=0. ! Evalua la precision
      do 11 i=1,n
        errmax=max(errmax,abs(yerr(i)/yscal(i)))
   11 enddo
      errmax=errmax/eps ! Escala relativo a la tolerancia
      if(errmax.gt.1.)then ! Error grande, reduce el paso
        htemp=SAFETY*h*(errmax**PSHRNK)
        h=sign(max(abs(htemp),0.1*abs(h)),h) ! No mas de factor de 10
        xnew=x+h
        if(xnew.eq.x) stop 'underflow en h'
        goto 1 ! Intenta de nuevo
      else ! Paso exitoso. Calcula el proximo paso
        if(errmax.gt.ERRCON)then
          hnext=SAFETY*h*(errmax**PGROW)
        else ! Incremento no mayor a un factor de 5
          hnext=5.*h
        endif
        hdid=h
        x=x+h
        do 12 i=1,n
          y(i)=ytemp(i)
   12   enddo
        return
      endif
      END

      !!====================SEGUNDA SUBRUTINA=================
      SUBROUTINE rkck(y,dydx,n,x,h,yout,yerr,derivs)
      INTEGER n,NMAX
      REAL h,x,dydx(n),y(n),yerr(n),yout(n)
      EXTERNAL derivs
      PARAMETER (NMAX=50)
!     Metodo Runge-Kutta de Cash-Karp de quinto orden.
!     Avanza la solucion un intervalo h y devuelve variables en yout.
!     Devuelve el error en yerr usando el metodo de 4to orden.
      INTEGER i
      REAL ak2(NMAX),ak3(NMAX),ak4(NMAX),ak5(NMAX),ak6(NMAX),
     * ytemp(NMAX),A2,A3,A4,A5,A6,B21,B31,B32,B41,B42,B43,B51,
     * B52,B53,B54,B61,B62,B63,B64,B65,C1,C3,C4,C6,DC1,DC3,
     * DC4,DC5,DC6
      PARAMETER (A2=.2,A3=.3,A4=.6,A5=1.,A6=.875,B21=.2,B31=3./40.,
     * B32=9./40.,B41=.3,B42=-.9,B43=1.2,B51=-11./54.,B52=2.5,
     * B53=-70./27.,B54=35./27.,B61=1631./55296.,B62=175./512.,
     * B63=575./13824.,B64=44275./110592.,B65=253./4096.,
     * C1=37./378.,C3=250./621.,C4=125./594.,C6=512./1771.,
     * DC1=C1-2825./27648.,DC3=C3-18575./48384.,
     * DC4=C4-13525./55296.,DC5=-277./14336.,DC6=C6-.25)

      do 11 i=1,n ! Primer paso
        ytemp(i)=y(i)+B21*h*dydx(i)
   11 enddo
      call derivs(x+A2*h,ytemp,ak2) ! Segundo paso
      do 12 i=1,n
        ytemp(i)=y(i)+h*(B31*dydx(i)+B32*ak2(i))
   12 enddo
      call derivs(x+A3*h,ytemp,ak3) ! Tercer paso
      do 13 i=1,n
        ytemp(i)=y(i)+h*(B41*dydx(i)+B42*ak2(i)+B43*ak3(i))
   13 enddo
      call derivs(x+A4*h,ytemp,ak4) ! Cuarto paso
      do 14 i=1,n
        ytemp(i)=y(i)+h*(B51*dydx(i)+B52*ak2(i)+B53*ak3(i)+
     * B54*ak4(i))
   14 enddo
      call derivs(x+A5*h,ytemp,ak5) ! Quinto paso
      do 15 i=1,n
        ytemp(i)=y(i)+h*(B61*dydx(i)+B62*ak2(i)+B63*ak3(i)+
     * B64*ak4(i)+B65*ak5(i))
   15 enddo
      call derivs(x+A6*h,ytemp,ak6) ! Sexto paso
      do 16 i=1,n ! Acumula incrementos con los pesos adecuados
        yout(i)=y(i)+h*(C1*dydx(i)+C3*ak3(i)+C4*ak4(i)+
     * C6*ak6(i))
   16 enddo
      do 17 i=1,n
!       Estima error como diferencia entre 4to y 5to orden
        yerr(i)=h*(DC1*dydx(i)+DC3*ak3(i)+DC4*ak4(i)+DC5*ak5(i)
     * +DC6*ak6(i))
   17 enddo
      return
      END
      !!====================TERCERA SUBRUTINA=================
      SUBROUTINE odeint(ystart,nvar,x1,x2,eps,h1,hmin,nok,nbad,
     * derivs,rkqs)
      INTEGER nbad,nok,nvar,MAXSTP,NMAX
      REAL eps,h1,hmin,x1,x2,ystart(nvar),TINY
      EXTERNAL derivs,rkqs
      PARAMETER (MAXSTP=10000,NMAX=50,TINY=1.e-30)
!     Driver Runge-Kutta adaptativo. Integra ystart desde x1 a x2
!     con precision eps.
      INTEGER i,nstp
      REAL h,hdid,hnext,x,dydx(NMAX),y(NMAX),yscal(NMAX)
      
      x=x1
      h=sign(h1,x2-x1)
      nok=0
      nbad=0
      do 11 i=1,nvar
        y(i)=ystart(i)
   11 enddo
   
      do 16 nstp=1,MAXSTP
        call derivs(x,y,dydx)
        do 12 i=1,nvar ! Escalado para monitorear la precision
          yscal(i)=abs(y(i))+abs(h*dydx(i))+TINY
   12   enddo
        
!       Si el paso sobrepasa el limite x2, lo ajusta
        if((x+h-x2)*(x+h-x1).gt.0.) h=x2-x
        call rkqs(y,dydx,nvar,x,h,eps,yscal,hdid,hnext,derivs)
        
        if(hdid.eq.h)then
          nok=nok+1
        else
          nbad=nbad+1
        endif
        
        if((x-x2)*(x2-x1).ge.0.)then ! Terminamos?
          do 14 i=1,nvar
            ystart(i)=y(i)
   14     enddo
          return
        endif
        
        if(abs(hnext).lt.hmin) stop 'paso menor que minimo'
        h=hnext
   16 enddo
      stop 'demasiados pasos en odeint'
      return
      END

      !!====================CUARTA SUBRUTINA=================
      SUBROUTINE derivs(x,y,dydt)
	      REAL x,y(1),dydt(1), a
	      a=3.3
	      dydt(1)=a*y(1) 
      return
      END

END RK5
