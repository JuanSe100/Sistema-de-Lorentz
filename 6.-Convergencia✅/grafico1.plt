# 0. ARCHIVO Y VARIABLES DEL ATRACTOR DE LORENZ
archivo = "caotico.dat"
sigma = 10.0
rho = 29.5
beta = 8.0/3.0
rho_crit = (sigma * (sigma + beta + 3.0)) / (sigma - beta - 1.0)   # ~24.7368

# 1. ESTÉTICA Y TÍTULOS
set title "Atractor de Lorenz" font ",16"
set xlabel "Eje X (x)"
set ylabel "Eje Y (y)"
set zlabel "Eje Z (z)"

# 2. CONFIGURACIÓN DEL ESPACIO 3D
set grid
set view 60, 45      
set ticslevel 0      
set key off          

# 3. GRAFICADO SEGÚN EL RÉGIMEN DE CONVERGENCIA
if (rho <= 1.0) {
    # Caso 1: Converge al origen (0,0,0)
    splot archivo using 2:3:4 every 3 w lp pt 7 ps 0.2 lw 0.3 lc rgb "blue", \
          '< echo "0 0 0"' using 1:2:3 w p pt 7 ps 1.5 lc rgb "red"

} else if (rho > 1.0 && rho < rho_crit) {
    # Caso 2: Converge a uno de los puntos fijos P+ o P- (verde) y el origen es inestable (rojo)
    xx = sqrt(beta * (rho - 1.0))
    zz = rho - 1.0
    splot archivo using 2:3:4 every 3 w lp pt 7 ps 0.2 lw 0.3 lc rgb "blue", \
          '< echo "0 0 0"' using 1:2:3 w p pt 7 ps 1.5 lc rgb "red", \
          '+' using (xx):(xx):(zz) every ::0::0 w p pt 7 ps 1.5 lc rgb "green", \
          '+' using (-xx):(-xx):(zz) every ::0::0 w p pt 7 ps 1.5 lc rgb "green"

} else {
    # Caso 3: Régimen caótico (no converge a ningún punto estable, solo se grafica la trayectoria)
    splot archivo using 2:3:4 every 3 w lp pt 7 ps 0.2 lw 0.3 lc rgb "blue"
}

pause -1

