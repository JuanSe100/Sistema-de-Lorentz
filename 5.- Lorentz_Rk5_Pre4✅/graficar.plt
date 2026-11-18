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

# 3. GRAFICADO  every grafica un parte de los puntos 2 la mitas 3 la tercera parte y asi
splot "datos.dat" using 2:3:4 every 1 w lp pt 7 ps 0.2 lw 0.3 lc rgb "red"

# 4. PAUSA (Mantiene la ventana abierta hasta que presiones Enter en la consola)
pause -1
