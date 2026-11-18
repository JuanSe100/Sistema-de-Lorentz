set xlabel "Tiempo"
set ylabel "y(t)"
set title "Gráfica de y1 e y2 en función del tiempo"
set xrange [900:1000]

plot "datos_Programa3.dat" using 1:2 with lines lc rgb "red" title "y1", \
     "datos_Programa3.dat" using 1:3 with lines lc rgb "blue" title "y2"

pause -1
