set xlabel "Tiempo"
set ylabel "|z1(t) - z2(t)|"
set title "Diferencia en la componente z "
set xrange [0:150]
set grid

set logscale y
plot "datos_Programa3.dat" using 1:4 with lines lc rgb "blue" title "|y1-y2|"

pause -1
