#!/opt/local/bin/gnuplot

set terminal png 
set grid x y
set output "ocl.png"
set yrange [0:50]
set xrange [0:35000]
set style data linespoints
set format y "%.0f ms"
set style fill solid border -1
set boxwidth 0.8
set size ratio 1

set origin 0.0,0.0

set multiplot
set lmargin 10
set rmargin 2
set ylabel "Temps d'exécution par itération (ms)"
set xlabel "Nombre d'atomes"
set key at screen 0.6,screen 0.9
plot "atoms_ompboxX.dat" using 1:($2/1000) title "OpenMP tri boites"
set key at screen 0.6,screen 0.86
plot "atoms_ocl.dat" using 1:($2/1000) lt 2 title "OpenCL"
