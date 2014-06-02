#!/opt/local/bin/gnuplot

set terminal png 
set output "omp3.png"
set grid x y
set yrange [0:80]
set xrange [0:20000]
set ylabel "Temps d'exécution par itération (ms)"
set xlabel "Nombre d'atomes"

set style data linespoints
set format y "%.0f ms"
set style fill solid border -1
set boxwidth 0.8
set size ratio 1
set origin 0.0,0.0
set multiplot
set lmargin 10
set rmargin 2

set key at screen 0.87,screen 0.9
plot "atoms_omp.dat" using 1:($2/1000) title "OpenMP"
set key at screen 0.87,screen 0.86
plot "atoms_ompZ.dat" using 1:($2/1000) lt 2 title "OpenMP tri Z"
set key at screen 0.87,screen 0.82
plot "atoms_ompbox.dat" using 1:($2/1000) lt 3 title "OpenMP tri boites"
