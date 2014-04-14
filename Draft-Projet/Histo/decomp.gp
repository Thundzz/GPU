#!/opt/local/bin/gnuplot

set term pdf color enhanced #monochrome dashed
set output "histo.pdf"
set yrange [0:100]
set format y "%.0f %%"
set style data histograms
set style histogram rowstacked
set style fill solid border -1
set boxwidth 0.8
set size ratio 1.2

set origin 0.0,0.0
set key invert invert reverse Left outside title "Kernel"
set key autotitle columnheader #box

set ylabel "% of total time"
set xlabel "Platform"
plot "decomp.data" using ((100.*$2)/($2+$3+$4+$5+$6+$7+$8)):xtic(1) title "reset",\
     for [i=3:8] '' using (100.*column(i)/($2+$3+$4+$5+$6+$7+$8)) title column(i)
