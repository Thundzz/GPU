#!/usr/bin/gnuplot

set logscale x 2
set terminal pdf
set output "fig.pdf"

set yrange [0 : 5]
plot "data" with linespoints
