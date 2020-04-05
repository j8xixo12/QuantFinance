set terminal png font " Times_New_Roman,12 "
set output "plot.png"

set multiplot layout 2,1 title "Greeks" font "14"
unset key
set title "Delta"
set xlabel "Time"
set ylabel "Delta"
set key right bottom font "Times_New_Roman, 7"
plot \
"out.csv" using 1:2 with lines linewidth 2 title "exact", \
"out.csv" using 1:3 with lines linewidth 2 title "CN", \
"out.csv" using 1:4 with lines linewidth 2 title "CubicSpline" \

unset key
set title "Gamma"
set xlabel "Time"
set ylabel "Gamma"
set key right top font "Times_New_Roman, 7"
plot \
"out.csv" using 1:5 with lines linewidth 2 title "exact", \
"out.csv" using 1:6 with lines linewidth 2 title "CN", \
"out.csv" using 1:7 with lines linewidth 2 title "CubicSpline" \



