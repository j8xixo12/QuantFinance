set terminal png size 1024, 628 font " Times_New_Roman,12 "
set output "plot.png"
set grid ytics xtics


set multiplot layout 3,1 title "Greeks" font "14"
unset key
set title "Delta"
set xlabel "x"
set ylabel "Delta"
set key right bottom font "Times_New_Roman, 7"
plot \
"out.csv" using 1:2 with lines linewidth 2 title "exact", \
"out.csv" using 1:3 with lines linewidth 2 title "Crank-Nicolson", \
"out.csv" using 1:4 with lines linewidth 2 title "Cubic-Spline" \

unset key
set title "Gamma"
set xlabel "x"
set ylabel "Gamma"
set yrange [0:0.06]
set key right top font "Times_New_Roman, 7"
plot \
"out.csv" using 1:5 with lines linewidth 2 title "exact", \
"out.csv" using 1:6 with lines linewidth 2 title "Crank-Nicolson", \
"out.csv" using 1:7 with lines linewidth 2 title "Cubic-Spline" \

unset key
unset yrange
set title "solution"
set xlabel "x"
set ylabel "price"
set key right top font "Times_New_Roman, 7"
plot \
"out.csv" using 1:8 with lines linewidth 2 \


