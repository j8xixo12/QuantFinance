set title "plot"
set xlabel "Time"
set ylabel "Price"
set terminal png font " Times_New_Roman,12 "
set output "plot.png"
set key right bottom 

plot \
"out.csv" using 1:2 with lines linewidth 2 title "exact", \
"out.csv" using 1:3 with lines linewidth 2 title "CN", \
"out.csv" using 1:4 with lines linewidth 2 title "CubicSpline" \