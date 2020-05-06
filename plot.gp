reset
set term png size 1024, 768
set output 'test.png'
set grid 
set xlabel 'T (year)'
set ylabel 'Price'

file = 'out.dat'
cols = int(system('head -1 '.file.' | wc -w'))
# cols = 21
plot for [i = 2:cols] 'out.dat' using 1:i with lines notitle

reset
set term png size 1024, 768
set output 'test1.png'
set grid 
set xlabel 'lag K'
set ylabel 'ACF'
plot 'out2.dat' with lines notitle