reset
set term png
set output 'test.png'
set grid 

file = 'out.dat'
# cols = int(system('head -1 '.file.' | wc -w'))
cols = 21
plot for [i = 2:cols] 'out.dat' using 1:i with lines notitle