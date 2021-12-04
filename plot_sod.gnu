reset
set term gif size 800, 800 animate delay 100

set out 'density.gif'

set xrange [0:4]
set yrange [0:1.5]
set xlabel "x"
set ylabel "density"
set title "2D Sod test - density"
set grid
# plot 'density.dat' i 0 u 1:2 w l lw 2 t ''

n=20-1

do for [k=0:n:2] {
    plot "density.dat" i k u 1:2 w l lw 2 t sprintf("it=%i",k)
} 


reset
set term gif size 800, 800 animate delay 100

set out 'bulk_velocity.gif'

set xrange [0:4]
set yrange [-0.04:0.1]
set xlabel "x"
set ylabel "bulk velocity"
set title "2D Sod test - bulk velocity"
set grid

do for [k=0:n:2] {
    plot "bulk_velocity.dat" i k u 1:2 w l lw 2 t sprintf("it=%i",k)
} 