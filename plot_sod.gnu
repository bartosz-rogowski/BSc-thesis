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

n=10-1

do for [k=0:n:1] {
    plot "density.dat" i k u 1:2 w l lw 2 t sprintf("it=%i",k)
} 

noh = 0
if (noh == 0){
    reset
    set term gif size 800, 800 animate delay 100

    set out 'bulk_velocity.gif'

    # set xrange [0:4]
    set yrange [-10:100]
    set xlabel "x"
    set ylabel "bulk velocity"
    set title "2D Sod test - bulk velocity"
    set grid

    do for [k=0:n:1] {
        plot "bulk_velocity.dat" i k u 1:2 w l lw 2 t sprintf("it=%i",k+1)
    } 
}


if (noh == 1){
    reset
    set term gif size 800, 800 animate delay 100
    set yrange [0:0.03]
    set out 'density_radial.gif'
    set grid
    set xlabel "radius"
    set ylabel "density"
    set title "2D Noh test - radial density"
    do for [k=0:n:5] {
        plot "density_radial.dat" i k u 1:2 w l lw 2 t sprintf("it=%i",k)
    } 
}