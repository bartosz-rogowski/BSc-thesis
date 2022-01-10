# if noh test is performed then parameter below should be set to 1
# otherwise to 0
noh = 1
reset
set term png size 1200, 1600
set out 'density.png'

# set yrange [0:14]
set xlabel "x"
set ylabel "gęstość"
if (noh == 0){
    set title "Test Soda - gęstość"
}
else{
    set title "Test Noha - gęstość"
}
set grid

n=121-1
set multiplot layout 4,2 rowsfirst

# plot "density.dat" i 0 u 1:2 w l lw 2 t sprintf("it=%i",0)
# plot "density.dat" i 50 u 1:2 w l lw 2 t sprintf("it=%i",50)
do for [k=0:n:20] {
    plot "density.dat" i k u 1:2 w l lw 2 t sprintf("it=%i",k)
} 
# plot "density.dat" i n u 1:2 w l lw 2 t sprintf("it=%i",n)
unset multiplot

if (noh == 0){
    reset
    set term png size 1200, 1600

    set out 'bulk_velocity.png'

    # set xrange [0:4]
    set yrange [-80:600]
    set xlabel "x"
    set ylabel "prędkość masowa"
    set title "Test Soda - prędkość masowa"
    set grid
    set multiplot layout 4,2 rowsfirst
    plot "bulk_velocity.dat" i 0 u 1:2 w l lw 2 t sprintf("it=%i",0)
    plot "bulk_velocity.dat" i 50 u 1:2 w l lw 2 t sprintf("it=%i",50)
    do for [k=100:n:100] {
        plot "bulk_velocity.dat" i k u 1:2 w l lw 2 t sprintf("it=%i",k)
    } 
    plot "bulk_velocity.dat" i n u 1:2 w l lw 2 t sprintf("it=%i",n)
    unset multiplot
}


if (noh == 1){
    reset
    set term png size 1200, 1600
    set yrange [0:12]
    set out 'density_radial.png'
    set grid
    set key left top
    set xlabel "promień"
    set ylabel "gęstość"
    set title "Test Noha - gęstość radialna"
    set multiplot layout 4,2 rowsfirst
    do for [k=0:n:20] {
        plot "density_radial.dat" i k u 1:2 w l lw 2 t sprintf("it=%i",k)
    } 
    unset multiplot

    reset
    set term png size 1200, 1600
    set out 'velocity_radial.png'
    set key left bottom
    set grid
    set xrange [0:0.56]
    set yrange [-0.15:0.05]
    set xlabel "promień"
    set ylabel "prędkość radialna"
    set title "Test Noha - prędkość radialna"
    set multiplot layout 4,2 rowsfirst
    do for [k=0:n:20] {
        plot "velocity_radial.dat" i k u 1:2 w l lw 2 t sprintf("it=%i",k)
    } 
    unset multiplot
}