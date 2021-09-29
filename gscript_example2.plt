# Gnuplot script file
set terminal postscript eps color enhanced "Arial" 20 size 7in,3.5in
set output "I_example2.eps"

set palette defined ( 0 '#000090',1 '#000fff',2 '#0090ff',3 '#0fffee',4 '#90ff70',5 '#ffee00',6 '#ff7000',7 '#ee0000',8 '#7f0000')

set pm3d map

set size square
set multiplot layout 1,2

set xrange [-1.5 : 1.5]
set yrange [-1.5 : 1.5]
set xtics -1.5, 0.5, 1.5
set ytics -1.5, 0.5, 1.5

set title "electric field intensity"
set xlabel "{/Arial-Italic x}"
set ylabel "{/Arial-Italic y}"
splot "Ie_xy.txt"

set title "magnetic field intensity"
set xlabel "{/Arial-Italic x}"
set ylabel "{/Arial-Italic y}"
splot "Ih_xy.txt"

unset multiplot
set terminal x11
