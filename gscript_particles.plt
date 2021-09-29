# Gnuplot script file for visualization of particles file.
# usage on command line : gnuplot -e 'nodes="particles_filename"' gscript_particles.plg

set terminal postscript eps color enhanced "Arial" 20 size 3.5in,3.5in
set output "particles.eps"

if ( !exists("nodes") ) nodes="ex.particles"

set size square

set xrange [-1.5 : 1.5]
set yrange [-1.5 : 1.5]
set xtics -1.5, 0.5, 1.5
set ytics -1.5, 0.5, 1.5
set xlabel "{/Arial-Italic x}"
set ylabel "{/Arial-Italic y}"

plot nodes with points pt 0 lc 3

