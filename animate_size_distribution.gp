reset

set pointsize 1
do for [j=1:1000] {
do for [i=1:9999:100] {
file=sprintf('temperature%d.dat', i)
temperature=system("cat " . file . " | awk '{print $1}'")+0

unset origin
unset size
set multiplot

set xrange [0:5]
set yrange [0:0.4]
plot sprintf('size_distribution%d.dat', i) using 1:2 with points pointtype 5 title sprintf('size_distribution%d.dat', i), "crystallsation_properties.dat" u 2:1 axis x1y2, temperature axis x1y2, sprintf('melt_temperature%d.dat', i) using 1:2 w linesp axis x1y2 title "melt/cryst size"

#unset object 1
set origin 0.5, 0.3
set size 0.5,0.5
#plot sin(x)
load "../trace.gp"
replot sprintf('time%d.dat', i) w lines lc "black"

unset multiplot


pause -1 # 0.5

}
}


#, temperature  axis x1y2

#plot sprintf('size_distribution%d.dat', i) using 1:2 with linespoints pointtype 5 smooth unique

#temperature="`cat six.dat | awk '{print $1}'`"+0
#print temperature


#stats "size_distribution.dat" nooutput
 
#set xrange [STATS_min_x:STATS_max_x]
#set yrange [STATS_min_y:STATS_max_y]

#do for [a = 1: int(STATS_blocks - 1)] {
# plot "size_distribution.dat" u 1:2 index a
# pause 0.1
#}

