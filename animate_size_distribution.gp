reset

nfile=system("ls size_distribution[0-9]*.dat | wc -l")+0
nskip=1; #nfile/1;

set pointsize 1
do for [j=1:1000] {
do for [i=1:nfile:nskip] {
file=sprintf('temperature%d.dat', i)
temperature=system("cat " . file . " | awk '{print $1}'")+0

l_min="`head -1 size_distribution_limits.dat | awk '{print $1}'`"
l_max="`head -1 size_distribution_limits.dat | awk '{print $2}'`"
n_max="`head -1 size_distribution_limits.dat | awk '{print $3}'`"
 
unset origin
unset size
set multiplot

set xrange [l_min:l_max]
set yrange [0:n_max]

#set xrange [0:5]
#set yrange [0:0.4]
set pointsize 1
plot sprintf('size_distribution%d.dat', i) using 1:2 with points pointtype 5 title sprintf('size_distribution%d.dat', i), "crystallsation_properties.dat" u 2:1 w lines axis x1y2, "crystallsation_properties.dat" u 2:4 w lines axis x1y2, temperature axis x1y2, sprintf('melt_temperature%d.dat', i) using 1:2 w linesp axis x1y2 title "melt/cryst size"


set origin 0.3, 0.3
set size 0.5,0.5
load "../trace.gp"
replot sprintf('time%d.dat', i) w lines lc "black"

unset multiplot

pause -1 # 0.5

}
}


#stats "size_distribution.dat" nooutput
#set xrange [STATS_min_x:STATS_max_x]
#set yrange [STATS_min_y:STATS_max_y]
#do for [a = 1: int(STATS_blocks - 1)] {
# plot "size_distribution.dat" u 1:2 index a
# pause 0.1
#}

