array var_name[8]
var_name[1]="time"
var_name[2]="temperature"
var_name[3]="l_c"
var_name[4]="cryst rate"
var_name[5]="dc primary"
var_name[6]="dc secondary"
var_name[7]="dc melt"
var_name[8]="c"


set xrange [*:*]
set yrange [*:*]

ix=1
iy1=2
iy2=8
iy3=5
iy4=6
iy5=7

plot "trace.dat" u ix:iy2 w lines lw 2 lc "blue" title var_name[iy2] axis x1y1, \
     "trace.dat" u ix:iy3 w lines lw 2 lc "black" title var_name[iy3]  axis x1y2, \
     "trace.dat" u ix:iy4 w lines lw 2 lc "gray" title var_name[iy4]  axis x1y2, \
     "trace.dat" u ix:iy5 w lines lw 2 lc "brown" title var_name[iy5]  axis x1y2



# "trace.dat" u ix:iy1 w lines lw 2 lc "red"  title var_name[iy1] axis x1y2, \
    