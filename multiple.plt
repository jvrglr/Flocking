set term png size 1600,1200
L=20
set xrange[0:L]
set yrange[0:L]
do for [t=1:300]{
set output "image.".t.".png"
plot "picture.".t.".dat" u 1:2:3:4 w vectors head filled lt rgb "black" lw 4 notitle
}
