set term png

set xl "x"
set yl "y"

set xrange [-2:27.6]
set yrange [-2:27.6]

set out "V_k16.png"
set title "V(x,y) k=16"
plot "V_k16.txt" u 1:2:3 notitle with image

set out "V_k8.png"
set title "V(x,y) k=8"
plot "V_k8.txt" u 1:2:3 notitle with image

set out "V_k4.png"
set title "V(x,y) k=4"
plot "V_k4.txt" u 1:2:3 notitle with image

set out "V_k2.png"
set title "V(x,y) k=2"
plot "V_k2.txt" u 1:2:3 notitle with image

set out "V_k1.png"
set title "V(x,y) k=1"
plot "V_k1.txt" u 1:2:3 notitle with image