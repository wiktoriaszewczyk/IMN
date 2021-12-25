set term png

set xl "x"
set yl "y"

# set xrange [-2:27.6]
# set yrange [-2:27.6]

set out "5a.png"
set title "V(x,y) 5a"
plot "5a.txt" u 1:2:3 notitle with image

set out "5b.png"
set title "V(x,y) 5b"
plot "5b.txt" u 1:2:3 notitle with image

set out "5c.png"
set title "V(x,y) 5c"
plot "5c.txt" u 1:2:3 notitle with image

set out "6a.png"
set title "V(x,y) 6a"
plot "6a.txt" u 1:2:3 notitle with image

set out "6b.png"
set title "V(x,y) 6b"
plot "6b.txt" u 1:2:3 notitle with image

set out "6c.png"
set title "V(x,y) 6c"
plot "6c.txt" u 1:2:3 notitle with image