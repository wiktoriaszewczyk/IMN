#!/usr/bin/gnuplot

set term png enhanced size 600,300

set size ratio -1
set contour
unset surface
set view map
unset key

#-1000 zeta
set cntrparam levels increment -200,10,350
set cbr [-200:350]
set o "zeta-1000.png"
set title "Q=-1000, zeta"
sp "Q-1000.txt" u 1:2:4:4 w l lw 2 palette

#-1000 psi
set cntrparam levels increment -55,0.1,-50
set cbr [-55:-50]
set o "psi-1000.png"
set title "Q=-1000, psi"
sp "Q-1000.txt" u 1:2:3:3 w l lw 2 palette

#-1000 u(x, y) 
set cbr [-2:16]
set o "u-1000.png"
set title "Q=-1000, u(x, y)"
sp "Q-1000.txt" u 1:2:5 with image 

#-1000 v(x, y) 
set cbr [-6:1]
set o "v-1000.png"
set title "Q=-1000, v(x, y)"
sp "Q-1000.txt" u 1:2:6 with image 

#-4000 psi
set cntrparam levels increment -218,1,-202
set cbr [-218:-202]
set o "psi-4000.png"
set title "Q=-4000, psi"
sp "Q-4000.txt" u 1:2:3:3 w l lw 2 palette

#-4000 zeta
set cntrparam levels increment -800,20,1200
set cbr [-800:1200]
set o "zeta-4000.png"
set title "Q=-4000, zeta"
sp "Q-4000.txt" u 1:2:4:4 w l lw 2 palette

#-4000 u(x, y) 
set cbr [-10:70]
set o "u-4000.png"
set title "Q=-4000, u(x, y)"
sp "Q-4000.txt" u 1:2:5 with image 

#-4000 v(x, y) 
set cbr [-14:4]
set o "v-4000.png"
set title "Q=-4000, v(x, y)"
sp "Q-4000.txt" u 1:2:6 with image 


#4000 psi
set cntrparam levels increment 202,1,218
set cbr [202:218]
set o "psi4000.png"
set title "Q=4000, psi"
sp "Q4000.txt" u 1:2:3:3 w l lw 2 palette