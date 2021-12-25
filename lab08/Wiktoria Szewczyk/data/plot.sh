set term png 
set xlabel "t"
set ylabel "c(t), x_s_r(t)"
set out "t_C_Xsr_D=0.000000.png"
plot "t_C_Xsr_D=0.000000.dat"  u 1:2 title "c(t)" with lines, "t_C_Xsr_D=0.000000.dat"  u 1:3 title "x_s_r(t)" with lines

set out "t_C_Xsr_D=0.100000.png"
plot "t_C_Xsr_D=0.100000.dat"  u 1:2 title "c(t)" with lines, "t_C_Xsr_D=0.100000.dat"  u 1:3 title "x_s_r(t)" with lines


reset
set term png
set size ratio -1
set xlabel "x"
set ylabel "y"
set pm3d map

set output "u_D=0.000000_1.png"
splot "u_D=0.000000_1.dat" u 1:2:3 w pm3d title "u(x,y,t1)"

set output "u_D=0.000000_2.png"
splot "u_D=0.000000_2.dat" u 1:2:3 w pm3d title "u(x,y,t2)"

set output "u_D=0.000000_3.png"
splot "u_D=0.000000_3.dat" u 1:2:3 w pm3d title "u(x,y,t3)"

set output "u_D=0.000000_4.png"
splot "u_D=0.000000_4.dat" u 1:2:3 w pm3d title "u(x,y,t4)"

set output "u_D=0.000000_5.png"
splot "u_D=0.000000_5.dat" u 1:2:3 w pm3d title "u(x,y,t5)"


set output "u_D=0.100000_1.png"
splot "u_D=0.100000_1.dat" u 1:2:3 w pm3d title "u(x,y,t1)"

set output "u_D=0.100000_2.png"
splot "u_D=0.100000_2.dat" u 1:2:3 w pm3d title "u(x,y,t2)"

set output "u_D=0.100000_3.png"
splot "u_D=0.100000_3.dat" u 1:2:3 w pm3d title "u(x,y,t3)"

set output "u_D=0.100000_4.png"
splot "u_D=0.100000_4.dat" u 1:2:3 w pm3d title "u(x,y,t4)"

set output "u_D=0.100000_5.png"
splot "u_D=0.100000_5.dat" u 1:2:3 w pm3d title "u(x,y,t5)"
