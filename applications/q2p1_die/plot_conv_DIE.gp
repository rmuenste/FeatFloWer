set terminal png size 1200,800 linewidth 2.0 font 'Arial,24'
set output 'convergence_M+S_flowrate.png'

set xlabel 'timesteps'
set ylabel 'throughput [kg/h]'
plot 'flow' u 4 w l notitle




