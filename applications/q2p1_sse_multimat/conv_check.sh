#!/usr/bin/env bash

j=0
t="plot "
s="plot "

mkdir _temp

for i in _data/prot_0*txt
do
 let j=j+1
 outf=_prot0/${j}.txt 
 echo $i, $j
 grep ressureDiff $i > _temp/${j}.pres
 grep acting $i > _temp/${j}.torque
 t="$t '_temp/${j}.pres' u 5 w l notitle," 
 s="$s '_temp/${j}.torque' u 5 w l notitle,"
# echo ${t}
done

t="$t 0.0 w l lc rgb 'black'  notitle" 
s="$s 0.0 w l lc rgb 'black'  notitle" 
echo ${t}
echo ${s}

echo 'set terminal png size 1600,1200 ' > _temp/plot_it.gp
echo 'set output "convergence.png" ' >> _temp/plot_it.gp

echo 'set multiplot layout 2, 1 title "Convergence" font ",14"'  >> _temp/plot_it.gp

echo 'set title "Pressure"'  >> _temp/plot_it.gp
echo ${t} >> _temp/plot_it.gp

echo " " >> _temp/plot_it.gp

echo 'set title "Power"'  >> _temp/plot_it.gp
echo ${s} >> _temp/plot_it.gp

echo 'unset multiplot ' >> _temp/plot_it.gp

gnuplot _temp/plot_it.gp
rm -fr _temp

exit 0
