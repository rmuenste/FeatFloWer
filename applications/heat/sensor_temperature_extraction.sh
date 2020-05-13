#!/usr/bin/env bash
usage()
{
cat << EOF
usage: $0 options
gathergmvsnaps [-nMin min -nMax max -File filename -Res "x,y"]

OPTIONS:
   -nMin   Min Sensor Index
   -nMax   Max Sensor Index
   -File   Filename for png output
   -Res    Resolution of the outputpicture
EOF
}

# if [ $# -eq 0 ]
# then
#   usage
#   exit 1
# fi

grep "Machine/Element" _data/heat.s3d | grep -n Element | tail -1 > temp
maxS=`awk -F ":" '{print  $1 }' temp`
rm -f temp
echo "Max files: "$maxS
# exit 1

res='1600,1000'
File='sensor'
minS=1

mkdir _tmp

while [ $# != 0 ]; do
 flag="$1"
 case "$flag" in
  -nMin) if [ $# -gt 1 ]; then
        minS="$2"
        shift
       else
        echo "You did not provide an argument for the -nMin flag"
        usage
        exit 1
       fi
   ;;
  -nMax) if [ $# -gt 1 ]; then
        maxS="$2"
        shift
       else
        echo "You did not provide an argument for the -nMax flag"
        usage
        exit 1
       fi
   ;;
  -File) if [ $# -gt 1 ]; then
        File="$2"
        shift
       else
        echo "You did not provide an argument for the -File flag"
        usage
        exit 1
       fi
   ;;
  -Res) if [ $# -gt 1 ]; then
        res="$2"
        shift
       else
        echo "You did not provide an argument for the -Res flag"
        usage
        exit 1
       fi
   ;;
   *) echo "Unrecognized flag or argument: $flag"
      usage
      exit 1
   ;;
 esac
 shift
done

if [[ $minS -lt 0 ]] || [[ $maxS -lt 0 ]]
then
 echo "No limit are defined -nMin -nMax "
 usage
 exit 0
fi

File=${File}".png"

echo 'set terminal png size '${res}' font "Arial,24" linewidth 2.0' > plot_sensor.gp
echo 'set output "'${File}'"' >> plot_sensor.gp
echo 'set ylabel "Temperature [C]"'  >> plot_sensor.gp
echo 'set xlabel "Time [s]"'  >> plot_sensor.gp
echo 'set key bottom'  >> plot_sensor.gp
echo 'set key left'  >> plot_sensor.gp
echo 'plot \'  >> plot_sensor.gp

i=0
for ((j=${minS}; j<=${maxS};j++))
do
 grep "Sensor\["$j"\]" _data/prot.txt > _tmp/s$j.txt
 
 NumOfLines=`wc -l _tmp/s$j.txt | awk '{ print $1 }'`
 echo 'NumOfLines='$NumOfLines
 
 if [[ $NumOfLines -ne 0 ]] 
 then
  let i=i+1
#  if [ $j -lt ${maxS} ]
#  then
   echo '"_tmp/s'${j}'.txt" u 2:3 w l t "sensor_'${i}'",\'  >> plot_sensor.gp
#  else
#   echo '"_tmp/s'${j}'.txt" u 2:3 w l t "sensor_'${i}  >> plot_sensor.gp
#  fi
 fi
done

grep "IntQuantMELT" _data/prot.txt > _tmp/m.txt
echo '"_tmp/m.txt" u 2:5 w l lc rgb "black" t "melt'  >> plot_sensor.gp

gnuplot plot_sensor.gp
#eog ${File}
rm -rf plot_sensor.gp _tmp

exit 0
