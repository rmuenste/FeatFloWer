#!/bin/bash

while [ $# != 0 ]; do
 flag="$1"
 case "$flag" in
  -n) if [ $# -gt 1 ]; then
       end_process="$2"
       shift
      else
       echo "You did not provide an argument for the -n flag"
        usage
        exit 1
      fi
  ;;
  -s) if [ $# -gt 1 ]; then
       n_slots="$2"
       shift
      else
       echo "You did not provide an argument for the -s flag"
        usage
        exit 1
      fi
   ;;
   *) echo "Unrecognized flag or argument: $flag"
      exit 1
   ;;
 esac
 shift
done

if [[ -z $n_slots ]]
then
 echo "no number of slots per nodes defined // exit"
 exit 0
fi

if [[ -z $end_process ]]
then
 echo "no number of processors defined // exit"
 exit 0
fi

node_list_array=($(scontrol show hostnames $SLURM_JOB_NODELIST))

# Print the character array
for element in "${node_list_array[@]}"; do
    echo "$element"
done
node_list_size=${#node_list_array[@]}

# echo node_list_size=$node_list_size

if [[ $((node_list_size * n_slots)) -lt $end_process ]]
then
 echo "too less nodes or slots are available // exit"
 exit 0
fi

end_process=$((end_process-1))

rankfile='myRankFile'

ihostname="${node_list_array[${node_list_size}-1]}"
line="${ihostname}:1"

echo "$line"
echo "$line" > ${rankfile}

# Loop over the process numbers
for ((iP = 0; iP <= node_list_size-2; iP++)); do
    ihostname="${node_list_array[${iP}]}"
    line="${ihostname}:${n_slots}" 
    echo "$line"
    echo "$line" >> ${rankfile}
done

ihostname="${node_list_array[${node_list_size}-1]}"
remainingnodes=$((end_process -(node_list_size-1) * n_slots))
if [[ ${remainingnodes} -gt 0 ]]
then
 line="${ihostname}:${remainingnodes}"
 echo "$line"
 echo "$line" >> ${rankfile}
fi

# # Loop over the process numbers
# for ((iP = 0; iP <= end_process-1; iP++)); do
#     jP=$((iP + 1))
#     islot=$((iP % n_slots))
#     index=$(( (iP - islot) / n_slots ))
#     ihostname="${node_list_array[${index}]}"
#     line="rank ${jP}=${ihostname} slot=${islot}" 
#      echo "$line"
#     echo "$line" >> ${rankfile}
# done

exit 0
