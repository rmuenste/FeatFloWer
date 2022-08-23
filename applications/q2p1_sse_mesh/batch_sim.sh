#!/usr/bin/env bash
usage()
{
cat << EOF
usage: $0 options
./batch_sim.sh [-f Folder]

OPTIONS:
   -f      Output folder
EOF
}

while [ $# != 0 ]; do
 flag="$1"
 case "$flag" in
   -f) if [ $# -eq 2 ]; then
       Folder="$2"
       shift
       shift
      else
       echo "You did not provide enough arguments for the -XY flag"
       usage
       exit 1
      fi
#       echo "You supplied the arguments for the -f flag: $xFields,$yFields"
   ;;
   *) echo "Unrecognized flag or argument: $flag"
      usage
      exit 1
   ;;
 esac
 shift
done

cp -r MASTER Case_${Folder}

command="s/305-P1/"${Folder}"/g"
echo ${command}
sed -i `echo ${command}` Case_${Folder}/slurm_batch.sh
cd Case_${Folder}

#less slurm_batch.sh

sbatch slurm_batch.sh
