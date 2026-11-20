#!/bin/bash
#SBATCH --job-name=fl11_movie
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=01:30:00
#SBATCH --exclude=82ndalltheway,christmastruce,gottmituns,redbaron,soldierofheaven,steelcommanders,stormtroopers,wolfpack
#SBATCH --output=/data/warehouse17/rmuenste/runs/fluid11/movie/render_R_proto_movie.log
# Particle movie of a case 11 run: check, render (4 workers), encode. Pattern of runs/batch10/movie/render_PL15_s1.sbatch.
# usage: sbatch render_R_proto_movie.sbatch [RUN [T_END [STRIDE]]]
#   RUN     run directory name under runs/fluid11 (default R_proto_movie)
#   T_END   the last particle frame must reach this time (default 7.5 s: refuses a run in progress)
#   STRIDE  every STRIDE-th log frame (default 1 = one image per logged frame, no interpolation)
# for another run also pass --output=.../movie/render_<RUN>.log on the sbatch command line
RUN=${1:-R_proto_movie}; TEND=${2:-7.5}; STRIDE=${3:-1}
F=/data/warehouse17/rmuenste/runs/fluid11
FR=$F/movie/frames_$RUN
MP4=$F/movie/fluid11_${RUN}_particles.mp4
TOOL="python3 /data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-fluid11/applications/q2p1_el_pipeflow/validation_cases/s11_min_fluidization/tools/fl11_movie.py $F/$RUN $FR --stride $STRIDE"
FFMPEG=/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/blender_viz/fac3d_movie/tools_bin/ffmpeg-static
hostname; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
mkdir -p $FR
# the check reads the logs once and writes the cache that the workers reuse
$TOOL --check > $F/movie/check_$RUN.log || exit 1
cat $F/movie/check_$RUN.log
read n tlast <<< "$(awk '/^check: frames/{print $3, $7}' $F/movie/check_$RUN.log)"
echo "images $n, last frame t = $tlast s (required $TEND)"
awk -v a="$tlast" -v b="$TEND" 'BEGIN{exit !(a >= b - 1e-6)}' || { echo "run $RUN has not reached $TEND s: not rendered"; exit 3; }
old=$(/bin/ls $FR | grep -c '^frame_.*\.png$')
[ "$old" -le "$n" ] || { echo "$old png frames in $FR, more than the $n of this render: move them away first"; exit 4; }
W=4; q=$(( (n + W - 1) / W ))
for w in $(seq 0 $((W - 1))); do
  $TOOL --frames $((w * q)) $(( (w + 1) * q - 1 )) &
done
wait
m=$(/bin/ls $FR/frame_*.png | wc -l); echo "png frames $m"
[ "$m" = "$n" ] || exit 2
$FFMPEG -y -loglevel error -framerate 30 -i $FR/frame_%05d.png -vf "tpad=stop_mode=clone:stop_duration=2" -c:v libx264 -crf 18 -profile:v high -pix_fmt yuv420p -movflags +faststart $MP4
status=$?; echo "ffmpeg exit=$status $MP4"; exit $status
