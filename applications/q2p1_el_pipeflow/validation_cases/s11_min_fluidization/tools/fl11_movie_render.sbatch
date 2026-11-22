#!/bin/bash
#SBATCH --job-name=fl11_movie
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=01:30:00
#SBATCH --prefer=nx
#SBATCH --exclude=worldgames,82ndalltheway,christmastruce,gottmituns,redbaron,soldierofheaven,steelcommanders,stormtroopers,wolfpack
#SBATCH --output=/data/warehouse17/rmuenste/runs/fluid11/movie/render-%j.log
# Particle movie of a case 11 run: check, render (4 workers), encode. Pattern of runs/batch10/movie/render_PL15_s1.sbatch.
# usage: sbatch render_R_proto_movie.sbatch [RUN [T_END [STRIDE]]]
#   RUN     run directory name under runs/fluid11 (default R_proto_movie)
#   T_END   the time the movie must reach (default: the protocol end of the evaluation, printed on the check
#           line). The end check uses the SAME --check invocation: the run finished (finish marker), its solver
#           end t0 + K dt reaches T_END, the last selected frame reaches T_END, and the last RAW frame of the
#           particle logs does not lie beyond the solver end (frames of an earlier, longer run; a stride could
#           hide them from the selected frames). A run in progress is refused.
#   STRIDE  every STRIDE-th log frame (default 1 = one image per logged frame, no interpolation)
# stdout goes to movie/render-<job id>.log (one file per job; pass --output=... to name it otherwise)
# Every render writes into a fresh frame directory frames_<RUN>.<job id> (never into the frames of an
# earlier render), waits for each worker and stops on the first failed one, and encodes only after
# all n frames of this render exist (frame_00000 .. frame_<n-1>). The check log and the temporary mp4
# live in the job-specific frame directory (the same filesystem as the movie, so the final mv is a
# rename): concurrent renders of one run share no file but the final movie, which is replaced only by
# a successful render; the status of that rename is the status of the job.
RUN=${1:-R_proto_movie}; TEND=${2:-}; STRIDE=${3:-1}
F=${FL11_RUNS:-/data/warehouse17/rmuenste/runs/fluid11}   # FL11_RUNS: another runs root (tests of this job)
FR=$F/movie/frames_$RUN.${SLURM_JOB_ID:-$$}
MP4=$F/movie/fluid11_${RUN}_particles.mp4
TOOL="python3 /data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-fluid11/applications/q2p1_el_pipeflow/validation_cases/s11_min_fluidization/tools/fl11_movie.py $F/$RUN $FR --stride $STRIDE"
FFMPEG=/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/blender_viz/fac3d_movie/tools_bin/ffmpeg-static
hostname; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
[ -e $FR ] && { echo "frame directory $FR exists: not rendered"; exit 4; }
mkdir -p $FR; echo "frames into $FR"
# the check reads the logs once and writes the cache that the workers reuse
CHECK=$FR/check.log
$TOOL --check > $CHECK || exit 1
cat $CHECK
read n tlast send pend fin traw <<< "$(awk '/^check: frames/{print $3, $7, $13, $15, $17, $19}' $CHECK)"
[ -n "$TEND" ] || TEND=$pend
echo "images $n, last selected frame t = $tlast s, last raw frame $traw s, solver end $send s, finished $fin, protocol end $pend s (required $TEND)"
[ "$fin" = 1 ] || { echo "run $RUN has not finished (no finish marker): not rendered"; exit 3; }
awk -v b="$TEND" 'BEGIN{exit !(b + 0 == b && b > 0)}' || { echo "required end time unknown ($TEND): not rendered"; exit 3; }
awk -v s="$send" -v b="$TEND" 'BEGIN{exit !(s + 0 == s && s >= b - 1e-6)}' || { echo "the solver clock of run $RUN ends at $send s, before $TEND s: not rendered"; exit 3; }
awk -v a="$tlast" -v b="$TEND" 'BEGIN{exit !(a >= b - 1e-6)}' || { echo "the particle frames of run $RUN end at $tlast s, before $TEND s: not rendered"; exit 3; }
awk -v a="$traw" -v s="$send" 'BEGIN{exit !(a + 0 == a && a <= s + 1e-6)}' || { echo "a particle frame at $traw s lies beyond the solver end $send s (frames of another run): not rendered"; exit 3; }
W=4; q=$(( (n + W - 1) / W )); pids=()
for w in $(seq 0 $((W - 1))); do
  $TOOL --frames $((w * q)) $(( (w + 1) * q - 1 )) > $FR/worker_$w.log 2>&1 &
  pids+=($!)
done
fail=0
for w in "${!pids[@]}"; do
  if wait "${pids[$w]}"; then echo "worker $w ok"; else rc=$?; echo "worker $w (pid ${pids[$w]}) FAILED, exit $rc:"; tail -5 $FR/worker_$w.log; fail=1; fi
done
[ $fail = 0 ] || { echo "a worker failed: not encoded"; exit 5; }
m=$(/bin/ls $FR | grep -c '^frame_.*\.png$'); echo "png frames $m (expected $n)"
for ((i = 0; i < n; i++)); do [ -s $(printf "$FR/frame_%05d.png" $i) ] || { echo "frame $i missing or empty: not encoded"; exit 2; }; done
[ "$m" = "$n" ] || { echo "$m frames in $FR, expected $n: not encoded"; exit 2; }
TMP=$FR/movie.part.mp4
$FFMPEG -y -loglevel error -framerate 30 -i $FR/frame_%05d.png -vf "tpad=stop_mode=clone:stop_duration=2" -c:v libx264 -crf 18 -profile:v high -pix_fmt yuv420p -movflags +faststart $TMP
status=$?; echo "ffmpeg exit=$status"
[ $status = 0 ] || exit $status
mv -f $TMP $MP4; status=$?
[ $status = 0 ] && echo "movie $MP4 (frames $FR)" || echo "rename of $TMP to $MP4 FAILED, exit $status"
exit $status
