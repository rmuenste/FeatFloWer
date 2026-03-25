mkdir -p results/cube64_force1e-1
cp build2/applications/q2p1_dns_drag/particle_force.log results/cube64_force1e-1/
cp build2/applications/q2p1_dns_drag/bulk_flow.log results/cube64_force1e-1/
cp applications/q2p1_dns_drag/_data/q2p1_param.dat results/cube64_force1e-1/
cp applications/q2p1_dns_drag/example.json results/cube64_force1e-1/
python3 applications/q2p1_dns_drag/summarize_dns_drag_case.py \
  --force-log results/cube64_force1e-1/particle_force.log \
  --bulk-log results/cube64_force1e-1/bulk_flow.log \
  --last 100 \
  --label cube64_force1e-1 | tee results/cube64_force1e-1/summary.txt
