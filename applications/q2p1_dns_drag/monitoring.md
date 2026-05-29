# q2p1_dns_drag Log Monitoring Helpers

This directory contains small helper scripts for inspecting logs while a
`q2p1_dns_drag` run is still active. They tolerate partially written lines and
print summaries for the last complete timesteps.

## Particle force log

Use `monitor_force_log.py` for `particle_force.log`:

```bash
python3 applications/q2p1_dns_drag/monitor_force_log.py \
  build2/applications/q2p1_dns_drag/particle_force.log
```

Optional: report more complete timesteps:

```bash
python3 applications/q2p1_dns_drag/monitor_force_log.py \
  build2/applications/q2p1_dns_drag/particle_force.log \
  --last 10
```

If no path is given, the script reads `particle_force.log` in the current
working directory.

## Bulk-flow log

Use `monitor_bulk_flow_log.py` for `bulk_flow.log`:

```bash
python3 applications/q2p1_dns_drag/monitor_bulk_flow_log.py \
  build2/applications/q2p1_dns_drag/bulk_flow.log
```

Optional: report more complete timesteps:

```bash
python3 applications/q2p1_dns_drag/monitor_bulk_flow_log.py \
  build2/applications/q2p1_dns_drag/bulk_flow.log \
  --last 10
```

If no path is given, the script reads `bulk_flow.log` in the current working
directory.
