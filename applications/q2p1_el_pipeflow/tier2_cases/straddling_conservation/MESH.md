# Mesh: shared QBOX9 periodic box

This case intentionally reuses the `QBOX9` 3x3x3 partitioned periodic box from
`tier2_cases/momentum_conservation/mesh`.

The particle center is placed just inside one PE-owned subdomain while its EL
kernel support crosses the x=1/3 partition face. This keeps PE ownership
unambiguous and exercises distributed EL deposit/feedback across rank-owned
mesh cells.
