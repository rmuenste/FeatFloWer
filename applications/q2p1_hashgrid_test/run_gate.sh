#!/usr/bin/env bash
# Run the Layer-2 FBM classification gate from a staged build directory.
#   cmake --build <tree> --target q2p1_hashgrid_test_stage
#   cd <tree>/applications/q2p1_hashgrid_test && ./run_gate.sh [np]
# Exit code is the solver's: 0 green; 62 a verdict stage failed; 63 the gate compared nothing;
# 64 pe's grid was inactive or mis-sized on some rank (the harness exercised nothing).
set -u
NP="${1:-4}"   # QBOX9Z3 has 3 subdomains + the CFD master
mpirun -np "$NP" ./q2p1_hashgrid_test 2>&1 | tee gate_run.log
rc=${PIPESTATUS[0]}
echo "---- gate summary (exit $rc) ----"
grep -E "DNS_QUERY_GATE_PREFLIGHT|DNS_QUERY_GATE time=|FAILED" gate_run.log
exit "$rc"
