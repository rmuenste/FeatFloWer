# Mesh for the long-horizon momentum-conservation case

This case reuses the committed `QBOX9` 3x3x3 periodic box from
`tier2_cases/momentum_conservation/mesh`.

The CMake Tier-2 staging target copies that shared mesh into the runtime
directory for `momentum_conservation_long`. Keep this case meshless in source so
there is one canonical copy of the QBOX9 project and pre-partitioned submeshes.

Acceptance uses the same element-integrated momentum audit as the short
momentum case, but over 10000 steps with `ELWriteDiagnostics=No` and
`ELMomentumAuditFreq=100`.
