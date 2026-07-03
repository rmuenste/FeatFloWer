# Euler-Lagrange Semi-Implicit Drag Coupling

`SimPar@ELDragCoupling = semi_implicit` moves the linear drag reaction on the fluid
from a fully explicit RHS into the fine-level Q2 velocity operator. The default remains
`explicit`.

For particle-side drag represented by PE's subcycled effective drag impulse
`F_d,eff` and particle-side lift `F_l`, the target fluid reaction is:

```text
-F_d,eff - F_l
```

The semi-implicit path linearizes this reaction around the sampled fluid velocity:
it deposits `B u_f - F_d,eff - F_l` on the explicit fluid RHS and adds `B u` to the
implicit fluid matrix diagonal. The old-time Crank-Nicolson contribution is subtracted
from the velocity defect with the same `theta` split used by the rest of the velocity
operator. The drag field is stored as a distributed Q2-nodal `drag_B_source` snapshot
from the pre-advance EL pass and is summed once by the normal velocity matrix
synchronization.

`EL_FEEDBACK_CONSERVATION mode=semi_implicit` is a deposit-fidelity check for the
explicit source split, not a full Newton-pair assertion. Total momentum drift from
`EL_MOMENTUM_ELEMINT` is the acceptance metric for the coupled semi-implicit path.

Initial np=28 acceptance data from the 60-step momentum case:
`drift_rel = 1.510406e-3`, `EL_FEEDBACK_CONSERVATION residual = 3.388132e-21`,
zero void-fraction clipping, and average MG-UVW iterations per nonlinear step = 1.0.
The opt-in harness tolerance is `3.1e-3`, about 2x the measured drift.
