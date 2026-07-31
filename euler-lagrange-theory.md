
## 1. Eulerian fluid phase with porosity

Let

[
\varepsilon_f(\mathbf x,t)
]

be the local fluid volume fraction, also called porosity, and

[
\alpha_p(\mathbf x,t)=1-\varepsilon_f(\mathbf x,t)
]

the local particle volume fraction.

For an incompressible carrier liquid, a volume-averaged fluid model can be written as

[
\frac{\partial \varepsilon_f}{\partial t}
+
\nabla \cdot (\varepsilon_f \mathbf u_f)
========================================

0,
]

or, if the liquid density is constant,

[
\nabla \cdot (\varepsilon_f \mathbf u_f)
========================================

-\frac{\partial \varepsilon_f}{\partial t}.
]

In many simpler dilute or moderately dense Euler–Lagrange implementations, one uses the approximation

[
\nabla \cdot \mathbf u_f = 0
]

and treats (\varepsilon_f) mainly as an input to drag and source-term corrections. But for a dense suspension model, the porosity-aware form is usually the more consistent starting point.

The fluid momentum equation can be written as

[
\frac{\partial}{\partial t}
\left(
\varepsilon_f \rho_f \mathbf u_f
\right)
+
\nabla \cdot
\left(
\varepsilon_f \rho_f \mathbf u_f \otimes \mathbf u_f
\right)
=

-\varepsilon_f \nabla p
+
\nabla \cdot
\left(
\varepsilon_f \boldsymbol \tau_f
\right)
+
\varepsilon_f \rho_f \mathbf g
+
\mathbf f_{p \to f}.
]

Here

[
\boldsymbol \tau_f
==================

\mu_f
\left(
\nabla \mathbf u_f + \nabla \mathbf u_f^T
\right)
]

for a Newtonian carrier liquid. If you want the Eulerian phase itself to use an effective suspension viscosity, then this would become something like

[
\boldsymbol \tau_{\mathrm{eff}}
===============================

\mu_{\mathrm{eff}}(\alpha_p,\dot\gamma,\ldots)
\left(
\nabla \mathbf u_f + \nabla \mathbf u_f^T
\right),
]

where (\mu_{\mathrm{eff}}) could be informed by your PR-DNS/numerical-viscometry data.

The term

[
\mathbf f_{p \to f}(\mathbf x,t)
]

is the crucial two-way coupling source: it is the force per unit volume exerted by the particles onto the fluid.

## 2. Lagrangian particle equations

Each particle (i) has position, velocity, angular velocity, mass, and inertia:

[
\mathbf X_i(t), \qquad
\mathbf U_i(t), \qquad
\boldsymbol\omega_i(t), \qquad
m_i, \qquad
\mathbf I_i.
]

The translational equation is

[
m_i \frac{d\mathbf U_i}{dt}
===========================

\mathbf F_i^{h}
+
\mathbf F_i^{g}
+
\mathbf F_i^{c}
+
\mathbf F_i^{lub}
+
\mathbf F_i^{other}.
]

A common decomposition would be

[
\mathbf F_i^{h}
===============

\mathbf F_i^{drag}
+
\mathbf F_i^{lift}
+
\mathbf F_i^{pressure}
+
\mathbf F_i^{added}
+
\mathbf F_i^{history}
+
\ldots
]

In many practical first implementations, you might begin with

[
m_i \frac{d\mathbf U_i}{dt}
===========================

\mathbf F_i^{drag}
+
m_i\mathbf g
------------

\rho_f V_i\mathbf g
+
\mathbf F_i^{c}
+
\mathbf F_i^{lub}.
]

The position update remains

[
\frac{d\mathbf X_i}{dt}
=======================

\mathbf U_i.
]

For rotation,

[
\mathbf I_i \frac{d\boldsymbol\omega_i}{dt}
+
\boldsymbol\omega_i \times
\left(
\mathbf I_i \boldsymbol\omega_i
\right)
=

\mathbf T_i^{h}
+
\mathbf T_i^{c}
+
\mathbf T_i^{lub}.
]

For spherical particles, the simplest hydrodynamic drag closure could be

[
\mathbf F_i^{drag}
==================

\frac{1}{2}
\rho_f C_D A_i
\left|
\mathbf u_f(\mathbf X_i)-\mathbf U_i
\right|
\left(
\mathbf u_f(\mathbf X_i)-\mathbf U_i
\right),
]

or, in the Stokes limit,

[
\mathbf F_i^{drag}
==================

3\pi \mu_f d_i
\left(
\mathbf u_f(\mathbf X_i)-\mathbf U_i
\right).
]

For denser unresolved suspensions, this should usually be corrected by porosity:

[
\mathbf F_i^{drag}
==================

3\pi \mu_f d_i
,\beta(\varepsilon_f,Re_p)
\left(
\mathbf u_f(\mathbf X_i)-\mathbf U_i
\right),
]

where (\beta) is a hindered-drag or dense-suspension correction. The exact closure is one of the places where your PR-DNS data could be valuable.

## 3. Particle-to-grid force spreading

For two-way coupling, Newton’s third law requires that the reaction of the force exerted by the fluid on particle (i) is applied back to the fluid.

If

[
\mathbf F_i^{h}
]

is the hydrodynamic force exerted by the fluid on the particle, then the force exerted by the particle on the fluid is

[
-\mathbf F_i^{h}.
]

This is distributed to the Eulerian mesh using a compact kernel (W_h):

[
\mathbf f_{p \to f}(\mathbf x,t)
================================

-\sum_{i=1}^{N_p}
\mathbf F_i^{h}
,W_h(\mathbf x-\mathbf X_i).
]

The kernel should satisfy the conservation condition

[
\int_{\Omega} W_h(\mathbf x-\mathbf X_i),d\mathbf x = 1.
]

Then the total momentum transferred to the fluid is exactly

[
\int_{\Omega}
\mathbf f_{p \to f}(\mathbf x,t)
,d\mathbf x
===========

-\sum_{i=1}^{N_p}
\mathbf F_i^{h}.
]

On a finite-volume mesh, this becomes

[
\mathbf f_{p \to f,K}
=====================

-\frac{1}{|K|}
\sum_{i=1}^{N_p}
w_{iK}\mathbf F_i^{h},
]

where (K) is a cell, (|K|) is the cell volume, and (w_{iK}) is the fraction of particle (i)'s force assigned to cell (K). Conservation requires

[
\sum_K w_{iK}=1.
]

For FEM, the equivalent nodal/source assembly would be

[
\mathbf b_a^{p}
===============

-\sum_{i=1}^{N_p}
\mathbf F_i^{h}
N_a(\mathbf X_i),
]

for a simple point-evaluation coupling, where (N_a) is the finite-element shape function associated with degree of freedom (a). A smoother version uses kernel-weighted projection:

[
\mathbf b_a^{p}
===============

-\sum_{i=1}^{N_p}
\mathbf F_i^{h}
\int_{\Omega}
N_a(\mathbf x)
W_h(\mathbf x-\mathbf X_i)
,d\mathbf x.
]

This is the FEM analogue of conservative force spreading.

## 4. Fluid-to-particle interpolation

The particle force closure needs the local fluid velocity at the particle position:

[
\mathbf u_{f,i}
===============

\mathbf u_f(\mathbf X_i).
]

In FEM form,

[
\mathbf u_{f,i}
===============

\sum_a
N_a(\mathbf X_i)\mathbf u_a.
]

For a kernel-smoothed version,

[
\mathbf u_{f,i}
===============

\int_{\Omega}
\mathbf u_f(\mathbf x)
W_h(\mathbf x-\mathbf X_i)
,d\mathbf x.
]

Similarly, you may need

[
\nabla \mathbf u_{f,i},
\qquad
p_i,
\qquad
\nabla p_i,
\qquad
\varepsilon_{f,i}.
]

For example,

[
\varepsilon_{f,i}
=================

\varepsilon_f(\mathbf X_i).
]

This gives the coupled loop:

[
\mathbf u_f, p, \varepsilon_f
\quad \longrightarrow \quad
\mathbf F_i^{h}, \mathbf T_i^{h}
\quad \longrightarrow \quad
\mathbf f_{p\to f}
\quad \longrightarrow \quad
\mathbf u_f, p.
]

That is the mathematical core of two-way Euler–Lagrange coupling.

## 5. Local void fraction / occupancy accumulation

For unresolved particles, the particle volume has to be represented statistically or through smoothed occupancy.

The particle volume fraction can be accumulated as

[
\alpha_p(\mathbf x,t)
=====================

\sum_{i=1}^{N_p}
V_i W_h(\mathbf x-\mathbf X_i),
]

and therefore

[
\varepsilon_f(\mathbf x,t)
==========================

1-
\alpha_p(\mathbf x,t).
]

On a cell (K),

[
\alpha_{p,K}
============

\frac{1}{|K|}
\sum_i
V_i w_{iK},
]

and

[
\varepsilon_{f,K}
=================

1-\alpha_{p,K}.
]

Again, the conservative condition is

[
\sum_K w_{iK}=1,
]

so that

[
\sum_K |K| \alpha_{p,K}
=======================

\sum_i V_i.
]

This is important. Otherwise the particles may inject the wrong displaced volume into the Eulerian model.

For dense Euler–Lagrange simulations, you may also need to limit or regularize the void fraction:

[
\varepsilon_f \geq \varepsilon_{\min},
]

because drag laws and pressure-gradient terms can become singular or unstable as (\varepsilon_f \to 0).

## 6. Torque coupling

If your Eulerian formulation only solves linear momentum, then most unresolved Euler–Lagrange solvers do not explicitly feed particle torque back into the fluid. They include torque only in the particle rotational equation.

But if you want angular momentum consistency, there are two main possibilities.

The first is to spread torque as a local couple source:

[
\mathbf c_{p\to f}(\mathbf x,t)
===============================

-\sum_i
\mathbf T_i^h
W_h(\mathbf x-\mathbf X_i).
]

This belongs naturally to a micropolar or Cosserat-type continuum formulation, where the fluid has an angular momentum equation.

The second possibility is to represent the torque through an equivalent force couple. That means constructing a local force distribution whose net force is zero but whose net moment is

[
\int_{\Omega}
(\mathbf x-\mathbf X_i)
\times
\mathbf f_i^{torque}(\mathbf x)
,d\mathbf x
===========

-\mathbf T_i^h.
]

For a standard Navier–Stokes solver, I would treat torque feedback as an advanced extension. The first two-way coupling milestone should be conservative linear momentum feedback.

## 7. Source-term assembly into Navier–Stokes

In weak FEM form, the fluid momentum equation with particle feedback would contain an additional right-hand-side term:

[
\int_{\Omega}
\mathbf f_{p\to f}
\cdot
\mathbf v_h
,d\mathbf x.
]

Using the force-spreading definition,

[
\int_{\Omega}
\mathbf f_{p\to f}
\cdot
\mathbf v_h
,d\mathbf x
===========

-\sum_i
\mathbf F_i^h
\cdot
\int_{\Omega}
W_h(\mathbf x-\mathbf X_i)
\mathbf v_h(\mathbf x)
,d\mathbf x.
]

With point-evaluation coupling,

[
\int_{\Omega}
\mathbf f_{p\to f}
\cdot
\mathbf v_h
,d\mathbf x
\approx
-\sum_i
\mathbf F_i^h
\cdot
\mathbf v_h(\mathbf X_i).
]

In nodal form, this contributes

[
\mathbf b_a^{p}
===============

-\sum_i
\mathbf F_i^h N_a(\mathbf X_i).
]

So the discrete fluid system changes from roughly

[
\mathbf M \dot{\mathbf u}
+
\mathbf N(\mathbf u)
+
\mathbf K \mathbf u
+
\mathbf G p
===========

\mathbf b
]

to

[
\mathbf M_\varepsilon \dot{\mathbf u}
+
\mathbf N_\varepsilon(\mathbf u)
+
\mathbf K_\varepsilon \mathbf u
+
\mathbf G_\varepsilon p
=======================

\mathbf b
+
\mathbf b^{p}.
]

Here the subscript (\varepsilon) indicates that mass, convection, viscosity, and pressure terms may be weighted by the fluid volume fraction (\varepsilon_f).

## 8. Practical minimal two-way coupled model

A good first version for implementation would be:

Fluid:

[
\rho_f
\left(
\frac{\partial \mathbf u_f}{\partial t}
+
\mathbf u_f \cdot \nabla \mathbf u_f
\right)
=

-\nabla p
+
\nabla \cdot \boldsymbol \tau_f
+
\rho_f \mathbf g
+
\mathbf f_{p\to f},
]

[
\nabla \cdot \mathbf u_f = 0.
]

Particles:

[
m_i \frac{d\mathbf U_i}{dt}
===========================

\mathbf F_i^{drag}
+
\mathbf F_i^{g/b}
+
\mathbf F_i^{c}
+
\mathbf F_i^{lub},
]

[
\frac{d\mathbf X_i}{dt}
=======================

\mathbf U_i.
]

Drag:

[
\mathbf F_i^{drag}
==================

3\pi\mu_f d_i
,\beta(\varepsilon_{f,i},Re_{p,i})
\left(
\mathbf u_f(\mathbf X_i)-\mathbf U_i
\right).
]

Feedback:

[
\mathbf f_{p\to f}(\mathbf x)
=============================

-\sum_i
\mathbf F_i^{drag}
W_h(\mathbf x-\mathbf X_i).
]

Occupancy:

[
\alpha_p(\mathbf x)
===================

\sum_i
V_i W_h(\mathbf x-\mathbf X_i),
\qquad
\varepsilon_f=1-\alpha_p.
]

This would already be a genuine **two-way coupled Euler–Lagrange model**.

A denser, more complete version would use

[
\mathbf f_{p\to f}
==================

-\sum_i
\left(
\mathbf F_i^{drag}
+
\mathbf F_i^{lift}
+
\mathbf F_i^{pressure}
+
\mathbf F_i^{added}
+
\mathbf F_i^{lub,h}
\right)
W_h(\mathbf x-\mathbf X_i),
]

but I would be careful with contact forces. Dry contact forces between particles are internal to the particle phase; they should not automatically be applied to the fluid. Hydrodynamic interaction forces should be fed back to the fluid. Contact forces should influence the particle motion, but they are not necessarily a direct fluid momentum source.

## 9. Relation to your DNS method

In your current particle-resolved DNS method, the force (\mathbf F_i^h) comes from the resolved stress integral,

[
\mathbf F_i^h
=============

-\int_{\partial \Omega_i}
\boldsymbol\sigma \cdot \mathbf n
,d\Gamma,
]

or, in FBM form, from a volume-integral equivalent. In unresolved Euler–Lagrange, this is replaced by a closure model:

[
\mathbf F_i^h
\approx
\mathcal F
\left(
\mathbf u_f(\mathbf X_i),
\mathbf U_i,
\nabla \mathbf u_f(\mathbf X_i),
p(\mathbf X_i),
\varepsilon_f(\mathbf X_i),
Re_p,
\alpha_p,
\ldots
\right).
]

That is the major modeling step.

Your DNS can provide the closure information for

[
\mathcal F,
\qquad
\beta(\varepsilon_f,Re_p),
\qquad
\mu_{\mathrm{eff}}(\alpha_p,\dot\gamma),
\qquad
\text{dispersion/stress closures},
]

while the Euler–Lagrange solver provides the scalable apparatus-level model.

So the concise mathematical description of the upgrade is:

[
\boxed{
\text{resolved particle boundary stresses}
\quad \longrightarrow \quad
\text{unresolved force closures plus conservative force spreading}
}
]

with

[
\boxed{
\mathbf f_{p\to f}(\mathbf x,t)
===============================

-\sum_i
\mathbf F_i^h
W_h(\mathbf x-\mathbf X_i)
}
]

and

[
\boxed{
\varepsilon_f(\mathbf x,t)
==========================

1-
\sum_i
V_i W_h(\mathbf x-\mathbf X_i)
}.
]

Those two boxes are essentially the mathematical heart of the two-way Euler–Lagrange extension.
