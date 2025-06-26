The **Flow Map** publication does *not* explicitly define quantitative metrics like the axial distribution, radial distribution index, or vertical asymmetry factor using formulas such as:

* $\overline{\phi}(z)$
* $I_r = 1 - \frac{\sum_i|\phi_i - \phi_{\text{bulk}}|}{2\phi_{\text{bulk}}}$
* $A_y = \frac{\phi_\text{upper} - \phi_\text{lower}}{\phi_\text{upper} + \phi_\text{lower}}$

However, it **does describe these effects qualitatively** using categorical zones (e.g., red, yellow, dotted green, green) based on the **Shield’s parameter** $\Theta$ and the **modified Dean number ratio** $Dn^*/Dn$ to distinguish different suspension behaviors. Specifically:

---

### 🌡 Suspension States in the Flow Map

The following **qualitative thresholds** are introduced for the Shield’s parameter $\Theta$:

* **Θ < 10** → Rear accumulation only (red)
* **10 ≤ Θ ≤ 30** → Horizontal spread, mostly in lower vortex (yellow)
* **30 ≤ Θ ≤ 180** → Vertical spread into lower vortex increases (dotted green)
* **Θ > 180** → Nearly uniform distribution (green)

These zones do correspond roughly to:

* **Axial distribution**: from rear-loaded slugs (Θ < 10) to more spread-out distributions.
* **Radial index**: related to whether particles pile on the walls (low $I_r$) or spread evenly (high $I_r$).
* **Vertical asymmetry**: becomes important in yellow → dotted green → green transitions.

The authors acknowledge that these zones are currently **qualitative and heuristic**, and that Θ and $Dn^*/Dn$ *do not account for solid content, particle shape, or slug length*. They suggest these limitations could be addressed in future work with **more quantitative descriptors**.

---

### 📏 Clarifying the Proposed Metrics

Here's a more detailed explanation of the **quantitative metrics** I suggested, which are meant to *augment* the flow map with DNS-level simulation data:

---

#### 1. **Axial Distribution** $\overline{\phi}(z)$

**Definition:**
Average particle volume fraction along the axial coordinate $z$ of the slug, integrated radially and azimuthally:

$$
\overline{\phi}(z) = \frac{1}{A_{\text{cross}}} \int_A \phi(x,y,z)\,dx\,dy
$$

**Use:**
A flat profile means particles are distributed evenly in the slug. A spike near the back indicates rear-loading.

---

#### 2. **Radial Distribution Index $I_r$**

**Definition:**
Measures deviation from uniform radial distribution in cylindrical coordinates:

$$
I_r = 1 - \frac{1}{2\phi_{\text{bulk}}} \sum_i \left|\phi_i - \phi_{\text{bulk}}\right|
$$

* $\phi_i$: average particle volume fraction in radial bin $i$
* $\phi_{\text{bulk}}$: expected value if distribution were uniform

**Interpretation:**

* $I_r = 1$: perfectly uniform
* $I_r \ll 1$: strong radial inhomogeneity (e.g., particles all near wall)

---

#### 3. **Vertical Asymmetry Factor $A_y$**

**Definition:**
Measures how asymmetric the particle distribution is with respect to the vertical mid-plane:

$$
A_y = \frac{\phi_\text{upper} - \phi_\text{lower}}{\phi_\text{upper} + \phi_\text{lower}}
$$

* $\phi_\text{upper}$: volume fraction in the top half of the tube
* $\phi_\text{lower}$: in the bottom half

**Interpretation:**

* $A_y = 0$: vertically symmetric
* $A_y < 0$: particles dominate bottom half (gravity)
* $A_y > 0$: particles migrate upward (Dean vortices strong)

---

### 🔎 Why this matters

While the **Flow Map paper** captures behavior using flow regime markers and Shield's thresholds, your **DNS or LPT simulation data can compute these metrics explicitly**, making it possible to:

* Quantify where a specific experiment or simulation lies within a color zone.
* Validate the flow map’s qualitative predictions against numerical metrics.
* Enhance the publication with precise numbers that capture particle dynamics.

---

