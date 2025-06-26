import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D   # keeps the 3-D backend active

# --- helix parameters ----------------------------------------------------
R       = 25.0        # mm – radius (½ of d_ct)
P       = 7.0         # mm – pitch
N_coil  = 45          # total turns
phi0    = 0.0         # rad – start phase
handed  = "right"     # "left" would flip the sign of z

# --- parametric curve ----------------------------------------------------
θ = np.linspace(0, 2*np.pi*N_coil, 2500)
x = R * np.cos(θ + phi0)
y = R * np.sin(θ + phi0)
z =  P * θ / (2*np.pi)
if handed.lower().startswith("l"):
    z = -z            # flip for left-handed helix

# --- plot ---------------------------------------------------------------
fig = plt.figure(figsize=(8, 6))
ax  = fig.add_subplot(111, projection="3d")
ax.plot(x, y, z, label="ATC centre-line helix")
ax.set_xlabel("x [mm]"); ax.set_ylabel("y [mm]"); ax.set_zlabel("z [mm]")
ax.set_title("Analytical helix: R = 25 mm, P = 7 mm, N = 45, right-handed")
ax.legend()
ax.set_box_aspect([1, 1, (z.max()-z.min())/(2*R)])  # keep aspect ratio sensible
plt.tight_layout(); plt.show()
