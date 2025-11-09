import numpy as np
import matplotlib.pyplot as plt

# ------------------------
# 1) Paramètres physiques
# ------------------------
alpha = 9e-5   # diffusivité thermique (m^2/s) pour l'aluminium (ordre de grandeur)
T_cast = 700.0   # température de la pièce juste coulée (°C)
T_mold = 250.0   # température imposée sur le bord (°C)

# ------------------------
# 2) Domaine et maillage
# ------------------------
Lx = 0.20  # longueur en x (m)
Ly = 0.10  # longueur en y (m)
Nx = 81    # nombre de points en x
Ny = 41    # nombre de points en y

dx = Lx / (Nx - 1)
dy = Ly / (Ny - 1)

# pour simplifier on impose dx = dy sinon il faut adapter la formule de stabilité
assert abs(dx - dy) < 1e-12, "Ici on suppose dx = dy pour le schéma explicite."

# ------------------------
# 3) Temps de simulation
# ------------------------
t_final = 6.0  # secondes de refroidissement à simuler
dt_stable = dx**2 / (4 * alpha)
dt = 0.2 * dt_stable   # on prend un peu en dessous pour être sûr
Nt = int(t_final / dt)

print(f"dx = {dx:.4e} m, dt max stable ~ {dt_stable:.4e} s, dt choisi = {dt:.4e} s, Nt = {Nt}")

# ------------------------
# 4) Initialisation du champ de température
# ------------------------
T = np.full((Ny, Nx), T_cast)

def impose_bords(T):
    """
    Imposer la température du moule sur les 4 côtés.
    Attention: on travaille avec T[j,i] = y,x
    """
    T[0, :]  = T_mold   # bord bas
    T[-1, :] = T_mold   # bord haut
    T[:, 0]  = T_mold   # bord gauche
    T[:, -1] = T_mold   # bord droit
    return T

T = impose_bords(T)

# ------------------------
# 5) Boucle en temps
# ------------------------
T_list = [T.copy()]  # pour stocker quelques états pour visualiser après
save_every = max(1, Nt // 5)  # on garde ~5 images

for n in range(Nt):
    Tn = T.copy()
    # noyau de diffusion 2D, sur les points intérieurs
    T[1:-1, 1:-1] = (
        Tn[1:-1, 1:-1]
        + alpha * dt * (
            (Tn[1:-1, 2:] - 2*Tn[1:-1, 1:-1] + Tn[1:-1, 0:-2]) / dx**2
            + (Tn[2:, 1:-1] - 2*Tn[1:-1, 1:-1] + Tn[0:-2, 1:-1]) / dy**2
        )
    )
    # réimposer les bords
    T = impose_bords(T)

    if (n+1) % save_every == 0:
        T_list.append(T.copy())

# ------------------------
# 6) Visualisation simple
# ------------------------
fig, axes = plt.subplots(1, len(T_list), figsize=(4*len(T_list), 3))
if len(T_list) == 1:
    axes = [axes]

for ax, Ti in zip(axes, T_list):
    im = ax.imshow(Ti, origin='lower', extent=[0, Lx, 0, Ly])
    ax.set_title("T (°C)")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    fig.colorbar(im, ax=ax)

plt.tight_layout()
plt.show()


from mpl_toolkits.mplot3d import Axes3D  # juste pour activer le 3D
X = np.linspace(0, Lx, Nx)
Y = np.linspace(0, Ly, Ny)
X, Y = np.meshgrid(X, Y)

fig = plt.figure(figsize=(5,4))
ax = fig.add_subplot(111, projection='3d')
surf = ax.plot_surface(X, Y, T, cmap='viridis', linewidth=0, antialiased=False)
fig.colorbar(surf, ax=ax, label="Température (°C)")
ax.set_xlabel("x (m)")
ax.set_ylabel("y (m)")
ax.set_zlabel("T (°C)")
ax.set_title("Température en surface 3D")
plt.show()

#%%

import numpy as np

# ------------------------
# 1) Données matériau
# ------------------------
alpha_diff = 7.4e-5     # diffusivité thermique (m^2/s)
k = 180.0               # conductivité W/mK
rho = 2700.0
cp = 900.0

E = 70e9                # Pa (70 GPa)
alpha_th = 2.3e-5       # dilatation thermique 1/K
T_ref = 250.0           # ref pour la contrainte

# ------------------------
# 2) Domaine
# ------------------------
Lx, Ly = 0.20, 0.10
Nx, Ny = 81, 41
dx = Lx / (Nx - 1)
dy = Ly / (Ny - 1)
assert abs(dx - dy) < 1e-12

# ------------------------
# 3) Temps
# ------------------------
t_final = 6.0
dt_stable = dx**2 / (4 * alpha_diff)
dt = 0.2 * dt_stable
Nt = int(t_final / dt)

# Paramètres convection
h = 200.0          # W/m2K (exagéré pour voir l'effet)
T_fluid = 250.0    # fluide/moule

def simulate(T_mold=None, use_convection=False):
    """
    Simule le refroidissement et renvoie le champ final.
    Si use_convection=True, on applique une loi de convection sur les bords.
    Sinon, Dirichlet T=T_mold.
    """
    T = np.full((Ny, Nx), 700.0)

    def impose_dirichlet(T):
        T[0, :]  = T_mold
        T[-1, :] = T_mold
        T[:, 0]  = T_mold
        T[:, -1] = T_mold
        return T

    # pré-calcul pour convection
    beta = h * dx / k  # terme sans dimension

    # initialisation bords
    if not use_convection:
        T = impose_dirichlet(T)

    for n in range(Nt):
        Tn = T.copy()
        T[1:-1, 1:-1] = (
            Tn[1:-1, 1:-1]
            + alpha_diff * dt * (
                (Tn[1:-1, 2:] - 2*Tn[1:-1, 1:-1] + Tn[1:-1, 0:-2]) / dx**2
                + (Tn[2:, 1:-1] - 2*Tn[1:-1, 1:-1] + Tn[0:-2, 1:-1]) / dy**2
            )
        )

        if use_convection:
            # bord bas (j=0)
            T[0, 1:-1] = (T[1, 1:-1] + beta * T_fluid) / (1 + beta)
            # bord haut
            T[-1, 1:-1] = (T[-2, 1:-1] + beta * T_fluid) / (1 + beta)
            # bord gauche
            T[1:-1, 0] = (T[1:-1, 1] + beta * T_fluid) / (1 + beta)
            # bord droit
            T[1:-1, -1] = (T[1:-1, -2] + beta * T_fluid) / (1 + beta)
            # coins (simple)
            T[0,0] = T_fluid
            T[0,-1] = T_fluid
            T[-1,0] = T_fluid
            T[-1,-1] = T_fluid
        else:
            T = impose_dirichlet(T)

    return T

def gradient_max(T):
    dTx = (T[:, 2:] - T[:, 0:-2]) / (2*dx)
    dTy = (T[2:, :] - T[0:-2, :]) / (2*dy)

    # remettre aux mêmes dimensions pour faire une norme
    gx = dTx[1:-1, :]
    gy = dTy[:, 1:-1]
    grad = np.sqrt(gx**2 + gy**2)
    return grad.max()

def sigma_th_max(T):
    # contrainte thermique "bloquée" (très simplifiée)
    sigma = E * alpha_th * (T - T_ref)
    return sigma.max()

# ------------------------
# 4) Étude paramétrique
# ------------------------
T_molds = [250.0, 275.0, 300.0]

results = []
for Tm in T_molds:
    T_final = simulate(T_mold=Tm, use_convection=False)
    gmax = gradient_max(T_final)
    smax = sigma_th_max(T_final)
    results.append((Tm, gmax, smax))

print("=== Dirichlet (T imposée) ===")
for Tm, gmax, smax in results:
    print(f"T_mold={Tm:.0f}°C -> grad max={gmax:.2f} °C/m, sigma_th_max={smax/1e6:.1f} MPa")

# test avec convection
T_final_conv = simulate(use_convection=True)
gmax_conv = gradient_max(T_final_conv)
smax_conv = sigma_th_max(T_final_conv)
print("\n=== Convection ===")
print(f"grad max={gmax_conv:.2f} °C/m, sigma_th_max={smax_conv/1e6:.1f} MPa")


from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

fig = plt.figure(figsize=(6,4))
ax = fig.add_subplot(111, projection='3d')

surf = ax.plot_surface(
    X, Y, T,
    cmap='viridis',
    vmin=T_mold, vmax=T_cast,
    linewidth=0, antialiased=False
)

ax.set_xlabel("x (m)")
ax.set_ylabel("y (m)")
ax.set_zlabel("T (°C)")
ax.set_title("Température – surface 3D")
fig.colorbar(surf, ax=ax, label="Température (°C)")
plt.show()
