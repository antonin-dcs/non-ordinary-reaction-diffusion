# A le mérite d'exister mais pour l'instant absolument impraticable même à l'ordre 1 ...

# --- 0. MODULES ---

import csv
import numpy as np
import matplotlib.pyplot as plt
# import matplotlib.animation as animation
from pathlib import Path

# --- 1. PARAMÈTRES DU MODÈLE ---

L = 400.0          # Taille du domaine (en mètres)
T = 20.            # Temps total de simulation (en jours)
Nt = 481           # Nombres total d'itérations temporelles (conforme à la simulation)
dt = T/Nt

sigma = 1.0        # Coefficient de diffusion (volatilité du brownien)
D = (sigma**2)/2   # Coefficient de diffusion effectif dans l'EDP


nu = 0.1           # Taux de mort naturel
gamma = 10.0       # Efficacité des pièges
R_T = 0.1          # Rayon d'action des pièges
x_T,y_T = 10,40    # Position du piège (unique pour l'instant)
lambda_ = 1

N0 = 10000
x0,y0 = 0,0

# Condition CFL pour la stabilité (pas nécessaire ici car on ne résout pas numériquement une EDP)

dx = 100 * (4*D*dt)**.5 # dx >= (4*D*dt)**.5
Nx = int(L/dx) # Nombre de pas d'espace


print(f"Simulation sur grille {Nx}x{Nx}. Pas de temps dt={dt:.5f}, Étapes={Nt}")


# --- 2. INITIALISATION ---


# Création de l'espace
x = np.linspace(-L, L, Nx)
y = np.linspace(-L, L, Nx)
X, Y = np.meshgrid(x, y)

# Récupération des coordonnées des pièges
data_path_coordinates = Path(__file__).parent / "data" / "trap_coordinates(in).csv"
trap_coordinates = []

with data_path_coordinates.open(newline="") as f:
    reader = csv.DictReader(f)
    for row in reader:
        try:
            trap_coordinates.append({
                "id": int(row["id"]),
                "x": float(row["x"]),
                "y": float(row["y"]),
            })
        except (KeyError, ValueError):
            # ignore malformed rows
            continue
    
trap_coordinates.sort(key=lambda t: t["id"])
ids = [t["id"] for t in trap_coordinates]
xs = [t["x"] for t in trap_coordinates]
ys = [t["y"] for t in trap_coordinates]



# --- 3. RESOLUTION ---



def rho_0(x, y, t):
    """
    Solution libre (ordre 0) : diffusion 2D à partir d'une source ponctuelle.
    """
    r2 = (x - x0)**2 + (y - y0)**2
    prefactor = N0 / (4.0 * np.pi * D * t)
    return prefactor * np.exp(-r2 / (4.0 * D * t))


def T_trap(x, y):
    """
    Potentiel de piège lisse (gaussien).
    """
    r2 = (x - x_T)**2 + (y - y_T)**2
    return np.exp(-r2 / R_T**2)


def rho(x, y, t):
    """
    Approximation au premier ordre de Duhamel.

    Paramètres
    ----------
    x, y : float
        Coordonnées spatiales.
    t : float
        Temps.
    n_time : int
        Nombre de points pour l'intégration temporelle.
    n_space : int
        Nombre de points par dimension pour l'intégration spatiale.
    L : float
        Demi-taille du domaine d'intégration spatiale autour du piège.

    Retour
    ------
    rho : float
        Densité approchée au point (x, y) et à l'instant t.
    """
    # terme d'ordre 0
    rho_free = rho_0(x, y, t)

    # grilles d'intégration
    s_vals = np.linspace(0.0, t, Nt)
    xs = np.linspace(x_T - L, x_T + L, Nx)
    ys = np.linspace(y_T - L, y_T + L, Nx)
    dx = xs[1] - xs[0]
    dy = ys[1] - ys[0]
    ds = s_vals[1] - s_vals[0]

    correction = 0.0

    for s in s_vals[1:]:
        tau = t - s
        for xi in xs:
            for yj in ys:
                G = np.exp(-((x - xi)**2 + (y - yj)**2) / (4.0 * D * tau)) \
                    / (4.0 * np.pi * D * tau)

                correction += (
                    G
                    * T_trap(xi, yj)
                    * rho_0(xi, yj, s)
                )

    correction *= dx * dy * ds

    return rho_free - lambda_ * correction






# --- 4. VISUALISATION PLANE ---


# Temps auxquels on évalue la solution
times = [1.0, 5.0, 10.0, 20.0]  # jours

# Création de la figure
fig, axes = plt.subplots(2, 2, figsize=(12, 10))
axes = axes.flatten()

for idx, t in enumerate(times):
    # Calcul de la densité sur toute la grille
    Z = rho(X, Y, t)
    
    # Affichage
    ax = axes[idx]
    im = ax.imshow(Z, extent=[-L, L, -L, L], origin='lower', 
                   cmap='hot', interpolation='bilinear')
    ax.set_title(f'Densité à t = {t:.1f} jours')
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')
    ax.scatter([x0], [y0], c='blue', marker='x', s=30, label='Source initiale')
    ax.legend()
    
    # Barre de couleur
    plt.colorbar(im, ax=ax, label='Densité (moustiques/m$^2$)')

plt.tight_layout()
plt.savefig('figures/analytic_solution_evolution.png', dpi=150)
print("Figure sauvegardée: figures/analytic_solution_evolution.png")
plt.show()

# --- 5. PROFIL RADIAL ---

# Calcul du profil radial à différents temps
r_values = np.linspace(0.1, L/4, 13)
plt.figure(figsize=(10, 6))

for t in times:
    densities = [rho(r, 0, t) for r in r_values]
    plt.plot(r_values, densities, label=f't = {int(t)} jour')

plt.xlabel('Distance radiale r (m)')
plt.ylabel(r'Densité $\rho(r, t)$')
plt.title('Profil radial de la densité de moustiques')
plt.legend()
plt.grid(True, alpha=0.3)
plt.savefig('figures/analytic_solution_radial.png', dpi=150)
print("Figure sauvegardée: figures/analytic_solution_radial.png")
plt.show()

# --- 6. MASSE TOTALE ---

# Vérification de la conservation de la masse (avec mortalité)
time_array = np.linspace(0, T, 10)
masses = []

for t in time_array:
    # Masse totale théorique avec mortalité
    mass = N0 * np.exp(-nu * t)
    masses.append(mass)

plt.figure(figsize=(10, 5))
plt.plot(time_array, masses, 'b-', linewidth=2)
plt.xlabel('$t$ (jours)')
plt.ylabel('$N(t)$ (nombre total de moustiques)')
plt.title(f'Évolution de la masse totale ($N_0$ = {N0}, $\\nu$ = {nu})')
plt.grid(True, alpha=0.3)
plt.savefig('figures/analytic_solution_mass.png', dpi=150)
print("Figure sauvegardée: figures/analytic_solution_mass.png")
plt.show()

print(f"\nMasse initiale : {N0}")
print(f"Masse finale (t={T}) : {N0 * np.exp(-nu * T):.2f}")
print(f"Réduction : {(1 - np.exp(-nu * T)) * 100:.1f}%")


