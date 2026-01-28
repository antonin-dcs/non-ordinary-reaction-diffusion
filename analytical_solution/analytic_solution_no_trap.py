# --- 0. MODULES ---

import numpy as np
import matplotlib.pyplot as plt
# import matplotlib.animation as animation


# --- 1. PARAMÈTRES DU MODÈLE ---

L = 400.0          # Taille du domaine (en mètres)
T = 20.            # Temps total de simulation (en jours)
Nt = 481           # Nombres total d'itérations temporelles (conforme à la simulation)
dt = T/Nt

sigma = 10.0        # Coefficient de diffusion (volatilité du brownien)
D = (sigma**2)/2   # Coefficient de diffusion effectif dans l'EDP
nu = 0.1           # Taux de mort naturel
(x0,y0) = (0. ,0.) # Origine de la source
N0 = 10000         # Nombre total de moustiques


print(f"Simulation sur grille {L}x{L}. Pas de temps dt={dt:.5f}, Étapes={Nt}")



# --- 2. SOLUTION ---

def rho(x, y, t:float):
    """
    Solution analytique de l'équation de diffusion avec mortalité constante
    en 2D pour une source ponctuelle initiale.

    Paramètres
    ----------
    x, y : float
        Coordonnées spatiales.
    t : float
        Temps (t > 0).

    Retour
    ------
    rho : float
        Densité au point (x, y) et à l'instant t.
    """
    r2 = (x - x0)**2 + (y - y0)**2
    prefactor = N0 / (4.0 * np.pi * D * t)
    return prefactor * np.exp(-r2 / (4.0 * D * t)) * np.exp(-nu * t)



# --- 3. VISUALISATION PLANE ---

# Création de la grille spatiale
nx = 200  # Nombre de points en x
ny = 200  # Nombre de points en y
x = np.linspace(-L, L, nx)
y = np.linspace(-L, L, ny)
X, Y = np.meshgrid(x, y)

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
    plt.colorbar(im, ax=ax, label='Densité (moustiques/m²)')

plt.tight_layout()
plt.savefig('figures/analytic_solution_evolution.png', dpi=150)
print("Figure sauvegardée: figures/analytic_solution_evolution.png")
plt.show()

# --- 4. PROFIL RADIAL ---

# Calcul du profil radial à différents temps
r_values = np.linspace(0.1, L/4, 200)
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

# --- 5. MASSE TOTALE ---

# Vérification de la conservation de la masse (avec mortalité)
time_array = np.linspace(0, T, 100)
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


