import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D # Nécessaire pour la 3D
import MDF.MDF_function as MDF_function 
import time

# --- 1. CHARGEMENT ET CONFIGURATION ---
df_real = pd.read_csv('data_csv/resultats_capture_réelles.csv')
y_real = df_real['Total_Reel'].values

# Paramètres fixes (Scénario SIT)
N = 100
initial_totalmosq = 50000
x0, y0 = 0, 0
R_0 = 1

# On fixe sigma et nu (soit aux valeurs priors, soit aux valeurs optimisées précédemment)
sigma_fixed = 17.0 
R_fixed=3.42

# alpha_reg n'est pas nécessaire ici car on veut voir l'erreur brute (MSE)

# --- 2. DÉFINITION DE LA GRILLE (GRID SEARCH) ---
# Attention : La simulation est une "boîte noire", donc chaque point prend du temps.
# On choisit une résolution raisonnable (ex: 20x20 = 400 simulations).

# Plage de test pour Gamma (Efficacité)
gamma_vals = np.linspace(0.1, 1.5, 20) 

# Plage de test pour R_trap (Rayon)
# Important : On teste R_trap indépendamment de gamma ici pour voir le paysage complet
R_trap_vals = np.linspace(1.0, 10.0, 20)

# Création des matrices pour stocker les résultats
Gamma_grid, R_grid = np.meshgrid(gamma_vals, R_trap_vals)
Cost_grid = np.zeros_like(Gamma_grid)

print(f"Lancement du balayage de l'espace des paramètres ({len(gamma_vals)}x{len(R_trap_vals)} simulations)...")
start_time = time.time()

# --- 3. BOUCLE DE CALCUL (BOÎTE NOIRE) ---
for i in range(len(R_trap_vals)):
    for j in range(len(gamma_vals)):
        
        # Paramètres courants
        r_current = R_trap_vals[i]
        g_current = gamma_vals[j]
        
        # Exécution de la simulation
        # On passe sigma_fixed et nu_fixed, mais on fait varier g et r
        try:
            _, _, _, daily_capture_matrix = MDF_function.simulation_MDF(
                sigma_fixed, nu_fixed, g_current, r_current, N, initial_totalmosq, x0, y0, R_0
            )
            
            # Calcul de l'erreur (MSE uniquement, pour voir la physique pure)
            y_sim = np.sum(daily_capture_matrix, axis=1)
            mse = np.mean((y_sim - y_real)**2)
            
            # On stocke le coût. On peut appliquer un log si les valeurs varient trop.
            Cost_grid[i, j] = mse 
            
        except Exception as e:
            Cost_grid[i, j] = np.nan # Gérer les erreurs de simu

end_time = time.time()
print(f"Calcul terminé en {end_time - start_time:.2f} secondes.")

# --- 4. VISUALISATION 3D DU PAYSAGE ---
fig = plt.figure(figsize=(14, 6))

# -- Plot 1 : Surface 3D --
ax1 = fig.add_subplot(1, 2, 1, projection='3d')
surf = ax1.plot_surface(Gamma_grid, R_grid, Cost_grid, cmap='viridis', edgecolor='none', alpha=0.9)
ax1.set_xlabel('Gamma (Efficacité)')
ax1.set_ylabel('R_trap (Rayon)')
ax1.set_zlabel('Coût (MSE)')
ax1.set_title('Paysage de la Fonction de Coût')
fig.colorbar(surf, ax=ax1, shrink=0.5, aspect=5)

# -- Plot 2 : Carte de chaleur (Contour Plot) --
# C'est souvent plus lisible pour voir les vallées d'identifiabilité
ax2 = fig.add_subplot(1, 2, 2)
# On utilise une échelle log pour les couleurs si les différences sont énormes
contour = ax2.contourf(Gamma_grid, R_grid, Cost_grid, levels=20, cmap='inferno')
ax2.set_xlabel('Gamma (Efficacité)')
ax2.set_ylabel('R_trap (Rayon)')
ax2.set_title('Vue de dessus (Isolignes d\'erreur)')
fig.colorbar(contour, ax=ax2)

# Ajout d'une ligne théorique d'iso-efficacité (Hypothèse gamma * R^2 ~ constante)
# Si vous voyez la vallée suivre cette forme, c'est la preuve du problème.
# Exemple hypothétique pour illustration (à adapter selon vos observations)
# k = 10 # Constante arbitraire pour l'exemple
# ax2.plot(gamma_vals, np.sqrt(k/gamma_vals), 'w--', label='Courbe iso-capture théorique')

plt.tight_layout()
plt.show()

# --- 5. ANALYSE RAPIDE ---
# Recherche du minimum sur la grille
min_idx = np.unravel_index(np.argmin(Cost_grid, axis=None), Cost_grid.shape)
best_R = R_grid[min_idx]
best_G = Gamma_grid[min_idx]
print(f"\nMinimum trouvé sur la grille :")
print(f"R_trap = {best_R:.2f}, Gamma = {best_G:.2f}, MSE = {Cost_grid[min_idx]:.2f}")