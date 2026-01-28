import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import MDF_function as MDF_function  # Votre module
import time 

# --- CHARGEMENT DES DONNÉES RÉELLES ---
# On suppose que le CSV a une colonne 'Piege_ID' et une colonne 'Total'
# Il faut s'assurer que l'ordre des pièges correspond à celui de la simulation
df_real = pd.read_csv('data_csv/resultats_capture_réelles.csv')

#print(df_real.head())
# Si nécessaire, trier par ID pour être sûr que ça matche la simulation
# (Adaptez 'Piege_ID' et 'Total' aux vrais noms de vos colonnes)
# df_real = df_real.sort_values('Piege_ID') 
y_real = df_real['Total_Reel'].values  # Le vecteur cible (numpy array)
#valeurs du problème sit scenario
N=100
initial_totalmosq=50000
x0,y0=0,0



#R_trap=3.46

# Valeurs initiales (Priors)
sigma_init = 17
gamma_init = 0.8
nu_init=.2

R_0 = 1 #rayon de libération des moustiques

initial_params = np.array([sigma_init, gamma_init,nu_init])

# --- DÉFINITION DE LA FONCTION DE COÛT ---
def objective_function(params, *args):
    """
    Cette fonction calcule l'erreur entre la simulation et la réalité.
    C'est elle que l'optimiseur va essayer de minimiser (descendre le gradient).
    """
    # 1. Déballer les paramètres courants proposés par l'optimiseur
    sigma, gamma,nu_init = params
    R_trap = np.sqrt(25/np.log(gamma/0.1) ) #relation fonctionnelle entre R_trap et gamma pour optimiser un seul paramutres 
    
    # Arguments fixes passés via *args
    N, init_mosq, x0, y0, y_target, init_params_ref, alpha_reg = args

    # 2. Exécuter la simulation avec ces paramètres
    try:
        _, _, _, daily_capture_matrix = MDF_function.simulation_MDF(
            sigma, nu_init, gamma, R_trap,N, init_mosq, x0, y0, R_0
        )
        
        # 3. Calculer le résultat simulé (Total par piège)
        y_sim = np.sum(daily_capture_matrix, axis=1)
        
        # 4. Calculer l'erreur quadratique (MSE - Mean Squared Error)
        # C'est la distance entre la courbe réelle et simulée
        mse_loss = np.mean((y_sim - y_target)**2)
        
        # 5. Calculer la Régularisation (Pénalité L2)
        # On pénalise si on s'éloigne trop des valeurs initiales (priors)
        # Formule : alpha * somme((param - init)^2)
        reg_loss = alpha_reg * np.sum((params - init_params_ref)**2)
        
        # Coût total
        total_loss = mse_loss + reg_loss
        
        # (Optionnel) Print pour suivre l'évolution en temps réel
        print(f"Params: {np.round(params, 3)} | MSE: {mse_loss:.2f} | Reg: {reg_loss:.4f} | Total: {total_loss:.2f}")
        
        return total_loss

    except Exception as e:
        print(f"Erreur simulation: {e}")
        return 1e9  # Retourner une erreur géante si la simu plante

# --- PARAMÈTRES DE L'OPTIMISATION ---

# 1. Hyperparamètre de régularisation (alpha)
# Plus il est grand, plus on force les paramètres à rester proches de l'origine.
# Plus il est petit, plus on laisse l'optimiseur libre de trouver n'importe quelle valeur.
alpha = 10  # À ajuster selon l'échelle de votre erreur MSE

# 2. Bornes (Bounds) pour empêcher des valeurs physiques impossibles
# sigma > 0, nu >= 0, gamma > 0, R_trap > 0
bounds = [
    (15, 20),   # sigma (diffusion)
    (0.1, 1),   # gamma (efficacité des pièges)
    (10**(-4),1) #nu
    
]

# Arguments supplémentaires à passer à la fonction objective
args = (N, initial_totalmosq, x0, y0, y_real, initial_params, alpha)

print("Début de l'optimisation...")
start=time.time()

# --- LANCEMENT DE L'OPTIMISATION (L-BFGS-B) ---
result = minimize(
    objective_function,       # La fonction à minimiser
    initial_params,           # Point de départ
    args=args,                # Les données fixes
    method='L-BFGS-B',        # Méthode quasi-Newton (gère les bornes)
    bounds=bounds,            # Les limites min/max des variables
    options={'disp': True, 'maxiter':40} # limite à 40 itérations 
)
end=time.time()
print(f"Temps d'execution : {round((end-start)):.2f} secondes")
# --- RÉSULTATS ---
print("\n--- OPTIMISATION TERMINÉE ---")
print(f"Nombre d'itérations : {result.nit}")
print(f"Succès : {result.success}")
print(f"Message : {result.message}")
print(f"Paramètres optimisés : {result.x}")

# Récupération des meilleurs paramètres
sigma_opt, gamma_opt, nu_opt = result.x
R_trap_opt = np.sqrt(25/np.log(gamma_opt/0.1) )
# --- VÉRIFICATION VISUELLE ---
# On relance une dernière fois avec les meilleurs paramètres pour tracer
_, _, _, best_matrix = MDF_function.simulation_MDF(
    sigma_opt, nu_opt, gamma_opt, R_trap_opt, N, initial_totalmosq, x0, y0, R_0
)
y_best_sim = np.sum(best_matrix, axis=1)

# Plot comparaison
plt.figure(figsize=(10, 6))
indices = np.arange(len(y_real))
width = 0.35

plt.bar(indices - width/2, y_real, width, label='Données Réelles', color='skyblue',edgecolor='black', alpha=0.7)
plt.bar(indices + width/2, y_best_sim, width, label='Simulation Optimisée', color='orange', alpha=0.7)

plt.xlabel('ID Piège')
plt.ylabel('Total Moustiques Capturés')
plt.title(f'Comparaison après optimisation\nD={(sigma_opt**2)/2:.2f}, γ={gamma_opt:.2f}, R={R_trap_opt:.2f},nu={nu_opt:.2f}')
plt.legend()
plt.grid(True, alpha=0.3)
plt.show()
