import numpy as np
import math
import matplotlib.pyplot as plt
import pandas as pd
import MDF.MDF_function as MDF_function  # Ton fichier contenant la fonction
import time
#paramaters of Facundo experience
D=150 #m²/day diffusivity
pds=0.8 # probability of daily survival
pc0 = 0.8      # Prob. catch at 0m
pc5 = 0.1        # Prob. catch at 5m




# --- 1. PARAMÈTRES DU MODÈLE --- 
sigma = math.sqrt(2*150)         # Coefficient de diffusion microscopique

nu = - math.log(pds) # Taux de mort naturel (paramater of the exponential diffusion )
gamma = pc0 # Efficacité des pièges,
b = (math.log(pc0) - math.log(pc5)) / 5**2 #spatial decay
R_trap =   math.sqrt(1/b)    # Rayon d'action des pièges
R_0 = 1 #rayon de libération des moustiques 



# Condition initiale
x_0, y_0 = 0,0 

# Paramètres grille
N = 200
initial_total_mosquitoes = 50000

# --- 2. EXÉCUTION DE LA SIMULATION ---
print("Lancement de la simulation via MDF...")
start=time.time()
h, capture_map, traps, daily_capture_matrix = MDF_function.simulation_MDF(
    sigma, nu, gamma, R_trap, N, initial_total_mosquitoes, x_0, y_0, R_0
)
end=time.time()
# --- 3. PRÉPARATION DES DONNÉES POUR L'HISTOGRAMME ---
# On somme la matrice temporelle (axe des jours) pour avoir le total par piège
# daily_capture_matrix est de forme (21, 20) -> on obtient un vecteur (21,)
total_captures_per_trap = np.sum(daily_capture_matrix, axis=1)
trap_ids = [t['id'] for t in traps] # Liste [1, 2, ..., 21]

# --- 4. VISUALISATION (3 GRAPHIQUES) ---
plt.figure(figsize=(18, 5)) # Figure plus large pour accueillir 3 plots

# --- PLOT 1 : Densité Finale ---
plt.subplot(1, 3, 1)
plt.title("Densité finale de moustiques")
im1 = plt.imshow(h, extent=[-400, 400, -400, 400], origin='lower', cmap='hot', vmin=0)
plt.colorbar(im1, label='Densité')

# Affichage des pièges avec NUMÉROS
x_traps = [t['q'][0] for t in traps]
y_traps = [t['q'][1] for t in traps]
plt.scatter(x_traps, y_traps, c='cyan', marker='x',label='Pièges', alpha=0.7)
plt.legend(loc='upper right')

# Boucle pour ajouter le texte (numéro) sur chaque piège
for t in traps:
    plt.text(t['q'][0], t['q'][1] + 0.3, str(t['id']), color='cyan', 
             fontsize=9, fontweight='bold', ha='center')

plt.scatter([x_0],[y_0], c='green', marker='o', label='Départ')
plt.legend(loc='upper right')

# --- PLOT 2 : Carte de Capture Spatiale ---
plt.subplot(1, 3, 2)
plt.title("Carte de capture cumulée (Spatial)")
im2 = plt.imshow(capture_map, extent=[-400, 400, -400,400], origin='lower', cmap='Blues', vmin=0)
plt.colorbar(im2, label='Total capturé')
plt.scatter(x_traps, y_traps, c='red', marker='x',label='Pièges', s=30, alpha=0.6)
plt.legend(loc='upper right')

# Ajout des numéros aussi ici pour repérage
for t in traps:
    plt.text(t['q'][0], t['q'][1] + 0.3, str(t['id']), color='black', 
             fontsize=8, ha='center')

# --- PLOT 3 : Histogramme des captures ---
plt.subplot(1, 3, 3)
plt.title("Nombre de moustiques capturés par piège")
# Création des barres
bars = plt.bar(trap_ids, total_captures_per_trap, color='skyblue', edgecolor='black')

plt.xlabel("Numéro du piège (ID)")
plt.ylabel("Nombre de moustiques")
plt.xticks(trap_ids) # Affiche tous les entiers 1 à 21 en X
plt.grid(axis='y', linestyle='--', alpha=0.7)

# Optionnel : Ajouter le chiffre exact au-dessus de chaque barre
for bar in bars:
    height = bar.get_height()
    plt.text(bar.get_x() + bar.get_width()/2., height,
             f'{height:.1f}',
             ha='center', va='bottom', fontsize=7, rotation=90)

plt.tight_layout()
#plt.savefig("simulation_sitscenarioN200.pdf")


# --- 5. EXPORT CSV ET BILAN ---
# Création du DataFrame
df_res = pd.DataFrame(daily_capture_matrix, 
                      index=trap_ids,
                      columns=[f"Jour_{j+1}" for j in range(daily_capture_matrix.shape[1])])
df_res.index.name = "Piege_ID"
df_res['Total'] = df_res.sum(axis=1)


# Sauvegarde
df_res.to_csv("resultats_captures.csv")

total_captured_traps=df_res['Total'].sum()
print("\n--- Bilan Rapide ---")
print(f"Total capturé globalement : {df_res['Total'].sum():.2f}")
print(f"Restants sur la grille   : {np.sum(h):.2f}")
print(f"Morts naturels (approx)  : {initial_total_mosquitoes - total_captured_traps - np.sum(h):.2f}")

print("Top 3 des pièges les plus efficaces :")
print(df_res['Total'].nlargest(3))
print(f"Temps d'execution : {round((end-start)):.2f} secondes")
plt.show()


