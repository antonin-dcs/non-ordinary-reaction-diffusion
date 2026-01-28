import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# --- CHARGEMENT DES DONNÉES RÉELLES ---
# On suppose que le CSV a une colonne 'Piege_ID' et une colonne 'Total'
# Il faut s'assurer que l'ordre des pièges correspond à celui de la simulation
df_real = pd.read_csv('data_csv/resultats_capture_réelles.csv')

#print(df_real.head())
# Si nécessaire, trier par ID pour être sûr que ça matche la simulation
# (Adaptez 'Piege_ID' et 'Total' aux vrais noms de vos colonnes)
# df_real = df_real.sort_values('Piege_ID') 
y_real = df_real['Total_Reel'].values  # Le vecteur cible (numpy array)

# Plot comparaison
plt.figure(figsize=(10, 6))
indices = np.arange(len(y_real))
width = 0.35
trap_ids=[i for i in range(1,22)]
plt.xticks(trap_ids) # Affiche tous les entiers 1 à 21 en X

plt.bar(indices - width/2, y_real, width, label='Données Réelles', color='forestgreen',edgecolor='black', alpha=0.7)
#plt.bar(indices + width/2, y_best_sim, width, label='Simulation Optimisée', color='orange', alpha=0.7)

plt.xlabel('ID Piège')
plt.ylabel('Total Moustiques Capturés')
plt.title("Histogramme des moustiques capturés pour chaque piège d'après les données réelles")
plt.legend()
plt.grid(True, alpha=0.3)
plt.show()
