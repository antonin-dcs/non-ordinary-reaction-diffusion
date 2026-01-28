import pandas as pd
import numpy as np
import matplotlib.pyplot as plt



# --- CHARGEMENT DES DONNÉES  ---
df_real = pd.read_csv('data_csv/resultats_capture_réelles.csv')
#print(df_real.head())
y_real = df_real['Total_Reel'].values  # Le vecteur cible (numpy array)

df_mdf = pd.read_csv('resultats_captures.csv')
y_mdf = df_mdf['Total'].values

#print(df_mdf.head())
# Plot comparaison
plt.figure(figsize=(10, 6))
indices = np.arange(len(y_real))
width = 0.35

plt.bar(indices - width/2, y_real, width, label='Modèle probabiliste (référence)', color='skyblue',edgecolor='black', alpha=0.7)
plt.bar(indices + width/2, y_mdf, width, label='Modèle déterministe macroscopique', color='forestgreen', alpha=0.7)

plt.xlabel('ID Piège')
plt.ylabel('Total Moustiques Capturés')
plt.title(f'Moustiques capturés par piège simulés par les modèles probabiliste ou déterministe ')
plt.legend()
plt.grid(True, alpha=0.3)
plt.savefig("comparaison_new_title ")
plt.show()



