import pandas as pd
import matplotlib.pyplot as plt


# 1. Charger le fichier
df = pd.read_csv('trap_counts(in).csv')

# 2. Nettoyer les noms de colonnes (retire les espaces avant/après : 't      ' devient 't')
df.columns = df.columns.str.strip()
# Crée la liste des jours souhaités
jours_souhaites = list(range(1, 21))

# Filtre le DataFrame
df_filtre = df[df['t'].isin(jours_souhaites)]
resultat=df[df['t']==20]

df_reel_pivote = df_filtre.pivot(index='id', columns='t', values='n_cum')

# 3. Renommer les colonnes pour correspondre exactement au format "Jour_x"
df_reel_pivote.columns = [f"Jour_{int(c)}" for c in df_reel_pivote.columns]

# 4. Ajout de la colonne Total (comme dans la simulation)
df_reel_pivote['Total_Reel'] = df_reel_pivote.iloc[:, -1] # Pour du cumulé, le total est la dernière colonne

print(df_reel_pivote.head())
csv_filename = "resultats_capture_réelles.csv"
df_reel_pivote.to_csv(csv_filename)
print(f"\nFichier sauvegardé : {csv_filename}")

if __name__ == '__main__':
    #print(resultat)
    plt.title("Nombre de moustiques capturés par piège (données réelles)")
    plt.xlabel("Numéro du piège")
    plt.ylabel("Moustiques capturés")
    plt.xticks(resultat['id']) # Affiche tous les entiers en X
    plt.bar(resultat['id'], resultat['n_cum'],color='skyblue', edgecolor='black')
    plt.grid(axis='y', linestyle='--', alpha=0.7)

    plt.tight_layout()

    plt.savefig("captures_reelles.png")
    plt.show()
