import numpy as np
from scipy.spatial import Delaunay
from scipy.sparse import lil_matrix, diags
from scipy.sparse.linalg import factorized
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation

# --- 1. Paramètres de Simulation ---
target_time = 20.0
L_min, L_max = -400.0, 400.0
nx, ny = 180, 180   ###################################################################################################### précision du maillage
N_points = nx * ny
D = 150 #############################################################################  coefficient de diffusion


# --- 2. Configuration des Pièges et Terme Source ---
dt = 0.1
steps_to_simulate = int(target_time / dt)
natural_death_rate = 0.2 ############################################################################  taux de mortalité naturelle
trap_amplitude = 0.7  #############################################################################  efficacité des pièges
trap_sigma = 2.5 #############################################################################  étendue d'action des pièges

TRAP_POS = np.array([
    (-131.68, 131.68), (-124.68, 51.64), (-108.33, -44.87), (42.25, 0.0),
    (49.17, 118.72), (-37.34, 90.16), (12.76, -30.80), (0.0, -76.87),
    (73.57, -177.61), (0.0, -362.84), (-130.79, -315.75), (-86.98, -210.00),
    (-50.62, -122.21), (335.92, 0.0), (233.86, -233.86), (200.78, -200.78),
    (68.89, -166.31), (76.60, 0.0), (-26.42, 26.42), (-245.38, 245.38),
    (-85.74, 207.00)
])

# --- 3. Génération du Maillage ---
def repartition(x,y,radius):
    val = 0.15
    for i in range(len(TRAP_POS)):
        dist = np.linalg.norm(np.array([x,y]) - TRAP_POS[i])
        if dist < radius:
            val = 1.0
            break
        elif radius <= dist < 2*radius:
            val1 = 1 - 0.7*(dist - radius)/radius
            if val1 > val:
                val = val1
    return val

radius = 40
maxim_repartition = 1.0
points = []
ctn = 0
while len(points) < N_points:
    ctn += 1
    x_rand = np.random.uniform(L_min, L_max)
    y_rand = np.random.uniform(L_min, L_max)
    prob = repartition(x_rand, y_rand, radius) / maxim_repartition
    if np.random.rand() < prob:
        points.append([x_rand, y_rand])



np.random.seed(100)
points = np.random.rand(N_points, 2) * (L_max - L_min) + L_min
corners = np.array([[L_min, L_min], [L_max, L_min], [L_max, L_max], [L_min, L_max]])
points = np.vstack([points, corners])

delaunay = Delaunay(points)
triangles = delaunay.simplices
n_triangles = len(triangles)

# --- 4. Matrices Géométriques ---
def compute_matrices(delaunay, points):
    pts = points
    tris = delaunay.simplices
    centers = pts[tris].mean(axis=1)
    
    p1 = pts[tris[:, 0]]
    p2 = pts[tris[:, 1]]
    p3 = pts[tris[:, 2]]
    areas = 0.5 * np.abs(p1[:,0]*(p2[:,1]-p3[:,1]) + 
                         p2[:,0]*(p3[:,1]-p1[:,1]) + 
                         p3[:,0]*(p1[:,1]-p2[:,1]))
    
    A = lil_matrix((len(tris), len(tris)))
    neighbors = delaunay.neighbors
    for i in range(len(tris)):
        for j in range(3):
            neigh_idx = neighbors[i, j]
            if neigh_idx != -1:
                v1 = tris[i, (j+1)%3]
                v2 = tris[i, (j+2)%3]
                edge_len = np.linalg.norm(pts[v1] - pts[v2])
                dist = np.linalg.norm(centers[i] - centers[neigh_idx])
                trans = edge_len / dist
                A[i, neigh_idx] += trans
                A[i, i] -= trans
    return A.tocsr(), centers, areas

A_diff, centers, areas = compute_matrices(delaunay, points)
A_diff = D * A_diff



# --- NOUVEAU : Pré-calcul des profils individuels des pièges ---
# trap_profiles[j, i] contient l'efficacité du piège j sur le triangle i
trap_profiles = np.zeros((len(TRAP_POS), n_triangles))
R_total = np.full(n_triangles, natural_death_rate) # Le vecteur total pour le solveur

for j, trap_loc in enumerate(TRAP_POS):
    # Distance au carré entre ce piège et tous les triangles
    dists_sq = np.sum((centers - trap_loc)**2, axis=1)
    # Gaussienne pour ce piège
    gaussian = trap_amplitude * np.exp(-dists_sq / (2 * trap_sigma**2))

    trap_profiles[j, :] = gaussian
    R_total += gaussian # On ajoute au total pour la matrice système

# Création du solveur
V = diags(areas)
M_reaction = diags(R_total * areas)
System_Matrix = (V - dt * A_diff + dt * M_reaction).tocsc()
solve_step = factorized(System_Matrix)

# --- 5. État Initial et Variables de Suivi ---
u = np.zeros(n_triangles)
for i, center in enumerate(centers):
    if np.linalg.norm(center) <= 5.0: ############################################################# initialisation circulaire autour de l'origine
        u[i] = 450.0 ############################################################################################ densité initiale proche du centre

# Initialisation des compteurs
initial_population = np.sum(u * areas)
trapped_counts = np.zeros(len(TRAP_POS)) # Cumul pour chaque piège
history_total_pop = []

# --- 6. Simulation avec Comptage ---
print(f"Population initiale : {initial_population:.2f}")

for step in range(steps_to_simulate):
    # 1. Calcul de la dynamique (solveur)
    rhs = V.dot(u)
    u_new = solve_step(rhs)
    
    # 2. Comptage des moustiques capturés durant ce pas de temps
    # Quantité de moustiques dans chaque triangle pondérée par l'aire
    mass_per_triangle = u_new * areas 
    
    # Pour chaque piège, capture = Somme(Efficacité_locale * Masse_locale) * dt
    # On utilise le produit matriciel pour aller vite : (N_pieges x N_tri) . (N_tri)
    captures_instant = trap_profiles.dot(mass_per_triangle)
    
    # Mise à jour des cumuls
    trapped_counts += captures_instant * dt
    
    u = u_new
    
    # Suivi population totale
    current_pop = np.sum(u * areas)
    history_total_pop.append(current_pop)

final_population = np.sum(u * areas)
total_trapped = np.sum(trapped_counts)
died_naturally = initial_population - final_population - total_trapped





# --- 7. Affichage et Sauvegarde des Résultats ---
plt.style.use('seaborn-v0_8-whitegrid') 

# === FIGURE 1 : Évolution Temporelle ===
fig1, ax1 = plt.subplots(figsize=(8, 6)) # On crée une nouvelle fenetre

time_axis = np.linspace(0, target_time, len(history_total_pop))
ax1.plot(time_axis, history_total_pop, linewidth=2, color='#2c3e50', label='Population Totale')
ax1.set_title("Dynamique de la Population", fontsize=14, fontweight='bold')
ax1.set_xlabel("Temps (s)", fontsize=12)
ax1.set_ylabel("Nombre d'individus", fontsize=12)
ax1.legend()
ax1.grid(True, linestyle='--', alpha=0.6)

# Pour sauvegarder cette image seule :
# fig1.savefig("figure_1_dynamique.png", dpi=300)


# === FIGURE 2 : Carte de Densité Finale ===
fig2, ax2 = plt.subplots(figsize=(8, 7)) # Format un peu plus carré pour la carte
ax2.set_aspect('equal')

triang = Triangulation(points[:, 0], points[:, 1], triangles)
# Note : on garde shading='flat' et facecolors comme corrigé précédemment
tpc = ax2.tripcolor(triang, facecolors=u, shading='flat', cmap='inferno', vmin=0, vmax=np.max(u))

# On affiche les pièges
ax2.scatter(TRAP_POS[:,0], TRAP_POS[:,1], c='cyan', marker='x', s=40, label='Pièges')
ax2.set_title(f"Densité Finale à t={target_time}", fontsize=14, fontweight='bold')
# Attention : on utilise fig2 pour la colorbar ici
cbar = fig2.colorbar(tpc, ax=ax2, label='Densité ($u$)')
ax2.legend(loc='upper right')

# Pour sauvegarder cette image seule :
# fig2.savefig("figure_2_carte.png", dpi=300)


# === FIGURE 3 : Efficacité par Piège ===
fig3, ax3 = plt.subplots(figsize=(10, 6)) # Format large pour l'histogramme

trap_indices = np.arange(len(TRAP_POS))
bars = ax3.bar(trap_indices, trapped_counts, color='#e74c3c', edgecolor='black', alpha=0.8)
donnees_reelles=np.array([0,1,9,260,7,26,447,67,0,0,0,0,4,0,0,0,0,53,387,0,0]) ############################################# données réelles fournies
bars_real = ax3.bar(trap_indices, donnees_reelles, color='#3498db', edgecolor='black', alpha=0.8, width=0.4, label='Données Réelles')


ax3.set_xlabel("Identifiant du Piège", fontsize=12)
ax3.set_ylabel("Nombre de captures cumulées", fontsize=12)
ax3.set_title("Efficacité cumulée des pièges", fontsize=14, fontweight='bold')
ax3.set_xticks(trap_indices)
ax3.grid(axis='y', linestyle='--', alpha=0.5)

# Annotation des valeurs max
for bar in bars:
    height = bar.get_height()
    if height > 5: 
        ax3.text(bar.get_x() + bar.get_width()/2., height,
                 f'{int(height)}', ha='center', va='bottom', fontsize=9)

# Pour sauvegarder cette image seule :
# fig3.savefig("figure_3_efficacite.png", dpi=300)

# Affiche toutes les fenêtres en même temps
plt.show()



