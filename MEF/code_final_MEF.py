import numpy as np
from scipy.sparse import lil_matrix
# -----------------------------
# Paramètres du domaine
# -----------------------------
Lx, Ly = 12.0, 12.0
Jx, Jy = 100,100
hx = Lx / Jx
hy = Ly / Jy
Jn = Jx * Jy  # nombre de nœuds globaux

# -----------------------------
# Fonctions utilitaires
# -----------------------------
def Coord_ji_to_n(j, i):
    """Indices nœud global avec périodicité"""
    j = (j-1) % Jx + 1
    i = (i-1) % Jy + 1
    return (i-1)*Jx + (j-1)

def phi_triangle(x, y, tri_nodes, local_index):
    """Fonction de base P1 sur triangle"""
    x1, y1 = tri_nodes[0]
    x2, y2 = tri_nodes[1]
    x3, y3 = tri_nodes[2]
    B = np.array([[1, x1, y1],
                  [1, x2, y2],
                  [1, x3, y3]])
    b = np.zeros(3)
    b[local_index] = 1
    coeffs = np.linalg.solve(B, b)
    return coeffs[0] + coeffs[1]*x + coeffs[2]*y

# -----------------------------
# Définition de K(x,y)
# -----------------------------
import numpy as np

# Paramètres des pièges
I = 0.003       # Intensité
r_trap = 0.7 # Rayon (écart-type de la gaussienne)
D = 0.02     # diffusion

# Positions originales des pièges
TRAP_POS_ORIG = np.array([
    (-131.68167734441144, 131.68167734441138),
    (-124.68511891759317,   51.6462672817693),
    (-108.33035524351197,  -44.87190235855795),
     ( 42.25450016994558,     0.0),
     ( 49.1781122967089,   118.7264656786217),
     (-37.346286147410396,   90.16191052134471),
     ( 12.761346832843234,  -30.808616597997073),
     (  0.0,               -76.8718590487307),
     ( 73.57099561684544, -177.61609547975),
     (  0.0,              -362.84280384397),
    (-130.79123869321234, -315.7579822927304),
    ( -86.98706983684231, -210.00536375120058),
    ( -50.62404310504939, -122.21725144637058),
    ( 335.9236323150811,     0.0),
    ( 233.86016009458638, -233.86016009458635),
    ( 200.78991294270418, -200.78991294270415),
    (  68.8919448551579, -166.31986760758153),
    (  76.60410820779381,     0.0),
    ( -26.42814033157892,   26.428140331578913),
    (-245.38328803096954,  245.38328803096945),
    ( -85.74349140349476,  207.0030998315381)
])

# Normalisation sur [0,10] et arrondi à 0.01
x_orig = TRAP_POS_ORIG[:,0]
y_orig = TRAP_POS_ORIG[:,1]

x_min, x_max = x_orig.min(), x_orig.max()
y_min, y_max = y_orig.min(), y_orig.max()

x_norm = np.round(10 * (x_orig - x_min) / (x_max - x_min), 2)
y_norm = np.round(10 * (y_orig - y_min) / (y_max - y_min), 2)

TRAP_POS = np.column_stack((x_norm, y_norm))

# -----------------------------
# Définition de K(x,y)
# -----------------------------
def K(x, y):
    """Somme de gaussiennes pour les pièges normalisés"""
    val = 0.0
    for x0, y0 in TRAP_POS:
        val += I * np.exp(-((x - x0)**2 + (y - y0)**2) / (2 * r_trap**2))
    return val


# -----------------------------
# Intégration numérique sur un triangle
# -----------------------------
def integrate_phi_i_phi_j_K(tri_nodes, i_local, j_local):
    """Quadrature 3 points de Gauss pour ∫ K φ_i φ_j"""
    bary = np.array([[1/6, 1/6, 2/3],
                     [1/6, 2/3, 1/6],
                     [2/3, 1/6, 1/6]])
    weights = np.array([1/3, 1/3, 1/3])

    x1, y1 = tri_nodes[0]
    x2, y2 = tri_nodes[1]
    x3, y3 = tri_nodes[2]

    val = 0.0
    for k in range(3):
        l1, l2, l3 = bary[k]
        x = l1*x1 + l2*x2 + l3*x3
        y = l1*y1 + l2*y2 + l3*y3
        val += weights[k] * K(x, y) * phi_triangle(x, y, tri_nodes, i_local) * phi_triangle(x, y, tri_nodes, j_local)

    area = 0.5 * abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1))
    return val * area

# -----------------------------
# Matrice R
# -----------------------------
R = lil_matrix((Jn, Jn))

for i in range(1, Jy+1):
    for j in range(1, Jx+1):
        x0, y0 = (j-1)*hx, (i-1)*hy
        x1, y1 = j*hx, i*hy

        # Triangle 1
        tri1_nodes = [(x0, y0), (x1, y0), (x1, y1)]
        nodes1 = [Coord_ji_to_n(j,i), Coord_ji_to_n(j+1,i), Coord_ji_to_n(j+1,i+1)]
        for a in range(3):
            for b in range(3):
                R[nodes1[a], nodes1[b]] += integrate_phi_i_phi_j_K(tri1_nodes, a, b)

        # Triangle 2
        tri2_nodes = [(x0, y0), (x1, y1), (x0, y1)]
        nodes2 = [Coord_ji_to_n(j,i), Coord_ji_to_n(j+1,i+1), Coord_ji_to_n(j,i+1)]
        for a in range(3):
            for b in range(3):
                R[nodes2[a], nodes2[b]] += integrate_phi_i_phi_j_K(tri2_nodes, a, b)

def integrate_phi_i_phi_j(tri_nodes, i_local, j_local):
    """Intégration numérique par quadrature 3 points de Gauss (exacte pour P1*P1)"""
    # Points de Gauss barycentriques et poids
    bary = np.array([[1/6, 1/6, 2/3],
                     [1/6, 2/3, 1/6],
                     [2/3, 1/6, 1/6]])
    weights = np.array([1/3, 1/3, 1/3])

    # Coordonnées cartésiennes des points
    x1, y1 = tri_nodes[0]
    x2, y2 = tri_nodes[1]
    x3, y3 = tri_nodes[2]

    val = 0.0
    for k in range(3):
        l1, l2, l3 = bary[k]
        x = l1*x1 + l2*x2 + l3*x3
        y = l1*y1 + l2*y2 + l3*y3
        val += weights[k] * phi_triangle(x, y, tri_nodes, i_local) * phi_triangle(x, y, tri_nodes, j_local)

    # Aire du triangle
    area = 0.5 * abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1))
    return val * area

# -----------------------------
# Matrice de masse M
# -----------------------------
M = lil_matrix((Jn, Jn))

for i in range(1, Jy+1):
    for j in range(1, Jx+1):
        x0, y0 = (j-1)*hx, (i-1)*hy
        x1, y1 = j*hx, i*hy

        # Triangle 1
        tri1_nodes = [(x0, y0), (x1, y0), (x1, y1)]
        nodes1 = [Coord_ji_to_n(j,i), Coord_ji_to_n(j+1,i), Coord_ji_to_n(j+1,i+1)]
        for a in range(3):
            for b in range(3):
                M[nodes1[a], nodes1[b]] += integrate_phi_i_phi_j(tri1_nodes, a, b)

        # Triangle 2
        tri2_nodes = [(x0, y0), (x1, y1), (x0, y1)]
        nodes2 = [Coord_ji_to_n(j,i), Coord_ji_to_n(j+1,i+1), Coord_ji_to_n(j,i+1)]
        for a in range(3):
            for b in range(3):
                M[nodes2[a], nodes2[b]] += integrate_phi_i_phi_j(tri2_nodes, a, b)

# Gradient de base pour un triangle donné
def grad_phi_triangle(tri_nodes):
    x1,y1 = tri_nodes[0]
    x2,y2 = tri_nodes[1]
    x3,y3 = tri_nodes[2]
    B = np.array([[1,x1,y1],[1,x2,y2],[1,x3,y3]])
    grads = []
    for k in range(3):
        b = np.zeros(3)
        b[k] = 1
        sol = np.linalg.solve(B,b)
        grads.append(sol[1:])
    return grads

# -----------------------------
# Matrice de rigidité
# -----------------------------
A = lil_matrix((Jn, Jn))

for i in range(1,Jy+1):
    for j in range(1,Jx+1):
        x0,y0 = (j-1)*hx, (i-1)*hy
        x1,y1 = j*hx, i*hy

        # Triangle bas-gauche, bas-droit, haut-droit
        tri1_nodes = [(x0,y0),(x1,y0),(x1,y1)]
        grads1 = grad_phi_triangle(tri1_nodes)
        nodes1 = [Coord_ji_to_n(j,i), Coord_ji_to_n(j+1,i), Coord_ji_to_n(j+1,i+1)]

        # Triangle bas-gauche, haut-droit, haut-gauche
        tri2_nodes = [(x0,y0),(x1,y1),(x0,y1)]
        grads2 = grad_phi_triangle(tri2_nodes)
        nodes2 = [Coord_ji_to_n(j,i), Coord_ji_to_n(j+1,i+1), Coord_ji_to_n(j,i+1)]

        area = hx*hy/2

        for a in range(3):
            for b in range(3):
                A[nodes1[a], nodes1[b]] += area * np.dot(grads1[a], grads1[b])
                A[nodes2[a], nodes2[b]] += area * np.dot(grads2[a], grads2[b])



# =========================================================
# Résolution en temps : Crank–Nicolson (FEM)
# =========================================================
M = M.tocsr()
A = A.tocsr()
R = R.tocsr()
from scipy.sparse.linalg import spsolve


import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# -----------------------------
# Paramètres physiques
# -----------------------------

dt = 0.05        # pas de temps
T = 20     # temps final
Nt = int(T / dt)
L_reel=[0,1,9,260,7,26,447,67,0,0,0,0,4,0,0,0,0,53,387,0,0]
print(len(L_reel))
# -----------------------------
# Matrices Crank–Nicolson
# -----------------------------
Kmat = D * A + R

LHS = M + 0.5 * dt * Kmat
RHS_mat = M - 0.5 * dt * Kmat



# -----------------------------
# Condition initiale : 10000 moustiques au centre
# -----------------------------
U = np.zeros(Jn)
x_center_orig = 0.0  # si tu veux vraiment le point (0,0) original
y_center_orig = 0.0

x_center = 10 * (x_center_orig - x_min) / (x_max - x_min)
y_center = 10 * (y_center_orig - y_min) / (y_max - y_min)

j_center = int(round(x_center / hx))
i_center = int(round(y_center / hy))

j_center = min(max(j_center, 0), Jx-1)
i_center = min(max(i_center, 0), Jy-1)

n_center = i_center * Jx + j_center
U = np.zeros(Jn)
U[n_center] = 10000

# -----------------------------
# Coordonnées des nœuds
# -----------------------------
x_nodes = np.array([ (n % Jx) * hx for n in range(Jn) ])
y_nodes = np.array([ (n // Jx) * hy for n in range(Jn) ])

# -----------------------------
# Fonction : captures par piège
# -----------------------------
def trap_contrib(U, k):
    x0, y0 = TRAP_POS[k]
    # Gaussienne du piège sur tous les nœuds
    Kk = I * np.exp(-((x_nodes - x0)**2 + (y_nodes - y0)**2)/(2*r_trap**2))
    return dt * (Kk * U).sum()

# -----------------------------
# Boucle en temps et sauvegarde
# -----------------------------
U_hist = [U.copy()]
time = [0.0]
cumul = np.zeros(len(TRAP_POS))
cumulative_captures = np.zeros((Nt+1, len(TRAP_POS)))

for n in range(Nt):
    # Crank–Nicolson
    RHS = RHS_mat @ U
    U = spsolve(LHS, RHS)

    # Captures par piège ce pas de temps
    captures_this_step = np.zeros(len(TRAP_POS))
    for k in range(len(TRAP_POS)):
        captures_this_step[k] = trap_contrib(U, k)

    cumul += captures_this_step
    cumulative_captures[n+1] = cumul.copy()  # stocker cumul

    # Sauvegarde de U pour animation
    U_hist.append(U.copy())
    time.append((n+1)*dt)

# -----------------------------
# Reconstruction 2D
# -----------------------------
def U_to_grid(U):
    grid = np.zeros((Jy, Jx))
    for n in range(Jn):
        j = n % Jx
        i = n // Jx
        grid[i, j] = U[n]
    return grid

# ----------------------------- Animation et captures par piège -----------------------------
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# ----------------- Préparer figure avec deux axes : densité et barres -----------------
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

# ----------------- Coordonnées des nœuds pour calcul des captures -----------------
x_nodes = np.array([ (n % Jx) * hx for n in range(Jn) ])
y_nodes = np.array([ (n // Jx) * hy for n in range(Jn) ])

# ----------------- Fonction pour calculer captures d’un piège k -----------------
def trap_contrib(U, k):
    x0, y0 = TRAP_POS[k]
    Kk = I * np.exp(-((x_nodes - x0)**2 + (y_nodes - y0)**2)/(2*r_trap**2))
    return dt * (Kk * U).sum()

# ----------------- Calcul cumulatif des captures -----------------
cumulative_captures = np.zeros((len(U_hist), len(TRAP_POS)))
cumul = np.zeros(len(TRAP_POS))

for t, U_t in enumerate(U_hist):
    captures_this_step = np.array([trap_contrib(U_t, k) for k in range(len(TRAP_POS))])
    cumul += captures_this_step
    cumulative_captures[t] = cumul.copy()  # Stocker le cumul par frame

# ----------------- Densité initiale -----------------
grid_init = U_to_grid(U_hist[0])
im = ax1.imshow(
    grid_init,
    origin="lower",
    extent=[0, Lx, 0, Ly],
    cmap="inferno"
)
ax1.set_title("Densité – t = 0")
'''
ax1.scatter(TRAP_POS[:,0], TRAP_POS[:,1], color='red', marker='x', s=100)
'''
# ----------------- Barres initiales : SIMULATION vs RÉEL -----------------
indices = np.arange(len(TRAP_POS))
width = 0.35

bars_sim = ax2.bar(
    indices - width/2,
    np.zeros(len(TRAP_POS)),
    width,
    label="Simulation",
    color="orange"
)

bars_real = ax2.bar(
    indices + width/2,
    L_reel,
    width,
    label="Observations réelles",
    color="steelblue"
)

# Échelle fixe
max_y = max(cumulative_captures[-1].max(), max(L_reel))
ax2.set_ylim(0, 1.1 * max_y)

ax2.set_xlabel("Piège #")
ax2.set_ylabel("Moustiques capturés")
ax2.set_xticks(indices)
ax2.set_xticklabels([f"P{i+1}" for i in range(len(TRAP_POS))], rotation=45)
ax2.legend()

# ----------------- Fonction d’update pour l’animation -----------------
def update(frame):
    # -------- Carte de densité --------
    grid = U_to_grid(U_hist[frame])
    im.set_array(grid)
    im.set_clim(0, 1.1 * grid.max())
    ax1.set_title(f"Densité – t = {time[frame]:.2f}")

    # -------- Barres SIMULATION --------
    for k, b in enumerate(bars_sim):
        b.set_height(cumulative_captures[frame, k])

    # -------- Titre stats --------
    ax2.set_title(
        f"Moustiques morts (sim) : {int(cumulative_captures[frame].sum())} | "
        f"Moustiques vivants : {int(U_hist[frame].sum())}"
    )

    return [im, *bars_sim, *bars_real]

ani = FuncAnimation(fig, update, frames=len(U_hist), interval=100, blit=False, repeat=False)
plt.show()

# a faire : nb de moustiques capturés dans chaque piège chaque jour
# =========================================================
# ERREUR QUADRATIQUE ENTRE SIMULATION ET DONNÉES RÉELLES
# =========================================================

def rmse(simulated, observed):
    simulated = np.array(simulated)
    observed = np.array(observed)
    return np.sqrt(np.mean((simulated - observed)**2))


# Captures simulées après 20 jours
sim_20_days = cumulative_captures[-1]

# Vérification cohérence
assert len(sim_20_days) == len(L_reel), "Mismatch between traps and real data"

# Calcul RMSE
error_rmse = rmse(sim_20_days, L_reel)

print("===================================")
print("CAPTURES APRÈS 20 JOURS (SIMULATION)")
print(sim_20_days.astype(int))
print("===================================")
print("DONNÉES RÉELLES")
print(L_reel)
print("===================================")
print(f"RMSE après 20 jours = {error_rmse:.2f}")
print("===================================")
