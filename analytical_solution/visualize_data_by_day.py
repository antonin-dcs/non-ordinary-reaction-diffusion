import csv
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

#!/usr/bin/env python3
# pieges.py
# Visualisation 3D de la capture des moustiques par les pièges au cours du temps.
# Lecture des fichiers relatif: "data/trap_coordinates(in).csv" et "data/trap_counts(in).csv"



# Coordonnées des pièges

data_path_coordinates = Path(__file__).parent / "data" / "trap_coordinates(in).csv"

trap_coordinates = []
with data_path_coordinates.open(newline="") as f:
    reader = csv.DictReader(f)
    for row in reader:
        try:
            trap_coordinates.append([float(row["x"]), 
                                     float(row["y"])]) # Attention : id = index + 1
        except (KeyError, ValueError):
            # ignore malformed rows
            continue

n = len(trap_coordinates)
ids = list(range(1,n+1))
xs = [trap[0] for trap in trap_coordinates]
ys = [trap[1] for trap in trap_coordinates]





# Donnée sur chaque piège (nombre de moustiques capturés à chaque instant et nombre cumulé de captures)

data_pat_counts = Path(__file__).parent / "data" / "trap_counts_by_day(in).csv"

trap_counts = []
ids = []

with data_pat_counts.open(newline="") as f:
    reader = csv.DictReader(f)
    for row in reader:
        try:
            L = list(row.values())
            L = [int(n_cum) for n_cum in L]
            ids.append(L[0])
            trap_counts.append(L[1:])
        except (KeyError, ValueError):
            # ignore malformed rows
            continue



N = len(trap_counts[0])


fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, projection="3d")

dx = dy = 10

def update(frame):
    L_n_cum  = [trap_counts[index][frame] for index in range(len(trap_counts))]
    ax.clear()
    for xi, yi, idi, n_cum_i in zip(xs, ys, ids, L_n_cum):
        ax.text(xi + dx, yi + dy, n_cum_i+10, str(idi), fontsize=8, ha="center", va="bottom", zorder=3)
    ax.set_xlabel("$x$")
    ax.set_ylabel("$y$")
    ax.set_zlabel("$N_\\text{cum}$")
    ax.set_zlim(0, 500)
    ax.set_title(f"Nombre de moustique attrapés par pièges (J={frame})")
    ax.grid(True)
    ax.scatter([0], [0], [0], color='red', s=150, edgecolor='black', linewidth=2, zorder=15, marker='o')
    bars = ax.bar3d(xs, ys, [0]*len(ids), dx, dy, L_n_cum, shade=True, color='blue', edgecolor='r', linewidth=3, alpha=0.85, zorder=2)
    plt.tight_layout()
    return [bars]


ani = FuncAnimation(fig, update, frames=N, interval=500)
plt.show()
plt.close()

# Calcul des statistiques
total_initial = 10000
total_final = sum([trap_counts[i][-1] for i in range(len(trap_counts))])
ratio = total_final / total_initial

print(f"\n=== Statistiques de capture ===")
print(f"Nombre total initial de moustiques : {total_initial}")
print(f"Nombre total final de moustiques capturés : {total_final}")
print(f"Ratio (capturés/initial) : {ratio*100:.1f}%")
