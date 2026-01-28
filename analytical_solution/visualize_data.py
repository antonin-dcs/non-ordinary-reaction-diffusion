import csv
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from matplotlib.colors import Normalize
from math import isclose
from mpl_toolkits.mplot3d import Axes3D  # Nécessaire pour enregistrer la projection 3D
from matplotlib.animation import FuncAnimation

#!/usr/bin/env python3
# pieges.py
# Visualisation 3D de la capture des moustiques par les pièges au cours du temps.
# Lecture des fichiers relatif: "data/trap_coordinates(in).csv" et "data/trap_counts(in).csv"

# TODO : mettre à jour en se basant sur visualize_data_by_day pour éviter le problème d'espaces dans le csv

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

data_pat_counts = Path(__file__).parent / "data" / "trap_counts(in).csv"

trap_counts = []
with data_pat_counts.open(newline="") as f:
    reader = csv.DictReader(f)
    for row in reader:
        try:
            trap_counts.append({
                "id": int(row["id"]),
                "t": float(row["t"]),
                "n": int(row["n"]),
                "n_cum": int(row["n_cum"]),
            })
        except (KeyError, ValueError):
            # ignore malformed rows
            continue


trap_counts_dict = {}

for d in trap_counts:
    id = d["id"]
    if id not in trap_counts_dict :
        trap_counts_dict[id] = {"t":[],"n":[],"n_cum":[]}
    trap_counts_dict[id]["t"].append(d["t"])
    trap_counts_dict[id]["n"].append(d["n"])
    trap_counts_dict[id]["n_cum"].append(d["n_cum"])



def n_cum_at_time_for_trap(trap_dict, t):
    """Retourne n_cum au temps le plus grand t <= T pour un piège donné.
    Si aucun temps <= T, retourne 0.
    Les listes dans trap_dict sont supposées triées par t croissant.
    """
    times = trap_dict.get("t", [])
    n_cum = trap_dict.get("n_cum", [])
    idx = None
    for i, tt in enumerate(times):
        if isclose(tt, t) or tt <= t:
            idx = i
        else:
            break
    if idx is None:
        return 0
    return n_cum[idx]


def get_n_cum_at_time(t):
    """Retourne une liste des n_cum pour tous les pièges au temps T, dans l'ordre des ids."""
    L_n_cum = []
    for tid in ids:
        trap_d = trap_counts_dict.get(tid, {})
        L_n_cum.append(n_cum_at_time_for_trap(trap_d, t))
    return L_n_cum


N_data = len(trap_counts_dict[21]["t"])
t_max = 20


fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, projection="3d")

dx = dy = 10

def update(frame):
    t = frame*t_max/N_data
    L_n_cum  = get_n_cum_at_time(t)
    
    ax.clear()
    for xi, yi, idi, n_cum_i in zip(xs, ys, ids, L_n_cum):
        ax.text(xi + dx, yi + dy, n_cum_i+10, str(idi), fontsize=8, ha="center", va="bottom", zorder=3)
    ax.set_xlabel("$x$")
    ax.set_ylabel("$y$")
    ax.set_zlabel("$N_\\text{cum}$")
    ax.set_zlim(0, 500)
    ax.set_title(f"Nombre de moustique attrapés par pièges (J={int(t)})")
    ax.grid(True)
    ax.scatter([0], [0], [0], color='red', s=150, edgecolor='black', linewidth=2, zorder=15, marker='o')
    bars = ax.bar3d(xs, ys, [0]*len(ids), dx, dy, L_n_cum, shade=True, color='blue', edgecolor='r', linewidth=3, alpha=0.85, zorder=2)
    plt.tight_layout()
    return [bars]


ani = FuncAnimation(fig, update, frames=N_data, interval=25)

# Frame à exporter (dernière frame)
frame_to_export = N_data - 1
update(frame_to_export)
out_path = Path(__file__).parent / "figures" / "moustiques_frame_last.pdf"
plt.savefig(out_path, format="pdf", bbox_inches="tight")  # ou "svg"

# ani.save("./figures/moustiques.gif", writer="pillow", fps=30)

plt.show()   # optionnel, après les exports
plt.close()
