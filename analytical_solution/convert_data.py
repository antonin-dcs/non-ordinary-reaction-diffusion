import csv
from pathlib import Path
from math import isclose


# Donnée sur chaque piège (nombre de moustiques capturés à chaque instant et nombre cumulé de captures)

data_path_counts = Path(__file__).parent / "data" / "trap_counts(in).csv"

Lt = []
trap_counts = []


with data_path_counts.open(newline="") as f:
    reader = csv.DictReader(f)
    for i, row in enumerate(reader):
        row = {name.strip(): value for name, value in row.items()}
        try:
            id = int(row["id"])
            index = id - 1
            if len(trap_counts) <= index:
                trap_counts.append([])
            trap_counts[index].append(int(row["n"]))
            if id == 1:
                t = round(float(row["t"]),5)
                Lt.append(t)
        except (KeyError, ValueError):
            # ignore malformed rows
            continue


trap_counts_day = [[0] for i in range(len(trap_counts))]

for idx in range(len(trap_counts)):
    n_cum_idx = 0
    for i in range(1,len(Lt)):
        n_cum_idx += trap_counts[idx][i]
        t = Lt[i]
        if isclose(t,round(t), abs_tol = .001):
            if round(t)==1:
                trap_counts_day[idx].append(n_cum_idx)
            else:
                trap_counts_day[idx].append(n_cum_idx+trap_counts_day[idx][-1])
            n_cum_idx = 0


with open(Path(__file__).parent / "data" / "trap_counts_by_day(in).csv", "w", newline="") as f:
    writer = csv.writer(f)
    headers = ["id"] + [f"day_{i}" for i in range(len(trap_counts_day[0]))]
    writer.writerow(headers)
    for idx, counts in enumerate(trap_counts_day):
        writer.writerow([idx + 1] + counts)
