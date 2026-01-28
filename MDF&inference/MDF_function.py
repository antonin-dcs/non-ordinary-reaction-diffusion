import numpy as np

def simulation_MDF(sigma,nu,gamma,R_trap,N,initial_total_mosquitoes,x0,y0,R0):
    # --- 1. CONFIGURATION ---
    # NEW: Define the physical boundaries
    L_min, L_max = -300, 300
    L_total = L_max - L_min  # Total length is 1000
    
    T_max = 20 #time of simulation in day
    dx = L_total / (N - 1)  # Modified to use total length (1000)
    D = (sigma**2) / 2     
    dt = 0.9 * (dx**2) / (4 * D) 
    steps = int(T_max / dt)

    print(f"Simulation MDF: Grid {N}x{N}, Range=[{L_min}, {L_max}], dt={dt:.5f}, Steps={steps}")

    # --- 2. INITIALISATION ---
    # NEW: Create grid centered on 0 (-500 to 500)
    x = np.linspace(L_min, L_max, N)
    y = np.linspace(L_min, L_max, N)
    X, Y = np.meshgrid(x, y)

    # Positions des pièges (Données brutes)
    raw_trap_positions = np.array([
        (-132, 132), (-125, 52), (-108, -45), (42, 0), (49, 119), (-37, 90), 
        (13, -31), (0, -77), (74,-178), (0, -363), (-131, -316), (-87, -210), 
        (-51, -122), (336, 0), (234, -234), (201, -201), (69, -166), (77, 0), 
        (-26, 26), (-245, 245), (-86, 207)
    ])

    # --- CHANGED: REMOVED NORMALIZATION ---
    # We now use the raw positions directly as they fit within -500 to 500
    traps = []
    for i in range(len(raw_trap_positions)):
        # raw_trap_positions is already x, y. We assign it directly.
        traps.append({'q': raw_trap_positions[i], 'id': i + 1})

    # --- PRÉ-CALCUL DES CARTES D'EFFICACITÉ ---
    V_natural = np.full((N, N), nu)
    V_traps_total = np.zeros((N, N)) 
    trap_maps = np.zeros((len(traps), N, N))

    for i, trap in enumerate(traps):
        # The distance calculation works automatically with negative coordinates
        dist_sq = (X - trap['q'][0])**2 + (Y - trap['q'][1])**2
        f_i = gamma * np.exp(-dist_sq / (R_trap**2)) 
        
        trap_maps[i, :, :] = f_i      
        V_traps_total += f_i          

    V = V_natural + V_traps_total 

    # --- VARIABLES DE SORTIE ---
    daily_captures_matrix = np.zeros((len(traps), T_max))
    capture_map_spatial = np.zeros((N, N))

    # Condition initiale h0
    # Ensure x0, y0 are set to 0 (or your desired start) in the function arguments
    h_0 = np.exp(-((X - x0)**2 + (Y - y0)**2) / (2 * R0**2)) 
    h = h_0 * initial_total_mosquitoes / np.sum(h_0)

    # --- 3. BOUCLE TEMPORELLE ---
    alpha = D * dt / (dx**2)

    for n in range(steps):
        h_inner = h[1:-1, 1:-1]
        
        # --- A. CALCUL DES CAPTURES ---
        spatial_loss = h_inner * V_traps_total[1:-1, 1:-1] * dt
        capture_map_spatial[1:-1, 1:-1] += spatial_loss
        
        current_day = int(n * dt)
        if current_day < T_max:
            captures_at_step = np.sum(trap_maps[:, 1:-1, 1:-1] * h_inner, axis=(1, 2)) * dt
            daily_captures_matrix[:, current_day] += captures_at_step

        # --- B. RÉACTION-DIFFUSION ---
        d2x = h[1:-1, 2:] - 2*h[1:-1, 1:-1] + h[1:-1, :-2]
        d2y = h[2:, 1:-1] - 2*h[1:-1, 1:-1] + h[:-2, 1:-1]
        laplacian = d2x + d2y
        
        reaction = V[1:-1, 1:-1] * h[1:-1, 1:-1]
        h[1:-1, 1:-1] = h[1:-1, 1:-1] + alpha * laplacian - dt * reaction
        
        # Conditions aux limites (Neumann)
        h[0, :] = h[1, :] 
        h[-1, :] = h[-2, :]
        h[:, 0] = h[:, 1]
        h[:, -1] = h[:, -2]

    return h, capture_map_spatial, traps, daily_captures_matrix