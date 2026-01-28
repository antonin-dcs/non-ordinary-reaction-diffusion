#include <iostream>  // Pour afficher du texte
#include <vector>    // Pour les tableaux dynamiques
#include <fstream>   // Pour écrire dans un fichier (CSV)
#include <cmath>     // Pour les maths

// Paramètres de simulation
const int N = 100;          // Nombre de points (taille grille)
const double L = 10.0;      // Longueur domaine
const double T_max = 1.0;   // Temps total
const double sigma = 1.0;   // Paramètre du modèle
const double D = (sigma * sigma) / 2.0; // Coefficient de diffusion

int main() {
    // 1. Initialisation des pas
    double dx = L / (N - 1);
    // Condition CFL stricte (avec marge de sécurité 0.9)
    double dt = 0.9 * (dx * dx) / (4.0 * D); 
    
    int n_steps = static_cast<int>(T_max / dt);

    std::cout << "Simulation C++ lancee..." << std::endl;
    std::cout << "Grille: " << N << "x" << N << std::endl;
    std::cout << "Pas de temps dt: " << dt << std::endl;
    std::cout << "Iterations totales: " << n_steps << std::endl;

    // 2. Création des vecteurs (Tableau 1D pour simuler la 2D)
    // h_old : état à t, h_new : état à t+1
    std::vector<double> h_old(N * N, 0.0);
    std::vector<double> h_new(N * N, 0.0);

    // 3. Conditions initiales (Un pic de moustiques au centre)
    int center = N / 2;
    h_old[center * N + center] = 100.0; 
    // Tu peux ajouter une boucle ici pour faire une tache gaussienne plus large

    // 4. Boucle Temporelle
    for (int t = 0; t < n_steps; t++) {
        
        // Boucle Spatiale (On évite les bords pour l'instant -> i de 1 à N-2)
        for (int i = 1; i < N - 1; i++) {
            for (int j = 1; j < N - 1; j++) {
                
                // Index actuel et voisins
                int idx = i * N + j;
                int idx_left = i * N + (j - 1);
                int idx_right = i * N + (j + 1);
                int idx_up = (i - 1) * N + j;
                int idx_down = (i + 1) * N + j;

                // Laplacien discret (Schéma différences finies)
                double laplacian = h_old[idx_left] + h_old[idx_right] + 
                                   h_old[idx_up] + h_old[idx_down] - 
                                   4.0 * h_old[idx];

                // Mise à jour explicite
                h_new[idx] = h_old[idx] + (D * dt / (dx * dx)) * laplacian;
            }
        }

        // Conditions aux limites (Neumann : réflexion)
        // Simplifié ici : on copie juste les valeurs voisines sur les bords
        // (À implémenter proprement si nécessaire, mais suffisant pour apprendre)

        // Mise à jour pour le pas suivant (échange des références ou copie)
        h_old = h_new;

        // Affichage progression tous les 1000 pas
        if (t % 1000 == 0) {
             std::cout << "Pas " << t << "/" << n_steps << " termine." << std::endl;
        }
    }

    // 5. Sauvegarde des résultats dans un fichier CSV
    // C++ ne sait pas tracer de courbes, on exporte les données pour Python
    std::ofstream file("resultats.csv");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            file << h_old[i * N + j];
            if (j < N - 1) file << ",";
        }
        file << "\n";
    }
    file.close();

    std::cout << "Termine ! Resultats dans 'resultats.csv'" << std::endl;
    return 0;
}