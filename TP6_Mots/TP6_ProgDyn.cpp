// Imagine++ project
// Project:  TP6_ProgDyn
// Author:   Dhia GARBAYA

using namespace std;
#include <vector>
#include <iostream>
#include <algorithm>
#include <time.h>



// COMPLEXITE : Temps: O(3^max(m,n))    (naive)
int minDistanceRecursive(const string& Mot1, const string& Mot2, int i=0, int j=0, bool damerau = false) {
    if (i >= Mot1.length()) return Mot2.length()-j;
    if (j >= Mot2.length()) return Mot1.length()-i;

    if (Mot1[i] == Mot2[j]) {
        return minDistanceRecursive(Mot1, Mot2, i + 1, j + 1);
    }
    else {
        return 1 + min({(damerau && i < Mot1.length() - 1 && j < Mot2.length() - 1 && Mot1[i] == Mot2[j + 1] && Mot1[i + 1] == Mot2[j]) ? minDistanceRecursive(Mot1, Mot2, i + 2, j + 2, damerau):999999,
                        minDistanceRecursive(Mot1, Mot2, i + 1, j),    //On supprime un caractère de s1
                        minDistanceRecursive(Mot1, Mot2, i, j + 1),    // On ajoute un caractère à s1
                        minDistanceRecursive(Mot1, Mot2, i + 1, j + 1)});  // On substitue un caractère de s1 par un de s2
    }
}


//Memoisée (sur 2 étapes)
//COMPLEXITE Temps et Espace : en O(m * n)
int minDistanceMemoisee(const string& Mot1, const string& Mot2, int i, int j, vector<vector<int>>& memo) {
    if (i == 0) return j;
    if (j == 0) return i;

    if (memo[i][j] != -1) {
        return memo[i][j];
    }

    if (Mot1[i - 1] == Mot2[j - 1]) {
        memo[i][j] = minDistanceMemoisee(Mot1, Mot2, i - 1, j - 1, memo);
    } else {
        memo[i][j] = 1 + min({ minDistanceMemoisee(Mot1, Mot2, i - 1, j, memo),
                                   minDistanceMemoisee(Mot1, Mot2, i, j - 1, memo),
                                   minDistanceMemoisee(Mot1, Mot2, i - 1, j - 1, memo) });
    }

    return memo[i][j];
}

int minDistance(const string& Mot1, const string& Mot2) {
    int m = Mot1.length();
    int n = Mot2.length();
    vector<vector<int>> memo(m + 1, vector<int>(n + 1, -1));
    return minDistanceMemoisee(Mot1, Mot2, m, n, memo);
}


//Ecrit la série des modifs effectué pour retrouver mot2 à partir de mot1 (pour la version iterative)
void printModifs(const vector<vector<int>>& DL, const string& Mot1, const string& Mot2) {
    int i = Mot1.length();
    int j = Mot2.length();

    vector<string> Modifs;

    while (i > 0 || j > 0) {
        if (i > 0 && j > 0 && Mot1[i - 1] == Mot2[j - 1]) {
            Modifs.push_back("RIEN");
            --i;
            --j;
        } else {
            if (i > 0 && DL[i][j] == DL[i - 1][j] + 1) {
                Modifs.push_back("SUPPRESSION");
                --i;
            } else if (j > 0 && DL[i][j] == DL[i][j - 1] + 1) {
                Modifs.push_back("INSERTION");
                --j;
            } else if (i > 0 && j > 0) {
                if ( i > 1 && j > 1 && Mot1[i - 2] == Mot2[j - 1] && Mot1[i - 1] == Mot2[j - 2]) {
                    Modifs.push_back("TRANSPOSITION");
                    i -= 2;
                    j -= 2;}
                else{
                    Modifs.push_back("SUBSTITUTION");
                    --i;
                    --j;}
            }
        }
    }
    cout << "Les modifications pour transformer '" << Mot1 << "' en '" << Mot2 << "' sont:" << endl;
    // Parcours dans le sens inverse
    for (int k = Modifs.size() - 1; k >= 0; --k) {
        cout << Modifs[k] << "   ";
    }
    cout << endl;
}


// COMPLEXITE Temps et Espace: en O(m*n)
int DL_iterative(const string& Mot1, const string& Mot2, bool damerau = false) {
    int m = Mot1.length();
    int n = Mot2.length();

    // Créer une matrice pour stocker les résultats intermédiaires
    vector<vector<int>> DL(m + 1, vector<int>(n + 1, 0));

    // Initialisation (1ère ligne et 1ère colonne)
    for (int i = 0; i <= m; ++i) {
        DL[i][0] = i;
    }
    for (int j = 0; j <= n; ++j) {
        DL[0][j] = j;
    }

    // Remplir le reste de la matrice
    for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= n; ++j) {
            if (Mot1[i - 1] == Mot2[j - 1]) {
                DL[i][j] = DL[i - 1][j - 1];
            } else {
                DL[i][j] = 1 + min({DL[i - 1][j], DL[i][j - 1], DL[i - 1][j - 1]});
                if (damerau && i>1 && j>1 && Mot1[i-1] == Mot2[j - 2] && Mot1[i - 2] == Mot2[j-1]){
                    DL[i][j] = min({DL[i][j], 1 + DL[i-2][j-2]});
                }
            }
        }
    }
    (damerau) ? cout<<"DAMERAU-LEVENSHTEIN : " : cout<<"LEVENSHTEIN : " ;
    printModifs(DL, Mot1, Mot2);
    // La valeur finale est
    return DL[m][n];
}


// COMPLEXITE Temps: O(m*n) / Espace: O( min (m,n) )
int DL_iterative_lineaire_espace(const string& Mot1, const string& Mot2, bool damerau = false) {
    int m = Mot1.length();
    int n = Mot2.length();

    // Utilisation de deux tableaux unidimensionnels (lignes)
    vector<int> lignePrecedente(n + 1, 0);
    vector<int> ligneActuelle(n + 1, 0);

    // Initialisation de la première ligne
    for (int j = 0; j <= n; ++j) {
        lignePrecedente[j] = j;
    }

    // Remplissage des tableaux
    for (int i = 1; i <= m; ++i) {
        // Initialisation de la première colonne
        ligneActuelle[0] = i;

        for (int j = 1; j <= n; ++j) {
            if (Mot1[i - 1] == Mot2[j - 1]) {
                ligneActuelle[j] = lignePrecedente[j - 1];
            } else {
                ligneActuelle[j] = 1 + min({lignePrecedente[j], ligneActuelle[j - 1], lignePrecedente[j - 1]});
                if (damerau && i > 1 && j > 1 && Mot1[i - 1] == Mot2[j - 2] && Mot1[i - 2] == Mot2[j - 1]) {
                    ligneActuelle[j] = min(ligneActuelle[j], 1 + lignePrecedente[j - 2]);
                }
            }
        }
        // Mise à jour de la ligne précédente (swap fait une permutation)
        swap(lignePrecedente, ligneActuelle);
    }
    // La valeur finale est dans la dernière case de la ligne précédente
    return lignePrecedente[n];
}








int main() {
    string s1 = "abcdefghijkl";
    string s2 = "khfiegfebbkl";
    // string s1 = "ecoles";
    // string s2 = "eclose";

    clock_t t1 = clock();
    int ddl_iter = DL_iterative(s1, s2, true);
    clock_t t2 = clock();
    float duration_ddl_iter = (float)(t2 - t1) / CLOCKS_PER_SEC;

    t1 = clock();
    int dl_iter_lin = DL_iterative_lineaire_espace(s1, s2);
    t2 = clock();
    float duration_dl_iter_lin = (float)(t2 - t1) / CLOCKS_PER_SEC;

    t1 = clock();
    int dl_iter = DL_iterative(s1, s2);
    t2 = clock();
    float duration_dl_iter = (float)(t2 - t1) / CLOCKS_PER_SEC;

    t1 = clock();
    int dl = minDistanceRecursive(s1, s2);
    t2 = clock();
    float duration_dl = (float)(t2 - t1) / CLOCKS_PER_SEC;

    t1 = clock();
    int dlmemo = minDistance(s1, s2);
    t2 = clock();
    float duration_dlmemo = (float)(t2 - t1) / CLOCKS_PER_SEC;

    cout << "Distance de Damerau Levenshtein (iterative): " << ddl_iter << " | Temps d'execution : " << duration_ddl_iter << " secondes" << endl;
    cout << "Distance de Levenshtein (iterative lineaire): " << dl_iter_lin << " | Temps d'execution : " << duration_dl_iter_lin << " secondes" << endl;
    cout << "Distance de Levenshtein (iterative): " << dl_iter << " | Temps d'execution : " << duration_dl_iter << " secondes" << endl;
    cout << "Distance de Levenshtein (recursive naive): " << dl << " | Temps d'execution : " << duration_dl << " secondes" << endl;
    cout << "Distance de Levenshtein (recursive Memoisee): " << dlmemo << " | Temps d'execution : " << duration_dlmemo << " secondes" << endl;

    return 0;
}


