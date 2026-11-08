# =====================================================
# MODULES
# =====================================================
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import random as rd
import heapq


# =====================================================
# 1. DONNÉES DU PROBLÈME (Ville et Clients)
# =====================================================
Graphe = {
    'A': ['B', 'D'],
    'B': ['A', 'C', 'E'],
    'C': ['B', 'F'],
    'D': ['A', 'E', 'G'],
    'E': ['B', 'D', 'F', 'H'],
    'F': ['C', 'E', 'I'],
    'G': ['D', 'H'],
    'H': ['G', 'E', 'I', 'J'],
    'I': ['F', 'H', 'K'],
    'J': ['H', 'L'],
    'K': ['I', 'L'],
    'L': ['J', 'K']
}

coords = {
    'A': (0, 6),
    'B': (2, 7),
    'C': (4, 6),
    'D': (0, 4.5),
    'E': (2, 5),
    'F': (4, 4),
    'G': (0.01, 2),
    'H': (2, 3),
    'I': (4, 2),
    'J': (2, 1),
    'K': (4, 0),
    'L': (6, 1)
}
LISTE_CLIENTS = list(Graphe.keys())
noms = LISTE_CLIENTS  # Variable globale utile pour l'affichage


# =====================================================
# 2. PRÉPARATION DU GRAPHE (Dijkstra & Utils)
# =====================================================


def distance_euclidienne(p1, p2):
    return np.hypot(p1[0] - p2[0], p1[1] - p2[1])


def creer_graphe_pondere(Graphe, coords):
    """Transforme le graphe simple en graphe avec poids (distances réelles)"""
    G_pondere = {}
    for u in Graphe:
        G_pondere[u] = []
        for v in Graphe[u]:
            dist = distance_euclidienne(coords[u], coords[v])
            G_pondere[u].append((v, dist))
    return G_pondere


def cle_min(d, L):
    """ renvoie la clé du minimum d'un dictionnaire"""
    cmin = None
    vmin = float("inf")
    for c in L:
        if d[c] < vmin:
            vmin = d[c]
            cmin = c
    return cmin


def dijkstra_vers_tous(G_pondere, depart):
    """Calcule la distance routière min du 'depart' vers tous les autres sommets
        Renvoie également le dictionnaire des parents du graphe pondere"""
    d = {x: float("inf") for x in G_pondere}
    parent = {}
    d[depart] = 0
    L = list(G_pondere.keys())


    while L:
        u = cle_min(d, L)
        if u is None: break
        L.remove(u)


        for v, w in G_pondere[u]:
            if d[u] + w < d[v]:
                d[v] = d[u] + w
                parent[v] = u
    return d, parent


def Reconstruire_Chemin(DictionnaireParents, Depart, Arrivee):
    """Renvoie le chemin du point de départ au point d'arrivée à partir du dicoParents """
    if Arrivee not in DictionnaireParents and Arrivee != Depart:
        return None
    Chemin_inverse = []
    Ville_actuelle = Arrivee
    while Ville_actuelle != Depart:
        Chemin_inverse.append(Ville_actuelle)
        if Ville_actuelle not in DictionnaireParents:
            return None
        Ville_actuelle = DictionnaireParents[Ville_actuelle]
    Chemin_inverse.append(Depart)
    Chemin_consruit = Chemin_inverse[::-1]
    # Note: On garde le chemin complet A -> ... -> Arrivée
    return Chemin_consruit


def construire_matrice_distances_completes(Graphe, coords, clients):
    """Construit la matrice NxN où A[i][j] est la distance routière réelle"""
    G_pondere = creer_graphe_pondere(Graphe, coords)
    n = len(clients)
    Matrice = np.zeros((n, n))


    for i in range(n):
        dists, _ = dijkstra_vers_tous(G_pondere, clients[i])
        for j in range(n):
            Matrice[i][j] = dists[clients[j]]
    return Matrice

def creer_matrice_euclidienne_adjacente(Graphe, coords, clients):
    """Crée une matrice NxN où:
    - distance euclidienne si connectés
    - np.inf sinon
    """
    n = len(clients)
    M = np.zeros((n, n))

    for i in range(n):
        for j in range(n):
            client_i = clients[i]
            client_j = clients[j]



            if client_j in Graphe[client_i] :
                M[i][j] = distance_euclidienne(coords[client_i], coords[client_j])
            else:
                M[i][j] = np.inf

    return M

def creer_dict_voisins_indices(Graphe, clients):
    """Crée un dictionnaire: indice -> liste des indices des voisins directs"""
    voisins = {}
    for i, client_i in enumerate(clients):
        voisins[i] = []
        for client_j in Graphe[client_i]:
            j = clients.index(client_j)
            voisins[i].append(j)
    return voisins

def cout(Chemin,n):
    """Chemin: une liste montrant le trajet des fourmis
    n: Indice de la derniere ville visitée """
    m=len(Chemin)
    if m<n:
        return None
    return m+1-n

# =====================================================
# 3. MOTEUR ACO (Outils communs)
# =====================================================


def PositionsInitiales(N, start_index=0):
    return [[start_index] for _ in range(N)]


def ActualisePheromones(L, A, B, n, N, r, mode):
    """
    mode=1 si cycle normal
    mode=2 si chemin (ACO Hybride)
    L: Liste de flottants
    A: Matrice d'adjacence du graphe (float)
    B: Matrice des phéromones associée au graphe (float)
    n: Nombre de sommets
    N: Nombre de fourmis
    r: Constante de calcul
    CORRIGÉ : Gère les chemins de longueurs variables (complets ou partiels)
    """
    # Évaporation
    for i in range(n):
        for j in range(n):
            B[i][j] *= (1 - r)


    # Dépôt
    for f in range(N):
        chemin = L[f]
        taille = len(chemin)


        # Sécurité pour éviter crash sur chemin vide
        if taille < 2: continue


        S = 0
        # Calcul de la longueur du chemin parcouru par la fourmi
        for i in range(taille - 1):
            S += A[chemin[i]][chemin[i + 1]]


        # Si boucle fermée (dernier == premier), on ajoute l'arête de retour si nécessaire
        if mode==1:
            S += A[chemin[-1]][chemin[0]]


        if S > 0:
            for i in range(taille - 1):
                u = chemin[i]
                v = chemin[i + 1]
                B[u][v] += 1.0 / S
                B[v][u] += 1.0 / S
            # Renforcement du retour (optionnel)
            if mode==1:
                B[chemin[-1]][chemin[0]] += 1.0/S
                B[chemin[0]][chemin[-1]] += 1.0/S


    return B


# =====================================================
# 4. ACO CLASSIQUE (Tournée complète)
# =====================================================
def TourComplet(A, B, M_euc_adj, voisins_indices, n, N, alpha, beta, depart_fixe=0):
    """Construit N parcours complets en respectant les voisins directs
    A: Matrice d'adjacence du graphe (float)
    B: Matrice des phéromones associée au graphe (float)
    M_euc_adj: Matrice euclidienne + adjacence
    voisins_indices: Dictionnaire des voisins directs par indice
    alpha,beta: Constantes de calcul
    n: Nombre de sommets
    N: Nombre de fourmis
    """
    L = []
    for _ in range(N):
        chemin = []
        non_visites = list(range(n))
        courant = depart_fixe
        chemin.append(courant)
        if courant in non_visites:
            non_visites.remove(courant)

        while non_visites:
            proba = []
            proba_candidates = []
            denom = 0

            for j in list(non_visites):
                # Vérifier si j est un voisin DIRECT de courant
                if j not in voisins_indices[courant]:
                    continue

                tau = B[courant][j] ** alpha
                dist = M_euc_adj[courant][j]

                if dist==np.inf:
                    non_visites.remove(j)
                    continue

                if dist == 0:
                    dist = 0.001

                eta = (1 / dist) ** beta
                val = tau * eta

                proba_candidates.append(j)
                proba.append(val)
                denom += val

            if denom == 0:
                voisins_revisitables = voisins_indices[courant]  # Tous les voisins (même visités)
                if voisins_revisitables:
                    suivant = rd.choice(voisins_revisitables)
                else:
                    # Aucun voisin du tout
                    break
            else:
                # Il y a des voisins nouveaux accessibles
                proba = np.array(proba)
                proba = proba / denom
                suivant = rd.choices(proba_candidates, weights=proba)[0]

            chemin.append(suivant)
            if suivant in non_visites:
                non_visites.remove(suivant)
            courant = suivant

        L.append(chemin)
    return L

def CheminLePlusCourt(A, M_euc_adj, voisins_indices, N, NV, r, alpha, beta):
    """Construit N parcours complets en respectant les voisins directs
    A: Matrice d'adjacence du graphe (float)
    M_euc_adj: Matrice euclidienne + adjacence
    voisins_indices: Dictionnaire des voisins directs par indice
    r: Constante de calcul
    alpha,beta: Constantes de calcul
    NV: Nombre de vagues de fourmis
    N: Nombre de fourmis
    """
    n = len(A)
    B = [[0.01 for _ in range(n)] for _ in range(n)]
    lmin = np.inf
    best_path_indices = []

    for _ in range(NV):
        L = TourComplet(A, B, M_euc_adj, voisins_indices, n, N, alpha, beta)
        B = ActualisePheromones(L, A, B, n, N, r, 1)

        for f in range(N):
            chemin = L[f]
            sommets_visites=set(chemin)
            print(f"Fourmi {f}: chemin_len={len(chemin)}, sommets={len(sommets_visites)}/{n}")
            if len(sommets_visites)<n:
                continue  # Ignorer les chemins incomplets

            # Vérifier que le retour (chemin[-1] -> chemin[0]) est une arête DIRECTE valide
            if chemin[-1] not in voisins_indices or chemin[0] not in voisins_indices[chemin[-1]]:
                # print(f"  ⚠️ Chemin invalide: {noms[chemin[-1]]} -> {noms[chemin[0]]} n'existe pas!")
                continue  # Chemin invalide, sauter

            S = 0
            for i in range(len(chemin) - 1):
                S += A[chemin[i]][chemin[i + 1]]
            S += A[chemin[-1]][chemin[0]] # Retour départ (validé)

            if S < lmin:
                lmin = S
                best_path_indices = list(chemin)

    # Conversion indices -> noms
    if best_path_indices:
        bestpathpoints = [noms[i] for i in best_path_indices]
        bestpathpoints.append(noms[best_path_indices[0]])
    else:
        bestpathpoints = []

    return lmin, bestpathpoints

def CheminLePlusCourt2(A,G_pondere, M_euc_adj, voisins_indices, N, NV, r, alpha, beta):
    """ACO Hybride - Construit N parcours de retour
    A: Matrice d'adjacence du graphe (float)
    M_euc_adj: Matrice euclidienne + adjacence
    voisins_indices: Dictionnaire des voisins directs par indice
    r: Constante de calcul
    alpha,beta: Constantes de calcul
    NV: Nombre de vagues de fourmis
    N: Nombre de fourmis
    ville_depart_index: Index de la ville de départ pour le retour
    Tabu_list_indices: Liste des indices à ne pas revisiter
    """
    n = len(A)
    B = [[0.01 for _ in range(n)] for _ in range(n)]

    lmin = np.inf
    best_aller_indices = []
    best_retour_indices = []

    for _ in range(NV):
        # === PHASE 1: ACO ALLER ===
        # Chaque fourmi construit un chemin sans revisites
        L = TourComplet(A, B, M_euc_adj, voisins_indices, n, N, alpha, beta)
        B = ActualisePheromones(L, A, B, n, N, r, 1)

        # === PHASE 2: DIJKSTRA RETOUR ===
        for f in range(N):
            chemin_aller = L[f]
            sommets_visites = set(chemin_aller)

            # Filtrer 1: chemin complet?
            if len(sommets_visites) < n:
                continue

            # Filtrer 2: première et dernière arête valides?
            depart = chemin_aller[0]   # A
            arrivee = chemin_aller[-1] # Dernière ville

            # Calcul distance aller
            dist_aller = 0
            for i in range(len(chemin_aller) - 1):
                dist_aller += A[chemin_aller[i]][chemin_aller[i + 1]]

            # Calcul retour via Dijkstra (de arrivee à depart)
            distances,parents = dijkstra_vers_tous(G_pondere, noms[arrivee])

            # Reconstruire le chemin retour
            chemin_retour_noms = Reconstruire_Chemin(parents, noms[arrivee],noms[depart])

            if not chemin_retour_noms:  # Pas de chemin!
                continue
            chemin_retour=[noms.index(nom) for nom in chemin_retour_noms]

            # Calcul distance retour
            dist_retour = distances[noms[depart]]

            # Distance totale
            dist_total = dist_aller + dist_retour

            if dist_total < lmin:
                lmin = dist_total
                best_aller_indices = list(chemin_aller)
                best_retour_indices = list(chemin_retour)

    # Conversion indices → noms
    if best_aller_indices:
        aller_points = [noms[i] for i in best_aller_indices]
        retour_points = [noms[i] for i in best_retour_indices]
    else:
        aller_points = []
        retour_points = []

    return lmin, aller_points, retour_points

# =====================================================
# 6. MODULE KRUSKAL (Borne Inférieure)
# =====================================================


def trouver_parent(parent, i):
    if parent[i] == i:
        return i
    parent[i] = trouver_parent(parent, parent[i])
    return parent[i]


def unir_ensembles(parent, rang, i, j):
    racine_i = trouver_parent(parent, i)
    racine_j = trouver_parent(parent, j)
    if racine_i != racine_j:
        if rang[racine_i] < rang[racine_j]:
            parent[racine_i] = racine_j
        elif rang[racine_i] > rang[racine_j]:
            parent[racine_j] = racine_i
        else:
            parent[racine_j] = racine_i
            rang[racine_i] += 1
        return True
    return False


def kruskal(Matrice_Distances, noms_clients):
    """Permet de déterminer la borne inférieure du parcours d'un graphe à partir de  la matrice des distances
    Ici noms_clients est en réalité la liste nom définie au début"""
    n = len(Matrice_Distances)
    aretes_possibles = []
    for i in range(n):
        for j in range(i + 1, n):
            dist = Matrice_Distances[i][j]
            if dist > 0 and dist != np.inf:
                aretes_possibles.append((i, j, dist))
    aretes_possibles.sort(key=lambda x: x[2])
    parent = {i: i for i in range(n)}
    rang = {i: 0 for i in range(n)}
    mst_poids = 0
    mst_aretes = []
    compteur_aretes = 0
    for u, v, poids in aretes_possibles:
        if unir_ensembles(parent, rang, u, v):
            mst_poids += poids
            mst_aretes.append((noms_clients[u], noms_clients[v]))
            compteur_aretes += 1
            if compteur_aretes == n - 1:
                break
    return mst_poids, mst_aretes

# =====================================================
# 7. EXÉCUTION PRINCIPALE (MAIN)
# =====================================================
def main(mode):
    print("--- 1. Calcul de la matrice des distances (Dijkstra) ---")


    A = construire_matrice_distances_completes(Graphe, coords, LISTE_CLIENTS)
    M_euc_adj = creer_matrice_euclidienne_adjacente(Graphe, coords, LISTE_CLIENTS)
    voisins_indices = creer_dict_voisins_indices(Graphe, LISTE_CLIENTS)
    G_pondere = creer_graphe_pondere(Graphe, coords)

    Distance_Kruskal = kruskal(A, LISTE_CLIENTS)[0]
    print(f" -> Distance Kruskal : {Distance_Kruskal:.2f}")

    # === ACO CLASSIQUE ===
    if mode == 1:
        print("--- 2. Lancement ACO CLASSIQUE ---")
        dist_classique, chemin_classique = CheminLePlusCourt(
            A, M_euc_adj, voisins_indices, 20, 50, 0.5, 1, 1
        )
        print(f" -> Distance Classique : {dist_classique:.2f}")
        print(f" -> Chemin : {' -> '.join(chemin_classique)}")

        # Visualisation mode 1
        plt.figure(figsize=(10, 8))
        # Graphe fond
        for u, voisins in Graphe.items():
            for v in voisins:
                plt.plot([coords[u][0], coords[v][0]], [coords[u][1], coords[v][1]], 'k-', alpha=0.1)

        # Trajet classique
        for k in range(len(chemin_classique)-1):
            u, v = chemin_classique[k], chemin_classique[k+1]
            plt.plot([coords[u][0], coords[v][0]], [coords[u][1], coords[v][1]],
                    'g:', linewidth=2, alpha=0.6, label='ACO Classique' if k==0 else "")

        plt.title(f"Trajet Classique ({dist_classique:.2f})")

    # === ACO HYBRIDE ===
    elif mode == 2:
        print("--- 4. Lancement ACO Hybride ---")
        #  CORRIGER: Ajouter G_pondere
        dist_hybride, aller, retour = CheminLePlusCourt2(A, G_pondere, M_euc_adj, voisins_indices, 20, 50, 0.5, 1, 1)


        print(f" -> Distance Hybride : {dist_hybride:.2f}")
        print(f" -> Chemin Aller: {' -> '.join(aller)}")
        print(f" -> Chemin Retour: {' -> '.join(retour)}")
        aller_def=aller
        aller_def.pop()
        traj_complet=aller_def+retour
        print(f" -> Chemin Complet: {' -> '.join(traj_complet)}")
        # Visualisation mode 2
        plt.figure(figsize=(10, 8))
        # Graphe fond
        for u, voisins in Graphe.items():
            for v in voisins:
                plt.plot([coords[u][0], coords[v][0]], [coords[u][1],coords[v][1]], 'k-', alpha=0.1)

        # Trajet ALLER (ACO)
        for k in range(len(aller) - 1):
            u, v = aller[k], aller[k + 1]
            plt.plot([coords[u][0], coords[v][0]], [coords[u][1], coords[v][1]],'g:', linewidth=2, alpha=0.6, label='Aller (ACO)' if k==0 else "")

        # Trajet RETOUR (Dijkstra)
        for k in range(len(retour) - 1):
            u, v = retour[k], retour[k + 1]
            plt.plot([coords[u][0], coords[v][0]], [coords[u][1], coords[v][1]],'b:', linewidth=2, alpha=0.6, label='Retour (Dijkstra)' if k==0 else "")

        plt.title(f"Trajet Hybride ({dist_hybride:.2f})")

    # Villes (pour les 2 modes)
    for sommet, (x, y) in coords.items():
        c = "orange" if sommet == 'A' else "white"
        plt.scatter(x, y, s=300, c=c, edgecolors="black", zorder=10)
        plt.text(x, y, sommet, fontsize=10, weight="bold", ha='center', va='center', zorder=11)

    plt.legend()
    plt.axis("equal")
    plt.show()


if __name__ == '__main__':
    mode = int(input("Mode (1=Classique, 2=Hybride): "))
    main(mode)
