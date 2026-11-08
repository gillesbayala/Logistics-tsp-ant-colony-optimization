# 🚚 Logistics TSP Optimization under Graph Topologic Constraints (ACO & Hybrid Dijkstra)

[![Python](https://img.shields.io/badge/Python-3.9%2B-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Optimization](https://img.shields.io/badge/Operations%20Research-NP--Hard-red.svg)](#theoretical-framework)

Implementation and comparative analysis of metaheuristic routing algorithms (Ant Colony Optimization, Dynamic ACO, and Hybrid ACO-Dijkstra) applied to logistics delivery networks with non-planar topologies, bottlenecks, and dead-ends.

---

## 📌 Problem Overview

Urban logistics and last-mile distribution face severe physical graph constraints: one-way roads, natural barriers, and cul-de-sacs that prevent simple cyclic returns.

* **Theoretical Formulation:** NP-Hard Traveling Salesperson Problem (TSP) on an incomplete non-oriented graph $G = (V, E, W)$.
* **Core Challenge:** Standard ACO fails on incomplete topology due to ant trapping in terminal dead-ends.
* **Proposed Solution:** A hybrid semi-stochastic, semi-deterministic framework coupling ant pheromone exploration with exact Dijkstra return pathfinding.

---

## ⚙️ Algorithmic Framework

### 1. Classical Ant Colony Optimization (Dorigo et al.)
Probabilistic transition rule for ant $k$ from node $i$ to unvisited neighbor $j \in \mathcal{N}_i$ :
$$P_{ij}^k = \frac{\left[\tau_{ij}\right]^\alpha \cdot \left[\eta_{ij}\right]^\beta}{\sum_{l \in \mathcal{N}_i} \left[\tau_{il}\right]^\alpha \cdot \left[\eta_{il}\right]^\beta}$$

where $\tau_{ij}$ is the pheromone density, $\eta_{ij} = \frac{1}{d_{ij}}$ represents the local visibility heuristic, and $\alpha, \beta$ are weighting hyper-parameters. Pheromone evaporation rule :
$$\tau_{ij}(t+1) = (1 - \rho) \cdot \tau_{ij}(t) + \sum_{k=1}^{m} \Delta \tau_{ij}^k, \quad \Delta \tau_{ij}^k = \frac{Q}{L_k}$$

### 2. Dynamic ACO (Vanishing Load Heuristic)
Adapts greedy sensitivity $\beta(c)$ as delivery parcels are delivered along the cycle :
$$\beta(c) = \beta_{\text{base}} + \gamma \left(\frac{Q - c}{c + 0.01}\right)$$

where $Q$ is the initial vehicle payload capacity, $c$ is the remaining parcels, and $\gamma$ sets the dynamic acceleration factor.

### 3. Hybrid ACO-Dijkstra Architecture
1. **Outward Leg (Stochastic ACO):** Explores delivery destinations using cooperative pheromone guidance until complete set coverage $|V|$ is achieved.
2. **Homeward Leg (Deterministic Dijkstra):** Directly computes the shortest path tree from the terminal node back to the depot $\mathbf{A}$, avoiding loop breakdown in dead-end streets.

---

## 🔬 Benchmark Topologies & Empirical Performance

Tested across diverse complex graph structures :
* **The Snail (Concentric Spiral)**: Severe topology with deep cul-de-sacs.
* **The Bottleneck**: Narrow isthmus linking two dense clusters.
* **Asymmetric Labyrinths**: High-degree dead-ends.
* **Complete Random Graphs ($N = 16$ to $32$)**.

| Benchmark Topology | Classical ACO Completion | Dynamic ACO Completion | Hybrid ACO-Dijkstra Completion |
| :--- | :---: | :---: | :---: |
| **Complete Metric Graph** | 100% | 100% (Faster Convergence) | **100% (Optimal)** |
| **The Snail (`L'Escargot`)** | 0% (Trapped / `inf`) | 0% (Trapped / `inf`) | **100% (Length: 188.00)** |
| **Bottleneck Isthmus** | Trapped frequently | Sensitive to $\gamma$ | **100% (Length: 42.12)** |
| **Asymmetric Maze** | Fails on dead ends | Traps in local optima | **100% (Length: 54.05)** |

---

## 📐 Theoretical Lower Bound (Kruskal's MST)

To quantify solution quality without known NP-Hard global minima, an absolute lower bound is established via Kruskal's Minimum Spanning Tree (MST) :
$$\text{Cost}(\text{TSP}^*) \ge \text{Weight}(\text{MST}(G))$$

Results confirm that the Hybrid ACO converges within a narrow gap ($\approx 15-18\%$) of the theoretical Kruskal lower bound on complex metric instances.

---

## 🚀 Quickstart

```bash
# Clone repository
git clone [https://github.com/gillesbayala/Circuits-logistiques.git](https://github.com/gillesbayala/Circuits-logistiques.git)
cd Circuits-logistiques

# Install dependencies
pip install -r requirements.txt

# Run solver
python tsp_solver.py
