# A Hashiwokakero solver using Integer Programming

### Introduction

[Hashiwokakero](https://en.wikipedia.org/wiki/Hashiwokakero) (or short Hashi) is a Japanese logic puzzle, where the aim is to find bridges on a grid of islands with the following rules:

- The number of bridges connected to each island must match the number on that island.
- Bridges can only run horizontally or vertically and cannot cross each other.
- Each bridge connects two islands, and there can be one or two bridges between any pair of islands.
- All islands must be connected in a single network, meaning every island must be reachable from every other island.

The problem was proven to be NP-complete using a reduction of Hamiltonian cycles in unit distance graphs [[1]](https://doi.org/10.1016/j.ipl.2009.07.017). It has already been encoded as a Integer Programming problem before [[2]](https://arxiv.org/abs/1905.00973), but the formulation there uses $O(2^n)$ inequalities to encode the connectivity of the $n$ islands and therefore requires a stepwise solution process with dynamically generated constraints. Here, a new model is presented that requires only $O(n^2)$ constraints in total. A Julia implementation using the JuMP framework is also provided, allowing the model to be solved either directly using the ILP solver [Gurobi](https://www.gurobi.com/) or as an the exported MPS file with SAT-based solvers such as [Exact](https://gitlab.com/nonfiction-software/exact) and Google's [OR-Tools](https://developers.google.com/optimization).

### ILP Model

We start by collecting all possible edges, so any direct orthogonal connection between two islands, in the set $E$ with $m = |E|$. Since each island can have at most 4 edges to other islands, $m = O(n)$. We now introduce the binary variables $x \in \\{0, 1\\}^{m \times 2}$:
```math
x_{e, k} = \begin{cases}
1, & \text{if at least $k$ bridges are drawn on edge $e$}\\ 0, & \text{otherwise}
\end{cases}
```

By our definition we know $x_{e, 2} \implies x_{e, 1}$, which we encode as $\\{\forall e \in E \\, | \\, x_{e, 2} \leq x_{e, 1}\\}$. And for all pairs $e, e^{\prime} \in E$ which intersect, we use $x_{e, 1} + x_{e^{\prime}, 1} \leq 1$ to enforce that only on one of the two edges bridges can be built.

To encode that the number of bridges connecting to island $v$ is equal to its number $k$, we simply sum over its adjacent edges $N_v$:
```math
\sum_{e \in N_v} x_{e, 1} + \sum_{e \in N_v} x_{e, 2} = k
```

This leaves only the most difficult constraint, which is that all islands must be connected to each other.


(c) Mia Müßig
