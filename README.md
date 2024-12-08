# A highly efficient Hashiwokakero solver

### Introduction

[Hashiwokakero](https://en.wikipedia.org/wiki/Hashiwokakero) (or short Hashi) is a Japanese logic puzzle, where the aim is to find bridges on a grid of islands with the following rules:

- The number of bridges connected to each island must match the number on that island.
- Bridges can only run horizontally or vertically and cannot cross each other.
- Each bridge connects two islands, and there can be one or two bridges between any pair of islands.
- All islands must be connected in a single network, meaning every island must be reachable from every other island.

The problem was proven to be NP-complete using a reduction of Hamiltonian cycles in unit distance graphs [[1]](https://doi.org/10.1016/j.ipl.2009.07.017). It has already been encoded as a Integer Programming problem before [[2]](https://arxiv.org/abs/1905.00973), but the formulation there uses $O(2^n)$ inequalities to encode the connectivity of the $n$ islands and therefore requires a stepwise solution process with dynamically generated constraints. Here, a new model is presented that requires only $O(n^2)$ constraints in total. A Julia implementation using the JuMP framework is also provided, allowing the model to be solved either directly using the ILP solver [Gurobi](https://www.gurobi.com/) or as an the exported MPS file with SAT-based solvers such as [Exact](https://gitlab.com/nonfiction-software/exact) and Google's [OR-Tools](https://developers.google.com/optimization). In addition, the combination with lazy constraints enables the solver to outperform currently existing solvers.

### ILP Model

We start by collecting all possible edges, so any direct orthogonal connection between two islands, in the set $E$ with $m = |E|$. Since each island can have at most 4 edges to other islands, $m = O(n)$. We now introduce the binary variables $x \in \\{0, 1\\}^{m \times 2}$:
```math
x_{e, k} = \begin{cases}
1, & \text{if at least $k$ bridges are drawn on edge $e$}\\ 0, & \text{otherwise}
\end{cases}
```

By our definition we know $x_{e, 2} \implies x_{e, 1}$, which we encode as $\\left\\{\forall e \in E \\, | \\, x_{e, 2} \leq x_{e, 1}\\right\\}$. And for all pairs $e, e^{\prime} \in E$ which intersect, we use $x_{e, 1} + x_{e^{\prime}, 1} \leq 1$ to enforce that only on one of the two edges bridges can be built.

To encode that the number of bridges connecting to island $v$ is equal to its number $k$, we simply sum over its adjacent edges $N_v$:
```math
\left\{\forall v \in I \, | \, \sum_{e \in N_v} x_{e, 1} + \sum_{e \in N_v} x_{e, 2} = k\right\}
```

This leaves only the most difficult constraint, namely that all islands must be connected. The first trick is to encode the constraint that all bridges are connected instead, which is equivalent since each island is connected to at least one bridge by our previous constraints.

For this we translate the BFS search procedure on edges into the model using the binary variable $y \in \\{0, 1\\}^{m \times m}$:
```math
y_{e, t} = \begin{cases}
1, & \text{if edge $e$ can be reached from source edge $s$ in at most $t - 1$ steps}\\ 0, & \text{otherwise}
\end{cases}
```

Source edges are the edges reachable in zero steps, so we ensure with $\sum_{e \in E} y_{e, 1} = 1$ that only one source edge exists. We also know that if an edge $e$ can be reached from $s$ in at most $t - 1$ steps, then it can also be reached in at most $t$ steps, which is reflected by:
```math
\{\forall e \in E \\, \forall t \in \{1, \ldots, m - 1\} \, | \, y_{e, t} \leq y_{e, t+1}\}
```

Inversily, if an edge $e = (v, w)$ can be reached from $s$ in at most $t$ steps, then either itself or one of its adjacent edges $\in N_v \cup N_w$ can be reached in at most $t - 1$ steps:
```math
\left\{\forall e \in E \\, \forall t \in \{1, \ldots, m - 1\} \, | \, y_{e, t + 1} \leq y_{e, t} + \sum_{e^{\prime} \in N_v, e^{\prime} \neq e} y_{e^{\prime}, t} + \sum_{e^{\prime} \in N_w, e^{\prime} \neq e} y_{e^{\prime}, t}\right\}
```

We have thus completed the BFS encoding and can now use $\\left\\{\forall e \in E \\, | \\, x_{e, 1} = y_{e, m}\\right\\}$ to ensure that all edges are reachable from the starting edge in at most $m - 1$ steps. Since with $m$ edges no shortest path can be longer than $m - 1$ steps, this is equivalent to all bridges being connected. This concludes the encoding of the Hashiwokakero puzzle, which is implemented as `solve(file::String)` in the Julia code.

### Lazy Constraints

In practice, however, the polynomial size connectivity encoding is slower than the dynamically constructed lazy connectivity constraints presented in [[2]](https://arxiv.org/abs/1905.00973). Altough we can obtain a highly efficient Hashiwokakero solver if we combine the based model presented here with lazy connectivity constraints, which is implemented as `solveLazy(file::String)` in the Julia code. To demonstrate that, we compare the average runtime of our lazy-solve with the runtime of the CLLV model from [[2]](https://arxiv.org/abs/1905.00973) on the common [dataset](https://w1.cirrelt.ca/~vidalt/resources/Hashi_Puzzles.zip) they provided.

<div align="center">

| **Method**       | **n = 100** | **n = 200** | **n = 300** | **n = 400** |
|:-----------------:|:-----------:|:-----------:|:-----------:|:-----------:|
| **CLLV**         |    0.16s    |    1.13s    |    6.32s    |   36.94s    |
| **lazy-solve**   |   <0.01s    |    0.04s    |    0.12s    |    1.48s    |

</div>

As the table shows, the new model represents a drastic improvement in performance for all instance sizes.

(c) Mia Müßig
