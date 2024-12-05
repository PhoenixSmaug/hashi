# A Hashiwokakero solver using Integer Programming

### Introduction

[Hashiwokakero](https://en.wikipedia.org/wiki/Hashiwokakero) (or short Hashi) is a Japanese logic puzzle, where the aim is to find bridges on a grid of islands with the following rules:

- The number of bridges connected to each island must match the number on that island.
- Bridges can only run horizontally or vertically and cannot cross each other.
- Each bridge connects two islands, and there can be one or two bridges between any pair of islands.
- All islands must be connected in a single network, meaning every island must be reachable from every other island.

The problem was proven to be NP-complete using a reduction of Hamiltonian cycles in unit distance graphs [[1]](https://doi.org/10.1016/j.ipl.2009.07.017). It has already been encoded as a Integer Programming problem before [[2]](https://arxiv.org/abs/1905.00973), but the formulation there uses $O(2^n)$ inequalities to encode the connectivity of the $n$ islands and therefore requires a stepwise solution process with dynamically generated constraints. Here, a new model is presented that requires only $O(n^2)$ constraints in total. A Julia implementation using the JuMP framework is also provided, allowing the model to be solved either directly using the ILP solver [Gurobi](https://www.gurobi.com/) or solving the exported MPS file using SAT-based solvers such as [Exact](https://gitlab.com/nonfiction-software/exact) and Google's [OR-Tools](https://developers.google.com/optimization).

### ILP Model

(c) Mia Müßig
