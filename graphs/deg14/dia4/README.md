# Degree 14, Diameter 4

A graph with 9,720 vertices, maximum degree 14 and diameter 4; the Moore bound is 33,321.

## Construction

A route-chart lift: a controller graph $`C`$ on 120 vertices whose darts carry the symbols $`a`$, $`a^{-1}`$, $`b`$, $`b^{-1}`$, $`p`$ and $`p^{-1}`$, a fibre $`\mathbb{F}_3^4`$ of 81 points over every controller vertex, and for every symbol a chart $`(M, d)`$, a $`4 \times 4`$ matrix over $`\mathbb{F}_3`$ and a direction. The lift has vertex set $`V(C) \times \mathbb{F}_3^4`$, and $`(v, x)`$ is joined to $`(t(v),\ M x + \lambda d)`$ for every symbol $`t`$, with the chart of $`t`$, and every $`\lambda \in \mathbb{F}_3`$: the three neighbours through $`t`$ form one line of the fibre over $`t(v)`$. The plain symbols $`p`$ and $`p^{-1}`$ give one neighbour each, $`(p(v),\ -x + c)`$ and $`(p^{-1}(v),\ -x + c)`$, so the degree is $`4 \cdot 3 + 2 = 14`$ and the order is $`120 \cdot 81 = 9720`$.

The controller is the Cayley graph of $`G = SL(2, 3) \times \mathbb{Z}_5`$ (order 120), whose elements are the pairs $`(X, i)`$ of a $`2 \times 2`$ matrix $`X`$ over $`\mathbb{F}_3`$ of determinant 1 and $`i \in \mathbb{Z}_5`$, with the symbols $`a`$, $`b`$, $`p`$ and their inverses, where every symbol $`t`$ sends $`g`$ to $`g t`$ and the three generators are

```math
a = \left( \begin{pmatrix} 0 & 1 \\ 2 & 2 \end{pmatrix}, 1 \right), \qquad b = \left( \begin{pmatrix} 0 & 1 \\ 2 & 1 \end{pmatrix}, 2 \right), \qquad p = (I, 1).
```

$`a`$ has order 15, $`b`$ has order 30 and $`p`$ is central of order 5; the controller is 6-regular on 120 vertices, and its four line symbols alone give a 4-regular graph of diameter 6.

The charts, with $`M`$ acting on column vectors and the 81 points of a fibre written as $`x = (x_0, x_1, x_2, x_3)`$, are

```math
M_{a} = \begin{pmatrix} 0 & 0 & 2 & 0 \\ 0 & 0 & 2 & 1 \\ 1 & 0 & 2 & 1 \\ 2 & 1 & 2 & 0 \end{pmatrix}, \quad d_{a} = \begin{pmatrix} 1 \\ 1 \\ 0 \\ 0 \end{pmatrix}, \qquad
M_{a^{-1}} = \begin{pmatrix} 0 & 2 & 1 & 0 \\ 2 & 2 & 1 & 1 \\ 2 & 0 & 0 & 0 \\ 2 & 1 & 0 & 0 \end{pmatrix}, \quad d_{a^{-1}} = \begin{pmatrix} 2 \\ 1 \\ 2 \\ 0 \end{pmatrix},
```

```math
M_{b} = \begin{pmatrix} 2 & 0 & 0 & 2 \\ 1 & 2 & 0 & 0 \\ 2 & 2 & 0 & 2 \\ 2 & 1 & 1 & 2 \end{pmatrix}, \quad d_{b} = \begin{pmatrix} 2 \\ 0 \\ 0 \\ 2 \end{pmatrix}, \qquad
M_{b^{-1}} = \begin{pmatrix} 1 & 1 & 2 & 0 \\ 1 & 0 & 2 & 0 \\ 1 & 0 & 1 & 1 \\ 1 & 2 & 1 & 0 \end{pmatrix}, \quad d_{b^{-1}} = \begin{pmatrix} 2 \\ 2 \\ 1 \\ 2 \end{pmatrix},
```

```math
M_{p} = -I, \quad c_{p} = \begin{pmatrix} 2 \\ 0 \\ 2 \\ 2 \end{pmatrix}.
```

The charts of $`a^{-1}`$ and $`b^{-1}`$ describe the same edges as those of $`a`$ and $`b`$ read from the other end; every edge is listed once in the edge file, under the symbol that gives it.

Vertex $`81 v + x_0 + 3 x_1 + 9 x_2 + 27 x_3`$ is the pair $`(v, x)`$ with $`v = 5 k + i`$ the group element $`(X_k, i)`$, where $`X_0 = I`$ and $`X_1, \dots, X_{23}`$ are the other matrices of determinant 1 in lexicographic order of their entries $`(a, b, c, d)`$, and $`x = (x_0, x_1, x_2, x_3)`$; the fibre over $`v`$ is the run of 81 numbers from $`81 v`$.

## How it was found

Mixed route-chart lift over $`\mathbb{F}_3^4`$ with four line symbols and a plain central pair $`p^{\pm 1}`$ of order 5 (chart $`M_p = -I`$) over $`\mathrm{Cay}(SL(2,3) \times \mathbb{Z}_5;\ a^{\pm 1}, b^{\pm 1}, p^{\pm 1})`$ (120 vertices); the controller was selected by the chart-independent hard-deficit screen of the catalogue groups of orders 103 to 193 (two groups reach deficit 0, of orders 105 and 120) and the chart by exact-lift annealing. The search code is in [search/route_chart_lifts](../../../search/route_chart_lifts).

## Files

- `deg14-dia4-n9720.edges.txt`: the 68,040 edges as `u v` lines, under a `#` heading per symbol.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

A breadth-first search from every vertex (`tools/verify/verify_bfs.c`) confirms the order, the degrees, that there are no repeated edges, and the diameter; the edge list is also checked against `construction.json`.

    python3 tools/ddgraph.py check graphs/deg14/dia4
