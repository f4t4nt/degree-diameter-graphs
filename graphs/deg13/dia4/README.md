# Degree 13, Diameter 4

A graph with 6,804 vertices, maximum degree 13 and diameter 4; the Moore bound is 24,506.

## Construction

A route-chart lift: a controller graph $`C`$ on 84 vertices whose darts carry the symbols $`a`$, $`a^{-1}`$, $`b`$, $`b^{-1}`$ and $`z`$, a fibre $`\mathbb{F}_3^4`$ of 81 points over every controller vertex, and for every symbol a chart $`(M, d)`$, a $`4 \times 4`$ matrix over $`\mathbb{F}_3`$ and a direction. The lift has vertex set $`V(C) \times \mathbb{F}_3^4`$, and $`(v, x)`$ is joined to $`(t(v),\ M x + \lambda d)`$ for every symbol $`t`$, with the chart of $`t`$, and every $`\lambda \in \mathbb{F}_3`$: the three neighbours through $`t`$ form one line of the fibre over $`t(v)`$. The plain symbol $`z`$ gives one neighbour, $`(z(v),\ -x + c)`$, so the degree is $`4 \cdot 3 + 1 = 13`$ and the order is $`84 \cdot 81 = 6804`$.

The controller is the Cayley graph of $`G = \mathbb{Z}_7 \rtimes_2 \mathbb{Z}_{12} = \langle x, y \mid x^7 = y^{12} = 1,\ y x y^{-1} = x^2 \rangle`$ (order 84) with the symbols $`a = x y`$ (order 12), $`b = y^4`$ (order 3), their inverses, and the central involution $`z = y^6`$: the symbol $`t`$ sends $`g`$ to $`g t`$. It is 5-regular on 84 vertices.

The charts, with $`M`$ acting on column vectors and the 81 points of a fibre written as $`x = (x_0, x_1, x_2, x_3)`$, are

```math
M_{a} = \begin{pmatrix} 0 & 1 & 0 & 1 \\ 2 & 0 & 2 & 0 \\ 2 & 1 & 1 & 1 \\ 1 & 0 & 1 & 2 \end{pmatrix}, \quad d_{a} = \begin{pmatrix} 0 \\ 0 \\ 1 \\ 1 \end{pmatrix}, \qquad
M_{a^{-1}} = \begin{pmatrix} 2 & 1 & 1 & 0 \\ 1 & 1 & 0 & 1 \\ 1 & 1 & 2 & 0 \\ 0 & 2 & 0 & 2 \end{pmatrix}, \quad d_{a^{-1}} = \begin{pmatrix} 1 \\ 1 \\ 2 \\ 2 \end{pmatrix},
```

```math
M_{b} = \begin{pmatrix} 0 & 0 & 2 & 0 \\ 2 & 1 & 2 & 2 \\ 2 & 1 & 2 & 0 \\ 0 & 1 & 2 & 2 \end{pmatrix}, \quad d_{b} = \begin{pmatrix} 1 \\ 2 \\ 2 \\ 2 \end{pmatrix}, \qquad
M_{b^{-1}} = \begin{pmatrix} 0 & 2 & 0 & 1 \\ 2 & 2 & 1 & 1 \\ 2 & 0 & 0 & 0 \\ 0 & 2 & 1 & 0 \end{pmatrix}, \quad d_{b^{-1}} = \begin{pmatrix} 0 \\ 1 \\ 2 \\ 0 \end{pmatrix},
```

```math
M_{z} = -I, \quad c_{z} = \begin{pmatrix} 0 \\ 0 \\ 1 \\ 0 \end{pmatrix}.
```

The charts of $`a^{-1}`$ and $`b^{-1}`$ describe the same edges as those of $`a`$ and $`b`$ read from the other end; every edge is listed once in the edge file, under the symbol that gives it.

Vertex $`81 v + x_0 + 3 x_1 + 9 x_2 + 27 x_3`$ is the pair $`(v, x)`$ with $`v = i + 7 j`$ the group element $`x^i y^j`$ and $`x = (x_0, x_1, x_2, x_3)`$; the fibre over $`v`$ is the run of 81 numbers from $`81 v`$.

## How it was found

Mixed route-chart lift over $`\mathbb{F}_3^4`$ with four line symbols and one plain central involution $`z`$ (chart $`M_z = -I`$, one neighbour) over the Cayley graph $`\mathrm{Cay}(\mathbb{Z}_7 \rtimes_2 \mathbb{Z}_{12};\ a^{\pm 1}, b^{\pm 1}, z)`$ with $`a = (1,1)`$, $`b = (0,4)`$, $`z = (0,6)`$ (84 vertices); the group and its involution were the only admissible pair (hard deficit 0: every needed fibre has a controllable line 4-word or a tiling) among the catalogue groups of orders 82 to 157, and the chart came from exact-lift annealing. The search code is in [search/route_chart_lifts](../../../search/route_chart_lifts).

## Files

- `deg13-dia4-n6804.edges.txt`: the 44,226 edges as `u v` lines, under a `#` heading per symbol.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

A breadth-first search from every vertex (`tools/verify/verify_bfs.c`) confirms the order, the degrees, that there are no repeated edges, and the diameter; the edge list is also checked against `construction.json`.

    python3 tools/ddgraph.py check graphs/deg13/dia4
