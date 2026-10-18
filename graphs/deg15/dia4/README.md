# Degree 15, Diameter 4

A graph with 13,608 vertices, maximum degree 15 and diameter 4; the Moore bound is 44,326.

## Construction

A route-chart lift: a controller graph $`C`$ on 168 vertices whose darts carry the symbols $`a`$, $`a^{-1}`$, $`b`$, $`b^{-1}`$ and $`s`$, a fibre $`\mathbb{F}_3^4`$ of 81 points over every controller vertex, and for every symbol a chart $`(M, d)`$, a $`4 \times 4`$ matrix over $`\mathbb{F}_3`$ and a direction. The lift has vertex set $`V(C) \times \mathbb{F}_3^4`$, and $`(v, x)`$ is joined to $`(t(v),\ M x + \lambda d)`$ for every symbol $`t`$, with the chart of $`t`$, and every $`\lambda \in \mathbb{F}_3`$: the three neighbours through $`t`$ form one line of the fibre over $`t(v)`$, so the degree is $`5 \cdot 3 = 15`$ and the order is $`168 \cdot 81 = 13608`$.

The controller is the regular $`\mathbb{Z}_{42}`$-lift of the 4-vertex Schreier base with the permutations $`\alpha_0 = \beta_0 = (0\,1\,2\,3)`$ and the matching $`\sigma_0 = (0\,2)(1\,3)`$: its vertices are the pairs $`(c, g)`$ with $`c \in \{0, 1, 2, 3\}`$ and $`g \in \mathbb{Z}_{42}`$, $`a`$ sends $`(c, g)`$ to $`(\alpha_0(c),\ g + v_a(c))`$ with $`v_a = (0, 9, 12, 0)`$, $`b`$ sends it to $`(\beta_0(c),\ g + v_b(c))`$ with $`v_b = (19, 38, 7, 41)`$, and $`s`$ swaps $`(0, g)`$ with $`(2, g)`$ and $`(1, g)`$ with $`(3, g + 21)`$, so that $`s`$ is its own inverse. It is 5-regular on 168 vertices and has diameter 4; the symbol $`s`$ alone is a perfect matching of the controller.

The charts, with $`M`$ acting on column vectors and the 81 points of a fibre written as $`x = (x_0, x_1, x_2, x_3)`$, are

```math
M_{a} = \begin{pmatrix} 1 & 1 & 2 & 2 \\ 0 & 1 & 0 & 1 \\ 0 & 0 & 0 & 2 \\ 0 & 1 & 2 & 0 \end{pmatrix}, \quad d_{a} = \begin{pmatrix} 0 \\ 2 \\ 0 \\ 0 \end{pmatrix}, \qquad
M_{a^{-1}} = \begin{pmatrix} 1 & 0 & 2 & 2 \\ 0 & 1 & 1 & 0 \\ 0 & 1 & 1 & 2 \\ 0 & 0 & 2 & 0 \end{pmatrix}, \quad d_{a^{-1}} = \begin{pmatrix} 0 \\ 2 \\ 2 \\ 0 \end{pmatrix},
```

```math
M_{b} = \begin{pmatrix} 1 & 2 & 0 & 1 \\ 1 & 1 & 0 & 0 \\ 2 & 1 & 1 & 2 \\ 0 & 2 & 1 & 1 \end{pmatrix}, \quad d_{b} = \begin{pmatrix} 2 \\ 0 \\ 0 \\ 1 \end{pmatrix}, \qquad
M_{b^{-1}} = \begin{pmatrix} 2 & 0 & 1 & 2 \\ 1 & 1 & 2 & 1 \\ 1 & 0 & 1 & 0 \\ 0 & 1 & 1 & 2 \end{pmatrix}, \quad d_{b^{-1}} = \begin{pmatrix} 0 \\ 0 \\ 2 \\ 2 \end{pmatrix},
```

```math
M_{s} = \begin{pmatrix} 1 & 1 & 0 & 2 \\ 0 & 0 & 0 & 1 \\ 0 & 1 & 1 & 2 \\ 0 & 1 & 0 & 0 \end{pmatrix}, \quad d_{s} = \begin{pmatrix} 0 \\ 1 \\ 2 \\ 1 \end{pmatrix}.
```

The charts of $`a^{-1}`$ and $`b^{-1}`$ describe the same edges as those of $`a`$ and $`b`$ read from the other end, and the chart of the self-inverse symbol $`s`$ satisfies $`M_s^2 = I`$ and $`M_s d_s = d_s`$, so its edges are symmetric too; every edge is listed once in the edge file, under the symbol that gives it.

![The 4-vertex Schreier base of the controller with its voltages](controller-base.svg)<br><sub>Every dart $`c \to t(c)`$ of the base carries a symbol $`t`$ and a voltage $`g`$, and the controller joins $`(c, h)`$ to $`(t(c), h + g)`$ for every $`h \in \mathbb{Z}_{42}`$. The dart of the inverse symbol runs the other way with the opposite voltage.</sub>

Vertex $`81 v + x_0 + 3 x_1 + 9 x_2 + 27 x_3`$ is the pair $`(v, x)`$ with $`v = 42 c + g`$ the controller vertex $`(c, g)`$ and $`x = (x_0, x_1, x_2, x_3)`$; the fibre over $`v`$ is the run of 81 numbers from $`81 v`$.

## How it was found

Route-chart lift over $`\mathbb{F}_3^4`$ with the five-symbol alphabet $`a^{\pm 1}, b^{\pm 1}, s`$ ($`M_s^2 = I`$, $`M_s d_s = d_s`$) over the regular $`\mathbb{Z}_{42}`$-lift of the 4-vertex base $`\alpha_0 = \beta_0 = (0\,1\,2\,3)`$, $`\sigma_0 = (0\,2)(1\,3)`$ with voltages $`v_a = (0,9,12,0)`$, $`v_b = (19,38,7,41)`$, $`v_s = 21`$ on the link $`1 - 3`$ (168 vertices, 5-regular); chart and controller co-designed: a pool of sampled controllers in which every pair at distance 1, 3 or 4 has a non-backtracking 4-walk, then chart annealing against the pool until one controller has a controllable word for every such pair and the chart-only conditions at distance 0 and 2 hold. The search code is in [search/route_chart_lifts](../../../search/route_chart_lifts).

## Files

- `deg15-dia4-n13608.edges.txt`: the 102,060 edges as `u v` lines, under a `#` heading per symbol.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

A breadth-first search from every vertex (`tools/verify/verify_bfs.c`) confirms the order, the degrees, that there are no repeated edges, and the diameter; the edge list is also checked against `construction.json`.

    python3 tools/ddgraph.py check graphs/deg15/dia4
