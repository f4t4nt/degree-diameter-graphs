# Degree 5, Diameter 7

A graph with 6,000 vertices, maximum degree 5 and diameter 7. The Moore bound for degree 5 and diameter 7 is 27,306.

## Construction

A regular voltage lift over the group $`G`$ described below, of order 3,000. The base graph has vertices 0 and 1, with three parallel links $`0 \to 1`$ (voltages $`1`$, $`((0, 1), c, z)\, y^{4}`$, $`((1, 4), c, z^{4})\, y^{7}`$), a loop at 0 (voltage $`y`$) and a loop at 1 (voltage $`((0, 2), c, z^{4})\, y^{4}`$).

The lift has vertex set $`\{0, 1\} \times G`$ and joins $`(u, g)`$ to $`(v, g a)`$ for every dart $`u \to v`$ with voltage $`a`$; a loop with voltage $`a`$ joins $`(v, g)`$ to $`(v, g a)`$ and $`(v, g a^{-1})`$. $`G`$ acts on each fibre $`\{v\} \times G`$ by left multiplication.

The group is $`G = ((\mathbb{F}_{25} \rtimes \langle c \rangle) \times \langle z \rangle) \rtimes \langle y \rangle`$ of order 3000: $`\mathbb{F}_{25} = \mathbb{F}_5^2`$ as column vectors, $`c`$ of order 3 acting on $`\mathbb{F}_{25}`$ by $`B`$ (so $`\mathbb{F}_{25} \rtimes \langle c \rangle`$ is a Frobenius group of order 75), $`z`$ central of order 5 in $`K = (\mathbb{F}_{25} \rtimes \langle c \rangle) \times \langle z \rangle`$, and $`y`$ of order 8 acting on $`K`$ by $`(w, c^i, z^j) \mapsto (A w, c^{2i}, z^{2j})`$, where

```math
B = \begin{pmatrix} 0 & 4 \\ 1 & 4 \end{pmatrix}, \qquad A = \begin{pmatrix} 1 & 1 \\ 2 & 4 \end{pmatrix}.
```

$`G`$ is not metacyclic; it has trivial centre and 25 involutions. Its elements are written $`(w, c^i, z^j)\, y^m`$.

Vertex $`3000 v + \mathrm{index}(g)`$ is the pair $`(v, g)`$, where $`g = ((a, b), c^i, z^j)\, y^m`$ has index $`(((5a + b) \cdot 3 + i) \cdot 5 + j) \cdot 8 + m`$, so vertex $`3000 v`$ is $`(v, e)`$.

## How it was found

Complete sweep of all shape-103 two-vertex lift assignments (one loop and three links per base vertex, one link trivial) on the non-metacyclic group $`((\mathbb{F}_{25} \rtimes C_3) \times \mathbb{Z}_5) \rtimes C_8`$ up to symmetry, the group selected by the conjugacy-class-count (commuting probability) rule. The search code is in [search/lift_nonmetacyclic_groups](../../../search/lift_nonmetacyclic_groups).

## Files

- `deg5-dia7-n6000.edges.txt`: the 15,000 edges, one `u v` per line in the numbering above, grouped by voltage under `#` headings.
- `construction.json`: the construction in the form `tools/ddgraph.py` reads.
- `verification.txt`: the verifier output.

## Verification

A breadth-first search from every vertex (`tools/verify/verify_bfs.c`) confirms the order, the degrees, that there are no repeated edges, and the diameter; the edge list is also checked to be exactly the construction in `construction.json`.

    python3 tools/ddgraph.py check graphs/deg5/dia7
