# Tools

`ddgraph.py` turns a graph's `construction.json` into its edge list and runs the verifiers in `verify/` on it. It
needs only the Python standard library; `check` compiles the verifiers with `cc` into `verify/bin/` the first time.

    python3 tools/ddgraph.py edges graphs/deg9/dia8 > deg9-dia8.edges.txt
    python3 tools/ddgraph.py check graphs/deg9/dia8
    python3 tools/ddgraph.py check --all --max-vertices 1000000000

## construction.json

- `degree`, `diameter`, `order`: the pair $`(d, k)`$ and the number of vertices.
- `construction.type`:
  - `cayley`: `group` and `connection_set`, a list of group elements $`S`$. Every $`g`$ is joined to $`g s`$ and
    $`g s^{-1}`$; an involution is listed once.
  - `lift`: `group`, `base_vertices`, and lists of `loops` ($`[v, a]`$), `links` ($`[u, v, a]`$) and `semi_edges`
    ($`[v, t]`$). Vertex $`(u, g)`$ is joined to $`(v, g a)`$ for every dart $`u \to v`$ with voltage $`a`$; a
    semi-edge carries an involution $`t`$ and joins $`(v, g)`$ to $`(v, g t)`$.
  - `edge_list`: the graph's edge file is the definition (the two hexagon graphs and the $`(4, 4)`$ graph).
- `group.type`:
  - `metacyclic`: `N`, `M`, `r` for
    $`\mathbb{Z}_N \rtimes_r \mathbb{Z}_M = \langle x, y \mid x^N = y^M = 1,\ y x y^{-1} = x^r \rangle`$. Elements are
    words such as `x^4406 y^18`, `y` or `1`, and $`x^i y^j`$ has index $`i + N j`$.
  - `kernel_by_cyclic`: the non-metacyclic group of the $`(5, 7)`$ graph,
    $`((\mathbb{F}_p^2 \rtimes_B \langle c \rangle) \times \langle z \rangle) \rtimes_A \langle y \rangle`$, with
    elements written `v(a,b) c^i z^j y^m`.

Vertex numbers are group indices: $`g`$ in a Cayley graph, $`v \cdot |G| + \operatorname{index}(g)`$ in a lift.

## Verifiers

The C programs in `verify/` share no code with the searches or with `ddgraph.py`'s group arithmetic. Each one checks
the group law it relies on ($`r`$ has order exactly $`M`$ modulo $`N`$) and that the generators give a simple graph
of the stated degree before searching.

| program | input | what it checks |
|---|---|---|
| `verify_bfs.c` | an edge list (lines starting with `#` are skipped, so it reads the files in `graphs/` directly) | BFS from every vertex: order, minimum and maximum degree, repeated edges, diameter, pairs farther apart than $`D`$ |
| `verify_cayley.c` | an edge list and $`(N, M, r, S)`$ | the edge set equals the Cayley graph exactly, then BFS from the identity and 50 sampled vertices |
| `verify_cayley_implicit.c` | $`(N, M, r, S)`$ | BFS from the identity and sampled vertices over the implicitly defined Cayley graph, 5 bytes per vertex |
| `verify_cayley_inv.c` | the same | the same, with one involution allowed in $`S`$ |
| `verify_lift2_implicit.c` | a 2-vertex lift | BFS from $`(0, e)`$, $`(1, e)`$ and sampled vertices |
| `verify_bits_implicit.c`, `verify_bits_inv.c` | a Cayley graph or 2-vertex lift | bitmap BFS, 3 bits per vertex, up to about $`4 \cdot 10^9`$ vertices |
| `verify_2plane_implicit.c`, `verify_2plane_inv.c` | the same | two bit planes, 2 bits per vertex, up to about $`3 \cdot 10^{10}`$ vertices |

The `_inv` variants differ from their partners only in allowing an involution in $`S`$ (it then contributes one
neighbour instead of two).

Why one search is enough: left multiplication by $`h`$ maps the Cayley graph to itself and the identity to $`h`$, so
every vertex has the eccentricity of the identity. In a voltage lift, left multiplication maps each fibre
$`\{v\} \times G`$ to itself transitively, so the eccentricities of $`(0, e)`$, $`(1, e)`$, ... are all there is to
compute.
