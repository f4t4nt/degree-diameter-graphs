#!/usr/bin/env python3
"""Rebuild and verify the graphs in this repository.

    python3 tools/ddgraph.py edges graphs/deg5-dia8 > out.edges.txt
    python3 tools/ddgraph.py check graphs/deg5-dia8
    python3 tools/ddgraph.py check --all

`edges` prints the edge list defined by a cell's construction.json: one "u v" per line, vertices
numbered as in the cell's README, grouped into "#" sections by generator or voltage. `check`
compiles the C verifiers in tools/verify/ with `cc` and runs the ones that fit the cell (see
tools/README.md). Only the Python standard library is needed.
"""

import argparse
import gzip
import itertools
import json
import pathlib
import re
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
VERIFY_SRC = ROOT / "tools" / "verify"
VERIFY_BIN = VERIFY_SRC / "bin"

ALL_PAIRS_LIMIT = 2 * 10**10  # n * m budget for the all-pairs BFS of verify_bfs
BYTE_BFS_LIMIT = 2 * 10**8  # largest n for the 5-bytes-per-vertex implicit verifiers
BIT_BFS_LIMIT = 4 * 10**9  # largest n for the 3-bits-per-vertex verifier


class Metacyclic:
    """Z_N x|_r Z_M = <x, y | x^N = y^M = 1, y x y^-1 = x^r>.

    The element x^i y^j is the pair (i, j) and has index i + N j, so
    (i1, j1)(i2, j2) = (i1 + r^j1 i2 mod N, j1 + j2 mod M).
    """

    def __init__(self, N, M, r):
        if pow(r, M, N) != 1 % N:
            raise ValueError(f"r = {r} does not satisfy r^M = 1 mod N for N = {N}, M = {M}")
        self.N, self.M, self.r = N, M, r
        self.order = N * M
        self.identity = (0, 0)
        self.rpow = [pow(r, j, N) for j in range(M)]

    def mul(self, a, b):
        return ((a[0] + self.rpow[a[1]] * b[0]) % self.N, (a[1] + b[1]) % self.M)

    def index(self, a):
        return a[0] + self.N * a[1]

    def elements(self):
        for j in range(self.M):
            for i in range(self.N):
                yield (i, j)

    def generator_power(self, name, k):
        if name == "x":
            return (k % self.N, 0)
        if name == "y":
            return (0, k % self.M)
        raise ValueError(f"unknown generator {name!r} for a metacyclic group")


def _mat_mul(a, b, p):
    return [[sum(a[i][k] * b[k][j] for k in range(2)) % p for j in range(2)] for i in range(2)]


def _mat_pow(a, e, p):
    out = [[1, 0], [0, 1]]
    for _ in range(e):
        out = _mat_mul(out, a, p)
    return out


class KernelByCyclic:
    """G = K x| <y> with K = (F_p^2 x|_B <c>) x <z>.

    c acts on the column vectors F_p^2 by the matrix B, z is central in K, and y acts on K by
    (v, c^i, z^j) -> (A v, c^(cp i), z^(zp j)). The element (v, c^i, z^j) y^m is the tuple
    (v1, v2, i, j, m) with index (((p v1 + v2) |c| + i) |z| + j) |y| + m, and
    (v, i, j, m)(v', i', j', m') = (v + B^i A^m v', i + cp^m i', j + zp^m j', m + m').
    """

    def __init__(self, p, B, A, c_order, c_power, z_order, z_power, y_order):
        self.p, self.c_order, self.z_order, self.y_order = p, c_order, z_order, y_order
        self.c_power, self.z_power = c_power, z_power
        ident = [[1, 0], [0, 1]]
        if _mat_pow(B, c_order, p) != ident or _mat_pow(A, y_order, p) != ident:
            raise ValueError("B must have order dividing |c| and A order dividing |y|")
        if _mat_mul(A, B, p) != _mat_mul(_mat_pow(B, c_power, p), A, p):
            raise ValueError("y does not act on F_p^2 x| <c> by an automorphism (need A B = B^cp A)")
        if pow(c_power, y_order, c_order) != 1 or pow(z_power, y_order, z_order) != 1:
            raise ValueError("the action of y on <c> and <z> must have order dividing |y|")
        self.order = p * p * c_order * z_order * y_order
        self.identity = (0, 0, 0, 0, 0)
        self.ba = [[_mat_mul(_mat_pow(B, i, p), _mat_pow(A, m, p), p) for m in range(y_order)]
                   for i in range(c_order)]

    def mul(self, a, b):
        v1, v2, i1, j1, m1 = a
        w1, w2, i2, j2, m2 = b
        t = self.ba[i1][m1]
        p = self.p
        return ((v1 + t[0][0] * w1 + t[0][1] * w2) % p,
                (v2 + t[1][0] * w1 + t[1][1] * w2) % p,
                (i1 + pow(self.c_power, m1, self.c_order) * i2) % self.c_order,
                (j1 + pow(self.z_power, m1, self.z_order) * j2) % self.z_order,
                (m1 + m2) % self.y_order)

    def index(self, a):
        v1, v2, i, j, m = a
        return (((self.p * v1 + v2) * self.c_order + i) * self.z_order + j) * self.y_order + m

    def elements(self):
        return itertools.product(range(self.p), range(self.p), range(self.c_order),
                                 range(self.z_order), range(self.y_order))

    def generator_power(self, name, k):
        if name == "c":
            return (0, 0, k % self.c_order, 0, 0)
        if name == "z":
            return (0, 0, 0, k % self.z_order, 0)
        if name == "y":
            return (0, 0, 0, 0, k % self.y_order)
        raise ValueError(f"unknown generator {name!r} for this group")


TOKEN = re.compile(r"v\((\d+),(\d+)\)|([a-z])(?:\^(-?\d+))?|1")


def parse_element(group, text):
    """Parse a word such as "x^4406 y^18" or "v(0,2) c z^4 y^4", multiplied left to right."""
    out = group.identity
    pos = 0
    text = text.strip()
    while pos < len(text):
        if text[pos].isspace():
            pos += 1
            continue
        m = TOKEN.match(text, pos)
        if not m:
            raise ValueError(f"cannot parse group element {text!r}")
        if m.group(1) is not None:
            factor = (int(m.group(1)), int(m.group(2)), 0, 0, 0)
        elif m.group(3) is not None:
            factor = group.generator_power(m.group(3), int(m.group(4) or 1))
        else:
            factor = group.identity
        out = group.mul(out, factor)
        pos = m.end()
    return out


def make_group(spec):
    if spec["type"] == "metacyclic":
        return Metacyclic(spec["N"], spec["M"], spec["r"])
    if spec["type"] == "kernel_by_cyclic":
        return KernelByCyclic(spec["p"], spec["B"], spec["A"], spec["c_order"], spec["c_power"],
                              spec["z_order"], spec["z_power"], spec["y_order"])
    raise ValueError(f"unknown group type {spec['type']!r}")


def is_involution(group, a):
    return a != group.identity and group.mul(a, a) == group.identity


def multiply_edges(group, a, src, dst, once):
    """(src, g) -- (dst, g a) for every g in G, vertex (v, g) numbered v |G| + index(g).

    With once=True (a semi-edge or an involution in a Cayley graph) each pair is listed once.
    """
    n = group.order
    for g in group.elements():
        u = group.index(g)
        w = group.index(group.mul(g, a))
        if not once or u < w:
            yield src * n + u, dst * n + w


def sections(spec):
    """[(title, edge iterator)] for an algebraic construction, one section per generator or dart.

    Returns None for a cell that is defined by its edge file.
    """
    c = spec["construction"]
    if c["type"] == "edge_list":
        return None
    group = make_group(c["group"])
    out = []
    if c["type"] == "cayley":
        for s in c["connection_set"]:
            a = parse_element(group, s)
            if is_involution(group, a):
                out.append((f"involution {s}: g -- g {s}, each pair once", multiply_edges(group, a, 0, 0, True)))
            else:
                out.append((f"generator {s}: g -- g {s}", multiply_edges(group, a, 0, 0, False)))
        return out
    if c["type"] == "lift":
        for v, s in c.get("loops", []):
            out.append((f"loop at base vertex {v}, voltage {s}: ({v}, g) -- ({v}, g {s})",
                        multiply_edges(group, parse_element(group, s), v, v, False)))
        for v, s in c.get("semi_edges", []):
            out.append((f"semi-edge at base vertex {v}, voltage {s}: ({v}, g) -- ({v}, g {s}), each pair once",
                        multiply_edges(group, parse_element(group, s), v, v, True)))
        for x, y, s in c.get("links", []):
            out.append((f"link {x} -> {y}, voltage {s}: ({x}, g) -- ({y}, g {s})",
                        multiply_edges(group, parse_element(group, s), x, y, False)))
        return out
    raise ValueError(f"unknown construction type {c['type']!r}")


def load(cell):
    return json.loads((pathlib.Path(cell) / "construction.json").read_text())


def edge_file(cell):
    files = sorted(pathlib.Path(cell).glob("*.edges.txt")) + sorted(pathlib.Path(cell).glob("*.edges.txt.gz"))
    return files[0] if files else None


def edge_pairs(path):
    """The edges of an edge file as integer pairs; "#" lines are comments.

    Vertex names that are not integers (such as A.L0) are numbered in order of first appearance.
    """
    names = {}

    def vid(token):
        return int(token) if token.isdigit() else names.setdefault(token, len(names))

    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as f:
        for line in f:
            if line.strip() and not line.startswith("#"):
                a, b = line.split()[:2]
                yield vid(a), vid(b)


def compiled(name):
    """Compile tools/verify/<name>.c once and return the binary path."""
    VERIFY_BIN.mkdir(exist_ok=True)
    binary = VERIFY_BIN / name
    source = VERIFY_SRC / f"{name}.c"
    if not binary.exists() or binary.stat().st_mtime < source.stat().st_mtime:
        subprocess.run(["cc", "-O2", "-o", str(binary), str(source)], check=True)
    return binary


def run_verifier(name, args, stdin=None):
    """Run a compiled verifier, echo its command and output, and report whether it passed."""
    shown = " ".join([f"tools/verify/{name}"] + [str(a) for a in args])
    print(f"$ {shown}" + (" < the edge list" if stdin is not None else ""))
    sys.stdout.flush()
    proc = subprocess.run([str(compiled(name))] + [str(a) for a in args], input=stdin,
                          stdout=subprocess.PIPE, text=True)
    print(proc.stdout.rstrip())
    return proc.returncode == 0 and "FAIL" not in proc.stdout and "EXCEEDS" not in proc.stdout


def element_ij(group, text):
    a = parse_element(group, text)
    return [a[0], a[1]]


def implicit_args(spec):
    """The verifier and arguments for a metacyclic Cayley graph or 2-vertex lift, else None."""
    c = spec["construction"]
    D, n = spec["diameter"], spec["order"]
    if c["type"] == "edge_list" or c["group"]["type"] != "metacyclic":
        return None
    group = make_group(c["group"])
    g = c["group"]
    if c["type"] == "cayley":
        gens = [v for s in c["connection_set"] for v in element_ij(group, s)]
        if n <= BYTE_BFS_LIMIT:
            return "verify_cayley_inv", [D, g["N"], g["M"], g["r"], 3] + gens
        if n <= BIT_BFS_LIMIT:
            return "verify_bits_inv", ["cay", D, g["N"], g["M"], g["r"]] + gens
        return "verify_2plane_inv", ["cay", D, g["N"], g["M"], g["r"]] + gens
    if c["type"] == "lift" and c["base_vertices"] == 2 and not c.get("semi_edges"):
        loops = c.get("loops", [])
        k0 = [s for v, s in loops if v == 0]
        k1 = [s for v, s in loops if v == 1]
        links = [link[2] for link in c["links"]]
        volts = [v for s in k0 + k1 + links for v in element_ij(group, s)]
        head = [D, g["N"], g["M"], g["r"], len(k0), len(k1), len(links)]
        if n <= BYTE_BFS_LIMIT:
            return "verify_lift2_implicit", head + volts
        return "verify_bits_inv", ["lift"] + head + volts
    return None


def check(cell, max_vertices):
    """Verify one cell and return "ok", "FAILED" or "skipped"."""
    cell = pathlib.Path(cell)
    spec = load(cell)
    deg, D, n = spec["degree"], spec["diameter"], spec["order"]
    print(f"== {cell.name}: degree {deg}, diameter {D}, {n} vertices")
    ok = True
    path = edge_file(cell)
    if path is not None:
        secs = sections(spec)
        lines = 0
        if secs is not None:
            file_edges = set()
            for u, v in edge_pairs(path):
                file_edges.add(min(u, v) * n + max(u, v))
                lines += 1
            built = set(min(u, v) * n + max(u, v) for _, it in secs for u, v in it)
            same = built == file_edges and lines == len(file_edges)
            print(f"{path.name} lists exactly the edges of the construction in construction.json: {same}"
                  f" ({lines} edges)")
            ok = ok and same
        else:
            lines = sum(1 for _ in edge_pairs(path))
        if n * lines <= ALL_PAIRS_LIMIT:
            text = "".join(f"{u} {v}\n" for u, v in edge_pairs(path))
            ok = run_verifier("verify_bfs", ["/dev/stdin", D], text) and ok
            return "ok" if ok else "FAILED"
    implicit = implicit_args(spec)
    if implicit is None:
        print("no verifier applies to this cell")
        return "FAILED"
    if n > max_vertices:
        print(f"skipped here (n > --max-vertices {max_vertices}); to run it:")
        print("$ tools/verify/" + " ".join(str(a) for a in [implicit[0]] + implicit[1]))
        return "skipped" if ok else "FAILED"
    ok = run_verifier(*implicit) and ok
    return "ok" if ok else "FAILED"


def write_edges(spec, out):
    for title, it in sections(spec) or []:
        out.write(f"\n# {title}\n")
        for u, v in it:
            out.write(f"{u} {v}\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="cmd", required=True)
    e = sub.add_parser("edges", help="print the edge list of a cell's construction")
    e.add_argument("cell")
    c = sub.add_parser("check", help="verify cells with the C verifiers")
    c.add_argument("cells", nargs="*")
    c.add_argument("--all", action="store_true", help="every cell under graphs/")
    c.add_argument("--max-vertices", type=int, default=10**8,
                   help="skip implicit BFS runs above this order (default 1e8)")
    args = parser.parse_args()
    if args.cmd == "edges":
        spec = load(args.cell)
        if sections(spec) is None:
            sys.exit("this cell is defined by its edge file, not by a construction")
        write_edges(spec, sys.stdout)
        return
    cells = sorted((ROOT / "graphs").iterdir()) if args.all else [pathlib.Path(p) for p in args.cells]
    results = [(cell.name, check(cell, args.max_vertices)) for cell in cells if (cell / "construction.json").exists()]
    if len(results) > 1:
        print()
        for name, status in results:
            print(f"{name}: {status}")
    sys.exit(1 if any(status == "FAILED" for _, status in results) else 0)


if __name__ == "__main__":
    main()
