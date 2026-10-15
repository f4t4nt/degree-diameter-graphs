"""Does any pair of connectors co-hold two ports in a known 104?

The tight equality matching of shapes106.py / shapes107.py (matching_ok exact=True) assumes
(c) "two connectors co-hold at most one port". Hub-to-port exactness does not imply it (the
second shared port is a block-j port, and the port lemma at a block-j port only concerns
foreign blocks), so the assumption is tested here on the 35 corpus classes and the new class 36.

Usage: check_coholding.py > runs/coholding104.json
"""
import collections
import json
import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "template_calculus"))
import verify  # noqa: E402


def coholding(edges):
    verts, adj, problems = verify.build(edges)
    assert not problems
    blocks = verify.k44_blocks(verts, adj)
    decs = verify.decompositions(blocks, len(verts))
    assert len(decs) == 1
    port_block = {}
    for b, key in enumerate(decs[0]):
        L, R, pmap = blocks[key]
        for q in pmap.values():
            port_block[q] = b
    conns = [x for x in verts if x not in port_block and all(y in port_block or y not in set().union(*[set(blocks[k][0]) | set(blocks[k][1]) for k in decs[0]]) for y in adj[x])]
    hubs = set().union(*[set(blocks[k][0]) | set(blocks[k][1]) for k in decs[0]])
    conns = [x for x in verts if x not in port_block and x not in hubs]
    holders = {q: [c for c in adj[q] if c in set(conns)] for q in port_block}
    pair_count = collections.Counter()
    for q, hs in holders.items():
        assert len(hs) == 2
        pair_count[frozenset(hs)] += 1
    multi = {tuple(sorted(k)): v for k, v in pair_count.items() if v > 1}
    profile_of = {c: tuple(sum(1 for y in adj[c] if port_block.get(y) == b) for b in range(3)) for c in conns}
    detail = []
    for (c1, c2), v in multi.items():
        shared = [q for q in port_block if set(holders[q]) == {c1, c2}]
        detail.append({"connectors": [c1, c2], "profiles": [profile_of[c1], profile_of[c2]], "shared_ports": shared,
                       "blocks_of_shared": [port_block[q] for q in shared], "adjacent": c2 in adj[c1]})
    return {"connector_pairs_coholding_2_or_more": len(multi), "detail": detail}


def main():
    corpus = verify.load_corpus()
    out = []
    for i, rep in enumerate(corpus["representatives"]):
        r = coholding([tuple(e) for e in rep["edges_sorted"]])
        r["class"] = i
        out.append(r)
    total = sum(1 for r in out if r["connector_pairs_coholding_2_or_more"])
    print(f"classes with a connector pair co-holding >= 2 ports: {total} of {len(out)}", file=sys.stderr)
    for r in out:
        if r["connector_pairs_coholding_2_or_more"]:
            print(json.dumps({k: r[k] for k in ("class", "connector_pairs_coholding_2_or_more", "detail")}), file=sys.stderr)
    json.dump(out, sys.stdout, indent=1)


if __name__ == "__main__":
    main()
