"""Map every known 104-vertex class (../classes/) to its skeleton among the 54 skeletons
(../sat_closure/data/skeletons54.json, the output of unit_designs2.py) by nauty certificate of the
unit-port-hub incidence graph built directly from the known graph (no label conventions involved).
Run from this directory."""
import collections
import json
import sys

import pynauty

import surgery
import unit_designs2 as U2


def certificate_from_template(T):
    units = T.pairs
    port_ids = sorted(T.port_block)
    hub_ids = sorted(T.hub_block)
    idx = {}
    for i, u in enumerate(units):
        idx[("u", u)] = i
    for j, q in enumerate(port_ids):
        idx[("p", q)] = 16 + j
    for k, h in enumerate(hub_ids):
        idx[("h", h)] = 64 + k
    adj = collections.defaultdict(set)
    for q in port_ids:
        l, r = T.port_hubs[q]
        for h in (l, r):
            adj[idx[("p", q)]].add(idx[("h", h)])
            adj[idx[("h", h)]].add(idx[("p", q)])
    for u in units:
        for b in range(3):
            for q in T.pair_ports[u][b]:
                adj[idx[("u", u)]].add(idx[("p", q)])
                adj[idx[("p", q)]].add(idx[("u", u)])
    g = pynauty.Graph(88, adjacency_dict={k: sorted(v) for k, v in adj.items()},
                      vertex_coloring=[set(range(16)), set(range(16, 64)), set(range(64, 88))])
    return pynauty.certificate(g)


def main():
    data = json.load(open("../sat_closure/data/skeletons54.json"))
    cert_to_class = {}
    for i, c in enumerate(data["classes"]):
        cert, _ = U2.incidence_certificate([tuple(t) for t in c["design"]])
        cert_to_class[cert] = i
    assert len(cert_to_class) == len(data["classes"])
    out = []
    for name, verts, adj, trip, blocks in surgery.load_classes():
        T = surgery.Template(verts, adj, trip, blocks)
        cert = certificate_from_template(T)
        out.append({"class": name, "skeleton": cert_to_class.get(cert)})
        print(name, "-> skeleton", cert_to_class.get(cert), file=sys.stderr, flush=True)
    print("skeletons used by known classes:", sorted(collections.Counter(o["skeleton"] for o in out).items()),
          file=sys.stderr)
    json.dump(out, sys.stdout, indent=1)


if __name__ == "__main__":
    main()
