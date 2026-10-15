"""Unit-level structure of the known 104s (prediction of the tightness analysis).

Units are matched connector pairs; each unit owns two ports per block. For block i, G_i is the
multigraph on units whose edges are the ports of block i (each port joins the two units at its slots).
Predictions from hub-to-port tightness (TEMPLATE_CALCULUS.md, section 7 follow-up):
  * G_i is a perfect matching of units with every edge doubled: 8 twin pairs sharing both ports;
  * the two 2-matchings of a twin pair ... are the two ports they share; the 8 twin pairs fall into
    two groups of four, each group using the 8 ports {l1,l2}x{r1,r2} u {l3,l4}x{r3,r4} of one
    L-pairing/R-pairing, the other group the complementary 8 ports;
  * G_j u G_k (j, k != i) has no edge between the two groups.
"""
import collections
import json

import surgery
import verify


def analyze(name, verts, adj, trip, blocks):
    T = surgery.Template(verts, adj, trip, blocks)
    units = T.pairs
    unit_of = {}
    for u in units:
        for c in u:
            unit_of[c] = u
    out = {"class": name}
    port_units = {}
    for q in T.port_block:
        us = [unit_of[c] for c in T.adj[q] if c in unit_of]
        port_units[q] = us
    same_unit_ports = sum(1 for q, us in port_units.items() if len(us) == 2 and us[0] == us[1])
    out["ports_with_both_slots_in_one_unit"] = same_unit_ports
    twins_ok = True
    groups_ok = True
    cross_ok = True
    per_block = []
    for i in range(3):
        G = collections.Counter()
        for q, us in port_units.items():
            if T.port_block[q] == i:
                G[tuple(sorted(us))] += 1
        twin = all(m == 2 for m in G.values()) and len(G) == 8
        twins_ok &= twin
        # groups: twin pairs sharing an L-pairing/R-pairing structure
        twin_ports = {}
        for q, us in port_units.items():
            if T.port_block[q] == i:
                twin_ports.setdefault(tuple(sorted(us)), []).append(T.port_hubs[q])
        # hubs covered by each twin pair's two ports
        info = []
        for key, ports in twin_ports.items():
            info.append({"units": [list(u) for u in key], "ports": ports})
        # cross structure: components of G_j u G_k restricted to units
        others = [b for b in range(3) if b != i]
        adjU = {u: set() for u in units}
        for q, us in port_units.items():
            if T.port_block[q] in others and len(us) == 2 and us[0] != us[1]:
                adjU[us[0]].add(us[1])
                adjU[us[1]].add(us[0])
        seen, comps = set(), []
        for u in units:
            if u in seen:
                continue
            stack, comp = [u], set()
            while stack:
                x = stack.pop()
                if x in comp:
                    continue
                comp.add(x)
                stack.extend(adjU[x] - comp)
            seen |= comp
            comps.append(sorted(len(c) for c in [comp])[0])
        per_block.append({"block": i, "twin_structure": twin, "G_i_edge_multiplicities": sorted(G.values()),
                          "components_of_G_j_union_G_k": sorted(comps)})
        groups_ok &= (sorted(comps) != [16])
    out["twins_in_every_block"] = twins_ok
    out["G_j_union_G_k_disconnected_in_every_block"] = groups_ok
    out["per_block"] = per_block
    return out


def main():
    classes = surgery.load_classes()
    res = [analyze(*c) for c in classes]
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
