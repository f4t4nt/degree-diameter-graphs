"""Port lemma check for the 105 gadget configurations that survived the counting screen.

Port lemma (hub-to-port tightness localised). In a tight member, every hub of block i reaches a
port q of block j != i exactly once along a-p-c-q or a-p-c-c'-q. The connectors that carry q within
one connector step are N[holders(q)] = the two slot connectors of q and their connector neighbours.
Each hub of block i reaches q once per block-i port carried by that set, so the block-i ports of
N[holders(q)] are exactly four distinct ports forming a perfect matching of K4,4(i). A normal pair
(two triples with two ports per block) contributes exactly its two block-i ports when one of its
members holds q.

For a gadget configuration (profiles of the special connectors, their connector edges, stub triples)
this is a finite check on how gadget slots share ports: each slot of a gadget connector is either
shared with a normal connector (contributing 2 to every other block) or with another gadget slot in
the same block. Enumerate all sharing patterns; a configuration survives only if some pattern gives
count exactly 4 for every gadget-held port and every other block. Count 4 is necessary, not
sufficient, so survivors would still need the exact completion search.
"""
import itertools
import json
import sys


def load_survivors():
    r = json.load(open("runs/shapes105_v3.json"))
    return r["survivors"]


def build(conf):
    """Return connectors dict id -> profile tuple, adjacency among gadget connectors."""
    prof = {}
    adj = {}
    units = conf["units"]
    for i, p in enumerate(conf["profiles"]):
        prof[i] = tuple(p)
        adj[i] = set()
    for a, b in conf["edges"]:
        adj[a].add(b)
        adj[b].add(a)
    nxt = len(units)
    for owner, stubs in enumerate(conf["stubs"]):
        for t in stubs:
            prof[nxt] = tuple(t)
            adj[nxt] = {owner}
            adj[owner].add(nxt)
            nxt += 1
    for q in conf.get("quads", []):
        prof[nxt] = tuple(q)
        adj[nxt] = set()
        nxt += 1
    return prof, adj


def slots(prof):
    """List of gadget slots: (connector, block, index within its ports of that block)."""
    out = []
    for c, p in prof.items():
        for b in range(3):
            for k in range(p[b]):
                out.append((c, b, k))
    return out


def matchings(items):
    """All partial matchings of a list (each item unmatched or paired with another)."""
    if not items:
        yield []
        return
    first, rest = items[0], items[1:]
    for m in matchings(rest):
        yield m  # first unmatched
    for i, other in enumerate(rest):
        rest2 = rest[:i] + rest[i + 1:]
        for m in matchings(rest2):
            yield [(first, other)] + m


def check(conf):
    prof, adj = build(conf)
    closed = {c: {c} | adj[c] for c in prof}
    carried = {c: [sum(prof[x][i] for x in closed[c]) for i in range(3)] for c in prof}
    per_block_slots = {b: [s for s in slots(prof) if s[1] == b] for b in range(3)}
    results = {}
    block_ok = []
    for b in range(3):
        ok_patterns = []
        for m in matchings(per_block_slots[b]):
            if any(x[0] == y[0] for x, y in m):
                continue  # a connector cannot hold the same port twice
            matched = {}
            for x, y in m:
                matched[x] = y
                matched[y] = x
            good = True
            for s in per_block_slots[b]:
                c = s[0]
                for i in range(3):
                    if i == b:
                        continue
                    cnt = carried[c][i]
                    if s in matched:
                        other = matched[s][0]
                        cnt = sum(prof[x][i] for x in (closed[c] | closed[other]))
                    else:
                        cnt += 2
                    if cnt != 4:
                        good = False
                        break
                if not good:
                    break
            if good:
                ok_patterns.append(m)
        block_ok.append(ok_patterns)
    results["patterns_passing_per_block"] = [len(p) for p in block_ok]
    if not all(block_ok):
        results["survives"] = False
        results["joint_patterns_passing"] = 0
        return results
    # joint patterns: the carried block-i slots of N[holders(q)] must be pairwise distinct ports
    # among the gadget slots (matched slots are one port); with one normal holder exactly 2 gadget
    # ports, with two gadget holders exactly 4.
    joint_ok = 0
    for combo in itertools.product(*block_ok):
        port_id = {}
        for b, m in enumerate(combo):
            for x, y in m:
                port_id[x] = port_id[y] = ("shared", x, y)
        for s in slots(prof):
            port_id.setdefault(s, ("own", s))
        good = True
        for s in slots(prof):
            c, b, _ = s
            holders = [c]
            pid = port_id[s]
            if pid[0] == "shared":
                holders.append(pid[2][0] if pid[1] == s else pid[1][0])
            group = set()
            for h in holders:
                group |= closed[h]
            for i in range(3):
                if i == b:
                    continue
                ports_i = {port_id[(x, i, k)] for x in group for k in range(prof[x][i])}
                need = 4 - 2 * (2 - len(holders))
                if len(ports_i) != need:
                    good = False
                    break
            if not good:
                break
        if good:
            joint_ok += 1
    results["joint_patterns_passing"] = joint_ok
    results["survives"] = joint_ok > 0
    return results


def main():
    survivors = load_survivors()
    out = []
    n_surv = 0
    for i, conf in enumerate(survivors):
        res = check(conf)
        res["index"] = i
        res["shape"] = conf["shape"]
        res["profiles"] = conf["profiles"]
        res["edges"] = conf["edges"]
        res["stubs"] = conf["stubs"]
        out.append(res)
        n_surv += res["survives"]
        kind = {3: "triangle", 2: "path", 1: "one-edge", 0: "no-edge"}[len(conf["edges"])]
        print(i, kind, conf["profiles"], conf["stubs"], "->", "SURVIVES" if res["survives"] else "dead",
              res["patterns_passing_per_block"], file=sys.stderr)
    print("configurations:", len(survivors), "surviving the port lemma:", n_surv, file=sys.stderr)
    json.dump(out, sys.stdout, indent=1)


if __name__ == "__main__":
    main()
