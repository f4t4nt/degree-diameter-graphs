# Degree-Diameter Graphs

My graphs for the degree-diameter problem: 50 of them have more vertices than the entries of the
[table of largest known graphs](https://web.mat.upc.edu/francesc.comellas/delta-d/taula_delta_d.html) on October 14, 2026, and two more, (3,7) and (4,4), are my entries
already in it.

The problem asks for the largest graph with maximum degree d and diameter k. Every graph here has a directory
`graphs/deg<d>-dia<k>/` with a README on its construction, the construction in `construction.json`, the edge list
when the graph has fewer than a million vertices, and the output of an independent verifier in `verification.txt`.

| (d, k) | vertices | construction |
|---|---:|---|
| [(3,7)](graphs/deg3-dia7) | 200 | 5-vertex lift over Z_5 x\|_4 Z_8 |
| [(4,4)](graphs/deg4-dia4) | 104 | three subdivided K_4,4 blocks and 32 connectors |
| [(4,6)](graphs/deg4-dia6) | 760 | generalized hexagon H(3) with 32 flags subdivided |
| [(4,9)](graphs/deg4-dia9) | 7,884 | 3-vertex lift over Z_73 x\|_18 Z_36 |
| [(4,10)](graphs/deg4-dia10) | 19,306 | 2-vertex lift over Z_197 x\|_16 Z_49 |
| [(5,6)](graphs/deg5-dia6) | 2,836 | generalized hexagon H(4) with 106 flags subdivided |
| [(5,7)](graphs/deg5-dia7) | 6,000 | 2-vertex lift over a non-metacyclic group of order 3000 |
| [(5,8)](graphs/deg5-dia8) | 18,810 | 2-vertex lift over Z_209 x\|_4 Z_45 |
| [(5,9)](graphs/deg5-dia9) | 60,686 | 2-vertex lift over Z_1597 x\|_3 Z_19 |
| [(5,10)](graphs/deg5-dia10) | 202,550 | 2-vertex lift over Z_4051 x\|_4 Z_25 |
| [(6,8)](graphs/deg6-dia8) | 80,050 | Cayley graph of Z_3202 x\|_1419 Z_25 |
| [(6,10)](graphs/deg6-dia10) | 1,394,189 | Cayley graph of Z_32423 x\|_1693 Z_43 |
| [(7,6)](graphs/deg7-dia6) | 12,654 | Cayley graph of Z_703 x\|_85 Z_18 |
| [(7,8)](graphs/deg7-dia8) | 264,440 | Cayley graph of Z_6611 x\|_193 Z_40 |
| [(7,9)](graphs/deg7-dia9) | 1,332,774 | Cayley graph of Z_24681 x\|_788 Z_54 |
| [(7,10)](graphs/deg7-dia10) | 6,641,680 | Cayley graph of Z_83021 x\|_4400 Z_80 |
| [(8,7)](graphs/deg8-dia7) | 138,050 | Cayley graph of Z_2761 x\|_63 Z_50 |
| [(8,8)](graphs/deg8-dia8) | 779,803 | Cayley graph of Z_13217 x\|_8875 Z_59 |
| [(8,9)](graphs/deg8-dia9) | 4,471,111 | Cayley graph of Z_66733 x\|_1441 Z_67 |
| [(8,10)](graphs/deg8-dia10) | 26,864,290 | Cayley graph of Z_141391 x\|_240 Z_190 |
| [(9,5)](graphs/deg9-dia5) | 9,020 | Cayley graph of Z_451 x\|_46 Z_20 |
| [(9,7)](graphs/deg9-dia7) | 302,082 | Cayley graph of Z_4577 x\|_174 Z_66 |
| [(9,8)](graphs/deg9-dia8) | 1,961,370 | Cayley graph of Z_21793 x\|_34 Z_90 |
| [(9,9)](graphs/deg9-dia9) | 13,219,318 | Cayley graph of Z_134891 x\|_4746 Z_98 |
| [(9,10)](graphs/deg9-dia10) | 89,543,818 | Cayley graph of Z_491999 x\|_477175 Z_182 |
| [(10,7)](graphs/deg10-dia7) | 616,915 | Cayley graph of Z_9491 x\|_86 Z_65 |
| [(10,8)](graphs/deg10-dia8) | 4,507,500 | Cayley graph of Z_45075 x\|_2 Z_100 |
| [(10,9)](graphs/deg10-dia9) | 34,182,150 | Cayley graph of Z_227881 x\|_1562 Z_150 |
| [(10,10)](graphs/deg10-dia10) | 265,612,889 | Cayley graph of Z_1111351 x\|_8147 Z_239 |
| [(11,7)](graphs/deg11-dia7) | 1,194,278 | Cayley graph of Z_20591 x\|_2507 Z_58 |
| [(11,8)](graphs/deg11-dia8) | 9,739,368 | Cayley graph of Z_135269 x\|_8594 Z_72 |
| [(11,9)](graphs/deg11-dia9) | 81,961,120 | Cayley graph of Z_512257 x\|_706 Z_160 |
| [(11,10)](graphs/deg11-dia10) | 697,280,200 | Cayley graph of Z_3486401 x\|_1834706 Z_200 |
| [(12,7)](graphs/deg12-dia7) | 2,114,927 | Cayley graph of Z_19403 x\|_6308 Z_109 |
| [(12,8)](graphs/deg12-dia8) | 19,144,241 | Cayley graph of Z_90731 x\|_61666 Z_211 |
| [(12,9)](graphs/deg12-dia9) | 177,775,277 | Cayley graph of Z_783151 x\|_1808 Z_227 |
| [(12,10)](graphs/deg12-dia10) | 1,692,658,087 | Cayley graph of Z_3655849 x\|_5568 Z_463 |
| [(13,7)](graphs/deg13-dia7) | 3,691,038 | Cayley graph of Z_47321 x\|_3981 Z_78 |
| [(13,8)](graphs/deg13-dia8) | 36,275,160 | Cayley graph of Z_302293 x\|_7818 Z_120 |
| [(13,9)](graphs/deg13-dia9) | 362,492,640 | Cayley graph of Z_1258655 x\|_263 Z_288 |
| [(13,10)](graphs/deg13-dia10) | 3,789,298,110 | Cayley graph of Z_9716149 x\|_514348 Z_390 |
| [(14,8)](graphs/deg14-dia8) | 64,298,753 | Cayley graph of Z_185299 x\|_47110 Z_347 |
| [(14,9)](graphs/deg14-dia9) | 709,125,173 | Cayley graph of Z_2104229 x\|_864 Z_337 |
| [(14,10)](graphs/deg14-dia10) | 8,052,187,457 | Cayley graph of Z_10837399 x\|_3457 Z_743 |
| [(15,7)](graphs/deg15-dia7) | 9,758,988 | Cayley graph of Z_90361 x\|_1443 Z_108 |
| [(15,8)](graphs/deg15-dia8) | 110,115,980 | Cayley graph of Z_423523 x\|_8 Z_260 |
| [(15,9)](graphs/deg15-dia9) | 1,315,107,690 | Cayley graph of Z_3372071 x\|_832 Z_390 |
| [(15,10)](graphs/deg15-dia10) | 15,810,645,560 | Cayley graph of Z_20803481 x\|_5726682 Z_760 |
| [(16,7)](graphs/deg16-dia7) | 14,995,367 | Cayley graph of Z_83773 x\|_72211 Z_179 |
| [(16,8)](graphs/deg16-dia8) | 185,394,793 | Cayley graph of Z_769273 x\|_381264 Z_241 |
| [(16,9)](graphs/deg16-dia9) | 2,343,053,633 | Cayley graph of Z_5060591 x\|_949 Z_463 |
| [(16,10)](graphs/deg16-dia10) | 29,899,503,497 | Cayley graph of Z_24170981 x\|_16573517 Z_1237 |

## Checking the graphs

Python 3 and a C compiler are all it takes:

    python3 tools/ddgraph.py check graphs/deg5-dia8     # one cell
    python3 tools/ddgraph.py check --all                # every cell up to 10^8 vertices
    python3 tools/ddgraph.py edges graphs/deg9-dia8 > deg9-dia8.edges.txt

`check` confirms that a stored edge list is exactly its construction and runs the breadth-first-search verifiers in
`tools/verify/`: from every vertex for the smaller graphs, and from one vertex per orbit for the Cayley graphs and
voltage lifts, whose groups act transitively on the vertices or on each fibre. The largest graphs need a lot of
memory for that search, so `check` skips cells above 10^8 vertices unless `--max-vertices` is raised; their
`verification.txt` holds the output of the original run.

## Layout

- `graphs/`: one directory per cell.
- `tools/ddgraph.py`: groups, constructions, edge lists and the `check` driver (see [tools/README.md](tools/README.md)).
- `tools/verify/`: the C verifiers.
- `search/`: the code that found the graphs.

## Notation

- `Z_N x|_r Z_M` is the metacyclic group <x, y | x^N = y^M = 1, y x y^-1 = x^r>; its element x^i y^j is vertex
  i + N j.
- The Cayley graph of G with connection set S joins g to g s and g s^-1 for every g in G and s in S.
- A voltage lift of a base graph B over G has vertex set V(B) x G and joins (u, g) to (v, g a) for every dart
  u -> v of B with voltage a.

Disclaimer: the searches and the code were developed with Claude Fable (Anthropic) under my direction.
