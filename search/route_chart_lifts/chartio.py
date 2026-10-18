#!/usr/bin/env python3
"""Convert charts between cdsa text format (A row-major, d) and alph json; verify a cdsa result exactly (rclift).
Usage: chartio.py tojson chart.txt out.json | totxt in.json out.txt | verify chart.txt k a0 b0 s0 group.tbl"""
import sys, json, re, subprocess, os
import alph, rcl, ctrl
HERE = os.path.dirname(os.path.abspath(__file__))

def read_txt(fn):
    A = []; d = []; ctrl_line = None
    for line in open(fn):
        if line.startswith('A'):
            nums = list(map(int, re.findall(r'-?\d+', line)))
            rows = [nums[4 * i:4 * i + 4] for i in range(4)]
            A.append([[rows[i][j] for i in range(4)] for j in range(4)])   # columns
            d.append(nums[16:20])
        elif line.startswith('controller'):
            ctrl_line = line[len('controller '):].strip()
    al = alph.Alphabet(3, 4, [1, 0, 3, 2, 4], A, d)
    return al, ctrl_line

def write_txt(al, fn):
    with open(fn, 'w') as f:
        for t in range(al.K):
            rows = [[al.A[t][j][i] for j in range(4)] for i in range(4)]
            f.write('A ' + ' '.join(str(x) for r in rows for x in r) + ' d ' + ' '.join(map(str, al.d[t])) + '\n')

def parse_ctrl(line):
    volts = {}
    for tok in line.split():
        key, val = tok.split('=')
        volts[key] = val if key == 'base' else list(map(int, val.split(',')))
    return volts

if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'tojson':
        al, cl = read_txt(sys.argv[2]); al.check()
        json.dump(dict(al.to_json(), controller=parse_ctrl(cl) if cl else None, score=al.score()), open(sys.argv[3], 'w'))
        print('score', al.score())
    elif cmd == 'totxt':
        al = alph.Alphabet.from_json(json.load(open(sys.argv[2]))); write_txt(al, sys.argv[3])
    elif cmd == 'verify':
        al, cl = read_txt(sys.argv[2]); al.check()
        k = int(sys.argv[3]); a0, b0, s0, grp = sys.argv[4:8]
        volts = parse_ctrl(cl)
        print('chart score', al.score(), 'controller', volts)
        T = ctrl.read_table(grp)
        base = ctrl.Base(k, list(map(int, a0.split(','))), list(map(int, b0.split(','))), list(map(int, s0.split(','))))
        N, out = ctrl.lift_controller(base, T, volts['va'], volts['vb'], volts['vs'])
        print('controller', ctrl.controller_stats(N, out))
        fn = sys.argv[2] + '.rcl'
        al.write_rcl(fn, N, out)
        res = subprocess.run([os.path.join(HERE, 'rclift'), fn, '--edges', fn + '.edges'], capture_output=True, text=True).stdout
        print(res)
        if 'ordered pairs at distance > 4: 0\n' in res:
            n = N * 81
            edges = sorted((min(int(a), int(b)), max(int(a), int(b))) for a, b in (l.split() for l in open(fn + '.edges')))
            outfn = 'verify_%d.edges' % n
            with open(outfn, 'w') as f:
                for a, b in edges: f.write('%d %d\n' % (a, b))
            print('DIAMETER 4 LIFT: n=%d written %s' % (n, outfn))
