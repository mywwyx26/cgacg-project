import sys
import matplotlib.pyplot as plt
from read_inputs import read_files

def seed_boundaries(n, parts):
    # split range(n) into contiguous pieces whose lengths differ by at most 1
    base, extra = divmod(n, parts)
    bounds, start = [], 0
    for i in range(parts):
        length = base + (1 if i < extra else 0)
        bounds.append((start, start + length))
        start += length
    return bounds

def build_index(t, k):
    index = {}
    for i in range(len(t) - k + 1):
        kmer = t[i:i + k]
        if kmer not in index:
            index[kmer] = [i]
        else:
            index[kmer].append(i)
    return index

def best_match(p, window, max_ed):
    # semi-global alignment: all of p must align, but it may start/end anywhere in window
    n, m = len(p), len(window)
    prev = [0] * (m + 1)
    for i in range(1, n + 1):
        curr = [i] + [0] * m
        row_min = i
        for j in range(1, m + 1):
            cost = 0 if p[i - 1] == window[j - 1] else 1
            curr[j] = min(prev[j] + 1, curr[j - 1] + 1, prev[j - 1] + cost)
            if curr[j] < row_min:
                row_min = curr[j]
        if row_min > max_ed:
            return None
        prev = curr
    dist = min(prev)
    if dist > max_ed:
        return None
    return dist, prev.index(dist)

def seed_and_extend(max_ed, p, t, indexes=None):
    # indexes: optional dict {k: k-mer index of t}, pass for every miRNA on the same UTR so each index is only made once
    if indexes is None:
        indexes = {}
    hits = {}

    for off, end in seed_boundaries(len(p), max_ed + 1):
        k = end - off
        if k not in indexes:
            indexes[k] = build_index(t, k)
        seed = p[off:end]
        for idx in indexes[k].get(seed, []):
            lo = max(0, idx - off - max_ed)
            hi = min(len(t), idx - off + len(p) + max_ed)
            result = best_match(p, t[lo:hi], max_ed)
            if result:
                dist, e = result
                pos = lo + e
                if pos not in hits or dist < hits[pos]:
                    hits[pos] = dist
    return dict(sorted(hits.items()))

def analyze_6mers(p, t, output, plot_output, histogram_output, coverage_output):
    # p: {name: (site, mature)}, site = reverse complement of Seed+m8
    # first 6 of the seed (miRNA positions 2-7) pair with site[1:]
    index = build_index(t, 6)
    rows = []          # (miRNA, 1-based offset)
    counts = {}        # miRNA -> number of matches
    matched = set()    # sites (= families) with at least one match
    coverage = [0] * len(t)   # position -> number of miRNAs whose 6-mer starts there
    for name, (site, mature) in p.items():
        positions = index.get(site[1:], [])
        counts[name] = len(positions)
        if positions:
            matched.add(site)
        for pos in positions:
            rows.append((name, pos + 1))
            coverage[pos] += 1
    offsets = [offset for name, offset in rows]

    with open(output, 'w') as file:
        file.write(f'Human miRNAs: {len(p)}\n\n')
        file.write(f'Human miRNA families: {len({site for site, mature in p.values()})}\n\n')
        file.write(f'Human miRNA families with a 6-mer match in the 3\' UTR: {len(matched)}\n\n')
        file.write('| miRNA | offset |\n')
        file.write('| --- | --- |\n')
        for name, offset in rows:
            file.write(f'| {name} | {offset} |\n')

    # matches along the 3' UTR
    fig, ax = plt.subplots(figsize=(14, 2.5))
    ax.hlines(0, 1, len(t), color='black', linewidth=2)
    ax.scatter(offsets, [0] * len(offsets), marker='|', s=400, color='tab:blue', alpha=0.3)
    ax.set_xlim(1, len(t))
    ax.set_ylim(-1, 1)
    ax.set_yticks([])
    ax.set_xlabel('position in 3\' UTR (nt)')
    fig.tight_layout()
    fig.savefig(plot_output, dpi=150)
    plt.close(fig)

    # number of miRNAs with 0, 1, 2, ... matches
    most = max(counts.values())
    num_mirnas = [0] * (most + 1)
    for c in counts.values():
        num_mirnas[c] += 1
    fig, ax = plt.subplots(figsize=(8, 5))
    bars = ax.bar(range(most + 1), num_mirnas, color='tab:blue')
    ax.bar_label(bars)
    ax.set_xticks(range(most + 1))
    ax.set_xlabel('number of matches in 3\' UTR')
    ax.set_ylabel('number of miRNAs')
    fig.tight_layout()
    fig.savefig(histogram_output, dpi=150)
    plt.close(fig)

    # number of miRNAs matching at each position
    fig, ax = plt.subplots(figsize=(14, 4))
    ax.bar(range(1, len(t) + 1), coverage, width=1, color='tab:blue')
    ax.set_xlim(1, len(t))
    ax.set_xlabel('position in 3\' UTR (nt)')
    ax.set_ylabel('number of miRNAs with a 6-mer match')
    fig.tight_layout()
    fig.savefig(coverage_output, dpi=150)
    plt.close(fig)

if __name__ == '__main__':
    if len(sys.argv) != 4:
        sys.exit('usage: python seed_and_extend.py [max edit distance] [miRNAs file] [sequence file]')

    try:
        max_ed = int(sys.argv[1])
    except ValueError:
        sys.exit(f'max edit distance must be an integer, got {sys.argv[1]!r}')
    
    try:
        max_ed, p, t = read_files(sys.argv[2], sys.argv[3], max_ed)
    except (ValueError, OSError) as e:
        sys.exit(f'Error: {e}')

    analyze_6mers(p, t, 'seed.md', 'plot.png', 'histogram.png', 'coverage.png')

    indexes = {}
    for name, (site, mature) in p.items():
        hits = seed_and_extend(max_ed, site, t, indexes)
        for end, dist in hits.items():
            print(f'{name}\tmatch ending at t[{end}] (exclusive)\tedit distance {dist}')

    print(len(p), 'length p')
    print(len(t), 'length t')