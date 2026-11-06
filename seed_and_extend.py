import sys
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

    indexes = {}
    for name, (site, mature) in p.items():
        hits = seed_and_extend(max_ed, site, t, indexes)
        for end, dist in hits.items():
            print(f'{name}\tmatch ending at t[{end}] (exclusive)\tedit distance {dist}')