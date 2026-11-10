import re
import sys

VALID = set('ACGU')
COMP = str.maketrans('ACGU', 'UGCA')

def prepare_inputs(max_ed, p, t):
    # p: {name: (seed, mature)}, t: 3' UTR. Normalizes, validates, and reverse complements each seed
    t = t.strip().upper().replace('T', 'U')
    if len(t) == 0:
        raise ValueError('reference (t) must be non-empty')
    bad = set(t) - VALID
    if bad:
        raise ValueError(f'reference (t) contains invalid characters: {sorted(bad)}')
    if max_ed < 0:
        raise ValueError(f'max edit distance ({max_ed}) must be nonnegative')

    prepared = {}
    for name, (seed, mature) in p.items():
        seed = seed.strip().upper().replace('T', 'U')
        mature = mature.strip().upper().replace('T', 'U')
        if len(seed) == 0:
            raise ValueError(f'{name}: query (p) must be non-empty')
        if max_ed >= len(seed):
            raise ValueError(f'{name}: max edit distance ({max_ed}) must be less than len(p) = {len(seed)}')
        if len(seed) > len(t):
            raise ValueError(f'{name}: query length ({len(seed)}) must not exceed reference length ({len(t)})')
        bad = set(seed) - VALID
        if bad:
            raise ValueError(f'{name}: query (p) contains invalid characters: {sorted(bad)}')
        prepared[name] = (seed.translate(COMP)[::-1], mature)  # reverse complement of the seed
    return max_ed, prepared, t

def read_mirnas(mirna):
    p = {}
    with open(mirna, 'r') as file:
        next(file) # skip header line
        for line in file:
            line = line.split()
            if line[2] != '9606': # human only
                continue
            p[line[3]] = (line[1], line[4])
    return p

def read_sequence(sequence):
    seq = []
    in_origin = False
    start = 0 # index of where 3' UTR starts
    with open(sequence, 'r') as file:
        for line in file:
            if line.startswith('//'):
                break
            if in_origin:
                seq.append(re.sub(r'[^A-Za-z]', '', line))
            elif line.startswith('ORIGIN'):
                in_origin = True
            elif line.strip().startswith('CDS '):
                start = int(line.split('..')[1])
    return ''.join(seq)[start:]

def read_files(mirna, sequence, max_ed):
    p = read_mirnas(mirna)
    t = read_sequence(sequence)
    return prepare_inputs(max_ed, p, t)

if __name__ == '__main__':
    if len(sys.argv) != 3:
        sys.exit('usage: python read_inputs.py [miRNAs file] [sequence file]')

    try:
        max_ed, p, t = read_files(sys.argv[1], sys.argv[2], 0)
    except (ValueError, OSError) as e:
        sys.exit(f'Error: {e}')
    print(len(p), len(t))