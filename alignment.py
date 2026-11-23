"""Smith-Waterman local DNA alignment implemented with linear gap penalties."""
from array import array


def clean_dna(sequence):
    sequence = ''.join(sequence.split()).upper()
    if set(sequence) - set('ACGTN'):
        raise ValueError('DNA must contain only A, C, G, T or N.')
    return sequence


def smith_waterman(query, target, match=2, mismatch=-1, gap=-2):
    """Return one optimal alignment. Coordinates are zero-based, half-open.

    Ties resolve diagonal, then up, then left. N always receives mismatch.
    Time O(mn), space O(mn), compact int32 scores and byte traceback.
    """
    query, target = clean_dna(query), clean_dna(target)
    if match <= 0 or mismatch > 0 or gap >= 0:
        raise ValueError('Use a positive match and negative gap; mismatch must be <= 0.')
    width = len(target) + 1
    cells = (len(query) + 1) * width
    if cells > 50_000_000:
        raise ValueError('This traceback prototype supports at most 50 million DP cells per pair.')
    scores = array('i', [0]) * cells
    trace = bytearray(cells)
    best_score = best_index = 0
    for i, a in enumerate(query, 1):
        row = i * width
        previous = row - width
        for j, b in enumerate(target, 1):
            diagonal = scores[previous + j - 1] + (match if a == b and a != 'N' else mismatch)
            up = scores[previous + j] + gap
            left = scores[row + j - 1] + gap
            value = max(0, diagonal, up, left)
            index = row + j
            scores[index] = value
            if value > 0:
                trace[index] = 1 if value == diagonal else (2 if value == up else 3)
            if value > best_score:
                best_score, best_index = value, index
    i, j = divmod(best_index, width)
    query_end, target_end = i, j
    aligned_query, aligned_target = [], []
    while trace[i * width + j]:
        direction = trace[i * width + j]
        if direction == 1:
            i -= 1
            j -= 1
            aligned_query.append(query[i])
            aligned_target.append(target[j])
        elif direction == 2:
            i -= 1
            aligned_query.append(query[i])
            aligned_target.append('-')
        else:
            j -= 1
            aligned_query.append('-')
            aligned_target.append(target[j])
    aligned_query = ''.join(reversed(aligned_query))
    aligned_target = ''.join(reversed(aligned_target))
    length = len(aligned_query)
    matches = sum(a == b and a != 'N' for a, b in zip(aligned_query, aligned_target))
    gaps = aligned_query.count('-') + aligned_target.count('-')
    return dict(score=best_score, alignment_length=length, matches=matches,
                mismatches=length - matches - gaps, gaps=gaps,
                percent_identity=100 * matches / length if length else 0,
                query_coverage=100 * (query_end - i) / len(query) if query else 0,
                target_coverage=100 * (target_end - j) / len(target) if target else 0,
                query_start=i, query_end=query_end, target_start=j, target_end=target_end,
                aligned_query=aligned_query, aligned_target=aligned_target,
                dp_cells=len(query) * len(target), dp_buffer_bytes=cells * (scores.itemsize + 1))


def read_fasta(path):
    records = []
    header, parts = None, []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    records.append((header, clean_dna(''.join(parts))))
                header, parts = line[1:], []
            else:
                if header is None:
                    raise ValueError('FASTA sequence appears before its header.')
                parts.append(line)
    if header is not None:
        records.append((header, clean_dna(''.join(parts))))
    return records
