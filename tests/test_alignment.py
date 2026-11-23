"""Independent score checks against Biopython plus hand-checkable edge cases."""
import random
import unittest
from Bio.Align import PairwiseAligner
from alignment import smith_waterman


class AlignmentTests(unittest.TestCase):
    def test_known_cases(self):
        self.assertEqual(smith_waterman('ACGT', 'ACGT')['score'], 8)
        self.assertEqual(smith_waterman('AAAA', 'CCCC')['score'], 0)
        self.assertEqual(smith_waterman('', 'ACGT')['alignment_length'], 0)
        self.assertEqual(smith_waterman('NN', 'NN')['score'], 0)
        result = smith_waterman('TTACGTAA', 'GGACGTCC')
        self.assertEqual((result['query_start'], result['query_end']), (2, 6))
        self.assertEqual(result['percent_identity'], 100)

    def test_independent_scores_and_traceback(self):
        oracle = PairwiseAligner(mode='local', match_score=2, mismatch_score=-1,
                                 open_gap_score=-2, extend_gap_score=-2)
        rng = random.Random(691)
        for _ in range(100):
            a = ''.join(rng.choices('ACGT', k=rng.randint(1, 35)))
            b = ''.join(rng.choices('ACGT', k=rng.randint(1, 35)))
            result = smith_waterman(a, b)
            self.assertEqual(result['score'], oracle.score(a, b))
            self.assertEqual(result['score'], 2 * result['matches'] - result['mismatches'] - 2 * result['gaps'])
            self.assertEqual(result['aligned_query'].replace('-', ''), a[result['query_start']:result['query_end']])
            self.assertEqual(result['aligned_target'].replace('-', ''), b[result['target_start']:result['target_end']])


if __name__ == '__main__':
    unittest.main()
