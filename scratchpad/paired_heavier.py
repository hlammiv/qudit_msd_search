"""Proved weight-9 and weight-10 bounds for a paired-cap code [[686 - 2t, 2t, 8]]_7 (d' = 4, d'' = 8; WEIGHTED_SUPPORTS.md,
"Single-species codes").  N = 343 - t columns per hyperplane.  A full-code word with off-Y weights (x, y) on the two
hyperplanes has x, y in {0} u [4, oo) (d' = 4); x = 0 or y = 0 forces the other side into ker H6 (weight >= d'' = 8).
  one-sided, weight w: on a w-set the ker-H6 subcode has distance >= 8, so dimension <= kmax(w) by the Griesmer bound
      (w = 9: [9,2,8]_7 needs n >= 8 + 2 = 10, so k <= 1; w = 10: [10,3,8]_7 needs 8 + 2 + 1 = 11, so k <= 2); at most
      (7^kmax - 1) vectors per w-set, both hyperplanes.
  cross (x, y): for a fixed side-1 vector the valid side-2 vectors form a coset of ker H6 (distance >= 8), so two of
      weight y share <= y - 4 support points and differ in support (same support => distance <= y < 8): at most
      C(N, y - 3) / C(y, y - 3) of them.  Side-1 vectors: the exact d'-word list for weight 4 (x (p - 1)), and for weight 5
      at most 48 C(N, 5) (Singleton: the ker-H5 subcode on a 5-set has distance >= 4, dimension <= 2).
Per-output failing bounds are taken equal to the word bounds.  All counts are vectors (not projective)."""
from math import comb


def heavier_bounds(N, dprime_words, p=7):
    a4 = (p - 1) * dprime_words                      # weight-4 side vectors (exact list)
    a5 = 48 * comb(N, 5)                             # weight-5 side vectors (bound)
    partners = lambda y: comb(N, y - 3) // comb(y, y - 3)
    N9 = 2 * (p ** 1 - 1) * comb(N, 9) + 2 * a4 * partners(5)
    N10 = 2 * (p ** 2 - 1) * comb(N, 10) + 2 * a4 * partners(6) + a5 * partners(5)
    return {9: (N9, N9), 10: (N10, N10)}


if __name__ == "__main__":
    import sys, json
    c = json.load(open(sys.argv[1])); N = 343 - c["k"] // 2
    print(heavier_bounds(N, c["dprime_words"]))
