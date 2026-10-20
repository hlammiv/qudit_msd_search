# Exact-vs-proxy validation of the Z_d colour-code T-check (D = 3, n = 7)

Dense data-level model: check C = A^dag X^{(x)n} A with A = prod_i G^{-tau_i}, G = exact T (escaping carrier)
or the Clifford proxy S = zeta_{2d}^{x^2}. Pauli faults at 4 slots per data qudit (before A, after A, after the
controlled X, after A^dag); post-selection on the check and on the trivial stabilizer syndrome.
u(E) = accepted-but-wrong weight.

Run from this directory with PYTHONPATH=<repo root>:
- `python proxy_validate.py 2`: single, pair and weight-3 logical-support faults at d=2, saved to proxy_validate_d2.json
- `python slots.py 2 7`: breakdown by slot pattern at d=2
- d=4 (about 1 h on lenore with 30 workers): `for K in $(seq 0 29); do python slots_w.py 4 $K 30 & done; wait; python slots_agg.py 4`
