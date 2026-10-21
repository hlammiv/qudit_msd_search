import networkx as nx

from primepower_msd.outer_code import biaffine, check_girth5, sensitivity_graph


def test_girth5_iff_order5_sensitivity():
    for t in range(30):
        G = nx.random_regular_graph(3, 12, seed=t)
        adj = [set(G[u]) for u in range(12)]
        assert check_girth5(adj) == sensitivity_graph(adj, 5)
    P = nx.petersen_graph()
    adj = [set(P[u]) for u in range(10)]
    assert check_girth5(adj) and sensitivity_graph(adj, 5)


def test_biaffine_regular_girth():
    A = biaffine(7, xs=range(5), ms=range(5))
    assert {len(a) for a in A} == {5} and check_girth5(A)
    G = nx.Graph(); G.add_edges_from((u, v) for u in range(len(A)) for v in A[u])
    assert nx.girth(G) >= 6
