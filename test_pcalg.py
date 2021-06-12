# -*- coding: utf-8 -*-
'''
Test suite for pcalg
'''
import networkx as nx
import numpy as np

from gsq.ci_tests import ci_test_bin, ci_test_dis
from gsq.gsq_testdata import bin_data, dis_data

import pytest

from pcalg import estimate_cpdag
from pcalg import estimate_skeleton

@pytest.mark.parametrize(('indep_test_func', 'data_matrix', 'g_answer'), [
    (ci_test_bin, np.array(bin_data).reshape((5000, 5)), nx.DiGraph({
        0: (1, ),
        1: (),
        2: (3, 4),
        3: (1, 2),
        4: (1, 2),
    })),
    (ci_test_dis, np.array(dis_data).reshape((10000, 5)), nx.DiGraph({
        0: (2, ),
        1: (2, 3),
        2: (),
        3: (),
        4: (3, ),
    })),
])
def test_estimate_cpdag(indep_test_func, data_matrix, g_answer, alpha=0.01):
    '''
    estimate_cpdag should reveal the answer
    '''
    (graph, sep_set) = estimate_skeleton(indep_test_func=indep_test_func,
                                         data_matrix=data_matrix,
                                         alpha=alpha)
    graph = estimate_cpdag(skel_graph=graph, sep_set=sep_set)
    error_msg = 'True edges should be: %s' % (g_answer.edges(), )
    assert nx.is_isomorphic(graph, g_answer), error_msg

def test_fixed_edges():
    '''
    The fixed edges shall appear in the skeleton
    '''
    data_matrix = np.array(bin_data).reshape((5000, 5))
    (graph, sep_set) = estimate_skeleton(indep_test_func=ci_test_bin,
                                         data_matrix=data_matrix,
                                         alpha=0.01)
    graph = estimate_cpdag(skel_graph=graph, sep_set=sep_set)
    assert not graph.has_edge(1, 2)

    fixed_edges = nx.DiGraph()
    fixed_edges.add_nodes_from(range(5))
    fixed_edges.add_edge(1, 2)
    with pytest.raises(ValueError):
        _ = estimate_skeleton(indep_test_func=ci_test_bin,
                              data_matrix=data_matrix,
                              alpha=0.01,
                              fixed_edges=((1,2), ))
    with pytest.raises(ValueError):
        _ = estimate_skeleton(indep_test_func=ci_test_bin,
                              data_matrix=data_matrix,
                              alpha=0.01,
                              fixed_edges=nx.DiGraph({0: (1, )}))
    (graph, _) = estimate_skeleton(indep_test_func=ci_test_bin,
                                   data_matrix=data_matrix,
                                   alpha=0.01,
                                   fixed_edges=fixed_edges)
    assert graph.has_edge(1, 2), graph.edges
