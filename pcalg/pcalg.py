#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""@package graph

A graph generator based on the PC algorithm [Kalisch2007].

[Kalisch2007] Markus Kalisch and Peter Bhlmann. Estimating
high-dimensional directed acyclic graphs with the pc-algorithm. In The
Journal of Machine Learning Research, Vol. 8, pp. 613-636, 2007.

License: GPLv2

"""

from itertools import combinations, permutations
import logging

import networkx as nx

_logger = logging.getLogger(__name__)

def _create_complete_graph(node_ids):
    """Create a complete graph from the list of node ids.

    @param node_ids: a list of node ids
    @return an undirected graph (as a networkx.Graph)
    """
    g = nx.Graph()
    g.add_nodes_from(node_ids)
    for (i, j) in combinations(node_ids, 2):
        g.add_edge(i, j)
        pass
    return g

def estimate_skeleton(suff_stat, indep_test_func, alpha, max_reach=None):
    """Estimate a skeleton graph from the statistis information.

    @param suff_stat: the sufficient statistics (as a numpy.ndarray).
    @param indep_test_func: the function name for a conditional
      independency test.
    @param alpha: the significance level.
    @return g: a skeleton graph (as a networkx.Graph).
    @return sep_set: a separation set (An 2D-array of set()).
    """
    node_ids = range(suff_stat.shape[1])
    g = _create_complete_graph(node_ids)

    node_size = suff_stat.shape[1]
    sep_set = [[set() for i in range(node_size)] for j in range(node_size)]

    l = 1
    while True:
        cont = False
        for (i, j) in permutations(node_ids, 2):
            adj_i = g.neighbors(i)
            if j not in adj_i:
                continue
            else:
                adj_i.remove(j)
                pass
            if len(adj_i) >= l:
                _logger.debug('testing %s and %s' % (i,j))
                _logger.debug('neighbors of %s are %s' % (i, str(adj_i)))
                if len(adj_i) < l:
                    continue
                for k in combinations(adj_i, l):
                    _logger.debug('indep prob of %s and %s with subset %s'
                                  % (i, j, str(k)))
                    prob = indep_test_func(i, j, set(k), suff_stat)
                    _logger.debug('prob is %s' % str(prob))
                    if prob > alpha:
                        if g.has_edge(i, j):
                            _logger.debug('remove edge (%s, %s)' % (i, j))
                            g.remove_edge(i, j)
                            pass
                        sep_set[i][j] |= set(k)
                        sep_set[j][i] |= set(k)
                        break
                    pass
                cont = True
                pass
            pass
        l += 1
        if cont is False:
            break
        if (max_reach is not None) and (l > max_reach):
            break
        pass

    return (g, sep_set)

def estimate_cpdag(skel_graph, sep_set):
    """Estimate a CPDAG from the skeleton graph and separation sets
    returned by the estimate_skeleton() function.

    @param skel_graph: A skeleton graph (an undirected networkx.Graph).
    @param sep_set: An 2D-array of separation set.
      The contents look like something sep_set[i][j] = set([k, l, m]).
    @return dag: An estimated DAG.
    """
    dag = skel_graph.to_directed()
    node_ids = skel_graph.nodes()
    for (i, j) in combinations(node_ids, 2):
        adj_i = set(skel_graph.neighbors(i))
        if j in adj_i:
            continue
        adj_j = set(skel_graph.neighbors(j))
        common_k = adj_i & adj_j
        for k in common_k:
            if k not in sep_set[i][j]:
                if dag.has_edge(k, i):
                    dag.remove_edge(k, i)
                    pass
                if dag.has_edge(k, j):
                    dag.remove_edge(k, j)
                    pass
                pass
            pass
        pass

    def _has_both_edges(dag, i, j):
        return dag.has_edge(i, j) and dag.has_edge(j, i)

    def _has_any_edge(dag, i, j):
        return dag.has_edge(i, j) or dag.has_edge(j, i)

    def _has_one_edge(dag, i, j):
        return ((dag.has_edge(i, j) and (not dag.has_edge(j, i))) or
                (not dag.has_edge(i, j)) and dag.has_edge(j, i))

    def has_no_edge(dag, i, j):
        return (not dag.has_edge(i, j)) and (not dag.has_edge(j, i))

    # For all the combination of nodes i and j, apply the following
    # rules.
    for (i, j) in combinations(node_ids, 2):
        # Rule 1: Orient i-j into i->j whenever there is an arrow k->i
        # such that k and j are nonadjacent.
        #
        # Check if i-j.
        if _has_both_edges(dag, i, j):
            # Look all the predecessors of i.
            for k in dag.predecessors(i):
                # Skip if there is an arrow i->k.
                if dag.has_edge(i, k):
                    continue
                # Skip if k and j are adjacent.
                if _has_any_edge(dag, k, j):
                    continue
                # Make i-j into i->j
                dag.remove_edge(j, i)
                pass
            pass

        # Rule 2: Orient i-j into i->j whenever there is a chain
        # i->k->j.
        #
        # Check if i-j.
        if _has_both_edges(dag, i, j):
            # Find nodes k where k is i->k.
            succs_i = set()
            for k in dag.successors(i):
                if not dag.has_edge(k, i):
                    succs_i.add(k)
                    pass
                pass
            # Find nodes j where j is k->j.
            preds_j = set()
            for k in dag.predecessors(j):
                if not dag.has_edge(j, k):
                    preds_j.add(k)
                    pass
                pass
            # Check if there is any node k where i->k->j.
            if len(succs_i & preds_j) > 0:
                # Make i-j into i->j
                dag.remove_edge(j, i)
                pass
            pass

        # Rule 3: Orient i-j into i->j whenever there are two chains
        # i-k->j and i-l->j such that k and l are nonadjacent.
        #
        # Check if i-j.
        if _has_both_edges(dag, i, j):
            # Find nodes k where i-k.
            adj_i = set()
            for k in dag.successors(i):
                if dag.has_edge(k, i):
                    adj_i.add(k)
                    pass
                pass
            # For all the pairs of nodes in adj_i,
            for (k, l) in combinations(adj_i, 2):
                # Skip if k and l are adjacent.
                if _has_any_edge(dag, k, l):
                    continue
                # Skip if not k->j.
                if dag.has_edge(j, k) or (not dag.has_edge(k, j)):
                    continue
                # Skip if not l->j.
                if dag.has_edge(j, l) or (not dag.has_edge(l, j)):
                    continue
                # Make i-j into i->j.
                dag.remove_edge(j, i)
                pass
            pass

        #R4
            
        pass

    return dag
