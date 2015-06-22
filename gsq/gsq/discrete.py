#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""@module discrete

A conditional independency test function for discrete data.

The code included in this package is logically copied and pasted from
the pcalg package for R developed by Markus Kalisch, Alain Hauser,
Martin Maechler, Diego Colombo, Doris Entner, Patrik Hoyer, Antti
Hyttinen, and Jonas Peters.

License: GPLv2
"""

from __future__ import print_function

import logging

from scipy.stats import chi2
import numpy as np

_logger = logging.getLogger(__name__)

def g_square_dis(x, y, s, levels, dm):
    def _calculate_tlog(x, y, s, dof, levels, dm):
        prod_levels = np.prod(map(lambda x: levels[x], s))
        nijk = np.zeros((levels[x], levels[y], prod_levels))
        s_size = len(s)
        for row_index in range(dm.shape[0]):
            i = dm[row_index, x]
            j = dm[row_index, y]
            k = []
            k_index = 0
            for s_index in range(s_size):
                if s_index == 0:
                    k_index += dm[row_index, s[s_index]]
                else:
                    lprod = np.prod(map(lambda x: levels[x], s[:s_index]))
                    k_index += (dm[row_index, s[s_index]] * lprod)
                    pass
                pass
            nijk[i, j, k_index] += 1
            pass
        nik = np.ndarray((levels[x], prod_levels))
        njk = np.ndarray((levels[y], prod_levels))
        for k_index in range(prod_levels):
            nik[:, k_index] = nijk[:, :, k_index].sum(axis = 1)
            njk[:, k_index] = nijk[:, :, k_index].sum(axis = 0)
            pass
        nk = njk.sum(axis = 0)
        tlog = np.zeros((levels[x], levels[y], prod_levels))
        tlog.fill(np.nan)
        for k in range(prod_levels):
            tx = np.array([nik[:, k]]).T
            ty = np.array([njk[:, k]])
            tdijk = tx.dot(ty)
            tlog[:, :, k] = nijk[:, :, k] * nk[k] / tdijk
            pass
        return (nijk, tlog)

    _logger.debug('Edge %d -- %d with subset: %s' % (x, y, s))
    row_size = dm.shape[0]
    s_size = len(s)
    dof = ((levels[x] - 1) * (levels[y] - 1)
           * np.prod(map(lambda x: levels[x], s)))
    _logger.debug('dof = %d' % dof)
    if row_size < 10 * dof:
        _logger.warning('Not enough samples.')
        return 1
    nijk = None
    if s_size < 5:
        if s_size == 0:
            nijk = np.zeros((levels[x], levels[y]))
            for row_index in range(row_size):
                i = dm[row_index, x]
                j = dm[row_index, y]
                nijk[i, j] += 1
                pass
            tx = np.array([nijk.sum(axis = 1)]).T
            ty = np.array([nijk.sum(axis = 0)])
            tdij = tx.dot(ty)
            tlog = nijk * row_size / tdij
            pass
        if s_size > 0:
            nijk, tlog = _calculate_tlog(x, y, s, dof, levels, dm)
            pass
        pass
    else:
        # s_size >= 5
        nijk = np.zeros((levels[x], levels[y], 1))
        i = dm[0, x]
        j = dm[0, y]
        k = []
        for z in s:
            k.append(dm[:, z])
            pass
        k = np.array(k).T
        parents_count = 1
        parents_val = np.array([k[0, :]])
        nijk[i, j, parents_count - 1] = 1
        for it_sample in range(1, row_size):
            is_new = True
            i = dm[it_sample, x]
            j = dm[it_sample, y]
            tcomp = parents_val[:parents_count, :] == k[it_sample, :]
            for it_parents in range(parents_count):
                if np.all(tcomp[it_parents, :]):
                    nijk[i, j, it_parents] += 1
                    is_new = False
                    break
                pass
            if is_new is True:
                parents_count += 1
                parents_val = np.r_[parents_val, [k[it_sample, :]]]
                nnijk = np.zeros((levels[x], levels[y], parents_count))
                for p in range(parents_count - 1):
                    nnijk[:, :, p] = nijk[:, :, p]
                    pass
                nnijk[i, j, parents_count - 1] = 1
                nijk = nnijk
                pass
            pass
        nik = np.ndarray((levels[x], parents_count))
        njk = np.ndarray((levels[y], parents_count))
        for k_index in range(parents_count):
            nik[:, k_index] = nijk[:, :, k_index].sum(axis = 1)
            njk[:, k_index] = nijk[:, :, k_index].sum(axis = 0)
            pass
        nk = njk.sum(axis = 0)
        tlog = np.zeros((levels[x], levels[y], parents_count))
        tlog.fill(np.nan)
        for k in range(parents_count):
            tx = np.array([nik[:, k]]).T
            ty = np.array([njk[:, k]])
            tdijk = tx.dot(ty)
            tlog[:, :, k] = nijk[:, :, k] * nk[k] / tdijk
            pass
        pass
    log_tlog = np.log(tlog)
    G2 = np.nansum(2 * nijk * log_tlog)
    _logger.debug('nijk = %s' % nijk)
    _logger.debug('tlog = %s' % tlog)
    _logger.debug('log(tlog) = %s' % log_tlog)
    _logger.debug('G2 = %f' % G2)
    ret = chi2.sf(G2, dof)
    _logger.info('probability = %s' % str(ret))
    return ret


if __name__ == '__main__':
    import gsq_testdata

    _logger.setLevel(logging.DEBUG)
    _console_handler = logging.StreamHandler()
    _logger.setLevel(logging.DEBUG)
    _logger.addHandler(_console_handler)

    dm = np.array([gsq_testdata.dis_data]).reshape((10000,5))

    g_square_dis(0, 1, [2,3], [3,2,3,4,2], dm)
