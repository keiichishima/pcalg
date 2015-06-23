#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""A conditional independency test function for binary data.

The code included in this package is logically copied and pasted from
the pcalg package for R developed by Markus Kalisch, Alain Hauser,
Martin Maechler, Diego Colombo, Doris Entner, Patrik Hoyer, Antti
Hyttinen, and Jonas Peters.

License: GPLv2
"""

from __future__ import print_function

import logging
from math import pow

from scipy.stats import chi2
import numpy as np

_logger = logging.getLogger(__name__)

def g_square_bin(dm, x, y, s):
    """G square test for a binary data.

    Args:
        dm: the data matrix to be used (as a numpy.ndarray).
        x: the first node (as an integer).
        y: the second node (as an integer).
        s: the set of neibouring nodes of x and y (as a set()).

    Returns:
        p_val: the p-value of conditional independence.
    """

    def _calculate_tlog(x, y, s, dof, dm):
        nijk = np.zeros((2, 2, dof))
        s_size = len(s)
        z = []
        for z_index in range(s_size):
            z.append(s.pop())
            pass
        for row_index in range(0, dm.shape[0]):
            i = dm[row_index, x]
            j = dm[row_index, y]
            k = []
            k_index = 0
            for z_index in range(s_size):
                k_index += dm[row_index, z[z_index]] * int(pow(2, z_index))
                pass
            nijk[i, j, k_index] += 1
            pass
        nik = np.ndarray((2, dof))
        njk = np.ndarray((2, dof))
        for k_index in range(dof):
            nik[:, k_index] = nijk[:, :, k_index].sum(axis = 1)
            njk[:, k_index] = nijk[:, :, k_index].sum(axis = 0)
            pass
        nk = njk.sum(axis = 0)
        tlog = np.zeros((2, 2 , dof))
        tlog.fill(np.nan)
        for k in range(dof):
            tx = np.array([nik[:,k]]).T
            ty = np.array([njk[:,k]])
            tdijk = tx.dot(ty)
            tlog[:,:,k] = nijk[:,:,k] * nk[k] / tdijk
            pass
        return (nijk, tlog)

    _logger.debug('Edge %d -- %d with subset: %s' % (x, y, s))
    row_size = dm.shape[0]
    s_size = len(s)
    dof = int(pow(2, s_size))
    _logger.debug('dof = %d' % dof)
    if row_size < 10 * dof:
        _logger.warning('Not enough samples.  %s is too small.'
                        % str(row_size))
        return 1
    nijk = None
    if s_size < 6:
        if s_size == 0:
            nijk = np.zeros((2, 2))
            for row_index in range(0, dm.shape[0]):
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
            nijk, tlog = _calculate_tlog(x, y, s, dof, dm)
            pass
        pass
    else:
        # s_size >= 6
        nijk = np.zeros((2, 2, 1))
        i = dm[0, x]
        j = dm[0, y]
        k = []
        for z in s:
            k.append(dm[:,z])
            pass
        k = np.array(k).T
        parents_count = 1
        parents_val = np.array([k[0,:]])
        nijk[i, j, parents_count - 1] = 1
        for it_sample in range(1, row_size):
            is_new = True
            i = dm[it_sample, x]
            j = dm[it_sample, y]
            tcomp = parents_val[:parents_count,:] == k[it_sample,:]
            for it_parents in range(parents_count):
                if np.all(tcomp[it_parents,:]):
                    nijk[i, j, it_parents] += 1
                    is_new = False
                    break
                pass
            if is_new is True:
                parents_count += 1
                parents_val = np.r_[parents_val, [k[it_sample,:]]]
                nnijk = np.zeros((2,2,parents_count))
                for p in range(parents_count - 1):
                    nnijk[:,:,p] = nijk[:,:,p]
                nnijk[i, j, parents_count - 1] = 1
                nijk = nnijk
                pass
            pass
        nik = np.ndarray((2, parents_count))
        njk = np.ndarray((2, parents_count))
        for k_index in range(parents_count):
            nik[:, k_index] = nijk[:, :, k_index].sum(axis = 1)
            njk[:, k_index] = nijk[:, :, k_index].sum(axis = 0)
            pass
        nk = njk.sum(axis = 0)
        tlog = np.zeros((2, 2 , parents_count))
        tlog.fill(np.nan)
        for k in range(parents_count):
            tX = np.array([nik[:,k]]).T
            tY = np.array([njk[:,k]])
            tdijk = tX.dot(tY)
            tlog[:,:,k] = nijk[:,:,k] * nk[k] / tdijk
            pass
        pass
    log_tlog = np.log(tlog)
    G2 = np.nansum(2 * nijk * log_tlog)
    _logger.debug('nijk = %s' % nijk)
    _logger.debug('tlog = %s' % tlog)
    _logger.debug('log(tlog) = %s' % log_tlog)
    _logger.debug('G2 = %f' % G2)
    p_val = chi2.sf(G2, dof)
    _logger.info('p_val = %s' % str(p_val))
    return p_val


if __name__ == '__main__':
    from math import frexp

    import gsq_testdata

    dm = np.array([gsq_testdata.bin_data]).reshape((5000,5))
    x = 0
    y = 1

    sets = [[], [2], [2, 3], [3, 4], [2, 3, 4]]
    for idx in range(len(sets)):
        print('x =', x, ', y =', y, ', s =', sets[idx], end='')
        p = g_square_bin(dm, 0, 1, set(sets[idx]))
        print(', p =', p, end='')
        fr_p = frexp(p)
        fr_a = frexp(gsq_testdata.bin_answer[idx])
        if (round(fr_p[0] - fr_a[0], 7) == 0
            and fr_p[1] == fr_a[1]):
            print(' => GOOD')
        else:
            print(' => WRONG')
            print('p =', fr_p)
            print('a =', fr_a)
