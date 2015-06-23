#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""A collection of wrapper functions for G-square conditional
independecy test functions.
"""

from __future__ import print_function

import logging

from binary import g_square_bin
from discrete import g_square_dis

_logger = logging.getLogger(__name__)

def ci_test_bin(data_matrix, x, y, s):
    return g_square_bin(data_matrix, x, y, s)

def ci_test_dis(data_matrix, x , y, s, **kwargs):
    levels = []
    if 'levels' in kwargs:
        levels = kwargs['levels']
    else:
        import numpy as np
        levels = np.amax(data_matrix, axis=0) + 1
    return g_square_dis(data_matrix, x, y, s, levels)
