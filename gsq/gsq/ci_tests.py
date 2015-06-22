#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""@module ci_tests

A collection of wrapper functions for G-square conditional independecy
test functions.

"""

from __future__ import print_function

import logging

from binary import g_square_bin
from discrete import g_square_dis

_logger = logging.getLogger(__name__)

def bin_ci_test(data_matrix, x, y, s):
    return g_square_bin(x, y, s, data_matrix)

def dis_ci_test(data_matrix, x , y, s, **kwargs):
    levels = kwargs['levels']
    return g_square_dis(x, y, s, levels, data_matrix)
