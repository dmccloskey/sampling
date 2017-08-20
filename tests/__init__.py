#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_sampling
----------------------------------

Tests for `sampling` module.
"""

from __future__ import absolute_import

from os.path import abspath, dirname, join

try:
    import pytest
    import pytest_benchmark
except ImportError:
    pytest = None

sampling_directory = abspath(join(dirname(abspath(__file__)), ".."))
sampling_location = abspath(join(sampling_directory, "sampling"))
data_dir = join(sampling_directory, "tests/data", "")

def test_all(args=None):
    """ alias for running all unit-tests on installed sampling
    """
    if pytest:
        args = args if args else []

        return pytest.main(
            ['--pyargs', 'sampling', '--benchmark-skip', '-v', '-rs'] + args
        )
    else:
        raise ImportError('missing package pytest and pytest_benchmark'
                          ' required for testing')