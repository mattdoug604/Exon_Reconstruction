#!/usr/bin/python3

import sys
sys.path.append("..")
import pytest
from exon_trap import frame_shift

def test_frame_shift_is_in_range_3():
    for frame in range(3):
        for i in range(1, 4):
            intron = '_', 1, i
            for dir in ('+<', '->', '+>', '-<'):
                new_frame = frame_shift(frame, intron, dir)
                assert new_frame in (0, 1, 2)

test_frame_shift_is_in_range_3()
