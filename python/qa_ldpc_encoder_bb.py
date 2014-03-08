#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 
# Copyright 2014 <+YOU OR YOUR COMPANY+>.
# 
# This is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3, or (at your option)
# any later version.
# 
# This software is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this software; see the file COPYING.  If not, write to
# the Free Software Foundation, Inc., 51 Franklin Street,
# Boston, MA 02110-1301, USA.
# 

from gnuradio import gr, gr_unittest
from gnuradio import blocks
import ldpc_ece535a_swig as ldpc_ece535a

class qa_ldpc_encoder_bb (gr_unittest.TestCase):

    def setUp (self):
        self.tb = gr.top_block ()

    def tearDown (self):
        self.tb = None

    def test_001_t (self):
        x = [ 0, 1, 0, 1, 0, 1 ]
        y_expected = [ 0, 1, 0, 1, 0, 1 ]

        # set up fg
        x_src = gr.vector_source_f(x, False)
        ldpc_encoder = ldpc_ece535a.ldpc_encoder_bb()
        y_dst = gr.vector_sink_f()

        self.tb.connect((x_src, 0), (ldpc_encoder, 0))
        self.tb.connect((ldpc_encoder, 0), (y_dst, 0))

        self.tb.run ()

        # check data
        y = dst.data()
        self.assertFloatTuplesAlmostEqual(y_expected, y, 6)


if __name__ == '__main__':
    gr_unittest.run(qa_ldpc_encoder_bb, "qa_ldpc_encoder_bb.xml")
