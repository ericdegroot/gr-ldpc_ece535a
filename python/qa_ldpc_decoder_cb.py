#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 
# ECE 535A Project
# LDPC Decoder Unit Test

from gnuradio import gr, gr_unittest
from gnuradio import blocks
import ldpc_ece535a_swig as ldpc_ece535a

class qa_ldpc_decoder_cb (gr_unittest.TestCase):

    def setUp (self):
        self.tb = gr.top_block ()

    def tearDown (self):
        self.tb = None

    def test_001_t (self):
        x = [ 112, 80, 93, 120, 20, 50 ]
        y_expected = [ 112, 80, 93, 120, 20, 50 ]

        # set up fg
        x_src = blocks.vector_source_c(x, False)
        ldpc_decoder = ldpc_ece535a.ldpc_decoder_cb()
        y_dst = blocks.vector_sink_b()

        self.tb.connect((x_src, 0), (ldpc_decoder, 0))
        self.tb.connect((ldpc_decoder, 0), (y_dst, 0))

        self.tb.run ()

        # check data
        y = y_dst.data()
        self.assertFloatTuplesAlmostEqual(y_expected, y, 6)


if __name__ == '__main__':
    gr_unittest.run(qa_ldpc_decoder_cb, "qa_ldpc_decoder_cb.xml")
