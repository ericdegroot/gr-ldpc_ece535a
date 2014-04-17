#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 
# ECE 535A Project
# LDPC encoder block unit test
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
        # data and check bits generated in MATLAB
        data  = ( 0b11101101, 0b00000010, 0b00110100, 0b10000000, 0b01010011, 0b00000110, 0b00110001, 0b11101000 )
        check = ( 0b10110100, 0b10101001, 0b11010101, 0b10001000, 0b01010000, 0b01011001, 0b10010000, 0b11110001 )

        # expected result is alternating check,data,check,data,... bytes
        expected = tuple()
        for pair in zip(check, data):
            expected += pair

        # set up fg
        src = blocks.vector_source_b(data, False)
        ldpc_encoder = ldpc_ece535a.ldpc_encoder_bb()
        dst = blocks.vector_sink_b()

        self.tb.connect((src, 0), (ldpc_encoder, 0))
        self.tb.connect((ldpc_encoder, 0), (dst, 0))

        self.tb.run ()

        # check data
        result = dst.data()

        self.assertTupleEqual(expected, result)


if __name__ == '__main__':
    gr_unittest.run(qa_ldpc_encoder_bb, "qa_ldpc_encoder_bb.xml")
