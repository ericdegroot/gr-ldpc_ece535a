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
        mod_data  = ( ( 1, 1, 1, -1, 1, 1, -1, 1 ),
                      ( -1, -1, -1, -1, -1, -1, 1, -1 ),
                      ( -1, -1, 1, 1, -1, 1, -1, -1 ),
                      ( 1, -1, -1, -1, -1, -1, -1, -1 ),
                      ( -1, 1, -1, 1, -1, -1, 1, 1 ),
                      ( -1, -1, -1, -1, -1, 1, 1, -1 ),
                      ( -1, -1, 1, 1, -1, -1, -1, 1 ),
                      ( 1, 1, 1, -1, 1, -1, -1, -1 ) )

        mod_check = ( ( 1, -1, 1, 1, -1, 1, -1, -1 ),
                      ( 1, -1, 1, -1, 1, -1, -1, 1 ),
                      ( 1, 1, -1, 1, -1, 1, -1, 1 ),
                      ( 1, -1, -1, -1, 1, -1, -1, -1 ),
                      ( -1, 1, -1, 1, -1, -1, -1, -1 ),
                      ( -1, 1, -1, 1, 1, -1, -1, 1 ),
                      ( 1, -1, -1, 1, -1, -1, -1, -1 ),
                      ( 1, 1, 1, 1, -1, -1, -1, 1 ) )

        mod_frames = tuple()
        for pair in zip(mod_check, mod_data):
            for item in pair:
                mod_frames += item

        expected_data  = ( 0b11101101, 0b00000010, 0b00110100, 0b10000000, 0b01010011, 0b00000110, 0b00110001, 0b11101000 )

        # set up fg
        src = blocks.vector_source_c(mod_frames, False)
        ldpc_decoder = ldpc_ece535a.ldpc_decoder_cb()
        dst = blocks.vector_sink_b()

        self.tb.connect((src, 0), (ldpc_decoder, 0))
        self.tb.connect((ldpc_decoder, 0), (dst, 0))

        self.tb.run ()

        # check data
        result = dst.data()
        self.assertTupleEqual(expected_data, result)


if __name__ == '__main__':
    gr_unittest.run(qa_ldpc_decoder_cb, "qa_ldpc_decoder_cb.xml")
