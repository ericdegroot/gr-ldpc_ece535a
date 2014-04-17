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
#>>> data  = ( 0b11101101, 0b00000010, 0b00110100, 0b10000000 )
#>>> check = ( 0b10110100, 0b10101001, 0b11010101, 0b10001000 )
#>>> print data
#(237, 2, 52, 128)
#>>> print check
#(180, 169, 213, 136)

        check = ( 0b10110100, 0b10101001, 0b11010101, 0b10001000 )
        data  = ( 0b11101101, 0b00000010, 0b00110100, 0b10000000 )

        print check
        print data

        expected = tuple()
        for pair in zip(check, data):
            expected += pair

        print expected

#check =
#     0     1     0     1     0     0     0     0
#data =
#     0     1     0     1     0     0     1     1
#check =
#     0     1     0     1     1     0     0     1
#data =
#     0     0     0     0     0     1     1     0
#check =
#     1     0     0     1     0     0     0     0
#data =
#     0     0     1     1     0     0     0     1
#check =
#     1     1     1     1     0     0     0     1
#data =
#     1     1     1     0     1     0     0     0

        # set up fg
        src = blocks.vector_source_b(data, False)
        ldpc_encoder = ldpc_ece535a.ldpc_encoder_bb()
        dst = blocks.vector_sink_b()

        self.tb.connect((src, 0), (ldpc_encoder, 0))
        self.tb.connect((ldpc_encoder, 0), (dst, 0))

        self.tb.run ()

        # check data
        result = dst.data()
        #print expected
        #print result
        self.assertTupleEqual(expected, result)


if __name__ == '__main__':
    gr_unittest.run(qa_ldpc_encoder_bb, "qa_ldpc_encoder_bb.xml")
