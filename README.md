Components
==========

LDPC Encoder
------------

lib/ldpc_encoder_bb_impl.cc

This is the C++ source code file for the LDPC Encoder block. DSP code will be added to the work method.

LDPC Decoder
------------

lib/ldpc_decoder_bb_impl.cc

This is the C++ source code file for the LDPC Decoder block. DSP code will be added to the work method.

GNU Radio Companion (GRC) Example
---------------------------------

examples/example1.grc

This is a GRC flow graph that includes the above LDPC blocks. When the example is run in GRC a histogram with the character distribution will be shown.

Application Example
-------------------

apps/ldpc_ece535a_dump

This is an example Python application that uses the LDPC blocks but does not require GRC to run. You can run it by:

./apps/ldpc_ece535a_dump

At the moment it will dump a repeating stream of random ASCII characters to your console. Stop it by using Ctrl-C.

Resources
=========
* [Installing GNU Radio From Source](http://gnuradio.org/redmine/projects/gnuradio/wiki/InstallingGRFromSource)
* [Developing Out-of-tree modules](http://gnuradio.org/redmine/projects/gnuradio/wiki/OutOfTreeModules)
