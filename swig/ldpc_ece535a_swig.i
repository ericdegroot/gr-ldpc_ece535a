/* -*- c++ -*- */

#define LDPC_ECE535A_API

%include "gnuradio.i"			// the common stuff

//load generated python docstrings
%include "ldpc_ece535a_swig_doc.i"

%{
#include "ldpc_ece535a/ldpc_encoder_bc.h"
#include "ldpc_ece535a/ldpc_decoder_cb.h"
#include "ldpc_ece535a/image_sink.h"
%}


%include "ldpc_ece535a/ldpc_encoder_bc.h"
GR_SWIG_BLOCK_MAGIC2(ldpc_ece535a, ldpc_encoder_bc);
%include "ldpc_ece535a/ldpc_decoder_cb.h"
GR_SWIG_BLOCK_MAGIC2(ldpc_ece535a, ldpc_decoder_cb);
%include "ldpc_ece535a/image_sink.h"
GR_SWIG_BLOCK_MAGIC2(ldpc_ece535a, image_sink);
