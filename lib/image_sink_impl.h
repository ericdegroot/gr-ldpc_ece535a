/* -*- c++ -*- */
/* 
 * ECE 535A Project
 * Image sink block
 */

#ifndef INCLUDED_LDPC_ECE535A_IMAGE_SINK_IMPL_H
#define INCLUDED_LDPC_ECE535A_IMAGE_SINK_IMPL_H

#include <ldpc_ece535a/image_sink.h>

#include <string>

namespace gr {
  namespace ldpc_ece535a {

    class image_sink_impl : public image_sink
    {
    private:
      std::string d_buffer;
      unsigned int d_fileSize;

    public:
      image_sink_impl();
      ~image_sink_impl();

      // Where all the action really happens
      int work(int noutput_items,
	       gr_vector_const_void_star &input_items,
	       gr_vector_void_star &output_items);
    };

  } // namespace ldpc_ece535a
} // namespace gr

#endif /* INCLUDED_LDPC_ECE535A_IMAGE_SINK_IMPL_H */

