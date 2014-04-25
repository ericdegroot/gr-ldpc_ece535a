/* -*- c++ -*- */
/* 
 * ECE 535A Project
 * Image sink block
 */

#ifndef INCLUDED_LDPC_ECE535A_IMAGE_SINK_H
#define INCLUDED_LDPC_ECE535A_IMAGE_SINK_H

#include <ldpc_ece535a/api.h>
#include <gnuradio/sync_block.h>

namespace gr {
  namespace ldpc_ece535a {

    /*!
     * \brief Image sink
     * \ingroup ldpc_ece535a
     *
     */
    class LDPC_ECE535A_API image_sink : virtual public gr::sync_block
    {
     public:
      typedef boost::shared_ptr<image_sink> sptr;

      /*!
       * \brief Return a shared_ptr to a new instance of ldpc_ece535a::image_sink.
       *
       * To avoid accidental use of raw pointers, ldpc_ece535a::image_sink's
       * constructor is in a private implementation
       * class. ldpc_ece535a::image_sink::make is the public interface for
       * creating new instances.
       */
      static sptr make();
    };

  } // namespace ldpc_ece535a
} // namespace gr

#endif /* INCLUDED_LDPC_ECE535A_IMAGE_SINK_H */

