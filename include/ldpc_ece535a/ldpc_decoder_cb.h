/* -*- c++ -*- */
/* 
 * ECE 535A Project
 * LDPC decoder block
 */


#ifndef INCLUDED_LDPC_ECE535A_LDPC_DECODER_CB_H
#define INCLUDED_LDPC_ECE535A_LDPC_DECODER_CB_H

#include <ldpc_ece535a/api.h>
#include <gnuradio/block.h>

namespace gr {
  namespace ldpc_ece535a {

    /*!
     * \brief LDPC decoder block
     * \ingroup ldpc_ece535a
     *
     */
    class LDPC_ECE535A_API ldpc_decoder_cb : virtual public gr::block
    {
     public:
      typedef boost::shared_ptr<ldpc_decoder_cb> sptr;

      /*!
       * \brief Return a shared_ptr to a new instance of ldpc_ece535a::ldpc_decoder_cb.
       *
       * To avoid accidental use of raw pointers, ldpc_ece535a::ldpc_decoder_cb's
       * constructor is in a private implementation
       * class. ldpc_ece535a::ldpc_decoder_cb::make is the public interface for
       * creating new instances.
       */
      static sptr make(const int method);
    };

  } // namespace ldpc_ece535a
} // namespace gr

#endif /* INCLUDED_LDPC_ECE535A_LDPC_DECODER_CB_H */

