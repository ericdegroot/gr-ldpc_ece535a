/* -*- c++ -*- */
/* 
 * ECE 535A Project
 * LDPC encoder block
 */


#ifndef INCLUDED_LDPC_ECE535A_LDPC_ENCODER_BC_H
#define INCLUDED_LDPC_ECE535A_LDPC_ENCODER_BC_H

#include <ldpc_ece535a/api.h>
#include <gnuradio/block.h>

namespace gr {
  namespace ldpc_ece535a {

    /*!
     * \brief LDPC encoder block
     * \ingroup ldpc_ece535a
     *
     */
    class LDPC_ECE535A_API ldpc_encoder_bc : virtual public gr::block
    {
     public:
      typedef boost::shared_ptr<ldpc_encoder_bc> sptr;

      /*!
       * \brief Return a shared_ptr to a new instance of ldpc_ece535a::ldpc_encoder_bc.
       *
       * To avoid accidental use of raw pointers, ldpc_ece535a::ldpc_encoder_bc's
       * constructor is in a private implementation
       * class. ldpc_ece535a::ldpc_encoder_bc::make is the public interface for
       * creating new instances.
       */
      static sptr make();
    };

  } // namespace ldpc_ece535a
} // namespace gr

#endif /* INCLUDED_LDPC_ECE535A_LDPC_ENCODER_BC_H */

