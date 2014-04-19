/* -*- c++ -*- */
/* 
 * Copyright 2014 <+YOU OR YOUR COMPANY+>.
 * 
 * This is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 * 
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this software; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street,
 * Boston, MA 02110-1301, USA.
 */

#ifndef INCLUDED_LDPC_ECE535A_LDPC_ENCODER_BB_IMPL_H
#define INCLUDED_LDPC_ECE535A_LDPC_ENCODER_BB_IMPL_H

#include <ldpc_ece535a/ldpc_encoder_bb.h>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <suitesparse/umfpack.h>

namespace gr {
  namespace ldpc_ece535a {

    using namespace boost::numeric;

    class ldpc_encoder_bb_impl : public ldpc_encoder_bb
    {
     private:
      unsigned int d_M;
      unsigned int d_N;

      ublas::matrix<int> d_H;
      ublas::matrix<int> d_L;
      ublas::matrix<int> d_U;

      ublas::vector<int> solve(const ublas::matrix<int> &A,
                               const ublas::vector<int> &B);
      void reorderHMatrix(ublas::matrix<int> &H,
                          ublas::matrix<int> &L,
                          ublas::matrix<int> &U);
      ublas::vector<int> makeParityCheck(const ublas::vector<int> &dSource,
                                         const ublas::matrix<int> &H,
                                         const ublas::matrix<int> &L,
                                         const ublas::matrix<int> &U);

      ublas::vector<int> mod2(const ublas::vector<int> &u);

     public:
      ldpc_encoder_bb_impl();
      ~ldpc_encoder_bb_impl();

      // Where all the action really happens
      void forecast (int noutput_items, gr_vector_int &ninput_items_required);

      int general_work(int noutput_items,
		       gr_vector_int &ninput_items,
		       gr_vector_const_void_star &input_items,
		       gr_vector_void_star &output_items);
    };

  } // namespace ldpc_ece535a
} // namespace gr

#endif /* INCLUDED_LDPC_ECE535A_LDPC_ENCODER_BB_IMPL_H */

