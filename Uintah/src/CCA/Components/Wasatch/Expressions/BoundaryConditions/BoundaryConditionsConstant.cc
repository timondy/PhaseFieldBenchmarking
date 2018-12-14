/*
 * The MIT License
 *
 * Copyright (c) 2012-2018 The University of Utah
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

/*
 *  NOTE: we have split the implementations into several files to reduce the size
 *        of object files when compiling with nvcc, since we were crashing the
 *        linker in some cases.
 */

#include <CCA/Components/Wasatch/Expressions/BoundaryConditions/BoundaryConditions.h>
#include <spatialops/structured/SpatialMask.h>

namespace WasatchCore{

  // ###################################################################
  // Generic constant BC
  template< typename FieldT >
  void
  ConstantBC<FieldT>::
  evaluate()
  {
    using namespace SpatialOps;
    if (this->spatialMask_) {
      FieldT& lhs = this->value();
      APPLY_CONSTANT_BC(lhs, bcValue_);
    }
  }

  // ###################################################################
  // ConstantBC Specialization for normal fluxes, X, Y, and Z. This works for DIRICHLET conditions ONLY on the normal fluxes
  template<>
  void
  ConstantBC<SpatialOps::SSurfXField>::
  evaluate()
  {
    typedef SpatialOps::SSurfXField FieldT;
    FieldT& lhs = this->value();
    masked_assign( SpatialOps::convert<FieldT>( *(this->svolSpatialMask_), SpatialOps::MINUS_SIDE, SpatialOps::PLUS_SIDE), lhs, bcValue_);
  }
  
  template<>
  void
  ConstantBC<SpatialOps::SSurfYField>::
  evaluate()
  {
    typedef SpatialOps::SSurfYField FieldT;
    FieldT& lhs = this->value();
    masked_assign( SpatialOps::convert<FieldT>( *(this->svolSpatialMask_), SpatialOps::MINUS_SIDE, SpatialOps::PLUS_SIDE), lhs, bcValue_);
  }

  template<>
  void
  ConstantBC<SpatialOps::SSurfZField>::
  evaluate()
  {
    typedef SpatialOps::SSurfZField FieldT;
    FieldT& lhs = this->value();
    masked_assign( SpatialOps::convert<FieldT>( *(this->svolSpatialMask_), SpatialOps::MINUS_SIDE, SpatialOps::PLUS_SIDE), lhs, bcValue_);
  }

  // ###################################################################
  // a necessary specialization for particle fields because the BCHelper automatically creates
  // ConstantBC for Dirichlet boundary conditions specified in the input file.
  template<>
  void
  ConstantBC<ParticleField>::
  evaluate()
  {}

  // ###################################################################
  // EXPLICIT INSTANTIATION
#include <CCA/Components/Wasatch/FieldTypes.h>

#define INSTANTIATE_BC_PROFILES(VOLT)   \
    template class ConstantBC<VOLT>;

  INSTANTIATE_BC_PROFILES(SVolField)
  INSTANTIATE_BC_PROFILES(XVolField)
  INSTANTIATE_BC_PROFILES(YVolField)
  INSTANTIATE_BC_PROFILES(ZVolField)

  template class ConstantBC<SpatialOps::SSurfXField>;
  template class ConstantBC<SpatialOps::SSurfYField>;
  template class ConstantBC<SpatialOps::SSurfZField>;

} // namespace WasatchCore
