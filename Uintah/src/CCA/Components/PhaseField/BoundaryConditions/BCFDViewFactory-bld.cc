/*
 * The MIT License
 *
 * Copyright (c) 1997-2018 The University of Utah
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

#include <CCA/Components/PhaseField/DataTypes/PureMetalProblem.h>
#include <CCA/Components/PhaseField/DataTypes/HeatProblem.h>
#include <CCA/Components/PhaseField/BoundaryConditions/BCFDView.h>
#include <CCA/Components/PhaseField/BoundaryConditions/BCFDViewFactory.h>

namespace Uintah {
namespace PhaseField {

template<> BCFactoryFDView < ScalarField<const double>, P3 >::FactoryMap BCFactoryFDView < ScalarField<const double>, P3 >::RegisteredNames = {};
template<> BCFactoryFDView < ScalarField<const double>, P5 >::FactoryMap BCFactoryFDView < ScalarField<const double>, P5 >::RegisteredNames = {};
template<> BCFactoryFDView < ScalarField<const double>, P7 >::FactoryMap BCFactoryFDView < ScalarField<const double>, P7 >::RegisteredNames = {};

template<> BCFactoryFDView < VectorField<const double, 1>, P5 >::FactoryMap BCFactoryFDView < VectorField<const double, 1>, P5 >::RegisteredNames = {};

template<> BCFactoryFDView < VectorField<const double, 3>, P7 >::FactoryMap BCFactoryFDView < VectorField<const double, 3>, P7 >::RegisteredNames = {};

} // namespace Uintah
} // namespace PhaseField
