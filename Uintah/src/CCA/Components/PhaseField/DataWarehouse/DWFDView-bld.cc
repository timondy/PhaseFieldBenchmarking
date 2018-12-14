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

#include <CCA/Components/PhaseField/DataWarehouse/DWFDView.h>
#include <CCA/Components/PhaseField/DataWarehouse/DWFDViewFactory.h>

namespace Uintah {
namespace PhaseField {

template<> DWFactoryFDView < ScalarField<const double>, P3 >::FactoryMap DWFactoryFDView < ScalarField<const double>, P3 >::RegisteredNames = {};
template<> DWFactoryFDView < ScalarField<const double>, P5 >::FactoryMap DWFactoryFDView < ScalarField<const double>, P5 >::RegisteredNames = {};
template<> DWFactoryFDView < ScalarField<const double>, P7 >::FactoryMap DWFactoryFDView < ScalarField<const double>, P7 >::RegisteredNames = {};

template<> DWFactoryFDView < VectorField<const double, 1>, P5 >::FactoryMap DWFactoryFDView < VectorField<const double, 1>, P5 >::RegisteredNames = {};

template<> DWFactoryFDView < VectorField<const double, 3>, P7 >::FactoryMap DWFactoryFDView < VectorField<const double, 3>, P7 >::RegisteredNames = {};

template<> const std::string DWFDView < ScalarField<const double>, P3, CC >::Name = "CC|P3";
template<> const std::string DWFDView < ScalarField<const double>, P3, NC >::Name = "NC|P3";
template<> const std::string DWFDView < ScalarField<const double>, P5, CC >::Name = "CC|P5";
template<> const std::string DWFDView < ScalarField<const double>, P5, NC >::Name = "NC|P5";
template<> const std::string DWFDView < ScalarField<const double>, P7, CC >::Name = "CC|P7";
template<> const std::string DWFDView < ScalarField<const double>, P7, NC >::Name = "NC|P7";

template<> const std::string DWFDView < VectorField<const double, 1>, P5, CC >::Name = "CC|P5";
template<> const std::string DWFDView < VectorField<const double, 1>, P5, NC >::Name = "NC|P5";

template<> const std::string DWFDView < VectorField<const double, 3>, P7, CC >::Name = "CC|P7";
template<> const std::string DWFDView < VectorField<const double, 3>, P7, NC >::Name = "NC|P7";

template class DWFDView < ScalarField<const double>, P3, CC >;
template class DWFDView < ScalarField<const double>, P3, NC >;
template class DWFDView < ScalarField<const double>, P5, CC >;
template class DWFDView < ScalarField<const double>, P5, NC >;
template class DWFDView < ScalarField<const double>, P7, CC >;
template class DWFDView < ScalarField<const double>, P7, NC >;

template class DWFDView < VectorField<const double, 1>, P5, CC >;
template class DWFDView < VectorField<const double, 1>, P5, NC >;

template class DWFDView < VectorField<const double, 3>, P7, CC >;
template class DWFDView < VectorField<const double, 3>, P7, NC >;

} // namespace Uintah
} // namespace PhaseField
