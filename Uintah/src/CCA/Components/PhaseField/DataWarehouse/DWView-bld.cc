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

#include <CCA/Components/PhaseField/DataWarehouse/DWView.h>
#include <CCA/Components/PhaseField/DataWarehouse/DWViewFactory.h>

namespace Uintah {
namespace PhaseField {

template<> DWFactoryView < ScalarField<int> >::FactoryMap DWFactoryView < ScalarField<int> >::RegisteredNames = {};
template<> DWFactoryView < ScalarField<double> >::FactoryMap DWFactoryView < ScalarField<double> >::RegisteredNames = {};
template<> DWFactoryView < ScalarField<const double> >::FactoryMap DWFactoryView < ScalarField<const double> >::RegisteredNames = {};
template<> DWFactoryView < ScalarField<Stencil7> >::FactoryMap DWFactoryView < ScalarField<Stencil7> >::RegisteredNames = {};
template<> DWFactoryView < ScalarField<const Stencil7> >::FactoryMap DWFactoryView < ScalarField<const Stencil7> >::RegisteredNames = {};
template<> DWFactoryView < VectorField<double, 1u> >::FactoryMap DWFactoryView < VectorField<double, 1u> >::RegisteredNames = {};
template<> DWFactoryView < VectorField<double, 2u> >::FactoryMap DWFactoryView < VectorField<double, 2u> >::RegisteredNames = {};
template<> DWFactoryView < VectorField<const double, 2u> >::FactoryMap DWFactoryView < VectorField<const double, 2u> >::RegisteredNames = {};
template<> DWFactoryView < VectorField<double, 3u> >::FactoryMap DWFactoryView < VectorField<double, 3u> >::RegisteredNames = {};
template<> DWFactoryView < VectorField<const double, 3u> >::FactoryMap DWFactoryView < VectorField<const double, 3u> >::RegisteredNames = {};

template<> const std::string DWView < ScalarField<int>, CC, D2 >::Name = "CC|D2|";
template<> const std::string DWView < ScalarField<int>, CC, D3 >::Name = "CC|D3|";
template<> const std::string DWView < ScalarField<int>, NC, D2 >::Name = "NC|D2|";
template<> const std::string DWView < ScalarField<int>, NC, D3 >::Name = "NC|D3|";

template<> const std::string DWView < ScalarField<double>, CC, D1 >::Name = "CC|D1|";
template<> const std::string DWView < ScalarField<double>, CC, D2 >::Name = "CC|D2|";
template<> const std::string DWView < ScalarField<double>, CC, D3 >::Name = "CC|D3|";
template<> const std::string DWView < ScalarField<double>, NC, D1 >::Name = "NC|D1|";
template<> const std::string DWView < ScalarField<double>, NC, D2 >::Name = "NC|D2|";
template<> const std::string DWView < ScalarField<double>, NC, D3 >::Name = "NC|D3|";

template<> const std::string DWView < ScalarField<const double>, CC, D1 >::Name = "CC|D1|";
template<> const std::string DWView < ScalarField<const double>, CC, D2 >::Name = "CC|D2|";
template<> const std::string DWView < ScalarField<const double>, CC, D3 >::Name = "CC|D3|";
template<> const std::string DWView < ScalarField<const double>, NC, D1 >::Name = "NC|D1|";
template<> const std::string DWView < ScalarField<const double>, NC, D2 >::Name = "NC|D2|";
template<> const std::string DWView < ScalarField<const double>, NC, D3 >::Name = "NC|D3|";

template<> const std::string DWView < ScalarField<Stencil7>, CC, D2 >::Name = "CC|D2|";
template<> const std::string DWView < ScalarField<Stencil7>, CC, D3 >::Name = "CC|D3|";
template<> const std::string DWView < ScalarField<Stencil7>, NC, D2 >::Name = "NC|D2|";
template<> const std::string DWView < ScalarField<Stencil7>, NC, D3 >::Name = "NC|D3|";

template<> const std::string DWView < ScalarField<const Stencil7>, CC, D2 >::Name = "CC|D2|";
template<> const std::string DWView < ScalarField<const Stencil7>, CC, D3 >::Name = "CC|D3|";
template<> const std::string DWView < ScalarField<const Stencil7>, NC, D2 >::Name = "NC|D2|";
template<> const std::string DWView < ScalarField<const Stencil7>, NC, D3 >::Name = "NC|D3|";

template<> const std::string DWView < VectorField<double, 1u>, CC, D2 >::Name = "CC|D2|";
template<> const std::string DWView < VectorField<double, 1u>, NC, D2 >::Name = "NC|D2|";

template<> const std::string DWView < VectorField<double, 2u>, CC, D2 >::Name = "CC|D2|";
template<> const std::string DWView < VectorField<double, 2u>, NC, D2 >::Name = "NC|D2|";

template<> const std::string DWView < VectorField<const double, 2u>, CC, D2 >::Name = "CC|D2|";
template<> const std::string DWView < VectorField<const double, 2u>, NC, D2 >::Name = "NC|D2|";

template<> const std::string DWView < VectorField<double, 3u>, CC, D3 >::Name = "CC|D3|";
template<> const std::string DWView < VectorField<double, 3u>, NC, D3 >::Name = "NC|D3|";

template<> const std::string DWView < VectorField<const double, 3u>, CC, D3 >::Name = "CC|D3|";
template<> const std::string DWView < VectorField<const double, 3u>, NC, D3 >::Name = "NC|D3|";

template class DWView < ScalarField<int>, CC, D2 >;
template class DWView < ScalarField<int>, CC, D3 >;
template class DWView < ScalarField<int>, NC, D2 >;
template class DWView < ScalarField<int>, NC, D3 >;

template class DWView < ScalarField<double>, CC, D1 >;
template class DWView < ScalarField<double>, CC, D2 >;
template class DWView < ScalarField<double>, CC, D3 >;
template class DWView < ScalarField<double>, NC, D1 >;
template class DWView < ScalarField<double>, NC, D2 >;
template class DWView < ScalarField<double>, NC, D3 >;

template class DWView < ScalarField<const double>, CC, D1 >;
template class DWView < ScalarField<const double>, CC, D2 >;
template class DWView < ScalarField<const double>, CC, D3 >;
template class DWView < ScalarField<const double>, NC, D1 >;
template class DWView < ScalarField<const double>, NC, D2 >;
template class DWView < ScalarField<const double>, NC, D3 >;

template class DWView < ScalarField<Stencil7>, CC, D2 >;
template class DWView < ScalarField<Stencil7>, CC, D3 >;
template class DWView < ScalarField<Stencil7>, NC, D2 >;
template class DWView < ScalarField<Stencil7>, NC, D3 >;

template class DWView < ScalarField<const Stencil7>, CC, D2 >;
template class DWView < ScalarField<const Stencil7>, CC, D3 >;
template class DWView < ScalarField<const Stencil7>, NC, D2 >;
template class DWView < ScalarField<const Stencil7>, NC, D3 >;

template class DWView < VectorField<double, 1u>, CC, D2 >;
template class DWView < VectorField<double, 1u>, NC, D2 >;

template class DWView < VectorField<double, 2u>, CC, D2 >;
template class DWView < VectorField<double, 2u>, NC, D2 >;

template class DWView < VectorField<const double, 2u>, CC, D2 >;
template class DWView < VectorField<const double, 2u>, NC, D2 >;

template class DWView < VectorField<double, 3u>, CC, D3 >;
template class DWView < VectorField<double, 3u>, NC, D3 >;

template class DWView < VectorField<const double, 3u>, CC, D3 >;
template class DWView < VectorField<const double, 3u>, NC, D3 >;

} // namespace Uintah
} // namespace PhaseField
