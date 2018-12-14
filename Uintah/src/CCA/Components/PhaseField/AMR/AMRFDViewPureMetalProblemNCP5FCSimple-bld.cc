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
#include <CCA/Components/PhaseField/BoundaryConditions/BCFDView.h>

namespace Uintah {
namespace PhaseField {

template<> const std::string BCFDView < PureMetalProblem<NC, P5>, 0, Patch::xminus | BC::FineCoarseInterface | FC::FCSimple >::Name = "PureMetalProblem|0|NC|xminus|FCSimple|";
template<> const std::string BCFDView < PureMetalProblem<NC, P5>, 0, Patch::xminus | BC::FineCoarseInterface | FC::FCSimple, Patch::yminus | BC::FineCoarseInterface | FC::FCSimple >::Name = "PureMetalProblem|0|NC|xminus|FCSimple|yminus|FCSimple|";
template<> const std::string BCFDView < PureMetalProblem<NC, P5>, 0, Patch::xminus | BC::FineCoarseInterface | FC::FCSimple, Patch::yplus | BC::FineCoarseInterface | FC::FCSimple >::Name = "PureMetalProblem|0|NC|xminus|FCSimple|yplus|FCSimple|";
template<> const std::string BCFDView < PureMetalProblem<NC, P5>, 0, Patch::xplus | BC::FineCoarseInterface | FC::FCSimple >::Name = "PureMetalProblem|0|NC|xplus|FCSimple|";
template<> const std::string BCFDView < PureMetalProblem<NC, P5>, 0, Patch::xplus | BC::FineCoarseInterface | FC::FCSimple, Patch::yminus | BC::FineCoarseInterface | FC::FCSimple >::Name = "PureMetalProblem|0|NC|xplus|FCSimple|yminus|FCSimple|";
template<> const std::string BCFDView < PureMetalProblem<NC, P5>, 0, Patch::xplus | BC::FineCoarseInterface | FC::FCSimple, Patch::yplus | BC::FineCoarseInterface | FC::FCSimple >::Name = "PureMetalProblem|0|NC|xplus|FCSimple|yplus|FCSimple|";
template<> const std::string BCFDView < PureMetalProblem<NC, P5>, 0, Patch::yminus | BC::FineCoarseInterface | FC::FCSimple >::Name = "PureMetalProblem|0|NC|yminus|FCSimple|";
template<> const std::string BCFDView < PureMetalProblem<NC, P5>, 0, Patch::yplus | BC::FineCoarseInterface | FC::FCSimple >::Name = "PureMetalProblem|0|NC|yplus|FCSimple|";
template<> const std::string BCFDView < PureMetalProblem<NC, P5>, 1, Patch::xminus | BC::FineCoarseInterface | FC::FCSimple >::Name = "PureMetalProblem|1|NC|xminus|FCSimple|";
template<> const std::string BCFDView < PureMetalProblem<NC, P5>, 1, Patch::xminus | BC::FineCoarseInterface | FC::FCSimple, Patch::yminus | BC::FineCoarseInterface | FC::FCSimple >::Name = "PureMetalProblem|1|NC|xminus|FCSimple|yminus|FCSimple|";
template<> const std::string BCFDView < PureMetalProblem<NC, P5>, 1, Patch::xminus | BC::FineCoarseInterface | FC::FCSimple, Patch::yplus | BC::FineCoarseInterface | FC::FCSimple >::Name = "PureMetalProblem|1|NC|xminus|FCSimple|yplus|FCSimple|";
template<> const std::string BCFDView < PureMetalProblem<NC, P5>, 1, Patch::xplus | BC::FineCoarseInterface | FC::FCSimple >::Name = "PureMetalProblem|1|NC|xplus|FCSimple|";
template<> const std::string BCFDView < PureMetalProblem<NC, P5>, 1, Patch::xplus | BC::FineCoarseInterface | FC::FCSimple, Patch::yminus | BC::FineCoarseInterface | FC::FCSimple >::Name = "PureMetalProblem|1|NC|xplus|FCSimple|yminus|FCSimple|";
template<> const std::string BCFDView < PureMetalProblem<NC, P5>, 1, Patch::xplus | BC::FineCoarseInterface | FC::FCSimple, Patch::yplus | BC::FineCoarseInterface | FC::FCSimple >::Name = "PureMetalProblem|1|NC|xplus|FCSimple|yplus|FCSimple|";
template<> const std::string BCFDView < PureMetalProblem<NC, P5>, 1, Patch::yminus | BC::FineCoarseInterface | FC::FCSimple >::Name = "PureMetalProblem|1|NC|yminus|FCSimple|";
template<> const std::string BCFDView < PureMetalProblem<NC, P5>, 1, Patch::yplus | BC::FineCoarseInterface | FC::FCSimple >::Name = "PureMetalProblem|1|NC|yplus|FCSimple|";
template<> const std::string BCFDView < PureMetalProblem<NC, P5>, 2, Patch::xminus | BC::FineCoarseInterface | FC::FCSimple >::Name = "PureMetalProblem|2|NC|xminus|FCSimple|";
template<> const std::string BCFDView < PureMetalProblem<NC, P5>, 2, Patch::xminus | BC::FineCoarseInterface | FC::FCSimple, Patch::yminus | BC::FineCoarseInterface | FC::FCSimple >::Name = "PureMetalProblem|2|NC|xminus|FCSimple|yminus|FCSimple|";
template<> const std::string BCFDView < PureMetalProblem<NC, P5>, 2, Patch::xminus | BC::FineCoarseInterface | FC::FCSimple, Patch::yplus | BC::FineCoarseInterface | FC::FCSimple >::Name = "PureMetalProblem|2|NC|xminus|FCSimple|yplus|FCSimple|";
template<> const std::string BCFDView < PureMetalProblem<NC, P5>, 2, Patch::xplus | BC::FineCoarseInterface | FC::FCSimple >::Name = "PureMetalProblem|2|NC|xplus|FCSimple|";
template<> const std::string BCFDView < PureMetalProblem<NC, P5>, 2, Patch::xplus | BC::FineCoarseInterface | FC::FCSimple, Patch::yminus | BC::FineCoarseInterface | FC::FCSimple >::Name = "PureMetalProblem|2|NC|xplus|FCSimple|yminus|FCSimple|";
template<> const std::string BCFDView < PureMetalProblem<NC, P5>, 2, Patch::xplus | BC::FineCoarseInterface | FC::FCSimple, Patch::yplus | BC::FineCoarseInterface | FC::FCSimple >::Name = "PureMetalProblem|2|NC|xplus|FCSimple|yplus|FCSimple|";
template<> const std::string BCFDView < PureMetalProblem<NC, P5>, 2, Patch::yminus | BC::FineCoarseInterface | FC::FCSimple >::Name = "PureMetalProblem|2|NC|yminus|FCSimple|";
template<> const std::string BCFDView < PureMetalProblem<NC, P5>, 2, Patch::yplus | BC::FineCoarseInterface | FC::FCSimple >::Name = "PureMetalProblem|2|NC|yplus|FCSimple|";
template<> const std::string BCFDView < PureMetalProblem<NC, P5>, 3, Patch::xminus | BC::FineCoarseInterface | FC::FCSimple >::Name = "PureMetalProblem|3|NC|xminus|FCSimple|";
template<> const std::string BCFDView < PureMetalProblem<NC, P5>, 3, Patch::xminus | BC::FineCoarseInterface | FC::FCSimple, Patch::yminus | BC::FineCoarseInterface | FC::FCSimple >::Name = "PureMetalProblem|3|NC|xminus|FCSimple|yminus|FCSimple|";
template<> const std::string BCFDView < PureMetalProblem<NC, P5>, 3, Patch::xminus | BC::FineCoarseInterface | FC::FCSimple, Patch::yplus | BC::FineCoarseInterface | FC::FCSimple >::Name = "PureMetalProblem|3|NC|xminus|FCSimple|yplus|FCSimple|";
template<> const std::string BCFDView < PureMetalProblem<NC, P5>, 3, Patch::xplus | BC::FineCoarseInterface | FC::FCSimple >::Name = "PureMetalProblem|3|NC|xplus|FCSimple|";
template<> const std::string BCFDView < PureMetalProblem<NC, P5>, 3, Patch::xplus | BC::FineCoarseInterface | FC::FCSimple, Patch::yminus | BC::FineCoarseInterface | FC::FCSimple >::Name = "PureMetalProblem|3|NC|xplus|FCSimple|yminus|FCSimple|";
template<> const std::string BCFDView < PureMetalProblem<NC, P5>, 3, Patch::xplus | BC::FineCoarseInterface | FC::FCSimple, Patch::yplus | BC::FineCoarseInterface | FC::FCSimple >::Name = "PureMetalProblem|3|NC|xplus|FCSimple|yplus|FCSimple|";
template<> const std::string BCFDView < PureMetalProblem<NC, P5>, 3, Patch::yminus | BC::FineCoarseInterface | FC::FCSimple >::Name = "PureMetalProblem|3|NC|yminus|FCSimple|";
template<> const std::string BCFDView < PureMetalProblem<NC, P5>, 3, Patch::yplus | BC::FineCoarseInterface | FC::FCSimple >::Name = "PureMetalProblem|3|NC|yplus|FCSimple|";

template class BCFDView < PureMetalProblem<NC, P5>, 0, Patch::xminus | BC::FineCoarseInterface | FC::FCSimple >;
template class BCFDView < PureMetalProblem<NC, P5>, 0, Patch::xminus | BC::FineCoarseInterface | FC::FCSimple, Patch::yminus | BC::FineCoarseInterface | FC::FCSimple >;
template class BCFDView < PureMetalProblem<NC, P5>, 0, Patch::xminus | BC::FineCoarseInterface | FC::FCSimple, Patch::yplus | BC::FineCoarseInterface | FC::FCSimple >;
template class BCFDView < PureMetalProblem<NC, P5>, 0, Patch::xplus | BC::FineCoarseInterface | FC::FCSimple >;
template class BCFDView < PureMetalProblem<NC, P5>, 0, Patch::xplus | BC::FineCoarseInterface | FC::FCSimple, Patch::yminus | BC::FineCoarseInterface | FC::FCSimple >;
template class BCFDView < PureMetalProblem<NC, P5>, 0, Patch::xplus | BC::FineCoarseInterface | FC::FCSimple, Patch::yplus | BC::FineCoarseInterface | FC::FCSimple >;
template class BCFDView < PureMetalProblem<NC, P5>, 0, Patch::yminus | BC::FineCoarseInterface | FC::FCSimple >;
template class BCFDView < PureMetalProblem<NC, P5>, 0, Patch::yplus | BC::FineCoarseInterface | FC::FCSimple >;
template class BCFDView < PureMetalProblem<NC, P5>, 1, Patch::xminus | BC::FineCoarseInterface | FC::FCSimple >;
template class BCFDView < PureMetalProblem<NC, P5>, 1, Patch::xminus | BC::FineCoarseInterface | FC::FCSimple, Patch::yminus | BC::FineCoarseInterface | FC::FCSimple >;
template class BCFDView < PureMetalProblem<NC, P5>, 1, Patch::xminus | BC::FineCoarseInterface | FC::FCSimple, Patch::yplus | BC::FineCoarseInterface | FC::FCSimple >;
template class BCFDView < PureMetalProblem<NC, P5>, 1, Patch::xplus | BC::FineCoarseInterface | FC::FCSimple >;
template class BCFDView < PureMetalProblem<NC, P5>, 1, Patch::xplus | BC::FineCoarseInterface | FC::FCSimple, Patch::yminus | BC::FineCoarseInterface | FC::FCSimple >;
template class BCFDView < PureMetalProblem<NC, P5>, 1, Patch::xplus | BC::FineCoarseInterface | FC::FCSimple, Patch::yplus | BC::FineCoarseInterface | FC::FCSimple >;
template class BCFDView < PureMetalProblem<NC, P5>, 1, Patch::yminus | BC::FineCoarseInterface | FC::FCSimple >;
template class BCFDView < PureMetalProblem<NC, P5>, 1, Patch::yplus | BC::FineCoarseInterface | FC::FCSimple >;
template class BCFDView < PureMetalProblem<NC, P5>, 2, Patch::xminus | BC::FineCoarseInterface | FC::FCSimple >;
template class BCFDView < PureMetalProblem<NC, P5>, 2, Patch::xminus | BC::FineCoarseInterface | FC::FCSimple, Patch::yminus | BC::FineCoarseInterface | FC::FCSimple >;
template class BCFDView < PureMetalProblem<NC, P5>, 2, Patch::xminus | BC::FineCoarseInterface | FC::FCSimple, Patch::yplus | BC::FineCoarseInterface | FC::FCSimple >;
template class BCFDView < PureMetalProblem<NC, P5>, 2, Patch::xplus | BC::FineCoarseInterface | FC::FCSimple >;
template class BCFDView < PureMetalProblem<NC, P5>, 2, Patch::xplus | BC::FineCoarseInterface | FC::FCSimple, Patch::yminus | BC::FineCoarseInterface | FC::FCSimple >;
template class BCFDView < PureMetalProblem<NC, P5>, 2, Patch::xplus | BC::FineCoarseInterface | FC::FCSimple, Patch::yplus | BC::FineCoarseInterface | FC::FCSimple >;
template class BCFDView < PureMetalProblem<NC, P5>, 2, Patch::yminus | BC::FineCoarseInterface | FC::FCSimple >;
template class BCFDView < PureMetalProblem<NC, P5>, 2, Patch::yplus | BC::FineCoarseInterface | FC::FCSimple >;
template class BCFDView < PureMetalProblem<NC, P5>, 3, Patch::xminus | BC::FineCoarseInterface | FC::FCSimple >;
template class BCFDView < PureMetalProblem<NC, P5>, 3, Patch::xminus | BC::FineCoarseInterface | FC::FCSimple, Patch::yminus | BC::FineCoarseInterface | FC::FCSimple >;
template class BCFDView < PureMetalProblem<NC, P5>, 3, Patch::xminus | BC::FineCoarseInterface | FC::FCSimple, Patch::yplus | BC::FineCoarseInterface | FC::FCSimple >;
template class BCFDView < PureMetalProblem<NC, P5>, 3, Patch::xplus | BC::FineCoarseInterface | FC::FCSimple >;
template class BCFDView < PureMetalProblem<NC, P5>, 3, Patch::xplus | BC::FineCoarseInterface | FC::FCSimple, Patch::yminus | BC::FineCoarseInterface | FC::FCSimple >;
template class BCFDView < PureMetalProblem<NC, P5>, 3, Patch::xplus | BC::FineCoarseInterface | FC::FCSimple, Patch::yplus | BC::FineCoarseInterface | FC::FCSimple >;
template class BCFDView < PureMetalProblem<NC, P5>, 3, Patch::yminus | BC::FineCoarseInterface | FC::FCSimple >;
template class BCFDView < PureMetalProblem<NC, P5>, 3, Patch::yplus | BC::FineCoarseInterface | FC::FCSimple >;

} // namespace Uintah
} // namespace PhaseField
