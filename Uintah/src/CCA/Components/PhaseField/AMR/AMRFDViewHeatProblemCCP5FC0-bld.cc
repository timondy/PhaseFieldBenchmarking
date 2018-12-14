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

#include <CCA/Components/PhaseField/DataTypes/HeatProblem.h>
#include <CCA/Components/PhaseField/BoundaryConditions/BCFDView.h>

namespace Uintah {
namespace PhaseField {

template<> const std::string BCFDView < HeatProblem<CC, P5>, 0, Patch::xminus | BC::FineCoarseInterface | FC::FC0 >::Name = "HeatProblem|0|CC|xminus|FC0|";
template<> const std::string BCFDView < HeatProblem<CC, P5>, 0, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::yminus | BC::FineCoarseInterface | FC::FC0 >::Name = "HeatProblem|0|CC|xminus|FC0|yminus|FC0|";
template<> const std::string BCFDView < HeatProblem<CC, P5>, 0, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::yminus | BC::Dirichlet >::Name = "HeatProblem|0|CC|xminus|FC0|yminus|Dirichlet|";
template<> const std::string BCFDView < HeatProblem<CC, P5>, 0, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::yminus | BC::Neumann >::Name = "HeatProblem|0|CC|xminus|FC0|yminus|Neumann|";
template<> const std::string BCFDView < HeatProblem<CC, P5>, 0, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::yplus | BC::FineCoarseInterface | FC::FC0 >::Name = "HeatProblem|0|CC|xminus|FC0|yplus|FC0|";
template<> const std::string BCFDView < HeatProblem<CC, P5>, 0, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::yplus | BC::Dirichlet >::Name = "HeatProblem|0|CC|xminus|FC0|yplus|Dirichlet|";
template<> const std::string BCFDView < HeatProblem<CC, P5>, 0, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::yplus | BC::Neumann >::Name = "HeatProblem|0|CC|xminus|FC0|yplus|Neumann|";
template<> const std::string BCFDView < HeatProblem<CC, P5>, 0, Patch::xminus | BC::Dirichlet, Patch::yminus | BC::FineCoarseInterface | FC::FC0 >::Name = "HeatProblem|0|CC|xminus|Dirichlet|yminus|FC0|";
template<> const std::string BCFDView < HeatProblem<CC, P5>, 0, Patch::xminus | BC::Dirichlet, Patch::yplus | BC::FineCoarseInterface | FC::FC0 >::Name = "HeatProblem|0|CC|xminus|Dirichlet|yplus|FC0|";
template<> const std::string BCFDView < HeatProblem<CC, P5>, 0, Patch::xminus | BC::Neumann, Patch::yminus | BC::FineCoarseInterface | FC::FC0 >::Name = "HeatProblem|0|CC|xminus|Neumann|yminus|FC0|";
template<> const std::string BCFDView < HeatProblem<CC, P5>, 0, Patch::xminus | BC::Neumann, Patch::yplus | BC::FineCoarseInterface | FC::FC0 >::Name = "HeatProblem|0|CC|xminus|Neumann|yplus|FC0|";
template<> const std::string BCFDView < HeatProblem<CC, P5>, 0, Patch::xplus | BC::FineCoarseInterface | FC::FC0 >::Name = "HeatProblem|0|CC|xplus|FC0|";
template<> const std::string BCFDView < HeatProblem<CC, P5>, 0, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::yminus | BC::FineCoarseInterface | FC::FC0 >::Name = "HeatProblem|0|CC|xplus|FC0|yminus|FC0|";
template<> const std::string BCFDView < HeatProblem<CC, P5>, 0, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::yminus | BC::Dirichlet >::Name = "HeatProblem|0|CC|xplus|FC0|yminus|Dirichlet|";
template<> const std::string BCFDView < HeatProblem<CC, P5>, 0, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::yminus | BC::Neumann >::Name = "HeatProblem|0|CC|xplus|FC0|yminus|Neumann|";
template<> const std::string BCFDView < HeatProblem<CC, P5>, 0, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::yplus | BC::FineCoarseInterface | FC::FC0 >::Name = "HeatProblem|0|CC|xplus|FC0|yplus|FC0|";
template<> const std::string BCFDView < HeatProblem<CC, P5>, 0, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::yplus | BC::Dirichlet >::Name = "HeatProblem|0|CC|xplus|FC0|yplus|Dirichlet|";
template<> const std::string BCFDView < HeatProblem<CC, P5>, 0, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::yplus | BC::Neumann >::Name = "HeatProblem|0|CC|xplus|FC0|yplus|Neumann|";
template<> const std::string BCFDView < HeatProblem<CC, P5>, 0, Patch::xplus | BC::Dirichlet, Patch::yminus | BC::FineCoarseInterface | FC::FC0 >::Name = "HeatProblem|0|CC|xplus|Dirichlet|yminus|FC0|";
template<> const std::string BCFDView < HeatProblem<CC, P5>, 0, Patch::xplus | BC::Dirichlet, Patch::yplus | BC::FineCoarseInterface | FC::FC0 >::Name = "HeatProblem|0|CC|xplus|Dirichlet|yplus|FC0|";
template<> const std::string BCFDView < HeatProblem<CC, P5>, 0, Patch::xplus | BC::Neumann, Patch::yminus | BC::FineCoarseInterface | FC::FC0 >::Name = "HeatProblem|0|CC|xplus|Neumann|yminus|FC0|";
template<> const std::string BCFDView < HeatProblem<CC, P5>, 0, Patch::xplus | BC::Neumann, Patch::yplus | BC::FineCoarseInterface | FC::FC0 >::Name = "HeatProblem|0|CC|xplus|Neumann|yplus|FC0|";
template<> const std::string BCFDView < HeatProblem<CC, P5>, 0, Patch::yminus | BC::FineCoarseInterface | FC::FC0 >::Name = "HeatProblem|0|CC|yminus|FC0|";
template<> const std::string BCFDView < HeatProblem<CC, P5>, 0, Patch::yplus | BC::FineCoarseInterface | FC::FC0 >::Name = "HeatProblem|0|CC|yplus|FC0|";

template class BCFDView < HeatProblem<CC, P5>, 0, Patch::xminus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < HeatProblem<CC, P5>, 0, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::yminus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < HeatProblem<CC, P5>, 0, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::yminus | BC::Dirichlet >;
template class BCFDView < HeatProblem<CC, P5>, 0, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::yminus | BC::Neumann >;
template class BCFDView < HeatProblem<CC, P5>, 0, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::yplus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < HeatProblem<CC, P5>, 0, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::yplus | BC::Dirichlet >;
template class BCFDView < HeatProblem<CC, P5>, 0, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::yplus | BC::Neumann >;
template class BCFDView < HeatProblem<CC, P5>, 0, Patch::xminus | BC::Dirichlet, Patch::yminus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < HeatProblem<CC, P5>, 0, Patch::xminus | BC::Dirichlet, Patch::yplus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < HeatProblem<CC, P5>, 0, Patch::xminus | BC::Neumann, Patch::yminus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < HeatProblem<CC, P5>, 0, Patch::xminus | BC::Neumann, Patch::yplus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < HeatProblem<CC, P5>, 0, Patch::xplus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < HeatProblem<CC, P5>, 0, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::yminus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < HeatProblem<CC, P5>, 0, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::yminus | BC::Dirichlet >;
template class BCFDView < HeatProblem<CC, P5>, 0, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::yminus | BC::Neumann >;
template class BCFDView < HeatProblem<CC, P5>, 0, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::yplus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < HeatProblem<CC, P5>, 0, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::yplus | BC::Dirichlet >;
template class BCFDView < HeatProblem<CC, P5>, 0, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::yplus | BC::Neumann >;
template class BCFDView < HeatProblem<CC, P5>, 0, Patch::xplus | BC::Dirichlet, Patch::yminus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < HeatProblem<CC, P5>, 0, Patch::xplus | BC::Dirichlet, Patch::yplus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < HeatProblem<CC, P5>, 0, Patch::xplus | BC::Neumann, Patch::yminus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < HeatProblem<CC, P5>, 0, Patch::xplus | BC::Neumann, Patch::yplus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < HeatProblem<CC, P5>, 0, Patch::yminus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < HeatProblem<CC, P5>, 0, Patch::yplus | BC::FineCoarseInterface | FC::FC0 >;

} // namespace Uintah
} // namespace PhaseField
