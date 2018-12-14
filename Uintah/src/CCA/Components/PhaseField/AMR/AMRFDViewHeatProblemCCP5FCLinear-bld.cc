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

template<> const std::string BCFDView < HeatProblem<CC, P5>, 0, Patch::xminus | BC::FineCoarseInterface | FC::FCLinear >::Name = "HeatProblem|0|CC|xminus|FCLinear|";
template<> const std::string BCFDView < HeatProblem<CC, P5>, 0, Patch::xminus | BC::FineCoarseInterface | FC::FCLinear, Patch::yminus | BC::FineCoarseInterface | FC::FCLinear >::Name = "HeatProblem|0|CC|xminus|FCLinear|yminus|FCLinear|";
template<> const std::string BCFDView < HeatProblem<CC, P5>, 0, Patch::xminus | BC::FineCoarseInterface | FC::FCLinear, Patch::yminus | BC::Dirichlet >::Name = "HeatProblem|0|CC|xminus|FCLinear|yminus|Dirichlet|";
template<> const std::string BCFDView < HeatProblem<CC, P5>, 0, Patch::xminus | BC::FineCoarseInterface | FC::FCLinear, Patch::yminus | BC::Neumann >::Name = "HeatProblem|0|CC|xminus|FCLinear|yminus|Neumann|";
template<> const std::string BCFDView < HeatProblem<CC, P5>, 0, Patch::xminus | BC::FineCoarseInterface | FC::FCLinear, Patch::yplus | BC::FineCoarseInterface | FC::FCLinear >::Name = "HeatProblem|0|CC|xminus|FCLinear|yplus|FCLinear|";
template<> const std::string BCFDView < HeatProblem<CC, P5>, 0, Patch::xminus | BC::FineCoarseInterface | FC::FCLinear, Patch::yplus | BC::Dirichlet >::Name = "HeatProblem|0|CC|xminus|FCLinear|yplus|Dirichlet|";
template<> const std::string BCFDView < HeatProblem<CC, P5>, 0, Patch::xminus | BC::FineCoarseInterface | FC::FCLinear, Patch::yplus | BC::Neumann >::Name = "HeatProblem|0|CC|xminus|FCLinear|yplus|Neumann|";
template<> const std::string BCFDView < HeatProblem<CC, P5>, 0, Patch::xminus | BC::Dirichlet, Patch::yminus | BC::FineCoarseInterface | FC::FCLinear >::Name = "HeatProblem|0|CC|xminus|Dirichlet|yminus|FCLinear|";
template<> const std::string BCFDView < HeatProblem<CC, P5>, 0, Patch::xminus | BC::Dirichlet, Patch::yplus | BC::FineCoarseInterface | FC::FCLinear >::Name = "HeatProblem|0|CC|xminus|Dirichlet|yplus|FCLinear|";
template<> const std::string BCFDView < HeatProblem<CC, P5>, 0, Patch::xminus | BC::Neumann, Patch::yminus | BC::FineCoarseInterface | FC::FCLinear >::Name = "HeatProblem|0|CC|xminus|Neumann|yminus|FCLinear|";
template<> const std::string BCFDView < HeatProblem<CC, P5>, 0, Patch::xminus | BC::Neumann, Patch::yplus | BC::FineCoarseInterface | FC::FCLinear >::Name = "HeatProblem|0|CC|xminus|Neumann|yplus|FCLinear|";
template<> const std::string BCFDView < HeatProblem<CC, P5>, 0, Patch::xplus | BC::FineCoarseInterface | FC::FCLinear >::Name = "HeatProblem|0|CC|xplus|FCLinear|";
template<> const std::string BCFDView < HeatProblem<CC, P5>, 0, Patch::xplus | BC::FineCoarseInterface | FC::FCLinear, Patch::yminus | BC::FineCoarseInterface | FC::FCLinear >::Name = "HeatProblem|0|CC|xplus|FCLinear|yminus|FCLinear|";
template<> const std::string BCFDView < HeatProblem<CC, P5>, 0, Patch::xplus | BC::FineCoarseInterface | FC::FCLinear, Patch::yminus | BC::Dirichlet >::Name = "HeatProblem|0|CC|xplus|FCLinear|yminus|Dirichlet|";
template<> const std::string BCFDView < HeatProblem<CC, P5>, 0, Patch::xplus | BC::FineCoarseInterface | FC::FCLinear, Patch::yminus | BC::Neumann >::Name = "HeatProblem|0|CC|xplus|FCLinear|yminus|Neumann|";
template<> const std::string BCFDView < HeatProblem<CC, P5>, 0, Patch::xplus | BC::FineCoarseInterface | FC::FCLinear, Patch::yplus | BC::FineCoarseInterface | FC::FCLinear >::Name = "HeatProblem|0|CC|xplus|FCLinear|yplus|FCLinear|";
template<> const std::string BCFDView < HeatProblem<CC, P5>, 0, Patch::xplus | BC::FineCoarseInterface | FC::FCLinear, Patch::yplus | BC::Dirichlet >::Name = "HeatProblem|0|CC|xplus|FCLinear|yplus|Dirichlet|";
template<> const std::string BCFDView < HeatProblem<CC, P5>, 0, Patch::xplus | BC::FineCoarseInterface | FC::FCLinear, Patch::yplus | BC::Neumann >::Name = "HeatProblem|0|CC|xplus|FCLinear|yplus|Neumann|";
template<> const std::string BCFDView < HeatProblem<CC, P5>, 0, Patch::xplus | BC::Dirichlet, Patch::yminus | BC::FineCoarseInterface | FC::FCLinear >::Name = "HeatProblem|0|CC|xplus|Dirichlet|yminus|FCLinear|";
template<> const std::string BCFDView < HeatProblem<CC, P5>, 0, Patch::xplus | BC::Dirichlet, Patch::yplus | BC::FineCoarseInterface | FC::FCLinear >::Name = "HeatProblem|0|CC|xplus|Dirichlet|yplus|FCLinear|";
template<> const std::string BCFDView < HeatProblem<CC, P5>, 0, Patch::xplus | BC::Neumann, Patch::yminus | BC::FineCoarseInterface | FC::FCLinear >::Name = "HeatProblem|0|CC|xplus|Neumann|yminus|FCLinear|";
template<> const std::string BCFDView < HeatProblem<CC, P5>, 0, Patch::xplus | BC::Neumann, Patch::yplus | BC::FineCoarseInterface | FC::FCLinear >::Name = "HeatProblem|0|CC|xplus|Neumann|yplus|FCLinear|";
template<> const std::string BCFDView < HeatProblem<CC, P5>, 0, Patch::yminus | BC::FineCoarseInterface | FC::FCLinear >::Name = "HeatProblem|0|CC|yminus|FCLinear|";
template<> const std::string BCFDView < HeatProblem<CC, P5>, 0, Patch::yplus | BC::FineCoarseInterface | FC::FCLinear >::Name = "HeatProblem|0|CC|yplus|FCLinear|";

template class BCFDView < HeatProblem<CC, P5>, 0, Patch::xminus | BC::FineCoarseInterface | FC::FCLinear >;
template class BCFDView < HeatProblem<CC, P5>, 0, Patch::xminus | BC::FineCoarseInterface | FC::FCLinear, Patch::yminus | BC::FineCoarseInterface | FC::FCLinear >;
template class BCFDView < HeatProblem<CC, P5>, 0, Patch::xminus | BC::FineCoarseInterface | FC::FCLinear, Patch::yminus | BC::Dirichlet >;
template class BCFDView < HeatProblem<CC, P5>, 0, Patch::xminus | BC::FineCoarseInterface | FC::FCLinear, Patch::yminus | BC::Neumann >;
template class BCFDView < HeatProblem<CC, P5>, 0, Patch::xminus | BC::FineCoarseInterface | FC::FCLinear, Patch::yplus | BC::FineCoarseInterface | FC::FCLinear >;
template class BCFDView < HeatProblem<CC, P5>, 0, Patch::xminus | BC::FineCoarseInterface | FC::FCLinear, Patch::yplus | BC::Dirichlet >;
template class BCFDView < HeatProblem<CC, P5>, 0, Patch::xminus | BC::FineCoarseInterface | FC::FCLinear, Patch::yplus | BC::Neumann >;
template class BCFDView < HeatProblem<CC, P5>, 0, Patch::xminus | BC::Dirichlet, Patch::yminus | BC::FineCoarseInterface | FC::FCLinear >;
template class BCFDView < HeatProblem<CC, P5>, 0, Patch::xminus | BC::Dirichlet, Patch::yplus | BC::FineCoarseInterface | FC::FCLinear >;
template class BCFDView < HeatProblem<CC, P5>, 0, Patch::xminus | BC::Neumann, Patch::yminus | BC::FineCoarseInterface | FC::FCLinear >;
template class BCFDView < HeatProblem<CC, P5>, 0, Patch::xminus | BC::Neumann, Patch::yplus | BC::FineCoarseInterface | FC::FCLinear >;
template class BCFDView < HeatProblem<CC, P5>, 0, Patch::xplus | BC::FineCoarseInterface | FC::FCLinear >;
template class BCFDView < HeatProblem<CC, P5>, 0, Patch::xplus | BC::FineCoarseInterface | FC::FCLinear, Patch::yminus | BC::FineCoarseInterface | FC::FCLinear >;
template class BCFDView < HeatProblem<CC, P5>, 0, Patch::xplus | BC::FineCoarseInterface | FC::FCLinear, Patch::yminus | BC::Dirichlet >;
template class BCFDView < HeatProblem<CC, P5>, 0, Patch::xplus | BC::FineCoarseInterface | FC::FCLinear, Patch::yminus | BC::Neumann >;
template class BCFDView < HeatProblem<CC, P5>, 0, Patch::xplus | BC::FineCoarseInterface | FC::FCLinear, Patch::yplus | BC::FineCoarseInterface | FC::FCLinear >;
template class BCFDView < HeatProblem<CC, P5>, 0, Patch::xplus | BC::FineCoarseInterface | FC::FCLinear, Patch::yplus | BC::Dirichlet >;
template class BCFDView < HeatProblem<CC, P5>, 0, Patch::xplus | BC::FineCoarseInterface | FC::FCLinear, Patch::yplus | BC::Neumann >;
template class BCFDView < HeatProblem<CC, P5>, 0, Patch::xplus | BC::Dirichlet, Patch::yminus | BC::FineCoarseInterface | FC::FCLinear >;
template class BCFDView < HeatProblem<CC, P5>, 0, Patch::xplus | BC::Dirichlet, Patch::yplus | BC::FineCoarseInterface | FC::FCLinear >;
template class BCFDView < HeatProblem<CC, P5>, 0, Patch::xplus | BC::Neumann, Patch::yminus | BC::FineCoarseInterface | FC::FCLinear >;
template class BCFDView < HeatProblem<CC, P5>, 0, Patch::xplus | BC::Neumann, Patch::yplus | BC::FineCoarseInterface | FC::FCLinear >;
template class BCFDView < HeatProblem<CC, P5>, 0, Patch::yminus | BC::FineCoarseInterface | FC::FCLinear >;
template class BCFDView < HeatProblem<CC, P5>, 0, Patch::yplus | BC::FineCoarseInterface | FC::FCLinear >;

} // namespace Uintah
} // namespace PhaseField
