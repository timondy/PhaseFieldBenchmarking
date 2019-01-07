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

template<> const std::string BCFDView < HeatProblem<NC, P5>, 0, Patch::xminus | BC::FineCoarseInterface | FC::FCBilinear >::Name = "HeatProblem|0|NC|xminus|FCBilinear|";
template<> const std::string BCFDView < HeatProblem<NC, P5>, 0, Patch::xminus | BC::FineCoarseInterface | FC::FCBilinear, Patch::yminus | BC::FineCoarseInterface | FC::FCBilinear >::Name = "HeatProblem|0|NC|xminus|FCBilinear|yminus|FCBilinear|";
template<> const std::string BCFDView < HeatProblem<NC, P5>, 0, Patch::xminus | BC::FineCoarseInterface | FC::FCBilinear, Patch::yminus | BC::Dirichlet >::Name = "HeatProblem|0|NC|xminus|FCBilinear|yminus|Dirichlet|";
template<> const std::string BCFDView < HeatProblem<NC, P5>, 0, Patch::xminus | BC::FineCoarseInterface | FC::FCBilinear, Patch::yminus | BC::Neumann >::Name = "HeatProblem|0|NC|xminus|FCBilinear|yminus|Neumann|";
template<> const std::string BCFDView < HeatProblem<NC, P5>, 0, Patch::xminus | BC::FineCoarseInterface | FC::FCBilinear, Patch::yplus | BC::FineCoarseInterface | FC::FCBilinear >::Name = "HeatProblem|0|NC|xminus|FCBilinear|yplus|FCBilinear|";
template<> const std::string BCFDView < HeatProblem<NC, P5>, 0, Patch::xminus | BC::FineCoarseInterface | FC::FCBilinear, Patch::yplus | BC::Dirichlet >::Name = "HeatProblem|0|NC|xminus|FCBilinear|yplus|Dirichlet|";
template<> const std::string BCFDView < HeatProblem<NC, P5>, 0, Patch::xminus | BC::FineCoarseInterface | FC::FCBilinear, Patch::yplus | BC::Neumann >::Name = "HeatProblem|0|NC|xminus|FCBilinear|yplus|Neumann|";
template<> const std::string BCFDView < HeatProblem<NC, P5>, 0, Patch::xminus | BC::Dirichlet, Patch::yminus | BC::FineCoarseInterface | FC::FCBilinear >::Name = "HeatProblem|0|NC|xminus|Dirichlet|yminus|FCBilinear|";
template<> const std::string BCFDView < HeatProblem<NC, P5>, 0, Patch::xminus | BC::Dirichlet, Patch::yplus | BC::FineCoarseInterface | FC::FCBilinear >::Name = "HeatProblem|0|NC|xminus|Dirichlet|yplus|FCBilinear|";
template<> const std::string BCFDView < HeatProblem<NC, P5>, 0, Patch::xminus | BC::Neumann, Patch::yminus | BC::FineCoarseInterface | FC::FCBilinear >::Name = "HeatProblem|0|NC|xminus|Neumann|yminus|FCBilinear|";
template<> const std::string BCFDView < HeatProblem<NC, P5>, 0, Patch::xminus | BC::Neumann, Patch::yplus | BC::FineCoarseInterface | FC::FCBilinear >::Name = "HeatProblem|0|NC|xminus|Neumann|yplus|FCBilinear|";
template<> const std::string BCFDView < HeatProblem<NC, P5>, 0, Patch::xplus | BC::FineCoarseInterface | FC::FCBilinear >::Name = "HeatProblem|0|NC|xplus|FCBilinear|";
template<> const std::string BCFDView < HeatProblem<NC, P5>, 0, Patch::xplus | BC::FineCoarseInterface | FC::FCBilinear, Patch::yminus | BC::FineCoarseInterface | FC::FCBilinear >::Name = "HeatProblem|0|NC|xplus|FCBilinear|yminus|FCBilinear|";
template<> const std::string BCFDView < HeatProblem<NC, P5>, 0, Patch::xplus | BC::FineCoarseInterface | FC::FCBilinear, Patch::yminus | BC::Dirichlet >::Name = "HeatProblem|0|NC|xplus|FCBilinear|yminus|Dirichlet|";
template<> const std::string BCFDView < HeatProblem<NC, P5>, 0, Patch::xplus | BC::FineCoarseInterface | FC::FCBilinear, Patch::yminus | BC::Neumann >::Name = "HeatProblem|0|NC|xplus|FCBilinear|yminus|Neumann|";
template<> const std::string BCFDView < HeatProblem<NC, P5>, 0, Patch::xplus | BC::FineCoarseInterface | FC::FCBilinear, Patch::yplus | BC::FineCoarseInterface | FC::FCBilinear >::Name = "HeatProblem|0|NC|xplus|FCBilinear|yplus|FCBilinear|";
template<> const std::string BCFDView < HeatProblem<NC, P5>, 0, Patch::xplus | BC::FineCoarseInterface | FC::FCBilinear, Patch::yplus | BC::Dirichlet >::Name = "HeatProblem|0|NC|xplus|FCBilinear|yplus|Dirichlet|";
template<> const std::string BCFDView < HeatProblem<NC, P5>, 0, Patch::xplus | BC::FineCoarseInterface | FC::FCBilinear, Patch::yplus | BC::Neumann >::Name = "HeatProblem|0|NC|xplus|FCBilinear|yplus|Neumann|";
template<> const std::string BCFDView < HeatProblem<NC, P5>, 0, Patch::xplus | BC::Dirichlet, Patch::yminus | BC::FineCoarseInterface | FC::FCBilinear >::Name = "HeatProblem|0|NC|xplus|Dirichlet|yminus|FCBilinear|";
template<> const std::string BCFDView < HeatProblem<NC, P5>, 0, Patch::xplus | BC::Dirichlet, Patch::yplus | BC::FineCoarseInterface | FC::FCBilinear >::Name = "HeatProblem|0|NC|xplus|Dirichlet|yplus|FCBilinear|";
template<> const std::string BCFDView < HeatProblem<NC, P5>, 0, Patch::xplus | BC::Neumann, Patch::yminus | BC::FineCoarseInterface | FC::FCBilinear >::Name = "HeatProblem|0|NC|xplus|Neumann|yminus|FCBilinear|";
template<> const std::string BCFDView < HeatProblem<NC, P5>, 0, Patch::xplus | BC::Neumann, Patch::yplus | BC::FineCoarseInterface | FC::FCBilinear >::Name = "HeatProblem|0|NC|xplus|Neumann|yplus|FCBilinear|";
template<> const std::string BCFDView < HeatProblem<NC, P5>, 0, Patch::yminus | BC::FineCoarseInterface | FC::FCBilinear >::Name = "HeatProblem|0|NC|yminus|FCBilinear|";
template<> const std::string BCFDView < HeatProblem<NC, P5>, 0, Patch::yplus | BC::FineCoarseInterface | FC::FCBilinear >::Name = "HeatProblem|0|NC|yplus|FCBilinear|";

template class BCFDView < HeatProblem<NC, P5>, 0, Patch::xminus | BC::FineCoarseInterface | FC::FCBilinear >;
template class BCFDView < HeatProblem<NC, P5>, 0, Patch::xminus | BC::FineCoarseInterface | FC::FCBilinear, Patch::yminus | BC::FineCoarseInterface | FC::FCBilinear >;
template class BCFDView < HeatProblem<NC, P5>, 0, Patch::xminus | BC::FineCoarseInterface | FC::FCBilinear, Patch::yminus | BC::Dirichlet >;
template class BCFDView < HeatProblem<NC, P5>, 0, Patch::xminus | BC::FineCoarseInterface | FC::FCBilinear, Patch::yminus | BC::Neumann >;
template class BCFDView < HeatProblem<NC, P5>, 0, Patch::xminus | BC::FineCoarseInterface | FC::FCBilinear, Patch::yplus | BC::FineCoarseInterface | FC::FCBilinear >;
template class BCFDView < HeatProblem<NC, P5>, 0, Patch::xminus | BC::FineCoarseInterface | FC::FCBilinear, Patch::yplus | BC::Dirichlet >;
template class BCFDView < HeatProblem<NC, P5>, 0, Patch::xminus | BC::FineCoarseInterface | FC::FCBilinear, Patch::yplus | BC::Neumann >;
template class BCFDView < HeatProblem<NC, P5>, 0, Patch::xminus | BC::Dirichlet, Patch::yminus | BC::FineCoarseInterface | FC::FCBilinear >;
template class BCFDView < HeatProblem<NC, P5>, 0, Patch::xminus | BC::Dirichlet, Patch::yplus | BC::FineCoarseInterface | FC::FCBilinear >;
template class BCFDView < HeatProblem<NC, P5>, 0, Patch::xminus | BC::Neumann, Patch::yminus | BC::FineCoarseInterface | FC::FCBilinear >;
template class BCFDView < HeatProblem<NC, P5>, 0, Patch::xminus | BC::Neumann, Patch::yplus | BC::FineCoarseInterface | FC::FCBilinear >;
template class BCFDView < HeatProblem<NC, P5>, 0, Patch::xplus | BC::FineCoarseInterface | FC::FCBilinear >;
template class BCFDView < HeatProblem<NC, P5>, 0, Patch::xplus | BC::FineCoarseInterface | FC::FCBilinear, Patch::yminus | BC::FineCoarseInterface | FC::FCBilinear >;
template class BCFDView < HeatProblem<NC, P5>, 0, Patch::xplus | BC::FineCoarseInterface | FC::FCBilinear, Patch::yminus | BC::Dirichlet >;
template class BCFDView < HeatProblem<NC, P5>, 0, Patch::xplus | BC::FineCoarseInterface | FC::FCBilinear, Patch::yminus | BC::Neumann >;
template class BCFDView < HeatProblem<NC, P5>, 0, Patch::xplus | BC::FineCoarseInterface | FC::FCBilinear, Patch::yplus | BC::FineCoarseInterface | FC::FCBilinear >;
template class BCFDView < HeatProblem<NC, P5>, 0, Patch::xplus | BC::FineCoarseInterface | FC::FCBilinear, Patch::yplus | BC::Dirichlet >;
template class BCFDView < HeatProblem<NC, P5>, 0, Patch::xplus | BC::FineCoarseInterface | FC::FCBilinear, Patch::yplus | BC::Neumann >;
template class BCFDView < HeatProblem<NC, P5>, 0, Patch::xplus | BC::Dirichlet, Patch::yminus | BC::FineCoarseInterface | FC::FCBilinear >;
template class BCFDView < HeatProblem<NC, P5>, 0, Patch::xplus | BC::Dirichlet, Patch::yplus | BC::FineCoarseInterface | FC::FCBilinear >;
template class BCFDView < HeatProblem<NC, P5>, 0, Patch::xplus | BC::Neumann, Patch::yminus | BC::FineCoarseInterface | FC::FCBilinear >;
template class BCFDView < HeatProblem<NC, P5>, 0, Patch::xplus | BC::Neumann, Patch::yplus | BC::FineCoarseInterface | FC::FCBilinear >;
template class BCFDView < HeatProblem<NC, P5>, 0, Patch::yminus | BC::FineCoarseInterface | FC::FCBilinear >;
template class BCFDView < HeatProblem<NC, P5>, 0, Patch::yplus | BC::FineCoarseInterface | FC::FCBilinear >;

} // namespace Uintah
} // namespace PhaseField
