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
#include <CCA/Components/PhaseField/BoundaryConditions/BCFDViewFactory.h>

namespace Uintah {
namespace PhaseField {

template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Dirichlet >::Name = "HeatProblem|0|CC|xminus|Dirichlet|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Dirichlet, Patch::yminus | BC::Dirichlet >::Name = "HeatProblem|0|CC|xminus|Dirichlet|yminus|Dirichlet|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Dirichlet, Patch::yminus | BC::Dirichlet, Patch::zminus | BC::Dirichlet >::Name = "HeatProblem|0|CC|xminus|Dirichlet|yminus|Dirichlet|zminus|Dirichlet|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Dirichlet, Patch::yminus | BC::Dirichlet, Patch::zminus | BC::Neumann >::Name = "HeatProblem|0|CC|xminus|Dirichlet|yminus|Dirichlet|zminus|Neumann|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Dirichlet, Patch::yminus | BC::Dirichlet, Patch::zplus | BC::Dirichlet >::Name = "HeatProblem|0|CC|xminus|Dirichlet|yminus|Dirichlet|zplus|Dirichlet|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Dirichlet, Patch::yminus | BC::Dirichlet, Patch::zplus | BC::Neumann >::Name = "HeatProblem|0|CC|xminus|Dirichlet|yminus|Dirichlet|zplus|Neumann|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Dirichlet, Patch::yminus | BC::Neumann >::Name = "HeatProblem|0|CC|xminus|Dirichlet|yminus|Neumann|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Dirichlet, Patch::yminus | BC::Neumann, Patch::zminus | BC::Dirichlet >::Name = "HeatProblem|0|CC|xminus|Dirichlet|yminus|Neumann|zminus|Dirichlet|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Dirichlet, Patch::yminus | BC::Neumann, Patch::zminus | BC::Neumann >::Name = "HeatProblem|0|CC|xminus|Dirichlet|yminus|Neumann|zminus|Neumann|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Dirichlet, Patch::yminus | BC::Neumann, Patch::zplus | BC::Dirichlet >::Name = "HeatProblem|0|CC|xminus|Dirichlet|yminus|Neumann|zplus|Dirichlet|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Dirichlet, Patch::yminus | BC::Neumann, Patch::zplus | BC::Neumann >::Name = "HeatProblem|0|CC|xminus|Dirichlet|yminus|Neumann|zplus|Neumann|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Dirichlet, Patch::yplus | BC::Dirichlet >::Name = "HeatProblem|0|CC|xminus|Dirichlet|yplus|Dirichlet|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Dirichlet, Patch::yplus | BC::Dirichlet, Patch::zminus | BC::Dirichlet >::Name = "HeatProblem|0|CC|xminus|Dirichlet|yplus|Dirichlet|zminus|Dirichlet|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Dirichlet, Patch::yplus | BC::Dirichlet, Patch::zminus | BC::Neumann >::Name = "HeatProblem|0|CC|xminus|Dirichlet|yplus|Dirichlet|zminus|Neumann|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Dirichlet, Patch::yplus | BC::Dirichlet, Patch::zplus | BC::Dirichlet >::Name = "HeatProblem|0|CC|xminus|Dirichlet|yplus|Dirichlet|zplus|Dirichlet|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Dirichlet, Patch::yplus | BC::Dirichlet, Patch::zplus | BC::Neumann >::Name = "HeatProblem|0|CC|xminus|Dirichlet|yplus|Dirichlet|zplus|Neumann|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Dirichlet, Patch::yplus | BC::Neumann >::Name = "HeatProblem|0|CC|xminus|Dirichlet|yplus|Neumann|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Dirichlet, Patch::yplus | BC::Neumann, Patch::zminus | BC::Dirichlet >::Name = "HeatProblem|0|CC|xminus|Dirichlet|yplus|Neumann|zminus|Dirichlet|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Dirichlet, Patch::yplus | BC::Neumann, Patch::zminus | BC::Neumann >::Name = "HeatProblem|0|CC|xminus|Dirichlet|yplus|Neumann|zminus|Neumann|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Dirichlet, Patch::yplus | BC::Neumann, Patch::zplus | BC::Dirichlet >::Name = "HeatProblem|0|CC|xminus|Dirichlet|yplus|Neumann|zplus|Dirichlet|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Dirichlet, Patch::yplus | BC::Neumann, Patch::zplus | BC::Neumann >::Name = "HeatProblem|0|CC|xminus|Dirichlet|yplus|Neumann|zplus|Neumann|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Dirichlet, Patch::zminus | BC::Dirichlet >::Name = "HeatProblem|0|CC|xminus|Dirichlet|zminus|Dirichlet|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Dirichlet, Patch::zminus | BC::Neumann >::Name = "HeatProblem|0|CC|xminus|Dirichlet|zminus|Neumann|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Dirichlet, Patch::zplus | BC::Dirichlet >::Name = "HeatProblem|0|CC|xminus|Dirichlet|zplus|Dirichlet|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Dirichlet, Patch::zplus | BC::Neumann >::Name = "HeatProblem|0|CC|xminus|Dirichlet|zplus|Neumann|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Neumann >::Name = "HeatProblem|0|CC|xminus|Neumann|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Neumann, Patch::yminus | BC::Dirichlet >::Name = "HeatProblem|0|CC|xminus|Neumann|yminus|Dirichlet|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Neumann, Patch::yminus | BC::Dirichlet, Patch::zminus | BC::Dirichlet >::Name = "HeatProblem|0|CC|xminus|Neumann|yminus|Dirichlet|zminus|Dirichlet|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Neumann, Patch::yminus | BC::Dirichlet, Patch::zminus | BC::Neumann >::Name = "HeatProblem|0|CC|xminus|Neumann|yminus|Dirichlet|zminus|Neumann|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Neumann, Patch::yminus | BC::Dirichlet, Patch::zplus | BC::Dirichlet >::Name = "HeatProblem|0|CC|xminus|Neumann|yminus|Dirichlet|zplus|Dirichlet|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Neumann, Patch::yminus | BC::Dirichlet, Patch::zplus | BC::Neumann >::Name = "HeatProblem|0|CC|xminus|Neumann|yminus|Dirichlet|zplus|Neumann|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Neumann, Patch::yminus | BC::Neumann >::Name = "HeatProblem|0|CC|xminus|Neumann|yminus|Neumann|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Neumann, Patch::yminus | BC::Neumann, Patch::zminus | BC::Dirichlet >::Name = "HeatProblem|0|CC|xminus|Neumann|yminus|Neumann|zminus|Dirichlet|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Neumann, Patch::yminus | BC::Neumann, Patch::zminus | BC::Neumann >::Name = "HeatProblem|0|CC|xminus|Neumann|yminus|Neumann|zminus|Neumann|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Neumann, Patch::yminus | BC::Neumann, Patch::zplus | BC::Dirichlet >::Name = "HeatProblem|0|CC|xminus|Neumann|yminus|Neumann|zplus|Dirichlet|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Neumann, Patch::yminus | BC::Neumann, Patch::zplus | BC::Neumann >::Name = "HeatProblem|0|CC|xminus|Neumann|yminus|Neumann|zplus|Neumann|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Neumann, Patch::yplus | BC::Dirichlet >::Name = "HeatProblem|0|CC|xminus|Neumann|yplus|Dirichlet|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Neumann, Patch::yplus | BC::Dirichlet, Patch::zminus | BC::Dirichlet >::Name = "HeatProblem|0|CC|xminus|Neumann|yplus|Dirichlet|zminus|Dirichlet|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Neumann, Patch::yplus | BC::Dirichlet, Patch::zminus | BC::Neumann >::Name = "HeatProblem|0|CC|xminus|Neumann|yplus|Dirichlet|zminus|Neumann|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Neumann, Patch::yplus | BC::Dirichlet, Patch::zplus | BC::Dirichlet >::Name = "HeatProblem|0|CC|xminus|Neumann|yplus|Dirichlet|zplus|Dirichlet|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Neumann, Patch::yplus | BC::Dirichlet, Patch::zplus | BC::Neumann >::Name = "HeatProblem|0|CC|xminus|Neumann|yplus|Dirichlet|zplus|Neumann|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Neumann, Patch::yplus | BC::Neumann >::Name = "HeatProblem|0|CC|xminus|Neumann|yplus|Neumann|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Neumann, Patch::yplus | BC::Neumann, Patch::zminus | BC::Dirichlet >::Name = "HeatProblem|0|CC|xminus|Neumann|yplus|Neumann|zminus|Dirichlet|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Neumann, Patch::yplus | BC::Neumann, Patch::zminus | BC::Neumann >::Name = "HeatProblem|0|CC|xminus|Neumann|yplus|Neumann|zminus|Neumann|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Neumann, Patch::yplus | BC::Neumann, Patch::zplus | BC::Dirichlet >::Name = "HeatProblem|0|CC|xminus|Neumann|yplus|Neumann|zplus|Dirichlet|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Neumann, Patch::yplus | BC::Neumann, Patch::zplus | BC::Neumann >::Name = "HeatProblem|0|CC|xminus|Neumann|yplus|Neumann|zplus|Neumann|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Neumann, Patch::zminus | BC::Dirichlet >::Name = "HeatProblem|0|CC|xminus|Neumann|zminus|Dirichlet|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Neumann, Patch::zminus | BC::Neumann >::Name = "HeatProblem|0|CC|xminus|Neumann|zminus|Neumann|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Neumann, Patch::zplus | BC::Dirichlet >::Name = "HeatProblem|0|CC|xminus|Neumann|zplus|Dirichlet|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Neumann, Patch::zplus | BC::Neumann >::Name = "HeatProblem|0|CC|xminus|Neumann|zplus|Neumann|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Dirichlet >::Name = "HeatProblem|0|CC|xplus|Dirichlet|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Dirichlet, Patch::yminus | BC::Dirichlet >::Name = "HeatProblem|0|CC|xplus|Dirichlet|yminus|Dirichlet|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Dirichlet, Patch::yminus | BC::Dirichlet, Patch::zminus | BC::Dirichlet >::Name = "HeatProblem|0|CC|xplus|Dirichlet|yminus|Dirichlet|zminus|Dirichlet|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Dirichlet, Patch::yminus | BC::Dirichlet, Patch::zminus | BC::Neumann >::Name = "HeatProblem|0|CC|xplus|Dirichlet|yminus|Dirichlet|zminus|Neumann|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Dirichlet, Patch::yminus | BC::Dirichlet, Patch::zplus | BC::Dirichlet >::Name = "HeatProblem|0|CC|xplus|Dirichlet|yminus|Dirichlet|zplus|Dirichlet|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Dirichlet, Patch::yminus | BC::Dirichlet, Patch::zplus | BC::Neumann >::Name = "HeatProblem|0|CC|xplus|Dirichlet|yminus|Dirichlet|zplus|Neumann|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Dirichlet, Patch::yminus | BC::Neumann >::Name = "HeatProblem|0|CC|xplus|Dirichlet|yminus|Neumann|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Dirichlet, Patch::yminus | BC::Neumann, Patch::zminus | BC::Dirichlet >::Name = "HeatProblem|0|CC|xplus|Dirichlet|yminus|Neumann|zminus|Dirichlet|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Dirichlet, Patch::yminus | BC::Neumann, Patch::zminus | BC::Neumann >::Name = "HeatProblem|0|CC|xplus|Dirichlet|yminus|Neumann|zminus|Neumann|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Dirichlet, Patch::yminus | BC::Neumann, Patch::zplus | BC::Dirichlet >::Name = "HeatProblem|0|CC|xplus|Dirichlet|yminus|Neumann|zplus|Dirichlet|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Dirichlet, Patch::yminus | BC::Neumann, Patch::zplus | BC::Neumann >::Name = "HeatProblem|0|CC|xplus|Dirichlet|yminus|Neumann|zplus|Neumann|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Dirichlet, Patch::yplus | BC::Dirichlet >::Name = "HeatProblem|0|CC|xplus|Dirichlet|yplus|Dirichlet|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Dirichlet, Patch::yplus | BC::Dirichlet, Patch::zminus | BC::Dirichlet >::Name = "HeatProblem|0|CC|xplus|Dirichlet|yplus|Dirichlet|zminus|Dirichlet|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Dirichlet, Patch::yplus | BC::Dirichlet, Patch::zminus | BC::Neumann >::Name = "HeatProblem|0|CC|xplus|Dirichlet|yplus|Dirichlet|zminus|Neumann|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Dirichlet, Patch::yplus | BC::Dirichlet, Patch::zplus | BC::Dirichlet >::Name = "HeatProblem|0|CC|xplus|Dirichlet|yplus|Dirichlet|zplus|Dirichlet|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Dirichlet, Patch::yplus | BC::Dirichlet, Patch::zplus | BC::Neumann >::Name = "HeatProblem|0|CC|xplus|Dirichlet|yplus|Dirichlet|zplus|Neumann|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Dirichlet, Patch::yplus | BC::Neumann >::Name = "HeatProblem|0|CC|xplus|Dirichlet|yplus|Neumann|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Dirichlet, Patch::yplus | BC::Neumann, Patch::zminus | BC::Dirichlet >::Name = "HeatProblem|0|CC|xplus|Dirichlet|yplus|Neumann|zminus|Dirichlet|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Dirichlet, Patch::yplus | BC::Neumann, Patch::zminus | BC::Neumann >::Name = "HeatProblem|0|CC|xplus|Dirichlet|yplus|Neumann|zminus|Neumann|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Dirichlet, Patch::yplus | BC::Neumann, Patch::zplus | BC::Dirichlet >::Name = "HeatProblem|0|CC|xplus|Dirichlet|yplus|Neumann|zplus|Dirichlet|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Dirichlet, Patch::yplus | BC::Neumann, Patch::zplus | BC::Neumann >::Name = "HeatProblem|0|CC|xplus|Dirichlet|yplus|Neumann|zplus|Neumann|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Dirichlet, Patch::zminus | BC::Dirichlet >::Name = "HeatProblem|0|CC|xplus|Dirichlet|zminus|Dirichlet|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Dirichlet, Patch::zminus | BC::Neumann >::Name = "HeatProblem|0|CC|xplus|Dirichlet|zminus|Neumann|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Dirichlet, Patch::zplus | BC::Dirichlet >::Name = "HeatProblem|0|CC|xplus|Dirichlet|zplus|Dirichlet|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Dirichlet, Patch::zplus | BC::Neumann >::Name = "HeatProblem|0|CC|xplus|Dirichlet|zplus|Neumann|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Neumann >::Name = "HeatProblem|0|CC|xplus|Neumann|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Neumann, Patch::yminus | BC::Dirichlet >::Name = "HeatProblem|0|CC|xplus|Neumann|yminus|Dirichlet|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Neumann, Patch::yminus | BC::Dirichlet, Patch::zminus | BC::Dirichlet >::Name = "HeatProblem|0|CC|xplus|Neumann|yminus|Dirichlet|zminus|Dirichlet|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Neumann, Patch::yminus | BC::Dirichlet, Patch::zminus | BC::Neumann >::Name = "HeatProblem|0|CC|xplus|Neumann|yminus|Dirichlet|zminus|Neumann|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Neumann, Patch::yminus | BC::Dirichlet, Patch::zplus | BC::Dirichlet >::Name = "HeatProblem|0|CC|xplus|Neumann|yminus|Dirichlet|zplus|Dirichlet|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Neumann, Patch::yminus | BC::Dirichlet, Patch::zplus | BC::Neumann >::Name = "HeatProblem|0|CC|xplus|Neumann|yminus|Dirichlet|zplus|Neumann|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Neumann, Patch::yminus | BC::Neumann >::Name = "HeatProblem|0|CC|xplus|Neumann|yminus|Neumann|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Neumann, Patch::yminus | BC::Neumann, Patch::zminus | BC::Dirichlet >::Name = "HeatProblem|0|CC|xplus|Neumann|yminus|Neumann|zminus|Dirichlet|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Neumann, Patch::yminus | BC::Neumann, Patch::zminus | BC::Neumann >::Name = "HeatProblem|0|CC|xplus|Neumann|yminus|Neumann|zminus|Neumann|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Neumann, Patch::yminus | BC::Neumann, Patch::zplus | BC::Dirichlet >::Name = "HeatProblem|0|CC|xplus|Neumann|yminus|Neumann|zplus|Dirichlet|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Neumann, Patch::yminus | BC::Neumann, Patch::zplus | BC::Neumann >::Name = "HeatProblem|0|CC|xplus|Neumann|yminus|Neumann|zplus|Neumann|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Neumann, Patch::yplus | BC::Dirichlet >::Name = "HeatProblem|0|CC|xplus|Neumann|yplus|Dirichlet|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Neumann, Patch::yplus | BC::Dirichlet, Patch::zminus | BC::Dirichlet >::Name = "HeatProblem|0|CC|xplus|Neumann|yplus|Dirichlet|zminus|Dirichlet|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Neumann, Patch::yplus | BC::Dirichlet, Patch::zminus | BC::Neumann >::Name = "HeatProblem|0|CC|xplus|Neumann|yplus|Dirichlet|zminus|Neumann|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Neumann, Patch::yplus | BC::Dirichlet, Patch::zplus | BC::Dirichlet >::Name = "HeatProblem|0|CC|xplus|Neumann|yplus|Dirichlet|zplus|Dirichlet|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Neumann, Patch::yplus | BC::Dirichlet, Patch::zplus | BC::Neumann >::Name = "HeatProblem|0|CC|xplus|Neumann|yplus|Dirichlet|zplus|Neumann|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Neumann, Patch::yplus | BC::Neumann >::Name = "HeatProblem|0|CC|xplus|Neumann|yplus|Neumann|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Neumann, Patch::yplus | BC::Neumann, Patch::zminus | BC::Dirichlet >::Name = "HeatProblem|0|CC|xplus|Neumann|yplus|Neumann|zminus|Dirichlet|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Neumann, Patch::yplus | BC::Neumann, Patch::zminus | BC::Neumann >::Name = "HeatProblem|0|CC|xplus|Neumann|yplus|Neumann|zminus|Neumann|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Neumann, Patch::yplus | BC::Neumann, Patch::zplus | BC::Dirichlet >::Name = "HeatProblem|0|CC|xplus|Neumann|yplus|Neumann|zplus|Dirichlet|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Neumann, Patch::yplus | BC::Neumann, Patch::zplus | BC::Neumann >::Name = "HeatProblem|0|CC|xplus|Neumann|yplus|Neumann|zplus|Neumann|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Neumann, Patch::zminus | BC::Dirichlet >::Name = "HeatProblem|0|CC|xplus|Neumann|zminus|Dirichlet|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Neumann, Patch::zminus | BC::Neumann >::Name = "HeatProblem|0|CC|xplus|Neumann|zminus|Neumann|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Neumann, Patch::zplus | BC::Dirichlet >::Name = "HeatProblem|0|CC|xplus|Neumann|zplus|Dirichlet|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Neumann, Patch::zplus | BC::Neumann >::Name = "HeatProblem|0|CC|xplus|Neumann|zplus|Neumann|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::yminus | BC::Dirichlet >::Name = "HeatProblem|0|CC|yminus|Dirichlet|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::yminus | BC::Dirichlet, Patch::zminus | BC::Dirichlet >::Name = "HeatProblem|0|CC|yminus|Dirichlet|zminus|Dirichlet|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::yminus | BC::Dirichlet, Patch::zminus | BC::Neumann >::Name = "HeatProblem|0|CC|yminus|Dirichlet|zminus|Neumann|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::yminus | BC::Dirichlet, Patch::zplus | BC::Dirichlet >::Name = "HeatProblem|0|CC|yminus|Dirichlet|zplus|Dirichlet|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::yminus | BC::Dirichlet, Patch::zplus | BC::Neumann >::Name = "HeatProblem|0|CC|yminus|Dirichlet|zplus|Neumann|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::yminus | BC::Neumann >::Name = "HeatProblem|0|CC|yminus|Neumann|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::yminus | BC::Neumann, Patch::zminus | BC::Dirichlet >::Name = "HeatProblem|0|CC|yminus|Neumann|zminus|Dirichlet|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::yminus | BC::Neumann, Patch::zminus | BC::Neumann >::Name = "HeatProblem|0|CC|yminus|Neumann|zminus|Neumann|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::yminus | BC::Neumann, Patch::zplus | BC::Dirichlet >::Name = "HeatProblem|0|CC|yminus|Neumann|zplus|Dirichlet|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::yminus | BC::Neumann, Patch::zplus | BC::Neumann >::Name = "HeatProblem|0|CC|yminus|Neumann|zplus|Neumann|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::yplus | BC::Dirichlet >::Name = "HeatProblem|0|CC|yplus|Dirichlet|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::yplus | BC::Dirichlet, Patch::zminus | BC::Dirichlet >::Name = "HeatProblem|0|CC|yplus|Dirichlet|zminus|Dirichlet|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::yplus | BC::Dirichlet, Patch::zminus | BC::Neumann >::Name = "HeatProblem|0|CC|yplus|Dirichlet|zminus|Neumann|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::yplus | BC::Dirichlet, Patch::zplus | BC::Dirichlet >::Name = "HeatProblem|0|CC|yplus|Dirichlet|zplus|Dirichlet|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::yplus | BC::Dirichlet, Patch::zplus | BC::Neumann >::Name = "HeatProblem|0|CC|yplus|Dirichlet|zplus|Neumann|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::yplus | BC::Neumann >::Name = "HeatProblem|0|CC|yplus|Neumann|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::yplus | BC::Neumann, Patch::zminus | BC::Dirichlet >::Name = "HeatProblem|0|CC|yplus|Neumann|zminus|Dirichlet|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::yplus | BC::Neumann, Patch::zminus | BC::Neumann >::Name = "HeatProblem|0|CC|yplus|Neumann|zminus|Neumann|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::yplus | BC::Neumann, Patch::zplus | BC::Dirichlet >::Name = "HeatProblem|0|CC|yplus|Neumann|zplus|Dirichlet|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::yplus | BC::Neumann, Patch::zplus | BC::Neumann >::Name = "HeatProblem|0|CC|yplus|Neumann|zplus|Neumann|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::zminus | BC::Dirichlet >::Name = "HeatProblem|0|CC|zminus|Dirichlet|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::zminus | BC::Neumann >::Name = "HeatProblem|0|CC|zminus|Neumann|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::zplus | BC::Dirichlet >::Name = "HeatProblem|0|CC|zplus|Dirichlet|";
template<> const std::string BCFDView < HeatProblem<CC, P7>, 0, Patch::zplus | BC::Neumann >::Name = "HeatProblem|0|CC|zplus|Neumann|";

template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Dirichlet >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Dirichlet, Patch::yminus | BC::Dirichlet >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Dirichlet, Patch::yminus | BC::Dirichlet, Patch::zminus | BC::Dirichlet >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Dirichlet, Patch::yminus | BC::Dirichlet, Patch::zminus | BC::Neumann >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Dirichlet, Patch::yminus | BC::Dirichlet, Patch::zplus | BC::Dirichlet >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Dirichlet, Patch::yminus | BC::Dirichlet, Patch::zplus | BC::Neumann >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Dirichlet, Patch::yminus | BC::Neumann >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Dirichlet, Patch::yminus | BC::Neumann, Patch::zminus | BC::Dirichlet >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Dirichlet, Patch::yminus | BC::Neumann, Patch::zminus | BC::Neumann >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Dirichlet, Patch::yminus | BC::Neumann, Patch::zplus | BC::Dirichlet >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Dirichlet, Patch::yminus | BC::Neumann, Patch::zplus | BC::Neumann >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Dirichlet, Patch::yplus | BC::Dirichlet >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Dirichlet, Patch::yplus | BC::Dirichlet, Patch::zminus | BC::Dirichlet >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Dirichlet, Patch::yplus | BC::Dirichlet, Patch::zminus | BC::Neumann >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Dirichlet, Patch::yplus | BC::Dirichlet, Patch::zplus | BC::Dirichlet >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Dirichlet, Patch::yplus | BC::Dirichlet, Patch::zplus | BC::Neumann >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Dirichlet, Patch::yplus | BC::Neumann >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Dirichlet, Patch::yplus | BC::Neumann, Patch::zminus | BC::Dirichlet >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Dirichlet, Patch::yplus | BC::Neumann, Patch::zminus | BC::Neumann >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Dirichlet, Patch::yplus | BC::Neumann, Patch::zplus | BC::Dirichlet >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Dirichlet, Patch::yplus | BC::Neumann, Patch::zplus | BC::Neumann >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Dirichlet, Patch::zminus | BC::Dirichlet >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Dirichlet, Patch::zminus | BC::Neumann >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Dirichlet, Patch::zplus | BC::Dirichlet >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Dirichlet, Patch::zplus | BC::Neumann >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Neumann >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Neumann, Patch::yminus | BC::Dirichlet >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Neumann, Patch::yminus | BC::Dirichlet, Patch::zminus | BC::Dirichlet >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Neumann, Patch::yminus | BC::Dirichlet, Patch::zminus | BC::Neumann >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Neumann, Patch::yminus | BC::Dirichlet, Patch::zplus | BC::Dirichlet >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Neumann, Patch::yminus | BC::Dirichlet, Patch::zplus | BC::Neumann >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Neumann, Patch::yminus | BC::Neumann >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Neumann, Patch::yminus | BC::Neumann, Patch::zminus | BC::Dirichlet >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Neumann, Patch::yminus | BC::Neumann, Patch::zminus | BC::Neumann >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Neumann, Patch::yminus | BC::Neumann, Patch::zplus | BC::Dirichlet >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Neumann, Patch::yminus | BC::Neumann, Patch::zplus | BC::Neumann >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Neumann, Patch::yplus | BC::Dirichlet >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Neumann, Patch::yplus | BC::Dirichlet, Patch::zminus | BC::Dirichlet >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Neumann, Patch::yplus | BC::Dirichlet, Patch::zminus | BC::Neumann >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Neumann, Patch::yplus | BC::Dirichlet, Patch::zplus | BC::Dirichlet >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Neumann, Patch::yplus | BC::Dirichlet, Patch::zplus | BC::Neumann >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Neumann, Patch::yplus | BC::Neumann >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Neumann, Patch::yplus | BC::Neumann, Patch::zminus | BC::Dirichlet >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Neumann, Patch::yplus | BC::Neumann, Patch::zminus | BC::Neumann >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Neumann, Patch::yplus | BC::Neumann, Patch::zplus | BC::Dirichlet >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Neumann, Patch::yplus | BC::Neumann, Patch::zplus | BC::Neumann >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Neumann, Patch::zminus | BC::Dirichlet >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Neumann, Patch::zminus | BC::Neumann >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Neumann, Patch::zplus | BC::Dirichlet >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xminus | BC::Neumann, Patch::zplus | BC::Neumann >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Dirichlet >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Dirichlet, Patch::yminus | BC::Dirichlet >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Dirichlet, Patch::yminus | BC::Dirichlet, Patch::zminus | BC::Dirichlet >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Dirichlet, Patch::yminus | BC::Dirichlet, Patch::zminus | BC::Neumann >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Dirichlet, Patch::yminus | BC::Dirichlet, Patch::zplus | BC::Dirichlet >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Dirichlet, Patch::yminus | BC::Dirichlet, Patch::zplus | BC::Neumann >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Dirichlet, Patch::yminus | BC::Neumann >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Dirichlet, Patch::yminus | BC::Neumann, Patch::zminus | BC::Dirichlet >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Dirichlet, Patch::yminus | BC::Neumann, Patch::zminus | BC::Neumann >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Dirichlet, Patch::yminus | BC::Neumann, Patch::zplus | BC::Dirichlet >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Dirichlet, Patch::yminus | BC::Neumann, Patch::zplus | BC::Neumann >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Dirichlet, Patch::yplus | BC::Dirichlet >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Dirichlet, Patch::yplus | BC::Dirichlet, Patch::zminus | BC::Dirichlet >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Dirichlet, Patch::yplus | BC::Dirichlet, Patch::zminus | BC::Neumann >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Dirichlet, Patch::yplus | BC::Dirichlet, Patch::zplus | BC::Dirichlet >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Dirichlet, Patch::yplus | BC::Dirichlet, Patch::zplus | BC::Neumann >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Dirichlet, Patch::yplus | BC::Neumann >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Dirichlet, Patch::yplus | BC::Neumann, Patch::zminus | BC::Dirichlet >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Dirichlet, Patch::yplus | BC::Neumann, Patch::zminus | BC::Neumann >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Dirichlet, Patch::yplus | BC::Neumann, Patch::zplus | BC::Dirichlet >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Dirichlet, Patch::yplus | BC::Neumann, Patch::zplus | BC::Neumann >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Dirichlet, Patch::zminus | BC::Dirichlet >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Dirichlet, Patch::zminus | BC::Neumann >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Dirichlet, Patch::zplus | BC::Dirichlet >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Dirichlet, Patch::zplus | BC::Neumann >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Neumann >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Neumann, Patch::yminus | BC::Dirichlet >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Neumann, Patch::yminus | BC::Dirichlet, Patch::zminus | BC::Dirichlet >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Neumann, Patch::yminus | BC::Dirichlet, Patch::zminus | BC::Neumann >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Neumann, Patch::yminus | BC::Dirichlet, Patch::zplus | BC::Dirichlet >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Neumann, Patch::yminus | BC::Dirichlet, Patch::zplus | BC::Neumann >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Neumann, Patch::yminus | BC::Neumann >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Neumann, Patch::yminus | BC::Neumann, Patch::zminus | BC::Dirichlet >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Neumann, Patch::yminus | BC::Neumann, Patch::zminus | BC::Neumann >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Neumann, Patch::yminus | BC::Neumann, Patch::zplus | BC::Dirichlet >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Neumann, Patch::yminus | BC::Neumann, Patch::zplus | BC::Neumann >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Neumann, Patch::yplus | BC::Dirichlet >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Neumann, Patch::yplus | BC::Dirichlet, Patch::zminus | BC::Dirichlet >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Neumann, Patch::yplus | BC::Dirichlet, Patch::zminus | BC::Neumann >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Neumann, Patch::yplus | BC::Dirichlet, Patch::zplus | BC::Dirichlet >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Neumann, Patch::yplus | BC::Dirichlet, Patch::zplus | BC::Neumann >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Neumann, Patch::yplus | BC::Neumann >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Neumann, Patch::yplus | BC::Neumann, Patch::zminus | BC::Dirichlet >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Neumann, Patch::yplus | BC::Neumann, Patch::zminus | BC::Neumann >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Neumann, Patch::yplus | BC::Neumann, Patch::zplus | BC::Dirichlet >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Neumann, Patch::yplus | BC::Neumann, Patch::zplus | BC::Neumann >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Neumann, Patch::zminus | BC::Dirichlet >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Neumann, Patch::zminus | BC::Neumann >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Neumann, Patch::zplus | BC::Dirichlet >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::xplus | BC::Neumann, Patch::zplus | BC::Neumann >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::yminus | BC::Dirichlet >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::yminus | BC::Dirichlet, Patch::zminus | BC::Dirichlet >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::yminus | BC::Dirichlet, Patch::zminus | BC::Neumann >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::yminus | BC::Dirichlet, Patch::zplus | BC::Dirichlet >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::yminus | BC::Dirichlet, Patch::zplus | BC::Neumann >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::yminus | BC::Neumann >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::yminus | BC::Neumann, Patch::zminus | BC::Dirichlet >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::yminus | BC::Neumann, Patch::zminus | BC::Neumann >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::yminus | BC::Neumann, Patch::zplus | BC::Dirichlet >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::yminus | BC::Neumann, Patch::zplus | BC::Neumann >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::yplus | BC::Dirichlet >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::yplus | BC::Dirichlet, Patch::zminus | BC::Dirichlet >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::yplus | BC::Dirichlet, Patch::zminus | BC::Neumann >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::yplus | BC::Dirichlet, Patch::zplus | BC::Dirichlet >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::yplus | BC::Dirichlet, Patch::zplus | BC::Neumann >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::yplus | BC::Neumann >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::yplus | BC::Neumann, Patch::zminus | BC::Dirichlet >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::yplus | BC::Neumann, Patch::zminus | BC::Neumann >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::yplus | BC::Neumann, Patch::zplus | BC::Dirichlet >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::yplus | BC::Neumann, Patch::zplus | BC::Neumann >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::zminus | BC::Dirichlet >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::zminus | BC::Neumann >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::zplus | BC::Dirichlet >;
template class BCFDView < HeatProblem<CC, P7>, 0, Patch::zplus | BC::Neumann >;

} // namespace Uintah
} // namespace PhaseField
