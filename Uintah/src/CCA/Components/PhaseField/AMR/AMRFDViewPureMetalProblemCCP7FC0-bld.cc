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

template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 0, Patch::xminus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|0|CC|xminus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 0, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::yminus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|0|CC|xminus|FC0|yminus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 0, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::yminus | BC::FineCoarseInterface | FC::FC0, Patch::zminus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|0|CC|xminus|FC0|yminus|FC0|zminus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 0, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::yminus | BC::FineCoarseInterface | FC::FC0, Patch::zplus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|0|CC|xminus|FC0|yminus|FC0|zplus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 0, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::yplus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|0|CC|xminus|FC0|yplus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 0, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::yplus | BC::FineCoarseInterface | FC::FC0, Patch::zminus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|0|CC|xminus|FC0|yplus|FC0|zminus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 0, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::yplus | BC::FineCoarseInterface | FC::FC0, Patch::zplus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|0|CC|xminus|FC0|yplus|FC0|zplus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 0, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::zminus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|0|CC|xminus|FC0|zminus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 0, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::zplus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|0|CC|xminus|FC0|zplus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 0, Patch::xplus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|0|CC|xplus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 0, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::yminus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|0|CC|xplus|FC0|yminus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 0, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::yminus | BC::FineCoarseInterface | FC::FC0, Patch::zminus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|0|CC|xplus|FC0|yminus|FC0|zminus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 0, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::yminus | BC::FineCoarseInterface | FC::FC0, Patch::zplus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|0|CC|xplus|FC0|yminus|FC0|zplus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 0, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::yplus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|0|CC|xplus|FC0|yplus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 0, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::yplus | BC::FineCoarseInterface | FC::FC0, Patch::zminus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|0|CC|xplus|FC0|yplus|FC0|zminus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 0, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::yplus | BC::FineCoarseInterface | FC::FC0, Patch::zplus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|0|CC|xplus|FC0|yplus|FC0|zplus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 0, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::zminus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|0|CC|xplus|FC0|zminus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 0, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::zplus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|0|CC|xplus|FC0|zplus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 0, Patch::yminus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|0|CC|yminus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 0, Patch::yminus | BC::FineCoarseInterface | FC::FC0, Patch::zminus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|0|CC|yminus|FC0|zminus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 0, Patch::yminus | BC::FineCoarseInterface | FC::FC0, Patch::zplus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|0|CC|yminus|FC0|zplus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 0, Patch::yplus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|0|CC|yplus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 0, Patch::yplus | BC::FineCoarseInterface | FC::FC0, Patch::zminus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|0|CC|yplus|FC0|zminus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 0, Patch::yplus | BC::FineCoarseInterface | FC::FC0, Patch::zplus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|0|CC|yplus|FC0|zplus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 0, Patch::zminus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|0|CC|zminus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 0, Patch::zplus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|0|CC|zplus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 1, Patch::xminus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|1|CC|xminus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 1, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::yminus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|1|CC|xminus|FC0|yminus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 1, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::yminus | BC::FineCoarseInterface | FC::FC0, Patch::zminus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|1|CC|xminus|FC0|yminus|FC0|zminus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 1, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::yminus | BC::FineCoarseInterface | FC::FC0, Patch::zplus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|1|CC|xminus|FC0|yminus|FC0|zplus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 1, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::yplus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|1|CC|xminus|FC0|yplus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 1, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::yplus | BC::FineCoarseInterface | FC::FC0, Patch::zminus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|1|CC|xminus|FC0|yplus|FC0|zminus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 1, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::yplus | BC::FineCoarseInterface | FC::FC0, Patch::zplus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|1|CC|xminus|FC0|yplus|FC0|zplus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 1, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::zminus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|1|CC|xminus|FC0|zminus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 1, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::zplus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|1|CC|xminus|FC0|zplus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 1, Patch::xplus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|1|CC|xplus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 1, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::yminus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|1|CC|xplus|FC0|yminus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 1, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::yminus | BC::FineCoarseInterface | FC::FC0, Patch::zminus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|1|CC|xplus|FC0|yminus|FC0|zminus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 1, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::yminus | BC::FineCoarseInterface | FC::FC0, Patch::zplus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|1|CC|xplus|FC0|yminus|FC0|zplus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 1, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::yplus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|1|CC|xplus|FC0|yplus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 1, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::yplus | BC::FineCoarseInterface | FC::FC0, Patch::zminus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|1|CC|xplus|FC0|yplus|FC0|zminus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 1, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::yplus | BC::FineCoarseInterface | FC::FC0, Patch::zplus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|1|CC|xplus|FC0|yplus|FC0|zplus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 1, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::zminus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|1|CC|xplus|FC0|zminus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 1, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::zplus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|1|CC|xplus|FC0|zplus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 1, Patch::yminus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|1|CC|yminus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 1, Patch::yminus | BC::FineCoarseInterface | FC::FC0, Patch::zminus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|1|CC|yminus|FC0|zminus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 1, Patch::yminus | BC::FineCoarseInterface | FC::FC0, Patch::zplus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|1|CC|yminus|FC0|zplus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 1, Patch::yplus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|1|CC|yplus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 1, Patch::yplus | BC::FineCoarseInterface | FC::FC0, Patch::zminus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|1|CC|yplus|FC0|zminus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 1, Patch::yplus | BC::FineCoarseInterface | FC::FC0, Patch::zplus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|1|CC|yplus|FC0|zplus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 1, Patch::zminus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|1|CC|zminus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 1, Patch::zplus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|1|CC|zplus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 2, Patch::xminus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|2|CC|xminus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 2, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::yminus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|2|CC|xminus|FC0|yminus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 2, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::yminus | BC::FineCoarseInterface | FC::FC0, Patch::zminus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|2|CC|xminus|FC0|yminus|FC0|zminus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 2, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::yminus | BC::FineCoarseInterface | FC::FC0, Patch::zplus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|2|CC|xminus|FC0|yminus|FC0|zplus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 2, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::yplus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|2|CC|xminus|FC0|yplus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 2, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::yplus | BC::FineCoarseInterface | FC::FC0, Patch::zminus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|2|CC|xminus|FC0|yplus|FC0|zminus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 2, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::yplus | BC::FineCoarseInterface | FC::FC0, Patch::zplus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|2|CC|xminus|FC0|yplus|FC0|zplus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 2, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::zminus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|2|CC|xminus|FC0|zminus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 2, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::zplus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|2|CC|xminus|FC0|zplus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 2, Patch::xplus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|2|CC|xplus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 2, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::yminus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|2|CC|xplus|FC0|yminus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 2, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::yminus | BC::FineCoarseInterface | FC::FC0, Patch::zminus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|2|CC|xplus|FC0|yminus|FC0|zminus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 2, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::yminus | BC::FineCoarseInterface | FC::FC0, Patch::zplus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|2|CC|xplus|FC0|yminus|FC0|zplus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 2, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::yplus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|2|CC|xplus|FC0|yplus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 2, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::yplus | BC::FineCoarseInterface | FC::FC0, Patch::zminus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|2|CC|xplus|FC0|yplus|FC0|zminus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 2, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::yplus | BC::FineCoarseInterface | FC::FC0, Patch::zplus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|2|CC|xplus|FC0|yplus|FC0|zplus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 2, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::zminus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|2|CC|xplus|FC0|zminus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 2, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::zplus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|2|CC|xplus|FC0|zplus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 2, Patch::yminus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|2|CC|yminus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 2, Patch::yminus | BC::FineCoarseInterface | FC::FC0, Patch::zminus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|2|CC|yminus|FC0|zminus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 2, Patch::yminus | BC::FineCoarseInterface | FC::FC0, Patch::zplus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|2|CC|yminus|FC0|zplus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 2, Patch::yplus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|2|CC|yplus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 2, Patch::yplus | BC::FineCoarseInterface | FC::FC0, Patch::zminus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|2|CC|yplus|FC0|zminus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 2, Patch::yplus | BC::FineCoarseInterface | FC::FC0, Patch::zplus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|2|CC|yplus|FC0|zplus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 2, Patch::zminus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|2|CC|zminus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 2, Patch::zplus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|2|CC|zplus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 3, Patch::xminus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|3|CC|xminus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 3, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::yminus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|3|CC|xminus|FC0|yminus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 3, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::yminus | BC::FineCoarseInterface | FC::FC0, Patch::zminus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|3|CC|xminus|FC0|yminus|FC0|zminus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 3, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::yminus | BC::FineCoarseInterface | FC::FC0, Patch::zplus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|3|CC|xminus|FC0|yminus|FC0|zplus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 3, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::yplus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|3|CC|xminus|FC0|yplus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 3, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::yplus | BC::FineCoarseInterface | FC::FC0, Patch::zminus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|3|CC|xminus|FC0|yplus|FC0|zminus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 3, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::yplus | BC::FineCoarseInterface | FC::FC0, Patch::zplus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|3|CC|xminus|FC0|yplus|FC0|zplus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 3, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::zminus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|3|CC|xminus|FC0|zminus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 3, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::zplus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|3|CC|xminus|FC0|zplus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 3, Patch::xplus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|3|CC|xplus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 3, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::yminus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|3|CC|xplus|FC0|yminus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 3, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::yminus | BC::FineCoarseInterface | FC::FC0, Patch::zminus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|3|CC|xplus|FC0|yminus|FC0|zminus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 3, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::yminus | BC::FineCoarseInterface | FC::FC0, Patch::zplus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|3|CC|xplus|FC0|yminus|FC0|zplus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 3, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::yplus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|3|CC|xplus|FC0|yplus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 3, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::yplus | BC::FineCoarseInterface | FC::FC0, Patch::zminus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|3|CC|xplus|FC0|yplus|FC0|zminus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 3, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::yplus | BC::FineCoarseInterface | FC::FC0, Patch::zplus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|3|CC|xplus|FC0|yplus|FC0|zplus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 3, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::zminus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|3|CC|xplus|FC0|zminus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 3, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::zplus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|3|CC|xplus|FC0|zplus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 3, Patch::yminus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|3|CC|yminus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 3, Patch::yminus | BC::FineCoarseInterface | FC::FC0, Patch::zminus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|3|CC|yminus|FC0|zminus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 3, Patch::yminus | BC::FineCoarseInterface | FC::FC0, Patch::zplus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|3|CC|yminus|FC0|zplus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 3, Patch::yplus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|3|CC|yplus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 3, Patch::yplus | BC::FineCoarseInterface | FC::FC0, Patch::zminus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|3|CC|yplus|FC0|zminus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 3, Patch::yplus | BC::FineCoarseInterface | FC::FC0, Patch::zplus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|3|CC|yplus|FC0|zplus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 3, Patch::zminus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|3|CC|zminus|FC0|";
template<> const std::string BCFDView < PureMetalProblem<CC, P7>, 3, Patch::zplus | BC::FineCoarseInterface | FC::FC0 >::Name = "PureMetalProblem|3|CC|zplus|FC0|";

template class BCFDView < PureMetalProblem<CC, P7>, 0, Patch::xminus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 0, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::yminus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 0, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::yminus | BC::FineCoarseInterface | FC::FC0, Patch::zminus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 0, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::yminus | BC::FineCoarseInterface | FC::FC0, Patch::zplus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 0, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::yplus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 0, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::yplus | BC::FineCoarseInterface | FC::FC0, Patch::zminus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 0, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::yplus | BC::FineCoarseInterface | FC::FC0, Patch::zplus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 0, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::zminus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 0, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::zplus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 0, Patch::xplus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 0, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::yminus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 0, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::yminus | BC::FineCoarseInterface | FC::FC0, Patch::zminus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 0, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::yminus | BC::FineCoarseInterface | FC::FC0, Patch::zplus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 0, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::yplus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 0, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::yplus | BC::FineCoarseInterface | FC::FC0, Patch::zminus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 0, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::yplus | BC::FineCoarseInterface | FC::FC0, Patch::zplus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 0, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::zminus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 0, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::zplus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 0, Patch::yminus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 0, Patch::yminus | BC::FineCoarseInterface | FC::FC0, Patch::zminus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 0, Patch::yminus | BC::FineCoarseInterface | FC::FC0, Patch::zplus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 0, Patch::yplus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 0, Patch::yplus | BC::FineCoarseInterface | FC::FC0, Patch::zminus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 0, Patch::yplus | BC::FineCoarseInterface | FC::FC0, Patch::zplus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 0, Patch::zminus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 0, Patch::zplus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 1, Patch::xminus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 1, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::yminus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 1, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::yminus | BC::FineCoarseInterface | FC::FC0, Patch::zminus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 1, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::yminus | BC::FineCoarseInterface | FC::FC0, Patch::zplus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 1, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::yplus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 1, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::yplus | BC::FineCoarseInterface | FC::FC0, Patch::zminus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 1, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::yplus | BC::FineCoarseInterface | FC::FC0, Patch::zplus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 1, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::zminus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 1, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::zplus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 1, Patch::xplus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 1, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::yminus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 1, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::yminus | BC::FineCoarseInterface | FC::FC0, Patch::zminus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 1, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::yminus | BC::FineCoarseInterface | FC::FC0, Patch::zplus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 1, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::yplus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 1, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::yplus | BC::FineCoarseInterface | FC::FC0, Patch::zminus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 1, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::yplus | BC::FineCoarseInterface | FC::FC0, Patch::zplus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 1, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::zminus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 1, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::zplus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 1, Patch::yminus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 1, Patch::yminus | BC::FineCoarseInterface | FC::FC0, Patch::zminus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 1, Patch::yminus | BC::FineCoarseInterface | FC::FC0, Patch::zplus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 1, Patch::yplus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 1, Patch::yplus | BC::FineCoarseInterface | FC::FC0, Patch::zminus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 1, Patch::yplus | BC::FineCoarseInterface | FC::FC0, Patch::zplus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 1, Patch::zminus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 1, Patch::zplus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 2, Patch::xminus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 2, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::yminus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 2, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::yminus | BC::FineCoarseInterface | FC::FC0, Patch::zminus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 2, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::yminus | BC::FineCoarseInterface | FC::FC0, Patch::zplus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 2, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::yplus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 2, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::yplus | BC::FineCoarseInterface | FC::FC0, Patch::zminus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 2, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::yplus | BC::FineCoarseInterface | FC::FC0, Patch::zplus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 2, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::zminus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 2, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::zplus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 2, Patch::xplus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 2, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::yminus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 2, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::yminus | BC::FineCoarseInterface | FC::FC0, Patch::zminus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 2, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::yminus | BC::FineCoarseInterface | FC::FC0, Patch::zplus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 2, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::yplus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 2, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::yplus | BC::FineCoarseInterface | FC::FC0, Patch::zminus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 2, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::yplus | BC::FineCoarseInterface | FC::FC0, Patch::zplus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 2, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::zminus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 2, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::zplus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 2, Patch::yminus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 2, Patch::yminus | BC::FineCoarseInterface | FC::FC0, Patch::zminus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 2, Patch::yminus | BC::FineCoarseInterface | FC::FC0, Patch::zplus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 2, Patch::yplus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 2, Patch::yplus | BC::FineCoarseInterface | FC::FC0, Patch::zminus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 2, Patch::yplus | BC::FineCoarseInterface | FC::FC0, Patch::zplus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 2, Patch::zminus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 2, Patch::zplus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 3, Patch::xminus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 3, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::yminus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 3, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::yminus | BC::FineCoarseInterface | FC::FC0, Patch::zminus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 3, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::yminus | BC::FineCoarseInterface | FC::FC0, Patch::zplus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 3, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::yplus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 3, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::yplus | BC::FineCoarseInterface | FC::FC0, Patch::zminus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 3, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::yplus | BC::FineCoarseInterface | FC::FC0, Patch::zplus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 3, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::zminus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 3, Patch::xminus | BC::FineCoarseInterface | FC::FC0, Patch::zplus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 3, Patch::xplus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 3, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::yminus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 3, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::yminus | BC::FineCoarseInterface | FC::FC0, Patch::zminus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 3, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::yminus | BC::FineCoarseInterface | FC::FC0, Patch::zplus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 3, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::yplus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 3, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::yplus | BC::FineCoarseInterface | FC::FC0, Patch::zminus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 3, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::yplus | BC::FineCoarseInterface | FC::FC0, Patch::zplus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 3, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::zminus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 3, Patch::xplus | BC::FineCoarseInterface | FC::FC0, Patch::zplus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 3, Patch::yminus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 3, Patch::yminus | BC::FineCoarseInterface | FC::FC0, Patch::zminus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 3, Patch::yminus | BC::FineCoarseInterface | FC::FC0, Patch::zplus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 3, Patch::yplus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 3, Patch::yplus | BC::FineCoarseInterface | FC::FC0, Patch::zminus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 3, Patch::yplus | BC::FineCoarseInterface | FC::FC0, Patch::zplus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 3, Patch::zminus | BC::FineCoarseInterface | FC::FC0 >;
template class BCFDView < PureMetalProblem<CC, P7>, 3, Patch::zplus | BC::FineCoarseInterface | FC::FC0 >;

} // namespace Uintah
} // namespace PhaseField
