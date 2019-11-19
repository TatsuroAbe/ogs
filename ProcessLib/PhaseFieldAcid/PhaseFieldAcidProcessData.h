/**
 * \file
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#pragma once

#include <Eigen/Eigen>
#include <memory>
#include <utility>

namespace MeshLib
{
template <typename T>
class PropertyVector;
}
namespace ParameterLib
{
template <typename T>
struct Parameter;
}
namespace ProcessLib::PhaseFieldAcid
{
struct PhaseFieldAcidProcessData
{

    MeshLib::PropertyVector<int> const* const material_ids = nullptr;
    Eigen::VectorXd const specific_body_force;

    ParameterLib::Parameter<double> const& chemical_diffusivity;
    ParameterLib::Parameter<double> const& alpha;
    ParameterLib::Parameter<double> const& rrc;
    ParameterLib::Parameter<double> const& epsilon;
    ParameterLib::Parameter<double> const& tau;

    static constexpr int concentration_process_id = 0;
    static constexpr int phasefield_process_id = 1;
};

}  // namespace PhaseFieldAcid

