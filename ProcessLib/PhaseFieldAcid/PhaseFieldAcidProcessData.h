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

#include "ParameterLib/Parameter.h"
#include "MaterialLib/Fluid/FluidType/FluidType.h"

#include <memory>
#include <utility>

#include <Eigen/Dense>


namespace ProcessLib
{
namespace PhaseFieldAcid
{
struct PhaseFieldAcidProcessData final
{
    MeshLib::PropertyVector<int> const* const material_ids = nullptr;
    ParameterLib::Parameter<double> const& intrinsic_permeability;
    ParameterLib::Parameter<double> const& fluid_viscosity;
    ParameterLib::Parameter<double> const& fluid_density;
    ParameterLib::Parameter<double> const& porosity;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

}  // namespace PhaseFieldAcid
}  // namespace ProcessLib
