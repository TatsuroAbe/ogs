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

#include <vector>

#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "ProcessLib/LocalAssemblerInterface.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "MeshLib/Elements/Elements.h"

namespace ProcessLib
{
namespace PhaseFieldAcid
{
struct PhaseFieldAcidLocalAssemblerInterface
    : public ProcessLib::LocalAssemblerInterface,
      public NumLib::ExtrapolatableElement
{
    Eigen::Vector3d getFlux(MathLib::Point3d const& pnt_local_coords,
                            double const t,
                            std::vector<double> const& local_x) const override =
        0;
};

}  // namespace PhaseFieldAcid
}  // namespace ProcessLib
