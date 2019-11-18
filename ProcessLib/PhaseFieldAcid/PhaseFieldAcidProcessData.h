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

#include <memory>
#include <utility>

#include <Eigen/Eigen>

namespace MaterialLib
{
}  // namespace MaterialLib
namespace ProcessLib
{
template <typename T>
struct Parameter;

namespace PhaseFieldAcid
{
struct PhaseFieldAcidProcessData
{

    Eigen::VectorXd const specific_body_force;

};

}  // namespace PhaseFieldAcid
}  // namespace ProcessLib
