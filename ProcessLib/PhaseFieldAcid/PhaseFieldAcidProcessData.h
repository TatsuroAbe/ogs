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
namespace ProcessLib
{
template <typename T>
struct Parameter;

namespace PhaseFieldAcid
{
struct PhaseFieldAcidProcessData
{
    PhaseFieldAcidProcessData(
        MeshLib::PropertyVector<int> const* const material_ids_,
        Eigen::VectorXd const& specific_body_force,
        Parameter<double> const& chemical_diffusivity_;
                              : specific_body_force(specific_body_force_),
                                  chemical_diffusivity(chemical_diffusivity_)
    )
    MeshLib::PropertyVector<int> const* const material_ids = nullptr;
    Eigen::VectorXd const specific_body_force;
    Parameter<double> const& chemical_diffusivity;
    static constexpr int concentration_process_id = 0;
    static constexpr int phasefield_process_id = 1;
};

}  // namespace PhaseFieldAcid
}  // namespace ProcessLib
