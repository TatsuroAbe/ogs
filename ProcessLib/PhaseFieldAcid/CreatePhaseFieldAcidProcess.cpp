/**
 * \file
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreatePhaseFieldAcidProcess.h"

#include <cassert>

#include "ParameterLib/Utils.h"
#include "PhaseFieldAcidProcess.h"
#include "PhaseFieldAcidProcessData.h"
#include "ProcessLib/Output/CreateSecondaryVariables.h"
#include "ProcessLib/Utils/ProcessUtils.h"

namespace ProcessLib
{
namespace PhaseField
{
std::unique_ptr<Process> createPhaseFieldAcidProcess(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    boost::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{prj__processes__process__type}
    config.checkConfigParameter("type", "PHASEFIELD_Acid");
    DBUG("Create PhaseFieldAcidProcess.");

    auto const staggered_scheme =
        //! \ogs_file_param{prj__processes__process__PHASEFIELD_ACID__coupling_scheme}
        config.getConfigParameterOptional<std::string>("coupling_scheme");
    const bool use_monolithic_scheme =
        !(staggered_scheme && (*staggered_scheme == "staggered"));

    // Process variable.

    //! \ogs_file_param{prj__processes__process__PHASEFIELD_ACID__process_variables}
    auto const pv_config = config.getConfigSubtree("process_variables");

    ProcessVariable* variable_ph;
    ProcessVariable* variable_c;
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>
        process_variables;
    if (use_monolithic_scheme)  // monolithic scheme.
    {
        OGS_FATAL("Monolithic implementation is not available.");
    }
    else  // staggered scheme.
    {
        using namespace std::string_literals;
        for (
            auto const& variable_name :
            {//! \ogs_file_param_special{prj__processes__process__PHASEFIELD_ACID__process_variables__concentration}
             "concentration"s,
             //! \ogs_file_param_special{prj__processes__process__PHASEFIELD_ACID__process_variables__phasefield}
             "phasefield"s})
        {
            auto per_process_variables =
                findProcessVariables(variables, pv_config, {variable_name});
            process_variables.push_back(std::move(per_process_variables));
        }
        variable_c = &process_variables[0][0].get();
        variable_ph = &process_variables[1][0].get();
    }

    DBUG("Associate concentration with process variable '%s'.",
         variable_c->getName().c_str());

    if (variable_c->getNumberOfComponents() != 1)
    {
        OGS_FATAL(
            "Phasefield process variable '%s' is not a scalar variable but has "
            "%d components.",
            variable_ph->getName().c_str(),
            variable_ph->getNumberOfComponents());
    }

    DBUG("Associate phase field with process variable '%s'.",
         variable_ph->getName().c_str());
    if (variable_ph->getNumberOfComponents() != 1)
    {
        OGS_FATAL(
            "Phasefield process variable '%s' is not a scalar variable but has "
            "%d components.",
            variable_ph->getName().c_str(),
            variable_ph->getNumberOfComponents());
    }

    auto const phasefield_parameters_config =
        //! \ogs_file_param{prj__processes__process__PHASEFIELD_ACID__phasefield_parameters}
        config.getConfigSubtree("phasefieldacid_parameters");

    // Specific body force parameter.
    Eigen::VectorXd specific_body_force;
    std::vector<double> const b =
        //! \ogs_file_param{prj__processes__process__HT__specific_body_force}
        config.getConfigParameter<std::vector<double>>("specific_body_force");
    assert(!b.empty() && b.size() < 4);
    if (b.size() < mesh.getDimension())
    {
        OGS_FATAL(
            "specific body force (gravity vector) has %d components, mesh "
            "dimension is %d",
            b.size(), mesh.getDimension());
    }
    bool const has_gravity = MathLib::toVector(b).norm() > 0;
    if (has_gravity)
    {
        specific_body_force.resize(b.size());
        std::copy_n(b.data(), b.size(), specific_body_force.data());
    }

    PhaseFieldAcidProcessData process_data{materialIDs(mesh),
                                           specific_body_force};

    SecondaryVariableCollection secondary_variables;

    NumLib::NamedFunctionCaller named_function_caller(
        {"PhaseField_displacement"});

    ProcessLib::createSecondaryVariables(config, secondary_variables,
                                         named_function_caller);

    return std::make_unique <
           PhaseFieldAcidProcess(
               std::move(name), mesh, std::move(jacobian_assembler), parameters,
               integration_order, std::move(process_variables),
               std::move(process_data), std::move(secondary_variables),
               std::move(named_function_caller), use_monolithic_scheme);
}

}  // namespace PhaseField
}  // namespace ProcessLib
