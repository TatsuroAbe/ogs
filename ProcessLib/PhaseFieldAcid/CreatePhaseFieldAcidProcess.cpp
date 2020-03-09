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

#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "ParameterLib/Parameter.h"
#include "ParameterLib/Utils.h"
#include "PhaseFieldAcidProcess.h"
#include "PhaseFieldAcidProcessData.h"
#include "ProcessLib/Output/CreateSecondaryVariables.h"
#include "ProcessLib/Utils/ProcessUtils.h"

namespace ProcessLib
{
namespace PhaseFieldAcid
{
std::unique_ptr<Process> createPhaseFieldAcidProcess(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    boost::optional<
        ParameterLib::CoordinateSystem> const& /*local_coordinate_system*/,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{prj__processes__process__type}
    config.checkConfigParameter("type", "PHASEFIELDACID");
    DBUG("Create PhaseFieldAcidProcess.");

    auto const staggered_scheme =
        //! \ogs_file_param{prj__processes__process__PHASEFIELD_ACID__coupling_scheme}
        config.getConfigParameterOptional<std::string>("coupling_scheme");
    const bool use_monolithic_scheme =
        !(staggered_scheme && (*staggered_scheme == "staggered"));

    // Process variable.

    //! \ogs_file_param{prj__processes__process__PHASEFIELD_ACID__process_variables}
    auto const pv_config = config.getConfigSubtree("process_variables");

    const int _concentration_process_id = 0;
    const int _phasefield_process_id = 1;
    ProcessVariable* variable_c;
    ProcessVariable* variable_ph;
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

        variable_c = &process_variables[_concentration_process_id][0].get();
        variable_ph = &process_variables[_phasefield_process_id][0].get();
    }


    DBUG("Associate concentration with process variable '%s'.",
         variable_c->getName().c_str());

    if (variable_c->getNumberOfComponents() != 1)
    {
        OGS_FATAL(
            "Phasefield process variable '%s' is not a scalar variable but has "
            "%d components.",
            variable_c->getName().c_str(),
            variable_c->getNumberOfComponents());
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

    auto const phasefieldacid_parameters_config =
        //! \ogs_file_param{prj__processes__process__PHASEFIELD_ACID__phasefieldacid_parameters}
        config.getConfigSubtree("phasefieldacid_parameters");

    // Chemical diffusivity
    auto& chemical_diffusivity = ParameterLib::findParameter<double>(
        phasefieldacid_parameters_config,
        //! \ogs_file_param_special{prj__processes__process__PHASEFIELD_ACID__phasefieldacid_parameters__checmial_diffusivity}
        "chemical_diffusivity", parameters, 1);
    DBUG("Use '%s' as chemical diffusivity.",
         chemical_diffusivity.name.c_str());

    // alpha
    auto& alpha = ParameterLib::findParameter<double>(
        phasefieldacid_parameters_config,
        //! \ogs_file_param_special{prj__processes__process__PHASEFIELD_ACID__phasefieldacid_parameters__alpha}
        "alpha", parameters, 1);
    DBUG("Use '%s' as alpha.", alpha.name.c_str());

    // rrc(reaction rate coefficient, k)
    auto& rrc = ParameterLib::findParameter<double>(
        phasefieldacid_parameters_config,
        //! \ogs_file_param_special{prj__processes__process__PHASEFIELD_ACID__phasefieldacid_parameters__rrc}
        "rrc", parameters, 1);
    DBUG("Use '%s' as rrc.", rrc.name.c_str());

    // epsilon
    auto& epsilon = ParameterLib::findParameter<double>(
        phasefieldacid_parameters_config,
        //! \ogs_file_param_special{prj__processes__process__PHASEFIELD_ACID__phasefieldacid_parameters__epsilon}
        "epsilon", parameters, 1);
    DBUG("Use '%s' as epsi.", epsilon.name.c_str());

    // tau
    auto& tau = ParameterLib::findParameter<double>(
        phasefieldacid_parameters_config,
        //! \ogs_file_param_special{prj__processes__process__PHASEFIELD_ACID__phasefieldacid_parameters__tau}
        "tau", parameters, 1);
    DBUG("Use '%s' as tau.", tau.name.c_str());

    // grad_phi_cutoff
    auto const grad_phi_cutoff =
        //! \ogs_file_param_special{prj__processes__process__PHASEFIELD_ACID__phasefieldacid_parameters__grad_phi_cutoff}
        phasefieldacid_parameters_config.getConfigParameter<double>(
            "grad_phi_cutoff");

    // Specific body force parameter.
    Eigen::VectorXd specific_body_force;
    std::vector<double> const b =
        //! \ogs_file_param{prj__processes__process__PHASEFIELD_ACID__specific_body_force}
        config.getConfigParameter<std::vector<double>>("specific_body_force");
    if (b.size() != mesh.getDimension())
    {
        OGS_FATAL(
            "specific body force (gravity vector) has %d components, mesh "
            "dimension is %d",
            b.size(), mesh.getDimension());
    }
    specific_body_force.resize(b.size());
    std::copy_n(b.data(), b.size(), specific_body_force.data());

    PhaseFieldAcidProcessData process_data{materialIDs(mesh),
                                           specific_body_force,
                                           chemical_diffusivity,
                                           alpha,
                                           rrc,
                                           epsilon,
                                           tau,
                                           grad_phi_cutoff};

    SecondaryVariableCollection secondary_variables;

    ProcessLib::createSecondaryVariables(config, secondary_variables);

    return std::make_unique<PhaseFieldAcidProcess>(
        std::move(name), mesh, std::move(jacobian_assembler), parameters,
        integration_order, std::move(process_variables),
        std::move(process_data), std::move(secondary_variables),
        use_monolithic_scheme,_concentration_process_id,_phasefield_process_id);
}

}  // namespace PhaseFieldAcid
}  // namespace ProcessLib
