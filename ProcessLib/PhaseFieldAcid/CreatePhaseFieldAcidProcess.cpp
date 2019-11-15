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
namespace PhaseFieldAcid
{
std::unique_ptr<Process> createPhaseFieldAcidProcess(
    std::string name, MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    unsigned const integration_order, BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<MeshLib::Mesh>> const& meshes)
{
    //! \ogs_file_param{prj__processes__process__type}
    config.checkConfigParameter("type", "PHASEFIELD_ACID");
    DBUG("Create PhaseFieldAcidProcess.");

    auto const staggered_scheme =
        //! \ogs_file_param{prj__processes__process__PHASEFIELD_ACID__coupling_scheme}
        config.getConfigParameterOptional<std::string>("coupling_scheme");
    const bool use_monolithic_scheme =
        !(staggered_scheme && (*staggered_scheme == "staggered"));

    // Process variable.

    //! \ogs_file_param{prj__processes__process__PHASEFIELD_ACID__process_variables}
    auto const pv_config = config.getConfigSubtree("process_variables");

    ProcessVariable* variable_c;
    ProcessVariable* variable_phi;
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>
        process_variables;
    if (use_monolithic_scheme)  // monolithic scheme.
    {
        auto per_process_variables = findProcessVariables(
            variables, pv_config,
            {//! \ogs_file_param_special{prj__processes__process__PHASEFIELD_ACID__process_variables__concentration}
             "concentration",
             //! \ogs_file_param_special{prj__processes__process__PHASEFIELD_ACID__process_variables__phasefield}
             "phasefield"});
        variable_c = &per_process_variables[0].get();
        variable_phi = &per_process_variables[1].get();
        process_variables.push_back(std::move(per_process_variables));
    }
    else  // staggered scheme.
    {
        using namespace std::string_literals;
        for (auto const& variable_name : {"concentration"s, "phasefield"s})
        {
            auto per_process_variables =
                findProcessVariables(variables, pv_config, {variable_name});
            process_variables.push_back(std::move(per_process_variables));
        }
        variable_c = &process_variables[0][0].get();
        variable_phi = &process_variables[1][0].get();
    }

    DBUG("Associate concentration with process variable '%s'.",
         variable_c->getName().c_str());
    if (variable_c->getNumberOfComponents() != 1)
    {
        OGS_FATAL(
            "Phasefield_Acid process variable '%s' is not a scalar variable "
            "but has "
            "%d components.",
            variable_c->getName().c_str(),
            variable_c->getNumberOfComponents());
    }

    DBUG("Associate phasefield with process variable '%s'.",
         variable_phi->getName().c_str());
    if (variable_phi->getNumberOfComponents() != 1)
    {
        OGS_FATAL(
            "Phasefield_Acid process variable '%s' is not a scalar variable "
            "but has "
            "%d components.",
            variable_phi->getName().c_str(),
            variable_phi->getNumberOfComponents());
    }

    // Intrinsic permeability
    auto& intrinsic_permeability = ParameterLib::findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__PHASEFIELD_ACID__intrinsic_permeability}
        "intrinsic_permeability", parameters, 1, &mesh);

    DBUG("Use '%s' as intrinsic conductivity parameter.",
         intrinsic_permeability.name.c_str());

    // Fluid viscosity
    auto& fluid_viscosity = ParameterLib::findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__PHASEFIELD_ACID__fluid_viscosity}
        "fluid_viscosity", parameters, 1, &mesh);
    DBUG("Use '%s' as fluid viscosity parameter.",
         fluid_viscosity.name.c_str());

    // Fluid density
    auto& fluid_density = ParameterLib::findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__PHASEFIELD_ACID___fluid_density}
        "fluid_density", parameters, 1, &mesh);
    DBUG("Use '%s' as fluid density parameter.", fluid_density.name.c_str());

    // Porosity
    auto& porosity = ParameterLib::findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__PHASEFIELD_ACID___porosity}
        "porosity", parameters, 1, &mesh);
    DBUG("Use '%s' as porosity parameter.", porosity.name.c_str());

    // Reference temperature
    double const reference_temperature =
        //! \ogs_file_param{prj__processes__process__PHASEFIELD_ACID___reference_temperature}
        config.getConfigParameter<double>(
            "reference_temperature", std::numeric_limits<double>::quiet_NaN());
    DBUG("Use 'reference_temperature' as reference temperature.");

    // Specific gas constant
    double const specific_gas_constant =
        //! \ogs_file_param{prj__processes__process__PHASEFIELD_ACID___specific_gas_constant}
        config.getConfigParameter<double>(
            "specific_gas_constant", std::numeric_limits<double>::quiet_NaN());
    DBUG("Use 'specific_gas_constant' as specific gas constant.");

    // Fluid compressibility
    double const fluid_compressibility =
        //! \ogs_file_param{prj__processes__process__PHASEFIELD_ACID___fluid_compressibility}
        config.getConfigParameter<double>(
            "fluid_compressibility", std::numeric_limits<double>::quiet_NaN());
    DBUG("Use 'fluid_compressibility' as fluid compressibility parameter.");

    auto const fluid_type = FluidType::strToFluidType(
        //! \ogs_file_param{prj__processes__process__PHASEFIELD_ACID___fluid_type}
        config.getConfigParameter<std::string>("fluid_type"));
    DBUG("Use 'fluid_type' as fluid type parameter.");

    if (!FluidType::checkRequiredParams(fluid_type, fluid_compressibility,
                                        reference_temperature,
                                        specific_gas_constant))
    {
        OGS_FATAL(FluidType::getErrorMsg(fluid_type));
    }

    PhaseFieldAcidProcessData process_data{materialIDs(mesh),
                                           intrinsic_permeability,
                                           fluid_viscosity,
                                           fluid_density,
                                           porosity};

    SecondaryVariableCollection secondary_variables;

    ProcessLib::createSecondaryVariables(config, secondary_variables);

    return std::make_unique<PhaseFieldAcidProcess>(
        std::move(name), mesh, std::move(jacobian_assembler), parameters,
        integration_order, std::move(process_variables),
        std::move(process_data), std::move(secondary_variables),
        use_monolithic_scheme);
}

}  // namespace PhaseFieldAcid
}  // namespace ProcessLib
