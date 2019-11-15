/**
 * \file
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "PhaseFieldAcidProcess.h"

#include <cassert>

#include "NumLib/DOF/DOFTableUtil.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "PhaseFieldAcidFEM.h"
#include "ProcessLib/SurfaceFlux/SurfaceFluxData.h"
#include "ProcessLib/Utils/CreateLocalAssemblers.h"

namespace ProcessLib
{
namespace PhaseFieldAcid
{
PhaseFieldAcidProcess::PhaseFieldAcidProcess(
    std::string name, MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    unsigned const integration_order,
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
        process_variables,
    PhaseFieldAcidProcessData&& process_data,
    SecondaryVariableCollection&& secondary_variables,
    bool const use_monolithic_scheme,
    std::unique_ptr<ProcessLib::SurfaceFluxData>&& surfaceflux,
    const int concentration_process_id, const int phasefield_process_id)
    : Process(std::move(name), mesh, std::move(jacobian_assembler), parameters,
              integration_order, std::move(process_variables),
              std::move(secondary_variables), use_monolithic_scheme),
      _process_data(std::move(process_data)),
      _surfaceflux(std::move(surfaceflux)),
      _concentration_process_id(concentration_process_id),
      _phasefield_process_id(phasefield_process_id)
{
    if (use_monolithic_scheme)
    {
        OGS_FATAL(
            "Monolithic scheme is not implemented for the PhaseField_Acid "
            "process.");
    }
}

void PhaseFieldAcidProcess::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
    const int process_id = 0;
    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    if (_use_monolithic_scheme)
    {
        OGS_FATAL(
            "Monolithic scheme is not implemented for the PhaseField_Acid "
            "process.");
    }
    else
    {
        ProcessLib::createLocalAssemblers<PhaseFieldAcidProcess>(
            mesh.getDimension(), mesh.getElements(), dof_table,
            pv.getShapeFunctionOrder(), _local_assemblers,
            mesh.isAxiallySymmetric(), integration_order, _process_data,
            _concentration_process_id, _phasefield_process_id);
    }

    _secondary_variables.addSecondaryVariable(
        "darcy_velocity",
        makeExtrapolator(
            mesh.getDimension(), getExtrapolator(), _local_assemblers,
            &PhaseFieldAcidLocalAssemblerInterface::getIntPtDarcyVelocity));
}

void PhaseFieldAcidProcess::assembleConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    int const process_id, GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b)
{
    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_tables;
    if (_use_monolithic_scheme)
    {
        OGS_FATAL(
            "Monolithic scheme is not implemented for the PhaseField_Acid "
            "process.");
    }
    else
    {
        if (process_id == _concentration_process_id)
        {
            DBUG(
                "Assemble the equations of concentration process within "
                "PhaseFieldAcidProcess.");
        }
        else
        {
            DBUG(
                "Assemble the equations of phase-field process "
                "within PhaseFieldAcidProcess.");
        }
        setCoupledSolutionsOfPreviousTimeStep();
        dof_tables.emplace_back(*_local_to_global_index_map);
        dof_tables.emplace_back(*_local_to_global_index_map);
    }

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];
    // Call global assembler for each local assembly item.
    GlobalExecutor::executeSelectedMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assemble, _local_assemblers,
        pv.getActiveElementIDs(), dof_tables, t, dt, x, process_id, M, K, b,
        _coupled_solutions);
}

void PhaseFieldAcidProcess::assembleWithJacobianConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    GlobalVector const& xdot, const double dxdot_dx, const double dx_dx,
    int const process_id, GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b,
    GlobalMatrix& Jac)
{
    DBUG("AssembleWithJacobian PhaseFieldAcidProcess.");

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_tables;
    if (!_use_monolithic_scheme)
    {
        OGS_FATAL(
            "Monolithic scheme is not implemented for the PhaseField_Acid "
            "process.");
    }
    else
    {
        dof_tables.emplace_back(std::ref(*_local_to_global_index_map));
        dof_tables.emplace_back(std::ref(*_local_to_global_index_map));
    }

    // Call global assembler for each local assembly item.
    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];
    GlobalExecutor::executeSelectedMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
        _local_assemblers, pv.getActiveElementIDs(), dof_tables, t, dt, x, xdot,
        dxdot_dx, dx_dx, process_id, M, K, b, Jac, _coupled_solutions);
}

void PhaseFieldAcidProcess::preTimestepConcreteProcess(
    std::vector<GlobalVector*> const& x, double const t, double const dt,
    const int process_id)
{
    assert(process_id < 2);

    if (_use_monolithic_scheme)
    {
        return;
    }

    if (!_xs_previous_timestep[process_id])
    {
        _xs_previous_timestep[process_id] =
            MathLib::MatrixVectorTraits<GlobalVector>::newInstance(
                *x[process_id]);
    }
    else
    {
        auto& x0 = *_xs_previous_timestep[process_id];
        MathLib::LinAlg::copy(*x[process_id], x0);
    }

    auto& x0 = *_xs_previous_timestep[process_id];
    MathLib::LinAlg::setLocalAccessibleVector(x0);
}


std::tuple<NumLib::LocalToGlobalIndexMap*, bool>
PhaseFieldAcidProcess::getDOFTableForExtrapolatorData() const
{
    if (!_use_monolithic_scheme)
    {
        // For single-variable-single-component processes reuse the existing DOF
        // table.
        const bool manage_storage = false;
        return std::make_tuple(_local_to_global_index_map.get(),
                               manage_storage);
    }

    // Otherwise construct a new DOF table.
    std::vector<MeshLib::MeshSubset> all_mesh_subsets_single_component{
        *_mesh_subset_all_nodes};

    const bool manage_storage = true;
    return std::make_tuple(new NumLib::LocalToGlobalIndexMap(
                               std::move(all_mesh_subsets_single_component),
                               // by location order is needed for output
                               NumLib::ComponentOrder::BY_LOCATION),
                           manage_storage);
}



void PhaseFieldAcidProcess::setCoupledSolutionsOfPreviousTimeStepPerProcess(
    const int process_id)
{
    const auto& x_t0 = _xs_previous_timestep[process_id];
    if (x_t0 == nullptr)
    {
        OGS_FATAL(
            "Memory is not allocated for the global vector of the solution of "
            "the previous time step for the staggered scheme.\n It can be done "
            "by overriding Process::preTimestepConcreteProcess (ref. "
            "HTProcess::preTimestepConcreteProcess) ");
    }

    _coupled_solutions->coupled_xs_t0[process_id] = x_t0.get();
}

void PhaseFieldAcidProcess::setCoupledSolutionsOfPreviousTimeStep()
{
    _coupled_solutions->coupled_xs_t0.resize(2);
    setCoupledSolutionsOfPreviousTimeStepPerProcess(_concentration_process_id);
    setCoupledSolutionsOfPreviousTimeStepPerProcess(_phasefield_process_id);
}

void PhaseFieldAcidProcess::postTimestepConcreteProcess(
    std::vector<GlobalVector*> const& x, double const t, double const dt,
    const int process_id)
{
    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];
}

}  // namespace PhaseFieldAcid
}  // namespace ProcessLib
