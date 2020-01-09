/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "PhaseFieldAcidProcess.h"

#include <cassert>

#include "NumLib/DOF/ComputeSparsityPattern.h"
#include "PhaseFieldAcidFEM.h"
#include "ProcessLib/Process.h"
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
    bool const use_monolithic_scheme, const int concentration_process_id,
    const int phasefield_process_id)
    : Process(std::move(name), mesh, std::move(jacobian_assembler), parameters,
              integration_order, std::move(process_variables),
              std::move(secondary_variables), use_monolithic_scheme),
      _process_data(std::move(process_data)),
      _concentration_process_id(concentration_process_id),
      _phasefield_process_id(phasefield_process_id)
{
    if (use_monolithic_scheme)
    {
        OGS_FATAL(
            "Monolithic scheme is not implemented for the PhaseField process.");
    }
}

bool PhaseFieldAcidProcess::isLinear() const
{
    return false;
}

MathLib::MatrixSpecifications PhaseFieldAcidProcess::getMatrixSpecifications(
    const int process_id) const
{
    auto const& l = *_local_to_global_index_map_single_component;
    return {l.dofSizeWithoutGhosts(), l.dofSizeWithoutGhosts(),
            &l.getGhostIndices(), &_sparsity_pattern_with_single_component};
}

NumLib::LocalToGlobalIndexMap const& PhaseFieldAcidProcess::getDOFTable(
    const int process_id) const
{
    // For the equation of phasefield
    return *_local_to_global_index_map_single_component;
}

void PhaseFieldAcidProcess::constructDofTable()
{
    const int concentration_process_id = 0;
    constructDofTableOfSpecifiedProsessStaggerdScheme(concentration_process_id);

    // TODO move the two data members somewhere else.
    // for extrapolation of secondary variables of stress or strain
    std::vector<MeshLib::MeshSubset> all_mesh_subsets_single_component{
        *_mesh_subset_all_nodes};
    _local_to_global_index_map_single_component =
        std::make_unique<NumLib::LocalToGlobalIndexMap>(
            std::move(all_mesh_subsets_single_component),
            // by location order is needed for output
            NumLib::ComponentOrder::BY_LOCATION);

    assert(_local_to_global_index_map_single_component);

    // For phase field equation.
    _sparsity_pattern_with_single_component = NumLib::computeSparsityPattern(
        *_local_to_global_index_map_single_component, _mesh);
}

void PhaseFieldAcidProcess::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
    const int process_id = 0;
    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];
    ProcessLib::createLocalAssemblers<PhaseFieldAcidLocalAssembler>(
        mesh.getDimension(), mesh.getElements(), dof_table,
        pv.getShapeFunctionOrder(), _local_assemblers,
        mesh.isAxiallySymmetric(), integration_order, _process_data,
        _concentration_process_id, _phasefield_process_id);

    _process_data.element_kappa = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "kappa_avg",
        MeshLib::MeshItemType::Cell, 1);

    // Initialize local assemblers after all variables have been set.
    GlobalExecutor::executeMemberOnDereferenced(
        &LocalAssemblerInterface::initialize, _local_assemblers,
        *_local_to_global_index_map);
}

void PhaseFieldAcidProcess::initializeBoundaryConditions()
{
    // Staggered scheme:
    // for the equations of deformation.
    const int concentration_process_id = 0;
    initializeProcessBoundaryConditionsAndSourceTerms(
        *_local_to_global_index_map, concentration_process_id);
    // for the phase field
    const int phasefield_process_id = 1;
    initializeProcessBoundaryConditionsAndSourceTerms(
        *_local_to_global_index_map_single_component, phasefield_process_id);
}

void PhaseFieldAcidProcess::assembleConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    int const process_id, GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b)
{
    DBUG("Assemble PhaseFieldAcidProcess.");

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_tables;
    switch (process_id)
    {
        case 1:
        {
            DBUG(
                "Assemble the equations of phase-field in "
                "PhaseFieldAcidProcess for the staggered scheme.");
            break;
        }
        case 0:
        {
            DBUG(
                "Assemble the equations of concentration in "
                "PhaseFieldAcidProcess for the staggered scheme.");
            break;
        }
    }

    setCoupledSolutionsOfPreviousTimeStep();
    dof_tables.emplace_back(*_local_to_global_index_map_single_component);
    dof_tables.emplace_back(*_local_to_global_index_map_single_component);

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
    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_tables;

    // For the staggered scheme
    if (process_id == 1)
    {
        DBUG(
            "Assemble the Jacobian equations of concentration in "
            "PhaseFieldAcidProcess for the staggered scheme.");
    }
    else
    {
        DBUG(
            "Assemble the Jacobian equations of phase-field in "
            "PhaseFieldAcidProcess for the staggered scheme.");
    }
    dof_tables.emplace_back(*_local_to_global_index_map_single_component);
    dof_tables.emplace_back(*_local_to_global_index_map_single_component);

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeSelectedMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
        _local_assemblers, pv.getActiveElementIDs(), dof_tables, t, dt, x, xdot,
        dxdot_dx, dx_dx, process_id, M, K, b, Jac, _coupled_solutions);
}

void PhaseFieldAcidProcess::preTimestepConcreteProcess(
    std::vector<GlobalVector*> const& x, double const t, double const dt,
    const int process_id)
{
    DBUG("PreTimestep PhaseFieldAcidProcess %d.", process_id);

    /*  TODO (yoshioka) Maybe needed for xdot.
    _x_previous_timestep =
        MathLib::MatrixVectorTraits<GlobalVector>::newInstance(x);
        */
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
    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &LocalAssemblerInterface::preTimestep, _local_assemblers,
        pv.getActiveElementIDs(), getDOFTable(process_id), *x[process_id], t,
        dt);
}

void PhaseFieldAcidProcess::postTimestepConcreteProcess(
    std::vector<GlobalVector*> const& x, const double t,
    const double dt, int const process_id)
{
    if (isPhaseFieldProcess(process_id))
    {
        DBUG("PostTimestep PhaseFieldProcess.");

        std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
            dof_tables;

        dof_tables.emplace_back(*_local_to_global_index_map_single_component);
        dof_tables.emplace_back(*_local_to_global_index_map_single_component);

        ProcessLib::ProcessVariable const& pv =
            getProcessVariables(process_id)[0];

        GlobalExecutor::executeSelectedMemberOnDereferenced(
            &LocalAssemblerInterface::postTimestep, _local_assemblers,
            pv.getActiveElementIDs(), getDOFTable(process_id), *x[process_id], t,
            dt);
    }
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

void PhaseFieldAcidProcess::postNonLinearSolverConcreteProcess(
    GlobalVector const& x, const double t, double const dt,
    const int process_id)
{
    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_tables;

    dof_tables.emplace_back(*_local_to_global_index_map_single_component);
    dof_tables.emplace_back(*_local_to_global_index_map_single_component);

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];
}

constexpr bool PhaseFieldAcidProcess::isPhaseFieldProcess(
    int const process_id) const
{
    return process_id == _phasefield_process_id;
}

}  // namespace PhaseFieldAcid
}  // namespace ProcessLib
