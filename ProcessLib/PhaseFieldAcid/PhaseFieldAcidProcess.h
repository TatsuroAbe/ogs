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

#include "LocalAssemblerInterface.h"
#include "PhaseFieldAcidProcessData.h"
#include "ProcessLib/Process.h"

namespace ProcessLib
{
namespace PhaseFieldAcid
{
class PhaseFieldAcidProcess final : public Process
{
public:
    PhaseFieldAcidProcess(
        std::string name,
        MeshLib::Mesh& mesh,
        std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&&
            jacobian_assembler,
        std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
            parameters,
        unsigned const integration_order,
        std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
            process_variables,
        PhaseFieldAcidProcessData&& process_data,
        SecondaryVariableCollection&& secondary_variables,
        bool const use_monolithic_scheme);

    //! \name ODESystem interface
    //! @{
    bool isLinear() const override;
    //! @}

    MathLib::MatrixSpecifications getMatrixSpecifications(
        const int process_id) const override;

    NumLib::LocalToGlobalIndexMap const& getDOFTable(
        const int process_id) const override;

private:
    using LocalAssemblerInterface = PhaseFieldAcidLocalAssemblerInterface;

    void constructDofTable() override;

    void initializeBoundaryConditions() override;

    void initializeConcreteProcess(
        NumLib::LocalToGlobalIndexMap const& dof_table,
        MeshLib::Mesh const& mesh,
        unsigned const integration_order) override;

    void assembleConcreteProcess(const double t, double const dt,
                                 std::vector<GlobalVector*> const& x,
                                 int const process_id, GlobalMatrix& M,
                                 GlobalMatrix& K, GlobalVector& b) override;

    void assembleWithJacobianConcreteProcess(
        const double t, double const dt, std::vector<GlobalVector*> const& x,
        GlobalVector const& xdot, const double dxdot_dx, const double dx_dx,
        int const process_id, GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b,
        GlobalMatrix& Jac) override;

    void preTimestepConcreteProcess(std::vector<GlobalVector*> const& x,
                                    double const t, double const dt,
                                    const int process_id) override;

    void postTimestepConcreteProcess(std::vector<GlobalVector*> const& x,
                                     const double t, const double delta_t,
                                     int const process_id) override;

    void postNonLinearSolverConcreteProcess(GlobalVector const& x,
                                            const double t, double const dt,
                                            int const process_id) override;

private:
    PhaseFieldAcidProcessData _process_data;

    std::vector<std::unique_ptr<LocalAssemblerInterface>> _local_assemblers;

    std::unique_ptr<NumLib::LocalToGlobalIndexMap>
        _local_to_global_index_map_single_component;

    /// Sparsity pattern for the phase field equation, and it is initialized
    ///  only if the staggered scheme is used.
    GlobalSparsityPattern _sparsity_pattern_with_single_component;

    /// Check whether the process represented by \c process_id is/has
    /// mechanical process. In the present implementation, the mechanical
    /// process has process_id == 0 in the staggered scheme.
    bool isPhaseFieldProcess(int const process_id) const;
};

}  // namespace PhaseFieldAcid
}  // namespace ProcessLib
