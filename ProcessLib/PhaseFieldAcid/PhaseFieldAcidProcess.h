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

#include "PhaseFieldAcidLocalAssemblerInterface.h"
#include "PhaseFieldAcidProcessData.h"
#include "ProcessLib/Process.h"

namespace NumLib
{
class LocalToGlobalIndexMap;
}

namespace ProcessLib
{
struct SurfaceFluxData;

namespace PhaseFieldAcid
{

class PhaseFieldAcidProcess final : public Process
{
public:
    PhaseFieldAcidProcess(
        std::string name, MeshLib::Mesh& mesh,
        std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&&
            jacobian_assembler,
        std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
            parameters,
        unsigned const integration_order,
        std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
            process_variables,
        PhaseFieldAcidProcessData&& process_data,
        SecondaryVariableCollection&& secondary_variables,
        bool const use_monolithic_scheme,
        std::unique_ptr<ProcessLib::SurfaceFluxData>&& surfaceflux,
        const int concentration_process_id, const int phasefield_process_id);

    //! \name ODESystem interface
    //! @{

    bool isLinear() const override;
    //! @}
    //!
    void setCoupledTermForTheStaggeredSchemeToLocalAssemblers(
        int const process_id) override;

    void postTimestepConcreteProcess(std::vector<GlobalVector*> const& x,
                                     const double t,
                                     const double delta_t,
                                     int const process_id) override;

private:
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

    void setCoupledSolutionsOfPreviousTimeStepPerProcess(const int process_id);

    void setCoupledSolutionsOfPreviousTimeStep();

    std::tuple<NumLib::LocalToGlobalIndexMap*, bool>
    getDOFTableForExtrapolatorData() const override;

    PhaseFieldAcidProcessData _process_data;

    void constructDofTable() override;

    void initializeBoundaryConditions() override;

    void postNonLinearSolverConcreteProcess(GlobalVector const& x,
                                            const double t, double const dt,
                                            int const process_id) override;

    std::vector<std::unique_ptr<PhaseFieldAcidLocalAssemblerInterface>> _local_assemblers;

    /// Solutions of the previous time step
    std::array<std::unique_ptr<GlobalVector>, 2> _xs_previous_timestep;

    std::unique_ptr<ProcessLib::SurfaceFluxData> _surfaceflux;

    const int _concentration_process_id;
    const int _phasefield_process_id;
};

}  // namespace PhaseFieldAcid
}  // namespace ProcessLib
