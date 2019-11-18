/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file
 *  Created on January 8, 2018, 3:00 PM
 */
#pragma once

#include "NumLib/DOF/DOFTableUtil.h"
#include "ProcessLib/CoupledSolutionsForStaggeredScheme.h"

namespace ProcessLib
{
namespace PhaseFieldAcid
{
template <typename ShapeFunction, typename IntegrationMethod, int GlobalDim>
void PhaseFieldAcidLocalAssembler<ShapeFunction, IntegrationMethod, GlobalDim>::
    assembleWithJacobianForStaggeredScheme(
        double const t, double const dt, std::vector<double> const& local_xdot,
        const double dxdot_dx, const double dx_dx, int const process_id,
        std::vector<double>& local_M_data, std::vector<double>& local_K_data,
        std::vector<double>& local_b_data, std::vector<double>& local_Jac_data,
        LocalCoupledSolutions const& local_coupled_solutions)
{
    // For the equations with phase field.
    if (process_id == 0)
    {
        assembleWithJacobianPhaseFiledEquations(
            t, dt, local_xdot, dxdot_dx, dx_dx, local_M_data, local_K_data,
            local_b_data, local_Jac_data, local_coupled_solutions);
        return;
    }

    // For the equations with concentration
    assembleWithJacobianForConcentrationEquations(
        t, dt, local_xdot, dxdot_dx, dx_dx, local_M_data, local_K_data,
        local_b_data, local_Jac_data, local_coupled_solutions);
}

template <typename ShapeFunction, typename IntegrationMethod, int GlobalDim>
void PhaseFieldAcidLocalAssembler<ShapeFunction, IntegrationMethod, GlobalDim>::
    assembleWithJacobianForConcentrationEquations(
        double const t, double const dt,
        std::vector<double> const& /*local_xdot*/, const double /*dxdot_dx*/,
        const double /*dx_dx*/, std::vector<double>& /*local_M_data*/,
        std::vector<double>& /*local_K_data*/,
        std::vector<double>& local_b_data, std::vector<double>& local_Jac_data,
        LocalCoupledSolutions const& local_coupled_solutions)
{

    auto const& local_phi =
        local_coupled_solutions.local_coupled_xs[_phase_field_process_id];
    auto const& local_c =
        local_coupled_solutions.local_coupled_xs[_concentration_process_id];
    assert(local_c.size() == concentration_size);
    assert(local_phi.size() == phasefield_size);

    auto phi = Eigen::Map<
        typename ShapeMatricesType::template VectorType<phasefield_size> const>(
        local_phi.data(), phasefield_size);

    auto c = Eigen::Map<
        typename ShapeMatricesType::template VectorType<concentration_size> const>(
        local_c.data(), concentration_size);

    auto c_dot = Eigen::Map<
        typename ShapeMatricesType::template VectorType<concentration_size> const>(
        local_xdot.data(), concentration_size);

    auto local_Jac = MathLib::createZeroedMatrix<
        typename ShapeMatricesType::template MatrixType<concentration_size,
                                                        concentration_size>>(
        local_Jac_data, concentration_size, concentration_size);

    auto local_rhs = MathLib::createZeroedVector<
        typename ShapeMatricesType::template VectorType<concentration_size>>(
        local_b_data, concentration_size);

    typename ShapeMatricesType::NodalMatrixType mass =
        ShapeMatricesType::NodalMatrixType::Zero(concentration_size, concentration_size);

    typename ShapeMatricesType::NodalMatrixType laplace =
        ShapeMatricesType::NodalMatrixType::Zero(concentration_size, concentration_size);



    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    int const n_integration_points = _integration_method.getNumberOfPoints();

    for (int ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);
        auto const& w = _ip_data[ip].integration_weight;
        auto const& N = _ip_data[ip].N;
        auto const& dNdx = _ip_data[ip].dNdx;
        double const phi_ip = N.dot(phi);

        auto const x_coord =
            interpolateXCoordinate<ShapeFunction, ShapeMatricesType>(_element,
                                                                     N);

        auto const& b = _process_data.specific_body_force;

        local_rhs.noalias() +=
            (1.0) * N * w;
        mass.noalias() += (1.0) *
                              N.transpose() * N * w;


        laplace.noalias()
          += ( dNdx.transpose() * dNdx) * w;
    }
    local_Jac.noalias() = laplace + mass / dt;

    local_rhs.noalias() -= laplace * c + mass * c_dot;

}

template <typename ShapeFunction, typename IntegrationMethod, int GlobaltDim>
void PhaseFieldLocalAssembler<ShapeFunction, IntegrationMethod,
                              DisplacementDim>::
    assembleWithJacobianPhaseFiledEquations(
        double const t, double const dt,
        std::vector<double> const& /*local_xdot*/, const double /*dxdot_dx*/,
        const double /*dx_dx*/, std::vector<double>& /*local_M_data*/,
        std::vector<double>& /*local_K_data*/,
        std::vector<double>& local_b_data, std::vector<double>& local_Jac_data,
        LocalCoupledSolutions const& local_coupled_solutions)
{

    auto const& local_phi =
        local_coupled_solutions.local_coupled_xs[_phase_field_process_id];
    auto const& local_c =
        local_coupled_solutions.local_coupled_xs[_concentration_process_id];
    assert(local_c.size() == concentration_size);
    assert(local_phi.size() == phasefield_size);

    auto phi = Eigen::Map<
        typename ShapeMatricesType::template VectorType<phasefield_size> const>(
        local_phi.data(), phasefield_size);

    auto c = Eigen::Map<
        typename ShapeMatricesType::template VectorType<concentration_size> const>(
        local_c.data(), concentration_size);

    auto phi_dot = Eigen::Map<
        typename ShapeMatricesType::template VectorType<phasefield_size> const>(
        local_xdot.data(), phasefield_size);

    auto local_Jac = MathLib::createZeroedMatrix<
        typename ShapeMatricesType::template MatrixType<phasefield_size,
                                                        concentration_size>>(
        local_Jac_data, phasefield_size, phasefield_size);

    auto local_rhs = MathLib::createZeroedVector<
        typename ShapeMatricesType::template VectorType<phasefield_size>>(
        local_b_data, phasefield_size);

    typename ShapeMatricesType::NodalMatrixType mass =
        ShapeMatricesType::NodalMatrixType::Zero(phasefield_size, phasefield_size);

    typename ShapeMatricesType::NodalMatrixType laplace =
        ShapeMatricesType::NodalMatrixType::Zero(phasefield_size, phasefield_size);



    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    int const n_integration_points = _integration_method.getNumberOfPoints();

    for (int ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);
        auto const& w = _ip_data[ip].integration_weight;
        auto const& N = _ip_data[ip].N;
        auto const& dNdx = _ip_data[ip].dNdx;
        double const phi_ip = N.dot(phi);

        auto const x_coord =
            interpolateXCoordinate<ShapeFunction, ShapeMatricesType>(_element,
                                                                     N);

        auto const& b = _process_data.specific_body_force;

        local_rhs.noalias() +=
            (1.0) * N * w;
        mass.noalias() += (1.0) *
                              N.transpose() * N * w;


        laplace.noalias()
          += ( dNdx.transpose() * dNdx) * w;
    }
    local_Jac.noalias() = laplace + mass / dt;

    local_rhs.noalias() -= laplace * c + mass * c_dot;
}

}  // namespace PhaseFieldAcid
}  // namespace ProcessLib
