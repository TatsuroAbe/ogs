/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file
 *  Created on November 29, 2017, 2:03 PM
 */

#pragma once

#include "MaterialLib/SolidModels/SelectSolidConstitutiveRelation.h"
#include "MathLib/KelvinVector.h"
#include "NumLib/Function/Interpolation.h"
#include "PhaseFieldAcidFEM.h"
#include "ProcessLib/CoupledSolutionsForStaggeredScheme.h"

namespace ProcessLib
{
namespace PhaseFieldAcid
{
template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
void PhaseFieldAcidFEM<ShapeFunction, IntegrationMethod, GlobalDim>::
    assembleForStaggeredScheme(double const t, double const dt,
                               int const process_id,
                               std::vector<double>& local_M_data,
                               std::vector<double>& local_K_data,
                               std::vector<double>& local_b_data,
                               LocalCoupledSolutions const& coupled_xs)
{
    if (process_id == _concentration_process_id)
    {
        assembleConcentrationEquation(t, local_M_data, local_K_data,
                                      local_b_data, coupled_xs);
        return;
    }

    assemblePhasefieldEquation(t, dt, local_M_data, local_K_data, local_b_data,
                               coupled_xs);
}

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
void PhaseFieldAcidLocalAssembler<ShapeFunction, IntegrationMethod, GlobalDim>::
    assembleWithJacobianForConcentrationEquations(
        const double t, double const dt, const std::vector<double>& local_xdot,
        const double /*dxdot_dx*/, const double /*dx_dx*/,
        std::vector<double>& /*local_M_data*/,
        std::vector<double>& /*local_K_data*/,
        std::vector<double>& local_b_data, std::vector<double>& local_Jac_data,
        const LocalCoupledSolutions& local_coupled_solutions)
{
    auto local_rhs =
        MathLib::createZeroedVector<typename ShapeMatricesTypeDisplacement::
                                        template VectorType<pressure_size>>(
            local_b_data, pressure_size);

    ParameterLib::SpatialPosition pos;
    pos.setElementID(this->_element.getID());

    auto const c =
        Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
            concentration_size> const>(
            &local_coupled_solutions.local_coupled_xs[concentration_index],
            concentration_size);

    auto const phi =
        Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
            phasefield_size> const>(
            &local_coupled_solutions.local_coupled_xs[phasefield_index],
            phasefield_size);

    auto c_dot =
        Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
            concentration_size> const>(local_xdot.data(), concentration_size);

    auto local_Jac = MathLib::createZeroedMatrix<
        typename ShapeMatricesTypeDisplacement::template MatrixType<
            concentration_size, concentration_size>>(
        local_Jac_data, concentration_size, concentration_size);
    typename ShapeMatricesTypePressure::NodalMatrixType mass =
        ShapeMatricesTypePressure::NodalMatrixType::Zero(concentration_size,
                                                         concentration_size);

    typename ShapeMatricesTypePressure::NodalMatrixType laplace =
        ShapeMatricesTypePressure::NodalMatrixType::Zero(concentration_size,
                                                         concentration_size);

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    int const n_integration_points = _integration_method.getNumberOfPoints();
    for (int ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);
        auto const& w = _ip_data[ip].integration_weight;

        auto const& N = _ip_data[ip].N;
        auto const& dNdx = _ip_data[ip].dNdx;

        double const K_over_mu =
            _process_data.intrinsic_permeability(t, x_position)[0] /
            _process_data.fluid_viscosity(t, x_position)[0];
        auto const alpha_b = _process_data.biot_coefficient(t, x_position)[0];

        laplace.noalias() += dNdx_p.transpose() * K_over_mu * dNdx_p * w;

        mass.noalias() +=
            N.transpose() * N * w *
            ((alpha_b - porosity) * (1.0 - alpha_b) / K_S + porosity * beta_p);

        auto const& b = _process_data.specific_body_force;
        local_rhs.noalias() += dNdx.transpose() * rho_fr * K_over_mu * b * w;

        const double dv_dt = 1.0 / dt;
        local_rhs.noalias() -= alpha_b * dv_dt * N * w;
    }
    local_Jac.noalias() = laplace + mass / dt;

    local_rhs.noalias() -= laplace * p + mass * p_dot;
}

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
void PhaseFieldAcidLocalAssembler<ShapeFunction, IntegrationMethod, GlobalDim>::
    assembleWithJacobianForConcentrationEquations(
        const double t, double const dt, const std::vector<double>& local_xdot,
        const double /*dxdot_dx*/, const double /*dx_dx*/,
        std::vector<double>& /*local_M_data*/,
        std::vector<double>& /*local_K_data*/,
        std::vector<double>& local_b_data, std::vector<double>& local_Jac_data,
        const LocalCoupledSolutions& local_coupled_solutions)
{
    auto local_rhs = MathLib::createZeroedVector<
        typename ShapeMatricesTypeDisplacement::template VectorType<
            concentration_size>>(local_b_data, concentration_size);

    ParameterLib::SpatialPosition pos;
    pos.setElementID(this->_element.getID());

    auto const c =
        Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
            concentration_size> const>(
            &local_coupled_solutions.local_coupled_xs[concentration_index],
            concentration_size);

    auto const phi =
        Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
            phasefield_size> const>(
            &local_coupled_solutions.local_coupled_xs[phasefield_index],
            phasefield_size);

    auto c_dot =
        Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
            concentration_size> const>(local_xdot.data(), concentration_size);

    auto local_Jac = MathLib::createZeroedMatrix<
        typename ShapeMatricesTypeDisplacement::template MatrixType<
            concentration_size, concentration_size>>(
        local_Jac_data, concentration_size, concentration_size);
    typename ShapeMatricesTypePressure::NodalMatrixType mass =
        ShapeMatricesTypePressure::NodalMatrixType::Zero(concentration_size,
                                                         concentration_size);

    typename ShapeMatricesTypePressure::NodalMatrixType laplace =
        ShapeMatricesTypePressure::NodalMatrixType::Zero(concentration_size,
                                                         concentration_size);

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    int const n_integration_points = _integration_method.getNumberOfPoints();
    for (int ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);
        auto const& w = _ip_data[ip].integration_weight;

        auto const& N = _ip_data[ip].N;
        auto const& dNdx = _ip_data[ip].dNdx;

        double const K_over_mu =
            _process_data.intrinsic_permeability(t, x_position)[0] /
            _process_data.fluid_viscosity(t, x_position)[0];
        auto const alpha_b = _process_data.biot_coefficient(t, x_position)[0];

        laplace.noalias() += dNdx_p.transpose() * K_over_mu * dNdx_p * w;

        mass.noalias() +=
            N.transpose() * N * w *
            ((alpha_b - porosity) * (1.0 - alpha_b) / K_S + porosity * beta_p);

        auto const& b = _process_data.specific_body_force;
        local_rhs.noalias() += dNdx.transpose() * rho_fr * K_over_mu * b * w;

        const double dv_dt = 1.0 / dt;
        local_rhs.noalias() -= alpha_b * dv_dt * N * w;
    }
    local_Jac.noalias() = laplace + mass / dt;

    local_rhs.noalias() -= laplace * p + mass * p_dot;
}

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
void PhaseFieldAcidLocalAssembler<ShapeFunction, IntegrationMethod, GlobalDim>::
    assembleWithJacobianForPhasefieldEquations(
        const double t, double const dt, const std::vector<double>& local_xdot,
        const double /*dxdot_dx*/, const double /*dx_dx*/,
        std::vector<double>& /*local_M_data*/,
        std::vector<double>& /*local_K_data*/,
        std::vector<double>& local_b_data, std::vector<double>& local_Jac_data,
        const LocalCoupledSolutions& local_coupled_solutions)
{
    auto local_rhs =
        MathLib::createZeroedVector<typename ShapeMatricesTypeDisplacement::
                                        template VectorType<phasefield_size>>(
            local_b_data, phasefield_size);

    ParameterLib::SpatialPosition pos;
    pos.setElementID(this->_element.getID());

    auto const phi =
        Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
            phasefield_size> const>(
            &local_coupled_solutions.local_coupled_xs[phasefield_index],
            phasefield_size);

    auto const c =
        Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
            concentration_size> const>(
            &local_coupled_solutions.local_coupled_xs[concentration_index],
            concentration_size);

    auto phi_dot =
        Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
            phasefield_size> const>(local_xdot.data(), phasefield_size);

    auto local_Jac = MathLib::createZeroedMatrix<
        typename ShapeMatricesTypeDisplacement::template MatrixType<
            phasefield_size, phasefield_size>>(
        local_Jac_data, phasefield_size, phasefield_size);

    typename ShapeMatricesTypePressure::NodalMatrixType mass =
        ShapeMatricesTypePressure::NodalMatrixType::Zero(phasefield_size,
                                                         phasefield_size);

    typename ShapeMatricesTypePressure::NodalMatrixType laplace =
        ShapeMatricesTypePressure::NodalMatrixType::Zero(phasefield_size,
                                                         phasefield_size);

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    int const n_integration_points = _integration_method.getNumberOfPoints();
    for (int ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);
        auto const& w = _ip_data[ip].integration_weight;

        auto const& N = _ip_data[ip].N;
        auto const& dNdx = _ip_data[ip].dNdx;

        double const K_over_mu =
            _process_data.intrinsic_permeability(t, x_position)[0] /
            _process_data.fluid_viscosity(t, x_position)[0];
        auto const alpha_b = _process_data.biot_coefficient(t, x_position)[0];

        laplace.noalias() += dNdx.transpose() * K_over_mu * dNdx * w;

        mass.noalias() +=
            N.transpose() * N * w *
            ((alpha_b - porosity) * (1.0 - alpha_b) / K_S + porosity * beta_p);

        auto const& b = _process_data.specific_body_force;
        local_rhs.noalias() += dNdx.transpose() * rho_fr * K_over_mu * b * w;

        const double dv_dt = 1.0 / dt;
        local_rhs.noalias() -= alpha_b * dv_dt * N * w;
    }
    local_Jac.noalias() = laplace + mass / dt;

    local_rhs.noalias() -= laplace * p + mass * p_dot;
}

}  // namespace PhaseFieldAcid
}  // namespace ProcessLib
