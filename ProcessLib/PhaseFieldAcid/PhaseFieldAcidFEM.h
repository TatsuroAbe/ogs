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

#include <memory>
#include <vector>

#include "IntegrationPointData.h"
#include "LocalAssemblerInterface.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ParameterLib/SpatialPosition.h"
#include "PhaseFieldAcidProcessData.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"
#include "SecondaryData.h"

namespace ProcessLib
{
namespace PhaseFieldAcid
{
template <typename ShapeFunction, typename IntegrationMethod, int GlobalDim>
class PhaseFieldAcidLocalAssembler
    : public PhaseFieldAcidLocalAssemblerInterface
{
public:
    using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction, GlobalDim>;

    using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;

    PhaseFieldAcidLocalAssembler(PhaseFieldAcidLocalAssembler const&) = delete;
    PhaseFieldAcidLocalAssembler(PhaseFieldAcidLocalAssembler&&) = delete;

    PhaseFieldAcidLocalAssembler(MeshLib::Element const& e,
                                 std::size_t const /*local_matrix_size*/,
                                 bool const is_axially_symmetric,
                                 unsigned const integration_order,
                                 PhaseFieldAcidProcessData& process_data)
        : _process_data(process_data),
          _integration_method(integration_order),
          _element(e),
          _is_axially_symmetric(is_axially_symmetric)
    {
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        _ip_data.reserve(n_integration_points);
        _secondary_data.N.resize(n_integration_points);

        auto const shape_matrices =
            initShapeMatrices<ShapeFunction, IntegrationMethod, GlobalDim>(
                e, is_axially_symmetric, _integration_method);

        ParameterLib::SpatialPosition x_position;
        x_position.setElementID(_element.getID());

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            auto& ip_data = _ip_data[ip];
            ip_data.integration_weight =
                _integration_method.getWeightedPoint(ip).getWeight() *
                shape_matrices[ip].integralMeasure * shape_matrices[ip].detJ;

            ip_data.dummy = 0.0;
            ip_data.dummy_prev = 0.0;

            ip_data.N = shape_matrices[ip].N;
            ip_data.dNdx = shape_matrices[ip].dNdx;

            _secondary_data.N[ip] = shape_matrices[ip].N;
        }
    }

    void assemble(double const /*t*/, double const /*dt*/,
                  std::vector<double> const& /*local_x*/,
                  std::vector<double>& /*local_M_data*/,
                  std::vector<double>& /*local_K_data*/,
                  std::vector<double>& /*local_rhs_data*/) override
    {
        OGS_FATAL(
            "PhaseFieldLocalAssembler: assembly without jacobian is not "
            "implemented.");
    }

    void assembleForStaggeredScheme(
        double const t, double const dt, int const process_id,
        std::vector<double>& local_M_data, std::vector<double>& local_K_data,
        std::vector<double>& local_b_data,
        LocalCoupledSolutions const& coupled_xs) override
    {
        switch (process_id)
        {
            case _process_data.phasefield_process_id:
            {
                assemblePhaseFiledEquations(t, dt, local_M_data, local_K_data,
                                            local_b_data, coupled_xs);
            }
            case _process_data.concentration_process_id:
            {
                // For the equations with concentration
                assembleConcentrationEquations(t, dt, local_M_data,
                                               local_K_data, local_b_data,
                                               coupled_xs);
            }
        }
    }

    void initializeConcrete() override
    {
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            _ip_data[ip].pushBackState();
        }
    }

    void postTimestepConcrete(std::vector<double> const& /*local_x*/,
                              double const /*t*/,
                              double const /*dt*/) override
    {
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            _ip_data[ip].pushBackState();
        }
    }

    Eigen::Map<const Eigen::RowVectorXd> getShapeMatrix(
        const unsigned integration_point) const override
    {
        auto const& N = _secondary_data.N[integration_point];

        // assumes N is stored contiguously in memory
        return Eigen::Map<const Eigen::RowVectorXd>(N.data(), N.size());
    }

private:
    void assemblePhaseFiledEquations(double const t,
                                     double const dt,
                                     std::vector<double>& local_M_data,
                                     std::vector<double>& local_K_data,
                                     std::vector<double>& local_b_data,
                                     LocalCoupledSolutions const& coupled_xs)
    {
        /*
        auto const& local_ph =
            local_coupled_solutions
                .local_coupled_xs[_process_data.phasefield_process_id];
        auto const& local_c =
            local_coupled_solutions
                .local_coupled_xs[_process_data.concentration_process_id];
        assert(local_c.size() == concentration_size);
        assert(local_ph.size() == phasefield_size);

        auto phi = Eigen::Map<typename ShapeMatricesType::template VectorType<
            phasefield_size> const>(local_ph.data(), phasefield_size);

        auto c = Eigen::Map<typename ShapeMatricesType::template VectorType<
            concentration_size> const>(local_c.data(), concentration_size);

        auto c_dot = Eigen::Map<typename ShapeMatricesType::template VectorType<
            concentration_size> const>(local_xdot.data(), concentration_size);

        auto local_Jac = MathLib::createZeroedMatrix<
            typename ShapeMatricesType::template MatrixType<
                concentration_size, concentration_size>>(
            local_Jac_data, concentration_size, concentration_size);

        auto local_rhs = MathLib::createZeroedVector<
            typename ShapeMatricesType::template VectorType<
                concentration_size>>(local_b_data, concentration_size);

        typename ShapeMatricesType::NodalMatrixType mass =
            ShapeMatricesType::NodalMatrixType::Zero(concentration_size,
                                                     concentration_size);

        typename ShapeMatricesType::NodalMatrixType laplace =
            ShapeMatricesType::NodalMatrixType::Zero(concentration_size,
                                                     concentration_size);

        ParameterLib::SpatialPosition x_position;
        x_position.setElementID(_element.getID());

        int const n_integration_points =
            _integration_method.getNumberOfPoints();

        for (int ip = 0; ip < n_integration_points; ip++)
        {
            x_position.setIntegrationPoint(ip);
            auto const& w = _ip_data[ip].integration_weight;
            auto const& N = _ip_data[ip].N;
            auto const& dNdx = _ip_data[ip].dNdx;
            double const phi_ip = N.dot(phi);

            auto const x_coord =
                interpolateXCoordinate<ShapeFunction, ShapeMatricesType>(
                    _element, N);

            auto const& b = _process_data.specific_body_force;

            local_rhs.noalias() += (1.0) * N * w;
            mass.noalias() += (1.0) * N.transpose() * N * w;

            laplace.noalias() += (dNdx.transpose() * dNdx) * w;
        }
        local_Jac.noalias() = laplace + mass / dt;

        local_rhs.noalias() -= laplace * c + mass * c_dot;
        */
    }

    void assembleConcentrationEquations(
        double const t,
        double const dt,
        std::vector<double>& local_M_data,
        std::vector<double>& local_K_data,
        std::vector<double>& local_b_data,
        LocalCoupledSolutions const& coupled_xs)
    {
        /*
        auto const& local_ph =
            local_coupled_solutions
                .local_coupled_xs[_process_data.phasefield_process_id];
        auto const& local_c =
            local_coupled_solutions
                .local_coupled_xs[_process_data.concentration_process_id];
        assert(local_c.size() == concentration_size);
        assert(local_ph.size() == phasefield_size);

        auto phi = Eigen::Map<typename ShapeMatricesType::template VectorType<
            phasefield_size> const>(local_ph.data(), phasefield_size);

        auto c = Eigen::Map<typename ShapeMatricesType::template VectorType<
            concentration_size> const>(local_c.data(), concentration_size);

        auto phi_dot =
            Eigen::Map<typename ShapeMatricesType::template VectorType<
                phasefield_size> const>(local_xdot.data(), phasefield_size);

        auto local_Jac = MathLib::createZeroedMatrix<
            typename ShapeMatricesType::template MatrixType<
                phasefield_size, concentration_size>>(
            local_Jac_data, phasefield_size, phasefield_size);

        auto local_rhs = MathLib::createZeroedVector<
            typename ShapeMatricesType::template VectorType<phasefield_size>>(
            local_b_data, phasefield_size);

        typename ShapeMatricesType::NodalMatrixType mass =
            ShapeMatricesType::NodalMatrixType::Zero(phasefield_size,
                                                     phasefield_size);

        typename ShapeMatricesType::NodalMatrixType laplace =
            ShapeMatricesType::NodalMatrixType::Zero(phasefield_size,
                                                     phasefield_size);

        ParameterLib::SpatialPosition x_position;
        x_position.setElementID(_element.getID());

        int const n_integration_points =
            _integration_method.getNumberOfPoints();

        for (int ip = 0; ip < n_integration_points; ip++)
        {
            x_position.setIntegrationPoint(ip);
            auto const& w = _ip_data[ip].integration_weight;
            auto const& N = _ip_data[ip].N;
            auto const& dNdx = _ip_data[ip].dNdx;
            double const phi_ip = N.dot(phi);

            auto const x_coord =
                interpolateXCoordinate<ShapeFunction, ShapeMatricesType>(
                    _element, N);

            auto const& b = _process_data.specific_body_force;

            local_rhs.noalias() += (1.0) * N * w;
            mass.noalias() += (1.0) * N.transpose() * N * w;

            laplace.noalias() += (dNdx.transpose() * dNdx) * w;
        }
        local_Jac.noalias() = laplace + mass / dt;

        local_rhs.noalias() -= laplace * phi + mass * phi_dot;
        */
    }

    PhaseFieldAcidProcessData& _process_data;

    std::vector<
        IntegrationPointData<ShapeMatricesType>,
        Eigen::aligned_allocator<IntegrationPointData<ShapeMatricesType>>>
        _ip_data;

    IntegrationMethod _integration_method;
    MeshLib::Element const& _element;
    SecondaryData<typename ShapeMatrices::ShapeType> _secondary_data;
    bool const _is_axially_symmetric;

    static const int concentration_index = 0;
    static const int concentration_size = ShapeFunction::NPOINTS;
    static const int phasefield_index = ShapeFunction::NPOINTS;
    static const int phasefield_size = ShapeFunction::NPOINTS;
};

}  // namespace PhaseFieldAcid
}  // namespace ProcessLib
