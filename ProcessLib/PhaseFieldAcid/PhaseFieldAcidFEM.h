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
#include "NumLib/Function/Interpolation.h"
#include "ParameterLib/SpatialPosition.h"
#include "PhaseFieldAcidProcessData.h"
#include "ProcessLib/CoupledSolutionsForStaggeredScheme.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"
#include "SecondaryData.h"

namespace ProcessLib
{
namespace PhaseFieldAcid
{
template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
class PhaseFieldAcidLocalAssembler
    : public PhaseFieldAcidLocalAssemblerInterface
{
    static const int concentration_index = 0;
    static const int concentration_size = ShapeFunction::NPOINTS;
    static const int phasefield_index = ShapeFunction::NPOINTS;
    static const int phasefield_size = ShapeFunction::NPOINTS;

    using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction, GlobalDim>;

    using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;

    using NodalVectorType = typename ShapeMatricesType::NodalVectorType;

    using GlobalDimVectorType = typename ShapeMatricesType::GlobalDimVectorType;
    using GlobalDimMatrixType = typename ShapeMatricesType::GlobalDimMatrixType;

    using GlobalDimNodalMatrixType =
        typename ShapeMatricesType::GlobalDimNodalMatrixType;

    using LocalBlockMatrixType =
        typename ShapeMatricesType::template MatrixType<phasefield_size,
                                                        phasefield_size>;
    using LocalSegmentVectorType =
        typename ShapeMatricesType::template VectorType<phasefield_size>;

public:
    PhaseFieldAcidLocalAssembler(PhaseFieldAcidLocalAssembler const&) = delete;
    PhaseFieldAcidLocalAssembler(PhaseFieldAcidLocalAssembler&&) = delete;

    PhaseFieldAcidLocalAssembler(MeshLib::Element const& e,
                                 std::size_t const /*local_matrix_size*/,
                                 bool const is_axially_symmetric,
                                 unsigned const integration_order,
                                 PhaseFieldAcidProcessData& process_data,
                                 const int concentration_process_id,
                                 const int phasefield_process_id)
        : _process_data(process_data),
          _integration_method(integration_order),
          _element(e),
          _is_axially_symmetric(is_axially_symmetric),
          _concentration_process_id(concentration_process_id),
          _phasefield_process_id(phasefield_process_id)
    {
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        _ip_data.reserve(n_integration_points);
        _secondary_data.N.resize(n_integration_points);

        auto const shape_matrices =
            initShapeMatrices<ShapeFunction, ShapeMatricesType,
                              IntegrationMethod, GlobalDim>(
                e, is_axially_symmetric, _integration_method);
        _shape_matrices_nodes =
            initShapeMatricesInNodes<ShapeFunction, ShapeMatricesType,
                                     GlobalDim>(e, is_axially_symmetric);

        ParameterLib::SpatialPosition x_position;
        x_position.setElementID(_element.getID());

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            _ip_data.emplace_back(
                shape_matrices[ip].N, shape_matrices[ip].dNdx,
                _integration_method.getWeightedPoint(ip).getWeight() *
                    shape_matrices[ip].integralMeasure *
                    shape_matrices[ip].detJ);
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
            "PhaseFieldLocalAssembler: assembly with jacobian is not "
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
            case 1:
            {
                assemblePhaseFiledEquations(t, dt, local_M_data, local_K_data,
                                            local_b_data, coupled_xs);
                break;
            }
            case 0:
            {
                // For the equations with concentration
                assembleConcentrationEquations(t, dt, local_M_data,
                                               local_K_data, local_b_data,
                                               coupled_xs);
                break;
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

        double kappa_avg = 0;
        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            _ip_data[ip].pushBackState();

            kappa_avg += _ip_data[ip].kappa;
        }

        (*_process_data.element_kappa)[_element.getID()] =
            kappa_avg / n_integration_points;
        DBUG("%d kappa_avg = %g", _element.getID(),
             (*_process_data.element_kappa)[_element.getID()]);
    }

    Eigen::Map<const Eigen::RowVectorXd> getShapeMatrix(
        const unsigned integration_point) const override
    {
        auto const& N = _secondary_data.N[integration_point];

        // assumes N is stored contiguously in memory
        return Eigen::Map<const Eigen::RowVectorXd>(N.data(), N.size());
    }

private:
    void assemblePhaseFiledEquations(
        double const t,
        double const /*dt*/,
        std::vector<double>& local_M_data,
        std::vector<double>& local_K_data,
        std::vector<double>& local_b_data,
        LocalCoupledSolutions const& local_coupled_solutions)
    {
        auto ph = Eigen::Map<const NodalVectorType>(
            &local_coupled_solutions.local_coupled_xs[phasefield_index],
            phasefield_size);
        auto ph0 = Eigen::Map<const NodalVectorType>(
            &local_coupled_solutions.local_coupled_xs0[phasefield_index],
            phasefield_size);

        auto c = Eigen::Map<const NodalVectorType>(
            &local_coupled_solutions.local_coupled_xs[concentration_index],
            concentration_size);
        auto c0 = Eigen::Map<const NodalVectorType>(
            &local_coupled_solutions.local_coupled_xs0[concentration_index],
            concentration_size);

        auto local_M = MathLib::createZeroedMatrix<LocalBlockMatrixType>(
            local_M_data, phasefield_size, phasefield_size);
        auto local_K = MathLib::createZeroedMatrix<LocalBlockMatrixType>(
            local_K_data, phasefield_size, phasefield_size);
        auto local_b = MathLib::createZeroedVector<LocalSegmentVectorType>(
            local_b_data, phasefield_size);

        GlobalDimNodalMatrixType v_at_nodes =
            GlobalDimNodalMatrixType::Zero(GlobalDim, phasefield_size);

        for (unsigned i = 0; i < _element.getNumberOfNodes(); i++)
        {
            v_at_nodes.col(i) = _shape_matrices_nodes[i].dNdx * ph0;
        }

        ParameterLib::SpatialPosition x_position;
        x_position.setElementID(_element.getID());

        int const n_integration_points =
            _integration_method.getNumberOfPoints();

        double const D = _process_data.chemical_diffusivity(t, x_position)[0];
        double const tau = _process_data.tau(t, x_position)[0];
        double const epsilon = _process_data.epsilon(t, x_position)[0];
        double const alpha = _process_data.alpha(t, x_position)[0];
        double const rrc = _process_data.rrc(t, x_position)[0];

        for (int ip = 0; ip < n_integration_points; ip++)
        {
            x_position.setIntegrationPoint(ip);
            auto const& w = _ip_data[ip].integration_weight;
            auto const& N = _ip_data[ip].N;
            auto const& dNdx = _ip_data[ip].dNdx;

            auto& kappa = _ip_data[ip].kappa;

            double c_ip = N.dot(c);
            double ph_ip = N.dot(ph0);

            GlobalDimMatrixType const grad_v_ip = dNdx * v_at_nodes.transpose();
            GlobalDimVectorType const v_ip = dNdx * ph0;
            double const squared_norm_v_ip = v_ip.squaredNorm();

            GlobalDimVectorType psi_ip =
                GlobalDimVectorType::Zero(GlobalDim, 1);
            if (squared_norm_v_ip > 1e-15)
            {
                psi_ip = grad_v_ip * v_ip / squared_norm_v_ip;
            }

            double const da = rrc * epsilon / D;
            double const lambda = D / epsilon / epsilon /
                                  (alpha * (5.0 / 3.0 + sqrt(2) / da));

            local_M.noalias() += w * tau * N.transpose() * N;

            local_K.noalias() +=
                epsilon * epsilon * dNdx.transpose() * dNdx * w;

            // f(phi) part
            local_b.noalias() +=
                (1 - ph_ip * ph_ip) * (ph_ip - lambda * c_ip) * N * w;
            // "kappa" part
            local_b.noalias() +=
                (N * psi_ip.dot(v_ip) + v_ip.transpose() * dNdx) * epsilon *
                epsilon * w;
        }
    }

    void assembleConcentrationEquations(
        double const t,
        double const dt,
        std::vector<double>& local_M_data,
        std::vector<double>& local_K_data,
        std::vector<double>& local_b_data,
        LocalCoupledSolutions const& local_coupled_solutions)
    {
        auto c = Eigen::Map<const NodalVectorType>(
            &local_coupled_solutions.local_coupled_xs[concentration_index],
            concentration_size);

        auto ph = Eigen::Map<const NodalVectorType>(
            &local_coupled_solutions.local_coupled_xs[phasefield_index],
            phasefield_size);
        auto ph0 = Eigen::Map<const NodalVectorType>(
            &local_coupled_solutions.local_coupled_xs0[phasefield_index],
            phasefield_size);

        auto local_M = MathLib::createZeroedMatrix<LocalBlockMatrixType>(
            local_M_data, concentration_size, concentration_size);
        auto local_K = MathLib::createZeroedMatrix<LocalBlockMatrixType>(
            local_K_data, concentration_size, concentration_size);
        auto local_b = MathLib::createZeroedVector<LocalSegmentVectorType>(
            local_b_data, concentration_size);

        ParameterLib::SpatialPosition x_position;
        x_position.setElementID(_element.getID());

        GlobalDimNodalMatrixType v_at_nodes =
            GlobalDimNodalMatrixType::Zero(GlobalDim, phasefield_size);

        for (unsigned i = 0; i < _element.getNumberOfNodes(); i++)
        {
            v_at_nodes.col(i) = _shape_matrices_nodes[i].dNdx * ph;
        }

        int const n_integration_points =
            _integration_method.getNumberOfPoints();

        double const D = _process_data.chemical_diffusivity(t, x_position)[0];
        double const alpha = _process_data.alpha(t, x_position)[0];
        double const rrc = _process_data.rrc(t, x_position)[0];

        for (int ip = 0; ip < n_integration_points; ip++)
        {
            auto const& w = _ip_data[ip].integration_weight;
            auto const& N = _ip_data[ip].N;
            auto const& dNdx = _ip_data[ip].dNdx;

            double ph_ip = N.dot(ph);
            double ph0_ip = N.dot(ph0);

            double const ph_dot = (ph_ip - ph0_ip) / dt;
            double const grad_ph_norm = (dNdx * ph).norm();

            GlobalDimMatrixType const grad_v_ip = dNdx * v_at_nodes.transpose();
            double const div_grad_phi_ip = grad_v_ip.diagonal().sum();

            double source = 1.;
            if (grad_ph_norm > std::numeric_limits<double>::epsilon())
            {
                source += (D * div_grad_phi_ip - ph_dot) / rrc / grad_ph_norm;
            }

            local_M.noalias() += w * N.transpose() * N;

            local_K.noalias() += D * dNdx.transpose() * dNdx * w;

            local_b.noalias() += alpha * ph_dot * source * N * w;
        }
    }

private:
    PhaseFieldAcidProcessData& _process_data;

    std::vector<
        IntegrationPointData<ShapeMatricesType>,
        Eigen::aligned_allocator<IntegrationPointData<ShapeMatricesType>>>
        _ip_data;
    std::vector<
        typename ShapeMatricesType::ShapeMatrices,
        Eigen::aligned_allocator<typename ShapeMatricesType::ShapeMatrices>>
        _shape_matrices_nodes;
    IntegrationMethod _integration_method;
    MeshLib::Element const& _element;
    SecondaryData<typename ShapeMatrices::ShapeType> _secondary_data;
    bool const _is_axially_symmetric;
    const int _concentration_process_id;
    const int _phasefield_process_id;
};

}  // namespace PhaseFieldAcid
}  // namespace ProcessLib
