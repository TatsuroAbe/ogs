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
            initShapeMatrices<ShapeFunction, ShapeMatricesType,
                              IntegrationMethod, GlobalDim>(
                e, is_axially_symmetric, _integration_method);

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
            case PhaseFieldAcidProcessData::phasefield_process_id:
            {
                assemblePhaseFiledEquations(t, dt, local_M_data, local_K_data,
                                            local_b_data, coupled_xs);
                break;
            }
            case PhaseFieldAcidProcessData::concentration_process_id:
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
    void assemblePhaseFiledEquations(
        double const t,
        double const dt,
        std::vector<double>& local_M_data,
        std::vector<double>& local_K_data,
        std::vector<double>& local_b_data,
        LocalCoupledSolutions const& local_coupled_solutions)
    {
        auto local_ph = Eigen::Map<const NodalVectorType>(
            &local_coupled_solutions.local_coupled_xs[phasefield_index],
            phasefield_size);

        auto local_c = Eigen::Map<const NodalVectorType>(
            &local_coupled_solutions.local_coupled_xs[concentration_index],
            concentration_size);
        auto local_c0 = Eigen::Map<const NodalVectorType>(
            &local_coupled_solutions.local_coupled_xs0[concentration_index],
            concentration_size);

        auto local_M = MathLib::createZeroedMatrix<LocalBlockMatrixType>(
            local_M_data, phasefield_size, phasefield_size);
        auto local_K = MathLib::createZeroedMatrix<LocalBlockMatrixType>(
            local_K_data, phasefield_size, phasefield_size);
        auto local_b = MathLib::createZeroedVector<LocalSegmentVectorType>(
            local_b_data, phasefield_size);

        /*        auto local_M = MathLib::createZeroedMatrix<
                    typename ShapeMatricesType::template
           MatrixType<phasefield_size, phasefield_size>>( local_M_data,
           phasefield_size, phasefield_size);

                auto local_K = MathLib::createZeroedMatrix<
                    typename ShapeMatricesType::template
           MatrixType<phasefield_size, phasefield_size>>( local_K_data,
           phasefield_size, phasefield_size);

                auto local_b = MathLib::createZeroedVector<
                    typename ShapeMatricesType::template
           VectorType<phasefield_size>>( local_b_data, phasefield_size);
        */
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

        double const tau = _process_data.tau(t, x_position)[0];
        double const epsilon = _process_data.epsilon(t, x_position)[0];
        auto const& b = _process_data.specific_body_force;
        for (int ip = 0; ip < n_integration_points; ip++)
        {
            x_position.setIntegrationPoint(ip);
            auto const& w = _ip_data[ip].integration_weight;
            auto const& N = _ip_data[ip].N;
            auto const& dNdx = _ip_data[ip].dNdx;

            double c_ip = 0.0;
            double ph_ip = 0.0;

            NumLib::shapeFunctionInterpolate(local_c, N, c_ip);
            NumLib::shapeFunctionInterpolate(local_ph, N, ph_ip);

            auto const x_coord =
                interpolateXCoordinate<ShapeFunction, ShapeMatricesType>(
                    _element, N);

            local_M.noalias() += w * N.transpose() * N;

            local_K.noalias() +=
                epsilon * epsilon / tau * dNdx.transpose() * dNdx * w;

            local_b.noalias() += 0.0 * w * dNdx.transpose() * b;
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
        auto local_c = Eigen::Map<const NodalVectorType>(
            &local_coupled_solutions.local_coupled_xs[concentration_index],
            concentration_size);

        auto local_ph = Eigen::Map<const NodalVectorType>(
            &local_coupled_solutions.local_coupled_xs[phasefield_index],
            phasefield_size);
        auto local_ph0 = Eigen::Map<const NodalVectorType>(
            &local_coupled_solutions.local_coupled_xs0[phasefield_index],
            phasefield_size);

        auto local_M = MathLib::createZeroedMatrix<LocalBlockMatrixType>(
            local_M_data, concentration_size, concentration_size);
        auto local_K = MathLib::createZeroedMatrix<LocalBlockMatrixType>(
            local_K_data, concentration_size, concentration_size);
        auto local_b = MathLib::createZeroedVector<LocalSegmentVectorType>(
            local_b_data, concentration_size);
        /*
        auto local_M = MathLib::createZeroedMatrix<
            typename ShapeMatricesType::template MatrixType<
                concentration_size, concentration_size>>(
            local_M_data, concentration_size, concentration_size);

        auto local_K = MathLib::createZeroedMatrix<
            typename ShapeMatricesType::template MatrixType<
                concentration_size, concentration_size>>(
            local_K_data, concentration_size, concentration_size);

        auto local_b = MathLib::createZeroedVector<
            typename ShapeMatricesType::template VectorType<
                concentration_size>>(local_b_data, concentration_size);
*/
        ParameterLib::SpatialPosition x_position;
        x_position.setElementID(_element.getID());

        int const n_integration_points =
            _integration_method.getNumberOfPoints();

        double const D = _process_data.chemical_diffusivity(t, x_position)[0];
        double const alpha = _process_data.alpha(t, x_position)[0];
        double const rrc = _process_data.rrc(t, x_position)[0];

        for (int ip = 0; ip < n_integration_points; ip++)
        {
            x_position.setIntegrationPoint(ip);
            auto const& w = _ip_data[ip].integration_weight;
            auto const& N = _ip_data[ip].N;
            auto const& dNdx = _ip_data[ip].dNdx;

            double c_ip = 0.0;
            double ph_ip = 0.0;

            NumLib::shapeFunctionInterpolate(local_c, N, c_ip);
            NumLib::shapeFunctionInterpolate(local_ph, N, ph_ip);

            auto const x_coord =
                interpolateXCoordinate<ShapeFunction, ShapeMatricesType>(
                    _element, N);

            auto const& b = _process_data.specific_body_force;

            local_M.noalias() += w * N.transpose() * N;

            local_K.noalias() += D * dNdx.transpose() * dNdx * w;

            local_b.noalias() += 0.0 * w * dNdx.transpose() * b;
        }
    }

private:
    PhaseFieldAcidProcessData& _process_data;

    std::vector<
        IntegrationPointData<ShapeMatricesType>,
        Eigen::aligned_allocator<IntegrationPointData<ShapeMatricesType>>>
        _ip_data;

    IntegrationMethod _integration_method;
    MeshLib::Element const& _element;
    SecondaryData<typename ShapeMatrices::ShapeType> _secondary_data;
    bool const _is_axially_symmetric;
};

}  // namespace PhaseFieldAcid
}  // namespace ProcessLib
