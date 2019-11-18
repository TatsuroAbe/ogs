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

#include "LocalAssemblerInterface.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ParameterLib/SpatialPosition.h"
#include "PhaseFieldAcidProcessData.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"

namespace ProcessLib
{
namespace PhaseFieldAcid
{
template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
struct IntegrationPointData final
{
    {
        IntegrationPointData(NodalRowVectorType const& N_,
                             GlobalDimNodalMatrixType const& dNdx_,
                             double const& integration_weight_,
                             NodalMatrixType const mass_operator_)
            : N(N_),
              dNdx(dNdx_),
              integration_weight(integration_weight_),
              mass_operator(mass_operator_)
        {
        }

        NodalRowVectorType const N;
        GlobalDimNodalMatrixType const dNdx;
        double const integration_weight;
        NodalMatrixType const mass_operator;

        double dummy;
        double dummy_prev;

        void pushBackState() { dummy_prev = dummy; }

        EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

        static constexpr int kelvin_vector_size =
            MathLib::KelvinVector::KelvinVectorDimensions<
                DisplacementDim>::value;
    };

    /// Used for the extrapolation of the integration point values. It is
    /// ordered (and stored) by integration points.
    template <typename ShapeMatrixType>
    struct SecondaryData
    {
        std::vector<ShapeMatrixType, Eigen::aligned_allocator<ShapeMatrixType>>
            N;
    };

    template <typename ShapeFunction, typename IntegrationMethod,
              unsigned GlobalDim>
    class PhaseFieldAcidLocalAssembler
        : public PhaseFieldAcidLocalAssemblerInterface
    {
    public:
        using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction>;

        using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;

        PhaseFieldAcidLocalAssembler(PhaseFieldAcidLocalAssembler const&) =
            delete;
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
                    shape_matrices[ip].integralMeasure *
                    shape_matrices[ip].detJ;

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

        void assembleWithJacobianForStaggeredScheme(
            double const t, double const dt,
            std::vector<double> const& local_xdot, const double dxdot_dx,
            const double dx_dx, int const process_id,
            std::vector<double>& local_M_data,
            std::vector<double>& local_K_data,
            std::vector<double>& local_b_data,
            std::vector<double>& local_Jac_data,
            LocalCoupledSolutions const& local_coupled_solutions) override;

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
        void assembleWithJacobianPhaseFiledEquations(
            double const t, double const dt,
            std::vector<double> const& local_xdot, const double dxdot_dx,
            const double dx_dx, std::vector<double>& local_M_data,
            std::vector<double>& local_K_data,
            std::vector<double>& local_b_data,
            std::vector<double>& local_Jac_data,
            LocalCoupledSolutions const& local_coupled_solutions);

        void assembleWithJacobianForConcentrationEquations(
            double const t, double const dt,
            std::vector<double> const& local_xdot, const double dxdot_dx,
            const double dx_dx, std::vector<double>& local_M_data,
            std::vector<double>& local_K_data,
            std::vector<double>& local_b_data,
            std::vector<double>& local_Jac_data,
            LocalCoupledSolutions const& local_coupled_solutions);

        PhaseFieldAcidProcessData& _process_data;

        std::vector<
            IntegrationPointData<ShapeFunction, ShapeMatricesType, GlobalDim>,
            Eigen::aligned_allocator<
                IntegrationPointData<ShapeFunction, ShapeMatricesType, Dim>>>
            _ip_data;

        IntegrationMethod _integration_method;
        MeshLib::Element const& _element;
        SecondaryData<typename ShapeMatrices::ShapeType> _secondary_data;
        bool const _is_axially_symmetric;

        static const int displacement_index = 0;
        static const int displacement_size = ShapeFunction::NPOINTS;
        static const int phasefield_index = ShapeFunction::NPOINTS;
        static const int phasefield_size = ShapeFunction::NPOINTS;
    };

}  // namespace PhaseFieldAcid
}  // namespace PhaseFieldAcid

#include "PhaseFieldAcidFEM-impl.h"
