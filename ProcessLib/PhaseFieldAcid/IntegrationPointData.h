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

namespace ProcessLib
{
namespace PhaseFieldAcid
{
template <typename ShapeMatrixType>
struct IntegrationPointData final
{
    typename ShapeMatrixType::NodalRowVectorType const N;
    typename ShapeMatrixType::GlobalDimNodalMatrixType const dNdx;
    double const integration_weight;
    typename ShapeMatrixType::NodalMatrixType const mass_operator;

    double dummy;
    double dummy_prev;

    void pushBackState() { dummy_prev = dummy; }

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};
}  // namespace PhaseFieldAcid
}  // namespace ProcessLib
