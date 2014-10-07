/**
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_GROUNDWATERFLOW_FEM_H_
#define PROCESS_LIB_GROUNDWATERFLOW_FEM_H_

namespace ProcessLib
{

namespace GroundwaterFlow
{

template <typename GlobalMatrix, typename GlobalVector>
class LocalAssemblerDataBase
{
public:
    virtual ~LocalAssemblerDataBase() = default;

    virtual void init(MeshLib::Element const& e,
            double const hydraulic_conductivity) = 0;

    virtual void assemble(std::size_t const rows, std::size_t const columns) = 0;

    virtual void addToGlobal(GlobalMatrix& A, GlobalVector& rhs,
            AssemblerLib::LocalToGlobalIndexMap::RowColumnIndices const&) const = 0;
};

template <typename ShapeFunction_,
         typename IntegrationMethod_,
         unsigned IntegrationOrder_,
         typename GlobalMatrix,
         typename GlobalVector>
class LocalAssemblerData : public LocalAssemblerDataBase<GlobalMatrix, GlobalVector>
{
public:
    using ShapeFunction = ShapeFunction_;
    using NodalMatrixType = typename ShapeMatrixPolicyType<ShapeFunction>::NodalMatrixType;
    using NodalVectorType = typename ShapeMatrixPolicyType<ShapeFunction>::NodalVectorType;

    using ShapeMatrices = typename ShapeMatrixPolicyType<ShapeFunction>::ShapeMatrices;


    static unsigned constexpr integration_order = IntegrationOrder_;
    static unsigned constexpr n_integration_points =
        MathLib::pow(integration_order, ShapeFunction::DIM);
    using IntegrationMethod = IntegrationMethod_;

    /// The hydraulic_conductivity factor is directly integrated into the local
    /// element matrix.
    void
    init(MeshLib::Element const& e, double const hydraulic_conductivity)
    {
    }

    void assemble(std::size_t const rows, std::size_t const columns)
    {
    }

    void addToGlobal(GlobalMatrix& A, GlobalVector& rhs,
            AssemblerLib::LocalToGlobalIndexMap::RowColumnIndices const& indices) const
    {
        A.add(indices, *localA);
        rhs.add(indices.rows, *localRhs);
    }

private:

    std::unique_ptr<NodalMatrixType> localA;
    std::unique_ptr<NodalVectorType> localRhs;
};

}   // namespace GroundwaterFlow
}   // namespace ProcessLib

#endif  // PROCESS_LIB_GROUNDWATERFLOW_FEM_H_
