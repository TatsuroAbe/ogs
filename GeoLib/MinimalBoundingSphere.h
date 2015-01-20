/**
 * \file   Calculation of a minimum bounding sphere for a vector of points
 * \author Karsten Rink
 * \date   2014-07-11
 * \brief  Definition of the MinimalBoundingSphere class.
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 */

#ifndef MINIMALBOUNDINGSPHERE_H_
#define MINIMALBOUNDINGSPHERE_H_

#include <vector>

#include "MathLib/Vector3.h"
#include "Point.h"

namespace GeoLib
{

/**
 * Calculated center and radius of a minimal bounding sphere for a given number of geometric points.
 */
class MinimalBoundingSphere
{
public:
    /// Copy constructor
    MinimalBoundingSphere(MinimalBoundingSphere const&) = default;
    /// Point-Sphere
    MinimalBoundingSphere(GeoLib::Point const& p, double radius = std::numeric_limits<double>::epsilon());
    /// Bounding sphere using two points
    MinimalBoundingSphere(GeoLib::Point const& p, GeoLib::Point const& q);
    /// Bounding sphere using three points
    MinimalBoundingSphere(GeoLib::Point const& p, GeoLib::Point const& q, GeoLib::Point const& r);
    /// Bounding sphere using four points
    MinimalBoundingSphere(GeoLib::Point const& p, GeoLib::Point const& q, GeoLib::Point const& r, GeoLib::Point const& s);
    /// Bounding sphere of n points
    MinimalBoundingSphere(std::vector<GeoLib::Point*> const& points);
    ~MinimalBoundingSphere() {}

    /// Returns the center point of the sphere
    GeoLib::Point getCenter() const { return GeoLib::Point(_center.getCoords()); }

    /// Returns the radius of the sphere
    double getRadius() const {return _radius; }

    /// Returns the squared euclidean distance of a point from the sphere (for points within the sphere distance is negative)
    double pointDistanceSquared(GeoLib::Point const& pnt) const;

    /// Creates n_points random points located on the surface of the bounding sphere (useful for visualisation)
    std::vector<GeoLib::Point*>* getRandomSpherePoints(std::size_t n_points) const;

private:
    /// Constructor using no points
    MinimalBoundingSphere();

    /**
     * Recursive method for calculating a minimal bounding sphere for an arbitrary number of points.
     * Note: This method will change the order of elements in the vector sphere_points.
     * \param sphere_points The vector of points for which the smallest enclosing sphere is calculated
     * \param start_idx Start index of the vector subrange analysed in the current recursive step
     * \param length Length of the vector subrange analysed in the current recursive step
     * \param n_boundary_points Number of found boundary points in the current recursive step
     *
     * Algorithm based the following two papers:
     *   Emo Welzl: Smallest enclosing disks (balls and ellipsoids). New Results and New Trends in Computer Science, pp. 359--370, 1991
     *   Bernd Gaertner: Fast and Robust Smallest Enclosing Balls. ESA99, pp. 325--338, 1999.
     * Code based on "Smallest Enclosing Spheres" implementation by Nicolas Capens on flipcode's Developer Toolbox (www.flipcode.com)
     */
    static MinimalBoundingSphere recurseCalculation(std::vector<GeoLib::Point*> sphere_points, std::size_t start_idx, std::size_t length, std::size_t n_boundary_points);

    double _radius;
    MathLib::Vector3 _center;
};

} // namespace

#endif /* MINIMALBOUNDINGSPHERE_H_ */
