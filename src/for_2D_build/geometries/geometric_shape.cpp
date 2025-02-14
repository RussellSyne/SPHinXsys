#include "geometric_shape.h"

namespace SPH
{
//=================================================================================================//
GeometricShapeBox::GeometricShapeBox(const Vec2d &halfsize, const std::string &shape_name)
    : GeometricBox(halfsize), Shape(shape_name){}
//=================================================================================================//
bool GeometricShapeBox::checkContain(const Vec2d &probe_point, bool BOUNDARY_INCLUDED)
{
    return GeometricBox::checkContain(probe_point);
}
//=================================================================================================//
Vec2d GeometricShapeBox::findClosestPoint(const Vec2d &probe_point)
{
    return GeometricBox::findClosestPoint(probe_point);
}
//=================================================================================================//
BoundingBox GeometricShapeBox::findBounds()
{
    return BoundingBox(-halfsize_, halfsize_);
}
//=================================================================================================//
GeometricShapeBall::GeometricShapeBall(const Vec2d &center, Real radius,
                                       const std::string &shape_name)
    : Shape(shape_name), center_(center), radius_(radius) {}
//=================================================================================================//
bool GeometricShapeBall::checkContain(const Vec2d &probe_point, bool BOUNDARY_INCLUDED)
{
    return (probe_point - center_).norm() < radius_;
}
//=================================================================================================//
Vec2d GeometricShapeBall::findClosestPoint(const Vec2d &probe_point)
{
    Vec2d displacement = probe_point - center_;
    Real distance = displacement.norm();
    Real cosine = (SGN(displacement[0]) * (ABS(displacement[0])) + TinyReal) / (distance + TinyReal);
    Real sine = displacement[1] / (distance + TinyReal);
    return probe_point + (radius_ - distance) * Vec2d(cosine, sine);
}
//=================================================================================================//
BoundingBox GeometricShapeBall::findBounds()
{
    Vec2d shift = Vec2d(radius_, radius_);
    return BoundingBox(center_ - shift, center_ + shift);
}
//=================================================================================================//
} // namespace SPH