#include "sphinxsys.h" //SPHinXsys Library.
using namespace SPH;   // Namespace cite here.
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 1.0; /**< Tank length. */
Real DH = 1.0; /**< Tank height. */
Real droplet_radius = 0.2;
Real particle_spacing_ref = droplet_radius / 10.0; /**< Initial reference particle spacing. */
Real BW = particle_spacing_ref * 4;                /**< Extending width for BCs. */
//----------------------------------------------------------------------
//	Material parameters.
//----------------------------------------------------------------------
Real rho0_f = 1.0;       /**< Reference density of water. */
Real rho0_a = 1.0;       /**< Reference density of air. */
Real U_ref = 1.0;        /**< Characteristic velocity. */
Real c_f = 10.0 * U_ref; /**< Reference sound speed. */
Real mu_f = 5.0e-2;      /**< Water viscosity. */
Real mu_a = 5.0e-2;      /**< Air viscosity. */
Real surface_tension_coef = 1.0;
Real transport_velocity_coef = 0.2;
Real U0 = 10;
Real Re_water = rho0_f * U_ref * 2 * droplet_radius / mu_f;
Real Re_air = rho0_a * U_ref * 2 * droplet_radius / mu_a;
//
// // Real rho0_f = 1.0;       /**< Reference density of water. */
// // Real rho0_a = 1.0e-3;       /**< Reference density of air. */
// Real rho0_f = 1000.0;       /**< Reference density of water. */
// Real rho0_a = 1.0;       /**< Reference density of air. */
// Real U_ref = 2.0;        /**< Characteristic velocity. */
// Real c_f = 10.0 * U_ref; /**< Reference sound speed. */
// // Real mu_f = 5.0e-2;        /**< Water viscosity. */
// // Real mu_a = 5.0e-4;        /**< Air viscosity. */
// Real mu_f = 5.0e-3;        /**< Water viscosity. */
// Real mu_a = 5.0e-5;        /**< Air viscosity. */
// Real surface_tension_coef = 1;
// Real transport_velocity_coef = 0.2;
// Real U0 = 0.0;
// Real Re_water = rho0_f *U_ref*2*droplet_radius/mu_f;
// Real Re_air = rho0_a *U_ref*2*droplet_radius/mu_a;
//----------------------------------------------------------------------
//	Geometric shapes used in this case.
//----------------------------------------------------------------------
Vec2d outer_wall_halfsize = Vec2d(0.5 * DL + BW, 0.5 * DH + BW);
Vec2d outer_wall_translation = Vec2d(0, 0);
Vec2d inner_wall_halfsize = Vec2d(0.5 * DL, 0.5 * DH);
Vec2d inner_wall_translation = Vec2d(0, 0);
Vecd air_halfsize = inner_wall_halfsize;       // local center at origin
Vecd air_translation = inner_wall_translation; // translation to global coordinates

Vecd droplet_center(0.0, 0.0);
Vecd droplet_halfsize = Vec2d(droplet_radius, droplet_radius); // local center at origin
Vecd droplet_translation = droplet_center;                     // translation to global coordinates
//----------------------------------------------------------------------
// Water body shape definition.
//----------------------------------------------------------------------
class WaterBlock : public MultiPolygonShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addACircle(droplet_center, droplet_radius, 100, ShapeBooleanOps::add);
    }
};
//----------------------------------------------------------------------
// Air body shape definition.
//----------------------------------------------------------------------/**
class AirBlock : public MultiPolygonShape
{
  public:
    explicit AirBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addABox(Transform(air_translation), air_halfsize, ShapeBooleanOps::add);
        multi_polygon_.addACircle(droplet_center, droplet_radius, 100, ShapeBooleanOps::sub);
    }
};
//----------------------------------------------------------------------
//	Complex shape for wall boundary, note that no partial overlap is allowed
//	for the shapes in a complex shape.
//----------------------------------------------------------------------
class WallBoundary : public ComplexShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TransformShape<GeometricShapeBox>>(Transform(outer_wall_translation), outer_wall_halfsize);
        subtract<TransformShape<GeometricShapeBox>>(Transform(inner_wall_translation), inner_wall_halfsize);
    }
};

using MultiPhaseLinearGradientCorrectionMatrixComplex = ComplexInteraction<LinearGradientCorrectionMatrix<Inner<>, Contact<>, Contact<>>>;
/**
 * application dependent initial velocity
 */
class InitialVelocity
    : public fluid_dynamics::FluidInitialCondition
{
  public:
    InitialVelocity(SPHBody &sph_body)
        : fluid_dynamics::FluidInitialCondition(sph_body) {};

    void update(size_t index_i, Real dt)
    {
        Real r0 = 0.05;
        Real r = pos_[index_i].norm();
        Real x = pos_[index_i][0];
        Real y = pos_[index_i][1];
        /** initial velocity profile */
        vel_[index_i][0] = U0 * x / r0 * (1 - y * y / (r0 * r)) * exp(-r / r0);
        vel_[index_i][1] = -U0 * y / r0 * (1 - x * x / (r0 * r)) * exp(-r / r0);
    }
};

class SetPositionCenterZero
    : public fluid_dynamics::FluidInitialCondition
{
  public:
    SetPositionCenterZero(SPHBody &sph_body)
        : fluid_dynamics::FluidInitialCondition(sph_body),
          pos_center_(particles_->registerSingularVariable<Vecd>("PositionCenter")->ValueAddress()) {};

    void update(size_t index_i, Real dt)
    {
        *pos_center_ = Vecd::Zero();
    }

  protected:
    Vecd *pos_center_;
};

class CalPositionCenter
    : public fluid_dynamics::FluidInitialCondition
{
  public:
    CalPositionCenter(SPHBody &sph_body)
        : fluid_dynamics::FluidInitialCondition(sph_body),
          pos_center_(particles_->getSingularVariableByName<Vecd>("PositionCenter")->ValueAddress()) {};

    void update(size_t index_i, Real dt)
    {
        *pos_center_ += pos_[index_i] / particles_->TotalRealParticles();
    }

  protected:
    Vecd *pos_center_;
};

class OutputPositionCenter
    : public fluid_dynamics::FluidInitialCondition
{
  public:
    OutputPositionCenter(SPHBody &sph_body)
        : fluid_dynamics::FluidInitialCondition(sph_body),
          pos_center_(particles_->getSingularVariableByName<Vecd>("PositionCenter")->ValueAddress()) {};

    void update(size_t index_i, Real dt)
    {
        if (abs((*pos_center_).norm()) > 1.0e-4)
            std::cout << "*pos_center_:" << *pos_center_ << std::endl;
    }

  protected:
    Vecd *pos_center_;
};

class GyrationTensorS1 : public LocalDynamicsReduce<ReduceSum<Vecd>>
{

  public:
    explicit GyrationTensorS1(SPHBody &sph_body) : LocalDynamicsReduce<ReduceSum<Vecd>>(sph_body), pos_(particles_->getVariableDataByName<Vecd>("Position")), pos_center_(particles_->getSingularVariableByName<Vecd>("PositionCenter")->ValueAddress())
    {
        quantity_name_ = "GyrationTensorS1";
    }
    virtual ~GyrationTensorS1() {};

    Vecd reduce(size_t index_i, Real dt = 0.0)
    {
        Vecd pos_diff = pos_[index_i] - *pos_center_;
        Matd ans = pos_diff * pos_diff.transpose() / particles_->TotalRealParticles();
        return ans * Vecd(1, 0);
    }

  protected:
    Vecd *pos_, *pos_center_;
};

class GyrationTensorS2 : public LocalDynamicsReduce<ReduceSum<Vecd>>
{

  public:
    explicit GyrationTensorS2(SPHBody &sph_body) : LocalDynamicsReduce<ReduceSum<Vecd>>(sph_body), pos_(particles_->getVariableDataByName<Vecd>("Position")), pos_center_(particles_->getSingularVariableByName<Vecd>("PositionCenter")->ValueAddress())
    {
        quantity_name_ = "GyrationTensorS2";
    }
    virtual ~GyrationTensorS2() {};

    Vecd reduce(size_t index_i, Real dt = 0.0)
    {
        Vecd pos_diff = pos_[index_i] - *pos_center_;
        Matd ans = pos_diff * pos_diff.transpose() / particles_->TotalRealParticles();
        return ans * Vecd(0, 1);
    }

  protected:
    Vecd *pos_, *pos_center_;
};

class PositionCenter : public LocalDynamicsReduce<ReduceSum<Vecd>>
{

  public:
    explicit PositionCenter(SPHBody &sph_body) : LocalDynamicsReduce<ReduceSum<Vecd>>(sph_body), pos_(particles_->getVariableDataByName<Vecd>("Position")), pos0_(particles_->registerStateVariableFrom<Vecd>("InitialPosition", "Position"))
    {
        quantity_name_ = "PositionCenter";
    }
    virtual ~PositionCenter() {};

    Vecd reduce(size_t index_i, Real dt = 0.0)
    {
        Vecd ans = Vecd::Zero();
        if (pos0_[index_i][0] > 0 && pos0_[index_i][1] > 0)
        {
            ans = 4 * pos_[index_i] / particles_->TotalRealParticles();
        }
        return ans;
    }

  protected:
    Vecd *pos_, *pos0_;
};

class PositionCenterV2 : public LocalDynamicsReduce<ReduceSum<Vecd>>
{

  public:
    explicit PositionCenterV2(SPHBody &sph_body) : LocalDynamicsReduce<ReduceSum<Vecd>>(sph_body), pos_(particles_->getVariableDataByName<Vecd>("Position")), pos0_(particles_->registerStateVariableFrom<Vecd>("InitialPosition", "Position"))
    {
        quantity_name_ = "PositionCenterV2";
    }
    virtual ~PositionCenterV2() {};

    Vecd reduce(size_t index_i, Real dt = 0.0)
    {
        Vecd ans = Vecd::Zero();
        // if (pos0_[index_i][0] > 0 && pos0_[index_i][1] > 0)
        if (pos_[index_i][0] > 0 && pos_[index_i][1] > 0)
        {
            ans = 4 * pos_[index_i] / particles_->TotalRealParticles();
        }
        return ans;
    }

  protected:
    Vecd *pos_, *pos0_;
};

class VelocityCenter : public LocalDynamicsReduce<ReduceSum<Vecd>>
{

  public:
    explicit VelocityCenter(SPHBody &sph_body) : LocalDynamicsReduce<ReduceSum<Vecd>>(sph_body), vel_(particles_->getVariableDataByName<Vecd>("Velocity")), pos0_(particles_->registerStateVariableFrom<Vecd>("InitialPosition", "Position"))
    {
        quantity_name_ = "VelocityCenter";
    }
    virtual ~VelocityCenter() {};

    Vecd reduce(size_t index_i, Real dt = 0.0)
    {
        Vecd ans = Vecd::Zero();
        if (pos0_[index_i][0] > 0 && pos0_[index_i][1] > 0)
        {
            ans = 4 * vel_[index_i] / particles_->TotalRealParticles();
        }
        return ans;
    }

  protected:
    Vecd *vel_, *pos0_;
};

class DeltaFunction : public LocalDynamicsReduce<ReduceSum<Real>>
{

  public:
    explicit DeltaFunction(SPHBody &sph_body) : LocalDynamicsReduce<ReduceSum<Real>>(sph_body), color_gradient_(particles_->getVariableDataByName<Vecd>("ColorGradient")), Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
                                                mass_(particles_->getVariableDataByName<Real>("Mass")), rho_(particles_->getVariableDataByName<Real>("Density"))
    {
        quantity_name_ = "DeltaFunction";
    }
    virtual ~DeltaFunction() {};

    Real reduce(size_t index_i, Real dt = 0.0)
    {
        return color_gradient_[index_i].norm() * Vol_[index_i];
        // return color_gradient_[index_i].norm() * mass_[index_i] / rho_[index_i];
        // return color_gradient_[index_i].norm() * particles_->getSPHBody().sph_adaptation_->ReferenceSpacing();
    }

  protected:
    Vecd *color_gradient_;
    Real *Vol_, *mass_, *rho_;
};

class ForceByHourglass : public LocalDynamicsReduce<ReduceSum<Real>>
{

  public:
    explicit ForceByHourglass(SPHBody &sph_body) : LocalDynamicsReduce<ReduceSum<Real>>(sph_body), force_by_hourglass_(particles_->getVariableDataByName<Vecd>("ForceByHourglass"))
    {
        quantity_name_ = "ForceByHourglass";
    }
    virtual ~ForceByHourglass() {};

    Real reduce(size_t index_i, Real dt = 0.0)
    {
        return force_by_hourglass_[index_i].norm();
    }

  protected:
    Vecd *force_by_hourglass_;
};