#include "general_reduce.h"
#include <limits>

namespace SPH
{
//=================================================================================================//
VelocityBoundCheck::
    VelocityBoundCheck(SPHBody &sph_body, Real velocity_bound)
    : LocalDynamicsReduce<ReduceOR>(sph_body),
      vel_(particles_->getVariableDataByName<Vecd>("Velocity")),
      velocity_bound_(velocity_bound) {}
//=================================================================================================//
bool VelocityBoundCheck::reduce(size_t index_i, Real dt)
{
    return vel_[index_i].norm() > velocity_bound_;
}
//=================================================================================================//
MaximumSpeed::MaximumSpeed(SPHBody &sph_body)
    : LocalDynamicsReduce<ReduceMax>(sph_body),
      vel_(particles_->getVariableDataByName<Vecd>("Velocity"))
{
    quantity_name_ = "MaximumSpeed";
}
//=================================================================================================//
Real MaximumSpeed::reduce(size_t index_i, Real dt)
{
    return vel_[index_i].norm();
}
//=================================================================================================//
PositionLowerBound::PositionLowerBound(SPHBody &sph_body)
    : LocalDynamicsReduce<ReduceLowerBound>(sph_body),
      pos_(particles_->getVariableDataByName<Vecd>("Position"))
{
    quantity_name_ = "PositionLowerBound";
}
//=================================================================================================//
Vecd PositionLowerBound::reduce(size_t index_i, Real dt)
{
    return pos_[index_i];
}
//=================================================================================================//
PositionUpperBound::PositionUpperBound(SPHBody &sph_body)
    : LocalDynamicsReduce<ReduceUpperBound>(sph_body),
      pos_(particles_->getVariableDataByName<Vecd>("Position"))
{
    quantity_name_ = "PositionUpperBound";
}
//=================================================================================================//
Vecd PositionUpperBound::reduce(size_t index_i, Real dt)
{
    return pos_[index_i];
}
//=================================================================================================//
TotalKineticEnergy::TotalKineticEnergy(SPHBody &sph_body)
    : LocalDynamicsReduce<ReduceSum<Real>>(sph_body),
      mass_(particles_->getVariableDataByName<Real>("Mass")),
      vel_(particles_->getVariableDataByName<Vecd>("Velocity"))
{
    quantity_name_ = "TotalKineticEnergy";
}
//=================================================================================================//
Real TotalKineticEnergy::reduce(size_t index_i, Real dt)
{
    return 0.5 * mass_[index_i] * vel_[index_i].squaredNorm();
}
//=================================================================================================//
TotalMechanicalEnergy::TotalMechanicalEnergy(SPHBody &sph_body, Gravity &gravity)
    : TotalKineticEnergy(sph_body),
      gravity_(gravity),
      pos_(particles_->getVariableDataByName<Vecd>("Position"))
{
    quantity_name_ = "TotalMechanicalEnergy";
}
//=================================================================================================//
Real TotalMechanicalEnergy::reduce(size_t index_i, Real dt)
{
    return TotalKineticEnergy::reduce(index_i, dt) +
           mass_[index_i] * gravity_.getPotential(pos_[index_i]);
}
//=================================================================================================//
// added //
//=================================================================================================//
ElasticStrainEnergy::ElasticStrainEnergy(SPHBody &sph_body)
    : LocalDynamicsReduce<ReduceSum<Real>>(sph_body),
      continuum_(DynamicCast<GeneralContinuum>(this, particles_->getBaseMaterial())),
      Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
      p_(particles_->getVariableDataByName<Real>("Pressure")),
      strain_tensor_(particles_->getVariableDataByName<Matd>("StrainTensor")),
      shear_stress_(particles_->getVariableDataByName<Matd>("ShearStress")),
      E_(continuum_.getYoungsModulus()),
      pos_(particles_->getVariableDataByName<Vecd>("Position"))
{
    quantity_name_ = "ElasticStrainEnergy";
}
//=================================================================================================//
Real ElasticStrainEnergy::reduce(size_t index_i, Real dt)
{
    Matd stress_tensor_i = shear_stress_[index_i] - p_[index_i] * Matd::Identity();
    //Matd stress_tensor_i = shear_stress_[index_i] + p_[index_i] * Matd::Identity();
    return 0.5 * Vol_[index_i] * (stress_tensor_i.cwiseProduct(strain_tensor_[index_i]).sum());
}
//=================================================================================================//
ElasticStrainEnergyV2::ElasticStrainEnergyV2(SPHBody &sph_body)
    : LocalDynamicsReduce<ReduceSum<Real>>(sph_body),
      continuum_(DynamicCast<GeneralContinuum>(this, particles_->getBaseMaterial())),
      Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
      p_(particles_->getVariableDataByName<Real>("Pressure")),
      strain_tensor_(particles_->getVariableDataByName<Matd>("StrainTensor")),
      shear_stress_(particles_->getVariableDataByName<Matd>("ShearStress")),
      E_(continuum_.getYoungsModulus())
{
    quantity_name_ = "ElasticStrainEnergyV2";
}
//=================================================================================================//
Real ElasticStrainEnergyV2::reduce(size_t index_i, Real dt)
{
    Matd stress_tensor_i = shear_stress_[index_i] - p_[index_i] * Matd::Identity();
    //return 0.5 * Vol_[index_i] * (stress_tensor_i.cwiseProduct(strain_tensor_[index_i]).sum());
    return 0.5 * Vol_[index_i] * (stress_tensor_i.cwiseProduct(stress_tensor_i).sum()) / E_;
}
//=================================================================================================//
ElasticStrainEnergyVersion3::ElasticStrainEnergyVersion3(SPHBody &sph_body)
    : LocalDynamicsReduce<ReduceSum<Real>>(sph_body),
      continuum_(DynamicCast<GeneralContinuum>(this, particles_->getBaseMaterial())),
      Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
      p_(particles_->getVariableDataByName<Real>("Pressure")),
      strain_tensor_(particles_->getVariableDataByName<Matd>("StrainTensor")),
      shear_stress_(particles_->getVariableDataByName<Matd>("ShearStress")),
      E_(continuum_.getYoungsModulus())
{
    quantity_name_ = "ElasticStrainEnergyV2";
}
//=================================================================================================//
Real ElasticStrainEnergyVersion3::reduce(size_t index_i, Real dt)
{
    //Matd stress_tensor_i = shear_stress_[index_i] - p_[index_i] * Matd::Identity();
    //return 0.5 * Vol_[index_i] * (stress_tensor_i.cwiseProduct(strain_tensor_[index_i]).sum());
    //return 0.5 * Vol_[index_i] * (stress_tensor_i.cwiseProduct(stress_tensor_i).sum()) / Youngs_modulus;
    return 0.5 * Vol_[index_i] * E_ * (strain_tensor_[index_i].cwiseProduct(strain_tensor_[index_i]).sum());
}
//=================================================================================================//
LinearMomentum::LinearMomentum(SPHBody &sph_body)
    : LocalDynamicsReduce<ReduceSum<Vecd>>(sph_body),
      Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
      rho_(particles_->getVariableDataByName<Real>("Density")),
      mass_(particles_->getVariableDataByName<Real>("Mass")),
      vel_(particles_->getVariableDataByName<Vecd>("Velocity")),
      pos_(particles_->getVariableDataByName<Vecd>("Position"))
{
    quantity_name_ = "LinearMomentumX";
}
//=================================================================================================//
Vecd LinearMomentum::reduce(size_t index_i, Real dt)
{
    return mass_[index_i]  * vel_[index_i];
}
//=================================================================================================//
TotalAngularMomentum::TotalAngularMomentum(SPHBody &sph_body)
    : LocalDynamicsReduce<ReduceSum<Real>>(sph_body),
      Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
      rho_(particles_->getVariableDataByName<Real>("Density")),
      mass_(particles_->getVariableDataByName<Real>("Mass")),
      vel_(particles_->getVariableDataByName<Vecd>("Velocity")),
      pos_(particles_->getVariableDataByName<Vecd>("Position"))
{
    quantity_name_ = "TotalAngularMomentum";
}
//=================================================================================================//
Real TotalAngularMomentum::reduce(size_t index_i, Real dt)
{
    return mass_[index_i] * pos_[index_i].x() * vel_[index_i].y() - pos_[index_i].y() * vel_[index_i].x();
}
} // namespace SPH
