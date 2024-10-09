#include "surface_tension.hpp"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
SurfaceTensionStress::
    SurfaceTensionStress(BaseContactRelation &contact_relation, Real surface_tension_coef)
    : LocalDynamics(contact_relation.getSPHBody()), DataDelegateContact(contact_relation),
      color_gradient_(particles_->registerStateVariable<Vecd>("ColorGradient")),
      norm_direction_(particles_->registerStateVariable<Vecd>("NormDirection")),
      surface_tension_stress_(particles_->registerStateVariable<Matd>("SurfaceTensionStress")),
      surface_tension_coef_(particles_->registerSingularVariable<Real>("SurfaceTensionCoef", surface_tension_coef)->ValueAddress()),
      B_(this->particles_->getVariableDataByName<Matd>("LinearGradientCorrectionMatrix"))
{
    particles_->addVariableToSort<Vecd>("ColorGradient");
    particles_->addVariableToWrite<Vecd>("ColorGradient");
    particles_->addVariableToSort<Matd>("SurfaceTensionStress");
    particles_->addVariableToWrite<Matd>("SurfaceTensionStress");
    particles_->addVariableToSort<Vecd>("NormDirection");
    particles_->addVariableToWrite<Vecd>("NormDirection");
    Real rho0 = getSPHBody().base_material_->ReferenceDensity();
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        contact_Vol_.push_back(contact_particles_[k]->getVariableDataByName<Real>("VolumetricMeasure"));
        Real rho0_k = contact_bodies_[k]->base_material_->ReferenceDensity();
        contact_fraction_.push_back(rho0 / (rho0 + rho0_k));
    }
}
//=================================================================================================//
void SurfaceTensionStress::interaction(size_t index_i, Real dt)
{
    color_gradient_[index_i] = ZeroData<Vecd>::value;
    surface_tension_stress_[index_i] = ZeroData<Matd>::value;
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        Vecd weighted_color_gradient = ZeroData<Vecd>::value;
        Real contact_fraction_k = contact_fraction_[k];
        Real *Vol_k = contact_Vol_[k];
        const Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            // weighted_color_gradient -= 2 * contact_fraction_k *
            //                            contact_neighborhood.dW_ij_[n] * Vol_k[index_j] * B_[index_i] * contact_neighborhood.e_ij_[n];
            weighted_color_gradient -= 2 * contact_fraction_k *
                                       contact_neighborhood.dW_ij_[n] * Vol_k[index_j] * contact_neighborhood.e_ij_[n];
        }
        color_gradient_[index_i] = weighted_color_gradient;
        norm_direction_[index_i] = weighted_color_gradient / (weighted_color_gradient.norm() + Eps);
        // norm_direction_[index_i] = B_[index_i] * weighted_color_gradient / ((B_[index_i]*weighted_color_gradient).norm() + Eps);
        surface_tension_stress_[index_i] += *surface_tension_coef_ *
                                            (Matd::Identity() - norm_direction_[index_i] * norm_direction_[index_i].transpose()) *
                                            weighted_color_gradient.norm();
        // Real norm = weighted_color_gradient.norm();
        // surface_tension_stress_[index_i] += *surface_tension_coef_ / (norm + Eps) *
        //                                     (norm * norm * Matd::Identity() -
        //                                      weighted_color_gradient * weighted_color_gradient.transpose());
    }
}
//=================================================================================================//
//=================================================================================================//
//=================================================================================================//
//=================================================================================================//
SurfaceTensionStressInner::
    SurfaceTensionStressInner(BaseInnerRelation &inner_relation)
    : LocalDynamics(inner_relation.getSPHBody()), DataDelegateInner(inner_relation),
      color_gradient_(particles_->registerStateVariable<Vecd>("ColorGradient")),
      norm_direction_(particles_->registerStateVariable<Vecd>("NormDirection")),
      surface_tension_stress_(particles_->registerStateVariable<Matd>("SurfaceTensionStress")),
      surface_tension_coef_(particles_->getSingularVariableByName<Real>("SurfaceTensionCoef")->ValueAddress()),
      Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
      B_(this->particles_->getVariableDataByName<Matd>("LinearGradientCorrectionMatrix")) {}
//=================================================================================================//
void SurfaceTensionStressInner::interaction(size_t index_i, Real dt)
{
    color_gradient_[index_i] = ZeroData<Vecd>::value;
    Vecd weighted_color_gradient = ZeroData<Vecd>::value;
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    if (norm_direction_[index_i].norm() > Eps)
    {
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            size_t index_j = inner_neighborhood.j_[n];
            // if (norm_direction_[index_j].norm() > Eps)
            // {
            //     // weighted_color_gradient -= 2*inner_neighborhood.dW_ij_[n] * Vol_[index_j] * B_[index_i] * inner_neighborhood.e_ij_[n];
            //     weighted_color_gradient += 2 * inner_neighborhood.dW_ij_[n] * Vol_[index_j] * inner_neighborhood.e_ij_[n];
            // }
            // weighted_color_gradient += 2*inner_neighborhood.dW_ij_[n] * Vol_[index_j] * B_[index_i] * inner_neighborhood.e_ij_[n];
            weighted_color_gradient += 2 * inner_neighborhood.dW_ij_[n] * Vol_[index_j] * inner_neighborhood.e_ij_[n];
        }
    }
    color_gradient_[index_i] = weighted_color_gradient;
    norm_direction_[index_i] = weighted_color_gradient / (weighted_color_gradient.norm() + Eps);
    // norm_direction_[index_i] = B_[index_i] * weighted_color_gradient / ((B_[index_i]*weighted_color_gradient).norm() + Eps);
    surface_tension_stress_[index_i] += *surface_tension_coef_ *
                                        (Matd::Identity() - norm_direction_[index_i] * norm_direction_[index_i].transpose()) *
                                        weighted_color_gradient.norm();
}
//=================================================================================================//
//=================================================================================================//
//=================================================================================================//
//=================================================================================================//
SurfaceStressForce<Inner<>>::SurfaceStressForce(BaseInnerRelation &inner_relation)
    : SurfaceStressForce<DataDelegateInner>(inner_relation),
      color_gradient_(particles_->getVariableDataByName<Vecd>("ColorGradient"))
{
    particles_->addVariableToSort<Vecd>("ForceByHourglass");
    particles_->addVariableToWrite<Vecd>("ForceByHourglass");
    particles_->addVariableToSort<Vecd>("ForceByHourglassContact");
    particles_->addVariableToWrite<Vecd>("ForceByHourglassContact");
}
//=================================================================================================//
void SurfaceStressForce<Inner<>>::interaction(size_t index_i, Real dt)
{
    Vecd summation = ZeroData<Vecd>::value;
    Vecd summation_hourglass = ZeroData<Vecd>::value;
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Real r_ij = inner_neighborhood.r_ij_[n];
        Vecd e_ij = inner_neighborhood.e_ij_[n];
        //===========================//
        /* hourglass control -- v1*/
        Real coef = 0.1;
        Vecd color_gradient_average = 0.5 * (color_gradient_[index_i] + color_gradient_[index_j]);
        Matd mismatch = Matd::Zero() - color_gradient_average * e_ij.transpose() * r_ij * (color_gradient_average * e_ij.transpose() * r_ij);
        Matd hourglass_correction = -coef * mismatch * *surface_tension_coef_ / r_ij;
        //===========================//
        /* hourglass control -- v2*/
        // Real coef = 0.1;
        // Vecd color_gradient_average = 0.5 * (color_gradient_[index_i] + color_gradient_[index_j]);
        // Matd mismatch = Matd::Zero() - color_gradient_average * e_ij.transpose() * r_ij * (color_gradient_average * e_ij.transpose() * r_ij) / ((color_gradient_average * e_ij.transpose() * r_ij).norm() + Eps);
        // Matd hourglass_correction = -coef * mismatch * *surface_tension_coef_ / r_ij;
        //===========================//
        summation += mass_[index_i] * inner_neighborhood.dW_ij_[n] * Vol_[index_j] *
                     (hourglass_correction +
                      surface_tension_stress_[index_i] + surface_tension_stress_[index_j]) *
                     inner_neighborhood.e_ij_[n];
        summation_hourglass += mass_[index_i] * inner_neighborhood.dW_ij_[n] * Vol_[index_j] *
                               hourglass_correction *
                               inner_neighborhood.e_ij_[n];
        //===========================//
        /* hourglass control -- v3*/
        // Real coef = 0.1;
        // Vecd color_gradient_average = 0.5 * (color_gradient_[index_i] + color_gradient_[index_j]);
        // Real hourglass_correction = 1.0 + coef * std::abs((0.0 - color_gradient_average.dot(e_ij) * r_ij)) / (r_ij * color_gradient_average.norm() + Eps);
        // // if (hourglass_correction >1)
        // //     std::cout << hourglass_correction << std::endl;
        // summation += mass_[index_i] * inner_neighborhood.dW_ij_[n] * Vol_[index_j] *
        //              hourglass_correction * (surface_tension_stress_[index_i] + surface_tension_stress_[index_j]) *
        //              inner_neighborhood.e_ij_[n];
        //===========================//
        /* without hourglass control*/
        // summation += mass_[index_i] * inner_neighborhood.dW_ij_[n] * Vol_[index_j] *
        //              (surface_tension_stress_[index_i] + surface_tension_stress_[index_j]) *
        //              inner_neighborhood.e_ij_[n];
    }
    surface_tension_force_[index_i] = summation / rho_[index_i];
    //
    force_by_hourglass_[index_i] = summation_hourglass / rho_[index_i];
}
//=================================================================================================//
SurfaceStressForce<Contact<>>::SurfaceStressForce(BaseContactRelation &contact_relation)
    : SurfaceStressForce<DataDelegateContact>(contact_relation)
{
    Real rho0 = getSPHBody().base_material_->ReferenceDensity();
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        Real rho0_k = contact_bodies_[k]->base_material_->ReferenceDensity();
        contact_fraction_.push_back(rho0 / (rho0 + rho0_k));
        contact_Vol_.push_back(contact_particles_[k]->getVariableDataByName<Real>("VolumetricMeasure"));
        contact_color_gradient_.push_back(
            contact_particles_[k]->getVariableDataByName<Vecd>("ColorGradient"));
        contact_surface_tension_stress_.push_back(
            contact_particles_[k]->getVariableDataByName<Matd>("SurfaceTensionStress"));
    }
}
//=================================================================================================//
void SurfaceStressForce<Contact<>>::interaction(size_t index_i, Real dt)
{
    Vecd summation = ZeroData<Vecd>::value;
    Vecd summation_hourglass = ZeroData<Vecd>::value;
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        Real contact_fraction_k = contact_fraction_[k];
        Real *Vol_k = contact_Vol_[k];
        Vecd *contact_color_gradient_k = contact_color_gradient_[k];
        Matd *contact_surface_tension_stress_k = contact_surface_tension_stress_[k];
        const Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            Real r_ij = contact_neighborhood.r_ij_[n];
            Vecd e_ij = contact_neighborhood.e_ij_[n];
            //===========================//
            /* hourglass control -- v-1*/
            // Real coef = 0.1; // for oscillating drop
            // Real mismatch = 1.0 - 0.5 * (color_gradient_[index_i] + contact_color_gradient_k[index_j]).dot(e_ij) * r_ij;
            // Matd hourglass_correction = -coef * mismatch * Matd::Identity();
            //===========================//
            /* hourglass control -- v0*/
            // Real coef = 0.1; // for oscillating drop
            // Real mismatch = 1.0 - 0.5 * (color_gradient_[index_i] + contact_color_gradient_k[index_j]).dot(e_ij) * r_ij;
            // Matd hourglass_correction = -coef * 2 * contact_fraction_k * mismatch * Matd::Identity() * *surface_tension_coef_ / r_ij;
            //===========================//
            /* hourglass control -- v1*/
            Real coef = 0.1; // for oscillating drop
            Vecd color_gradient_average = 0.5 * (color_gradient_[index_i] + contact_color_gradient_k[index_j]);
            Matd mismatch = Matd::Identity() - color_gradient_average * e_ij.transpose() * r_ij * (color_gradient_average * e_ij.transpose() * r_ij);
            Matd hourglass_correction = -coef * 2 * contact_fraction_k * mismatch * *surface_tension_coef_ / r_ij;
            //===========================//
            /* hourglass control -- v2*/
            // Real coef = 0.1; // for oscillating drop
            // Vecd color_gradient_average = 0.5 * (color_gradient_[index_i] + contact_color_gradient_k[index_j]);
            // Matd mismatch = Matd::Identity() - color_gradient_average * e_ij.transpose() * r_ij * (color_gradient_average * e_ij.transpose() * r_ij) / (color_gradient_average * e_ij.transpose() * r_ij).norm();
            // Matd hourglass_correction = -coef * 2 * contact_fraction_k * mismatch * *surface_tension_coef_ / r_ij;
            //===========================//
            summation += mass_[index_i] *
                         (hourglass_correction +
                          2 * (Real(1) - contact_fraction_k) * surface_tension_stress_[index_i] +
                          2 * contact_fraction_k * contact_surface_tension_stress_k[index_j]) *
                         contact_neighborhood.dW_ij_[n] * contact_neighborhood.e_ij_[n] * Vol_k[index_j];
            summation_hourglass += mass_[index_i] *
                                   hourglass_correction *
                                   contact_neighborhood.dW_ij_[n] * contact_neighborhood.e_ij_[n] * Vol_k[index_j];
            //===========================//
            /* hourglass control -- v3*/
            // Real coef = 0.1;
            // Vecd color_gradient_average = 0.5 * (color_gradient_[index_i] + contact_color_gradient_k[index_j]);
            // Real hourglass_correction = 1.0 + coef * std::abs((1.0 - color_gradient_average.dot(e_ij) * r_ij)) / (r_ij * color_gradient_average.norm() + Eps);
            // //std::cout << hourglass_correction << std::endl;
            // //
            // summation += mass_[index_i] * hourglass_correction *
            //              (2 * (Real(1) - contact_fraction_k) * surface_tension_stress_[index_i] + 2 * contact_fraction_k * contact_surface_tension_stress_k[index_j]) *
            //              contact_neighborhood.dW_ij_[n] * contact_neighborhood.e_ij_[n] * Vol_k[index_j];
        }
    }
    surface_tension_force_[index_i] += summation / rho_[index_i];
    //
    force_by_hourglass_contact_[index_i] = summation_hourglass / rho_[index_i];
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
