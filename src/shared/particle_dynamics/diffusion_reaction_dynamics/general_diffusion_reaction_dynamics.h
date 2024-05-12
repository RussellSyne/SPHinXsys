/* ------------------------------------------------------------------------- *
 *                                SPHinXsys                                  *
 * ------------------------------------------------------------------------- *
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for    *
 * physical accurate simulation and aims to model coupled industrial dynamic *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH   *
 * (smoothed particle hydrodynamics), a meshless computational method using  *
 * particle discretization.                                                  *
 *                                                                           *
 * SPHinXsys is partially funded by German Research Foundation               *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,            *
 *  HU1527/12-1 and HU1527/12-4.                                             *
 *                                                                           *
 * Portions copyright (c) 2017-2023 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file    general_diffusion_reaction_dynamics.h
 * @brief   This is the particle dynamics applicable for all type of particles.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef GENERAL_DIFFUSION_REACTION_DYNAMICS_H
#define GENERAL_DIFFUSION_REACTION_DYNAMICS_H

#include "all_particle_dynamics.h"
#include "diffusion_reaction.h"
#include "diffusion_reaction_particles.h"

namespace SPH
{
template <class ParticlesType>
using DiffusionReactionSimpleData = DataDelegateSimple<ParticlesType>;

template <class ParticlesType>
using DiffusionReactionInnerData = DataDelegateInner<ParticlesType>;

template <class ParticlesType, class ContactParticlesType>
using DiffusionReactionContactData =
    DataDelegateContact<ParticlesType, ContactParticlesType>;

/**
 * @class DiffusionReactionInitialCondition
 * @brief Pure abstract class for initial conditions
 */
template <class ParticlesType>
class DiffusionReactionInitialCondition
    : public LocalDynamics,
      public DiffusionReactionSimpleData<ParticlesType>
{
  public:
    explicit DiffusionReactionInitialCondition(SPHBody &sph_body);
    virtual ~DiffusionReactionInitialCondition(){};

  protected:
    StdLargeVec<Vecd> &pos_;
    StdVec<StdLargeVec<Real>> &all_species_;
};

/**
 * @class DiffusionBasedMapping
 * @brief Mapping inside of body according to diffusion.
 * This is a abstract class to be override for case specific implementation
 */
template <class ParticlesType>
class DiffusionBasedMapping
    : public LocalDynamics,
      public DiffusionReactionSimpleData<ParticlesType>
{
  public:
    explicit DiffusionBasedMapping(SPHBody &sph_body)
        : LocalDynamics(sph_body),
          DiffusionReactionSimpleData<ParticlesType>(sph_body),
          pos_(*this->base_particles_.getVariableByName<Vecd>("Position")),
          all_species_(this->particles_->all_species_){};
    virtual ~DiffusionBasedMapping(){};

  protected:
    StdLargeVec<Vecd> &pos_;
    StdVec<StdLargeVec<Real>> &all_species_;
};

} // namespace SPH
#endif // GENERAL_DIFFUSION_REACTION_DYNAMICS_H