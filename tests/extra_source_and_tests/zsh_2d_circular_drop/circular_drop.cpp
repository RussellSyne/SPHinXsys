/**
 * @file	rotation_patch.cpp
 * @brief	2D rotation patch example.
 * @author	Yaru Ren, Chi Zhang and Xiangyu Hu
 */

#include "sphinxsys.h" //SPHinXsys Library.
using namespace SPH;   // Namespace cite here.
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real droplet_radius = 1.0;
Real particle_spacing_ref = droplet_radius / 100; /**< Initial reference particle spacing. */
BoundingBox system_domain_bounds(Vec2d(-droplet_radius, -3 * droplet_radius), Vec2d(droplet_radius, 3 * droplet_radius));
//----------------------------------------------------------------------
//	Material parameters.
//----------------------------------------------------------------------
Real rho0_f = 1000.0;    /**< Reference density of fluid. */
Real U_max = 141.0;      /**< Characteristic velocity. */
Real c_f = 10.0 * U_max; /**< Reference sound speed. */

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
/**
 * application dependent initial velocity
 */
class InitialVelocity
    : public fluid_dynamics::FluidInitialCondition
{
  public:
    InitialVelocity(SPHBody &sph_body)
        : fluid_dynamics::FluidInitialCondition(sph_body){};

    void update(size_t index_i, Real dt)
    {
        /** initial velocity profile */
        vel_[index_i][0] = -100 * pos_[index_i][0];
        vel_[index_i][1] = 100 * pos_[index_i][1];
    }
};
/**
 * SymmetryCheck
 */
class SymmetryCheck
    : public fluid_dynamics::BaseIntegration<DataDelegateInner>
{
  public:
    SymmetryCheck(BaseInnerRelation &inner_relation)
        : fluid_dynamics::BaseIntegration<DataDelegateInner>(inner_relation){};

    void interaction(size_t index_i, Real dt)
    {
        Neighborhood &inner_neighborhood_ij = inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood_ij.current_size_; ++n)
        {
            size_t index_j = inner_neighborhood_ij.j_[n];
            Neighborhood &inner_neighborhood_ji = inner_configuration_[index_j];
            size_t is_symmetry = 0;
            for (size_t m = 0; m != inner_neighborhood_ji.current_size_; ++m)
            {
                size_t index_j_i = inner_neighborhood_ji.j_[m];
                if (index_j_i == index_i)
                    is_symmetry = 1;
            }
            if (!is_symmetry)
            {
                std::cout << "Warning! Asymmetric neighbour!" << std::endl;
                std::cout << "Particle_i_index:" << index_i << std::endl;
                std::cout << "Neighbours of particle i:" << std::endl;
                for (size_t n = 0; n != inner_neighborhood_ij.current_size_; ++n)
                {
                    size_t index_ij = inner_neighborhood_ij.j_[n];
                    std::cout << index_ij << std::endl;
                }
                std::cout << "Particle_j_index:" << index_j << std::endl;
                std::cout << "Neighbours of particle j:" << std::endl;
                for (size_t m = 0; m != inner_neighborhood_ji.current_size_; ++m)
                {
                    size_t index_ji = inner_neighborhood_ji.j_[m];
                    std::cout << index_ji << std::endl;
                }
                exit(0);
            }
        }
    }
};

//
// class TotalMomentumSetZero
//     : public fluid_dynamics::FluidInitialCondition
// {
//   public:
//     TotalMomentumSetZero(SPHBody &sph_body)
//         : fluid_dynamics::FluidInitialCondition(sph_body),
//         total_momentum_(particles->registerSingularVariable<Vecd>("TotalMomentum")->ValueAddress()) {};

//     void update(size_t index_i, Real dt)
//     {
//         total_momentum_ = Vecd::Zero();
//     }
//   protected:
//     Vecd &total_momentum_;
// };
//
class TotalMomentumUpdate
    : public fluid_dynamics::FluidInitialCondition
{
  public:
    TotalMomentumUpdate(SPHBody &sph_body)
        : fluid_dynamics::FluidInitialCondition(sph_body),
          total_momentum_(particles_->registerSingularVariable<Vecd>("TotalMomentum")->ValueAddress()),
          total_momentum_previous_(particles_->registerSingularVariable<Vecd>("TotalMomentumPrevious")->ValueAddress()),
          mass_(particles_->getVariableDataByName<Real>("Mass"))
    {
        total_momentum_previous_ = total_momentum_;
        *total_momentum_ = Vecd::Zero();
        std::cout << "Constructor of TotalMomentumUpdate." << std::endl;
    };

    void update(size_t index_i, Real dt)
    {
        *total_momentum_ += mass_[index_i] * vel_[index_i];
    }
    
    /* Destructor is not called! */
    virtual ~TotalMomentumUpdate()
    {
        std::cout << "Destructor of TotalMomentumUpdate." << std::endl;
        if ((*total_momentum_ - *total_momentum_previous_).norm() > 1.0e-8)
        {
            std::cout << "Destructor of TotalMomentumUpdate." << std::endl;
            std::cout << "Warning! Total momentum changed a lot!" << std::endl;
            std::cout << "Total momentum:" << *total_momentum_ << std::endl;
            std::cout << "Previous total momentum:" << *total_momentum_previous_ << std::endl;
            exit(0);
        }
    };
  protected:
    Vecd *total_momentum_, *total_momentum_previous_;
    Real *mass_;
};

class TotalMomentumCheck
    : public fluid_dynamics::FluidInitialCondition
{
  public:
    TotalMomentumCheck(SPHBody &sph_body)
        : fluid_dynamics::FluidInitialCondition(sph_body),
          total_momentum_(*(this->particles_->template getSingularVariableByName<Vecd>("TotalMomentum")->ValueAddress())),
          total_momentum_previous_(*(this->particles_->template getSingularVariableByName<Vecd>("TotalMomentumPrevious")->ValueAddress()))
          {
            if ((total_momentum_ - total_momentum_previous_).norm() > 1.0e-8)
            {
                std::cout << "TotalMomentumCheck." << std::endl;
                std::cout << "Warning! Total momentum changed a lot!" << std::endl;
                std::cout << "Previous total momentum:" << total_momentum_previous_.transpose() << std::endl;
                std::cout << "Current total momentum:" << total_momentum_.transpose() << std::endl;
                exit(0);
            }
          };

    void update(size_t index_i, Real dt)
    {
        // if ((total_momentum_ - total_momentum_previous_).norm() > 1.0e-7)
        // {
        //     std::cout << "TotalMomentumCheck." << std::endl;
        //     std::cout << "Warning! Total momentum changed a lot!" << std::endl;
        //     std::cout << "Previous total momentum:" << total_momentum_previous_.transpose() << std::endl;
        //     std::cout << "Current total momentum:" << total_momentum_.transpose() << std::endl;
        //     exit(0);
        // }
    }

  protected:
    Vecd &total_momentum_, &total_momentum_previous_;
};
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up an SPHSystem.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
    sph_system.handleCommandlineOptions(ac, av)->setIOEnvironment();
    //----------------------------------------------------------------------
    //	Creating bodies with corresponding materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineAdaptation<SPHAdaptation>(1.3, 1.0);
    water_block.defineMaterial<WeaklyCompressibleFluid>(rho0_f, c_f);
    // Using relaxed particle distribution if needed
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? water_block.generateParticles<BaseParticles, Reload>(water_block.getName())
        : water_block.generateParticles<BaseParticles, Lattice>();

    ObserverBody fluid_observer(sph_system, "FluidObserver");
    StdVec<Vecd> observation_location = {Vecd(0.0, 0.0)};
    fluid_observer.generateParticles<ObserverParticles>(observation_location);
    //----------------------------------------------------------------------
    // Combined relations built from basic relations
    // which is only used for update configuration.
    //----------------------------------------------------------------------
    InnerRelation water_body_inner(water_block);
    ContactRelation fluid_observer_contact(fluid_observer, {&water_block});
    //----------------------------------------------------------------------
    // Define the numerical methods used in the simulation.
    // Note that there may be data dependence on the sequence of constructions.
    // Generally, the geometric models or simple objects without data dependencies,
    // such as gravity, should be initiated first.
    // Then the major physical particle dynamics model should be introduced.
    // Finally, the auxillary models such as time step estimator, initial condition,
    // boundary condition and other constraints should be defined.
    //----------------------------------------------------------------------
    InteractionWithUpdate<SpatialTemporalFreeSurfaceIndicationInner> free_surface_indicator(water_body_inner);
    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> corrected_configuration_fluid(water_body_inner, 0.5);
    Dynamics1Level<fluid_dynamics::Integration1stHalfCorrectionInnerRiemann> fluid_pressure_relaxation_correct(water_body_inner);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfInnerRiemann> fluid_density_relaxation(water_body_inner);
    InteractionWithUpdate<fluid_dynamics::DensitySummationFreeSurfaceInner> update_density_by_summation(water_body_inner);
    ReduceDynamics<fluid_dynamics::AdvectionTimeStep> fluid_advection_time_step(water_block, U_max);
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> fluid_acoustic_time_step(water_block);
    SimpleDynamics<InitialVelocity> initial_condition(water_block);
    InteractionDynamics<SymmetryCheck> symmetry_check(water_body_inner);
    SimpleDynamics<TotalMomentumUpdate, SequencedPolicy> update_total_momentum(water_block);
    SimpleDynamics<TotalMomentumCheck, SequencedPolicy> check_total_momentum(water_block);
    // InteractionWithUpdate<fluid_dynamics::TransportVelocityCorrectionInner<NoLimiter, BulkParticles>, SequencedPolicy> transport_velocity_correction(water_body_inner);
    InteractionWithUpdate<fluid_dynamics::TransportVelocityCorrectionInner<NoLimiter, BulkParticles>> transport_velocity_correction(water_body_inner);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp body_states_recording(sph_system);
    body_states_recording.addToWrite<Real>(water_block, "DensitySummation");
    body_states_recording.addToWrite<int>(water_block, "Indicator");
    body_states_recording.addToWrite<Real>(water_block, "Density");
    body_states_recording.addToWrite<Real>(water_block, "Pressure");
    body_states_recording.addToWrite<Real>(water_block, "VolumetricMeasure");
    RestartIO restart_io(sph_system);
    RegressionTestDynamicTimeWarping<ReducedQuantityRecording<TotalKineticEnergy>> write_water_kinetic_energy(water_block);
    ReducedQuantityRecording<LinearMomentum> write_linear_momentum(water_block);
    ReducedQuantityRecording<TotalAngularMomentum> write_angular_momentum(water_block);
    ObservedQuantityRecording<Real> write_recorded_water_pressure("Pressure", fluid_observer_contact);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    initial_condition.exec();
    //----------------------------------------------------------------------
    //	Load restart file if necessary.
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    if (sph_system.RestartStep() != 0)
    {
        physical_time = restart_io.readRestartFiles(sph_system.RestartStep());
        water_block.updateCellLinkedList();
        water_body_inner.updateConfiguration();
    }
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    size_t number_of_iterations = sph_system.RestartStep();
    int screen_output_interval = 100;
    int restart_output_interval = screen_output_interval * 10;
    Real end_time = 0.01;
    Real output_interval = end_time / 20;
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    TimeInterval interval_computing_time_step;
    TimeInterval interval_computing_fluid_pressure_relaxation;
    TimeInterval interval_updating_configuration;
    TickCount time_instance;
    //----------------------------------------------------------------------
    //	First state recording before the main loop.
    //----------------------------------------------------------------------
    body_states_recording.writeToFile();
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    // while (physical_time < end_time)
    // {
    //     Real integration_time = 0.0;
    //     /** Integrate time (loop) until the next output time. */
    //     while (integration_time < output_interval)
    //     {
    //         /** outer loop for dual-time criteria time-stepping. */
    //         time_instance = TickCount::now();
    //         free_surface_indicator.exec();
    //         corner_indicator.exec();
    //         all_indicator.exec();
    //         corrected_configuration_fluid.exec();
    //         //update_density_by_summation.exec();
    //         free_surface_normal.exec();
    //         transport_velocity_correction.exec();
    //         interval_computing_time_step += TickCount::now() - time_instance;

    //         time_instance = TickCount::now();

    //         /** inner loop for dual-time criteria time-stepping.  */
    //         Real acoustic_dt = fluid_acoustic_time_step.exec();
    //         fluid_pressure_relaxation_correct.exec(acoustic_dt);
    //         fluid_density_relaxation.exec(acoustic_dt);
    //         integration_time += acoustic_dt;
    //         physical_time += acoustic_dt;

    //         interval_computing_fluid_pressure_relaxation += TickCount::now() - time_instance;

    //         /** screen output, write body reduced values and restart files  */
    //         if (number_of_iterations % screen_output_interval == 0)
    //         {
    //             std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
    //                       << physical_time
    //                       << "	 acoustic_dt = " << acoustic_dt << "\n";
    //             if (number_of_iterations % restart_output_interval == 0)
    //                 restart_io.writeToFile(number_of_iterations);
    //         }
    //         number_of_iterations++;

    //         /** Update cell linked list and configuration. */
    //         time_instance = TickCount::now();
    //         water_block.updateCellLinkedListWithParticleSort(100);
    //         water_body_inner.updateConfiguration();
    //         fluid_observer_contact.updateConfiguration();
    //         interval_updating_configuration += TickCount::now() - time_instance;
    //     }
    //     write_water_kinetic_energy.writeToFile(number_of_iterations);
    //     write_linear_momentum_x.writeToFile(number_of_iterations);
    //     write_linear_momentum_y.writeToFile(number_of_iterations);
    //     write_angular_momentum.writeToFile(number_of_iterations);
    //     write_recorded_water_pressure.writeToFile(number_of_iterations);
    //     body_states_recording.writeToFile();

    //     TickCount t2 = TickCount::now();
    //     TickCount t3 = TickCount::now();
    //     interval += t3 - t2;
    // }
    // dual time step
    while (physical_time < end_time)
    {
        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */
        while (integration_time < output_interval)
        {
            /** outer loop for dual-time criteria time-stepping. */
            time_instance = TickCount::now();
            Real Dt = fluid_advection_time_step.exec();

            free_surface_indicator.exec();
            corrected_configuration_fluid.exec();
            update_density_by_summation.exec();
            transport_velocity_correction.exec();
            interval_computing_time_step += TickCount::now() - time_instance;

            time_instance = TickCount::now();
            size_t inner_ite_dt = 0;
            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                /** inner loop for dual-time criteria time-stepping.  */
                // std::cout<<"number_of_iterations:"<<number_of_iterations<<std::endl;
                //symmetry_check.exec();
                //update_total_momentum.exec();
                //check_total_momentum.exec();
                Real acoustic_dt = fluid_acoustic_time_step.exec();
                fluid_pressure_relaxation_correct.exec(acoustic_dt);
                fluid_density_relaxation.exec(acoustic_dt);
                relaxation_time += acoustic_dt;
                integration_time += acoustic_dt;
                inner_ite_dt++;
                physical_time += acoustic_dt;
            }
            interval_computing_fluid_pressure_relaxation += TickCount::now() - time_instance;

            /** screen output, write body reduced values and restart files  */
            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << physical_time
                          << "	Dt = " << Dt << "	Dt / dt = " << inner_ite_dt << "\n";

                if (number_of_iterations % restart_output_interval == 0)
                    restart_io.writeToFile(number_of_iterations);

                write_water_kinetic_energy.writeToFile(number_of_iterations);
                write_recorded_water_pressure.writeToFile(number_of_iterations);
            }
            number_of_iterations++;

            /** Update cell linked list and configuration. */
            time_instance = TickCount::now();
            water_block.updateCellLinkedList();
            water_body_inner.updateConfiguration();
            fluid_observer_contact.updateConfiguration();
            interval_updating_configuration += TickCount::now() - time_instance;
        }
        write_linear_momentum.writeToFile(number_of_iterations);
        write_angular_momentum.writeToFile(number_of_iterations);
        body_states_recording.writeToFile();

        TickCount t2 = TickCount::now();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }

    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds()
              << " seconds." << std::endl;
    std::cout << std::fixed << std::setprecision(9) << "interval_computing_time_step ="
              << interval_computing_time_step.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_computing_fluid_pressure_relaxation = "
              << interval_computing_fluid_pressure_relaxation.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_updating_configuration = "
              << interval_updating_configuration.seconds() << "\n";

    return 0;
};