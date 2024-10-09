#include "oscillating_drop.h" //SPHinXsys Library.
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    std::cout << "Reynolds number of water:" << Re_water << std::endl;
    std::cout << "Reynolds number of air:" << Re_air << std::endl;
    //----------------------------------------------------------------------
    //	Build up an SPHSystem.
    //----------------------------------------------------------------------
    BoundingBox system_domain_bounds(Vec2d(-BW - DL / 2, -BW - DH / 2), Vec2d(BW + DL / 2, BW + DH / 2));
    SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
    sph_system.handleCommandlineOptions(ac, av)->setIOEnvironment();
    //----------------------------------------------------------------------
    //	Creating bodies with corresponding materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineMaterial<WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
    water_block.generateParticles<BaseParticles, Lattice>();

    FluidBody air_block(sph_system, makeShared<AirBlock>("AirBody"));
    air_block.defineMaterial<WeaklyCompressibleFluid>(rho0_a, c_f, mu_a);
    air_block.generateParticles<BaseParticles, Lattice>();

    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary"));
    wall_boundary.defineMaterial<Solid>();
    wall_boundary.generateParticles<BaseParticles, Lattice>();
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //----------------------------------------------------------------------
    InnerRelation water_inner(water_block);
    ContactRelation water_air_contact(water_block, {&air_block});
    ContactRelation water_wall_contact(water_block, {&wall_boundary});
    InnerRelation air_inner(air_block);
    ContactRelation air_water_contact(air_block, {&water_block});
    ContactRelation air_wall_contact(air_block, {&wall_boundary});
    //----------------------------------------------------------------------
    // Combined relations built from basic relations
    // which is only used for update configuration.
    //----------------------------------------------------------------------
    ComplexRelation water_air_complex(water_inner, {&water_air_contact, &water_wall_contact});
    ComplexRelation air_water_complex(air_inner, {&air_water_contact, &air_wall_contact});
    //----------------------------------------------------------------------
    //	Define the numerical methods used in the simulation.
    //	Note that there may be data dependence on the sequence of constructions.
    //----------------------------------------------------------------------
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
    InteractionWithUpdate<MultiPhaseLinearGradientCorrectionMatrixComplex> corrected_configuration_water(water_inner, water_air_contact, water_wall_contact);
    InteractionWithUpdate<MultiPhaseLinearGradientCorrectionMatrixComplex> corrected_configuration_air(air_inner, air_water_contact, air_wall_contact);

    Dynamics1Level<fluid_dynamics::MultiPhaseIntegration1stHalfWithWallRiemann> water_pressure_relaxation(water_inner, water_air_contact, water_wall_contact);
    Dynamics1Level<fluid_dynamics::MultiPhaseIntegration2ndHalfWithWallRiemann> water_density_relaxation(water_inner, water_air_contact, water_wall_contact);
    Dynamics1Level<fluid_dynamics::MultiPhaseIntegration1stHalfWithWallRiemann> air_pressure_relaxation(air_inner, air_water_contact, air_wall_contact);
    Dynamics1Level<fluid_dynamics::MultiPhaseIntegration2ndHalfWithWallRiemann> air_density_relaxation(air_inner, air_water_contact, air_wall_contact);

    InteractionWithUpdate<fluid_dynamics::BaseDensitySummationComplex<Inner<>, Contact<>, Contact<>>>
        update_air_density_by_summation(air_inner, air_water_contact, air_wall_contact);
    InteractionWithUpdate<fluid_dynamics::BaseDensitySummationComplex<Inner<>, Contact<>, Contact<>>>
        update_water_density_by_summation(water_inner, water_air_contact, water_wall_contact);
    InteractionWithUpdate<fluid_dynamics::MultiPhaseTransportVelocityCorrectionComplex<AllParticles>>
        air_transport_correction(ConstructorArgs(air_inner, 0.2), air_water_contact, air_wall_contact);
    InteractionWithUpdate<fluid_dynamics::MultiPhaseTransportVelocityCorrectionComplex<AllParticles>>
        water_transport_correction(ConstructorArgs(water_inner, 0.2), water_air_contact, water_wall_contact);

    InteractionWithUpdate<fluid_dynamics::MultiPhaseViscousForceWithWall> water_viscous_force(water_inner, water_air_contact, water_wall_contact);
    InteractionWithUpdate<fluid_dynamics::MultiPhaseViscousForceWithWall> air_viscous_force(air_inner, air_water_contact, air_wall_contact);

    InteractionDynamics<fluid_dynamics::SurfaceTensionStress> water_surface_tension_stress(water_air_contact, surface_tension_coef);
    InteractionDynamics<fluid_dynamics::SurfaceTensionStress> air_surface_tension_stress(air_water_contact, surface_tension_coef);
    InteractionWithUpdate<fluid_dynamics::SurfaceStressForceComplex> water_surface_tension_force(water_inner, water_air_contact);
    InteractionWithUpdate<fluid_dynamics::SurfaceStressForceComplex> air_surface_tension_force(air_inner, air_water_contact);

    ReduceDynamics<fluid_dynamics::AdvectionViscousTimeStep> get_water_advection_time_step_size(water_block, U_ref);
    ReduceDynamics<fluid_dynamics::AdvectionViscousTimeStep> get_air_advection_time_step_size(air_block, U_ref);
    // ReduceDynamics<fluid_dynamics::AdvectionTimeStep> get_water_advection_time_step_size(water_block, U_ref);
    // ReduceDynamics<fluid_dynamics::AdvectionTimeStep> get_air_advection_time_step_size(air_block, U_ref);
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> get_water_time_step_size(water_block, 0.2);
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> get_air_time_step_size(air_block, 0.2);
    //
    SimpleDynamics<InitialVelocity> initial_condition(water_block);
    SimpleDynamics<SetPositionCenterZero> set_position_center_zero(water_block);
    SimpleDynamics<CalPositionCenter, SequencedPolicy> calculate_position_center(water_block);
    SimpleDynamics<OutputPositionCenter> output_position_center(water_block);
    //
    SimpleDynamics<fluid_dynamics::ContinuumVolumeUpdate> water_volume_update(water_block);
    SimpleDynamics<fluid_dynamics::ContinuumVolumeUpdate> air_volume_update(air_block);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp body_states_recording(sph_system);
    body_states_recording.addToWrite<Real>(water_block, "Density");
    body_states_recording.addToWrite<Real>(water_block, "Pressure");
    body_states_recording.addToWrite<Matd>(water_block, "SurfaceTensionStress");
    body_states_recording.addToWrite<Vecd>(air_block, "ForcePrior");
    body_states_recording.addToWrite<Real>(air_block, "Density");
    body_states_recording.addToWrite<Real>(air_block, "Pressure");
    RegressionTestDynamicTimeWarping<ReducedQuantityRecording<TotalKineticEnergy>> write_water_kinetic_energy(water_block);
    ReducedQuantityRecording<LinearMomentum> write_water_linear_momentum(water_block);
    ReducedQuantityRecording<TotalAngularMomentum> write_water_angular_momentum(water_block);
    ReducedQuantityRecording<LinearMomentum> write_air_linear_momentum(air_block);
    ReducedQuantityRecording<TotalAngularMomentum> write_air_angular_momentum(air_block);
    ReducedQuantityRecording<PositionCenter> write_position_center(water_block);
    ReducedQuantityRecording<PositionCenterV2> write_position_center_v2(water_block);
    ReducedQuantityRecording<VelocityCenter> write_velocity_center(water_block);
    ReducedQuantityRecording<GyrationTensorS1> write_gyration_tensor_S1(water_block);
    ReducedQuantityRecording<GyrationTensorS2> write_gyration_tensor_S2(water_block);
    ReducedQuantityRecording<PositionLowerBound> write_position_lower_bond(water_block);
    ReducedQuantityRecording<PositionUpperBound> write_position_upper_bond(water_block);
    ReducedQuantityRecording<DeltaFunction> write_water_delta_function(water_block);
    ReducedQuantityRecording<DeltaFunction> write_air_delta_function(air_block);
    ReducedQuantityRecording<ForceByHourglass> write_water_hourglass_force(water_block);
    ReducedQuantityRecording<ForceByHourglass> write_air_hourglass_force(air_block);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    wall_boundary_normal_direction.exec();
    initial_condition.exec();
    calculate_position_center.exec();
    corrected_configuration_water.exec();
    corrected_configuration_air.exec();
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    size_t number_of_iterations = 0;
    int screen_output_interval = 500;
    Real end_time = 0.5;
    Real output_interval = end_time / 50; /**< Time stamps for output of body states. */
    Real dt = 0.0;                        /**< Default acoustic time step sizes. */
    /** statistics for computing CPU time. */
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    TimeInterval interval_computing_time_step;
    TimeInterval interval_computing_pressure_relaxation;
    TimeInterval interval_updating_configuration;
    TickCount time_instance;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    body_states_recording.writeToFile(0);
    write_water_kinetic_energy.writeToFile(number_of_iterations);
    write_water_linear_momentum.writeToFile(number_of_iterations);
    write_water_angular_momentum.writeToFile(number_of_iterations);
    write_air_linear_momentum.writeToFile(number_of_iterations);
    write_air_angular_momentum.writeToFile(number_of_iterations);
    write_position_center.writeToFile(number_of_iterations);
    write_position_center_v2.writeToFile(number_of_iterations);
    write_velocity_center.writeToFile(number_of_iterations);
    write_gyration_tensor_S1.writeToFile(number_of_iterations);
    write_gyration_tensor_S2.writeToFile(number_of_iterations);
    write_position_lower_bond.writeToFile(number_of_iterations);
    write_position_upper_bond.writeToFile(number_of_iterations);
    write_water_hourglass_force.writeToFile(number_of_iterations);
    write_air_hourglass_force.writeToFile(number_of_iterations);
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    /* dual time step*/
    while (physical_time < end_time)
    {
        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */
        while (integration_time < output_interval)
        {
            /** Force Prior due to viscous force and gravity. */
            time_instance = TickCount::now();

            Real Dt_f = get_water_advection_time_step_size.exec();
            Real Dt_a = get_air_advection_time_step_size.exec();
            Real Dt = 0.2 * SMIN(Dt_f, Dt_a);

            update_air_density_by_summation.exec();
            update_water_density_by_summation.exec();
            // water_volume_update.exec();
            // air_volume_update.exec();
            air_transport_correction.exec();
            water_transport_correction.exec();

            // air_viscous_force.exec();
            // water_viscous_force.exec();

            // water_surface_tension_stress.exec();
            // air_surface_tension_stress.exec();
            // water_surface_tension_force.exec();
            // air_surface_tension_force.exec();

            interval_computing_time_step += TickCount::now() - time_instance;

            /** Dynamics including pressure relaxation. */
            time_instance = TickCount::now();
            Real relaxation_time = 0.0;
            Real dt_f = 0.0;
            Real dt_a = 0.0;
            while (relaxation_time < Dt)
            {
                air_viscous_force.exec();
                water_viscous_force.exec();

                water_surface_tension_stress.exec();
                air_surface_tension_stress.exec();
                water_surface_tension_force.exec();
                air_surface_tension_force.exec();

                dt_f = get_water_time_step_size.exec();
                dt_a = get_air_time_step_size.exec();
                dt = SMIN(SMIN(dt_f, dt_a), Dt);
                water_pressure_relaxation.exec(dt);
                air_pressure_relaxation.exec(dt);

                water_density_relaxation.exec(dt);
                air_density_relaxation.exec(dt);

                relaxation_time += dt;
                integration_time += dt;
                physical_time += dt;
            }
            interval_computing_pressure_relaxation += TickCount::now() - time_instance;

            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << physical_time
                          << "	Dt = " << Dt << "	dt = " << dt << "	dt_f = " << dt_f << "	dt_a = " << dt_a << "\n";
            }
            number_of_iterations++;

            /** Update cell linked list and configuration. */
            time_instance = TickCount::now();

            water_block.updateCellLinkedList();
            water_air_complex.updateConfiguration();
            water_wall_contact.updateConfiguration();

            air_block.updateCellLinkedList();
            air_water_complex.updateConfiguration();
            air_wall_contact.updateConfiguration();

            // corrected_configuration_water.exec();
            // corrected_configuration_air.exec();

            interval_updating_configuration += TickCount::now() - time_instance;
        }
        write_water_kinetic_energy.writeToFile(number_of_iterations);
        write_water_linear_momentum.writeToFile(number_of_iterations);
        write_water_angular_momentum.writeToFile(number_of_iterations);
        write_air_linear_momentum.writeToFile(number_of_iterations);
        write_air_angular_momentum.writeToFile(number_of_iterations);
        write_position_center.writeToFile(number_of_iterations);
        write_position_center_v2.writeToFile(number_of_iterations);
        write_velocity_center.writeToFile(number_of_iterations);
        set_position_center_zero.exec();
        calculate_position_center.exec();
        // output_position_center.exec();
        write_gyration_tensor_S1.writeToFile(number_of_iterations);
        write_gyration_tensor_S2.writeToFile(number_of_iterations);
        write_position_lower_bond.writeToFile(number_of_iterations);
        write_position_upper_bond.writeToFile(number_of_iterations);
        write_water_delta_function.writeToFile(number_of_iterations);
        write_air_delta_function.writeToFile(number_of_iterations);
        write_water_hourglass_force.writeToFile(number_of_iterations);
        write_air_hourglass_force.writeToFile(number_of_iterations);
        TickCount t2 = TickCount::now();
        body_states_recording.writeToFile();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    /* single time step*/
    // while (physical_time < end_time)
    // {
    //     Real integration_time = 0.0;
    //     /** Integrate time (loop) until the next output time. */
    //     while (integration_time < output_interval)
    //     {
    //         /** Force Prior due to viscous force and gravity. */
    //         time_instance = TickCount::now();

    //         // update_air_density_by_summation.exec();
    //         // update_water_density_by_summation.exec();
    //         air_transport_correction.exec();
    //         water_transport_correction.exec();

    //         air_viscous_force.exec();
    //         water_viscous_force.exec();

    //         water_surface_tension_stress.exec();
    //         air_surface_tension_stress.exec();
    //         water_surface_tension_force.exec();
    //         air_surface_tension_force.exec();

    //         interval_computing_time_step += TickCount::now() - time_instance;

    //         /** Dynamics including pressure relaxation. */
    //         time_instance = TickCount::now();
    //         Real dt_f = get_water_time_step_size.exec();
    //         Real dt_a = get_air_time_step_size.exec();
    //         dt = SMIN(dt_f, dt_a);

    //         water_pressure_relaxation.exec(dt);
    //         air_pressure_relaxation.exec(dt);

    //         water_density_relaxation.exec(dt);
    //         air_density_relaxation.exec(dt);

    //         integration_time += dt;
    //         physical_time += dt;
    //         interval_computing_pressure_relaxation += TickCount::now() - time_instance;

    //         if (number_of_iterations % screen_output_interval == 0)
    //         {
    //             std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
    //                       << physical_time
    //                       << "	dt = " << dt << "\n";
    //         }
    //         number_of_iterations++;

    //         /** Update cell linked list and configuration. */
    //         time_instance = TickCount::now();

    //         water_block.updateCellLinkedList();
    //         water_air_complex.updateConfiguration();
    //         water_wall_contact.updateConfiguration();

    //         air_block.updateCellLinkedList();
    //         air_water_complex.updateConfiguration();
    //         air_wall_contact.updateConfiguration();

    //         interval_updating_configuration += TickCount::now() - time_instance;
    //     }
    //     write_water_kinetic_energy.writeToFile(number_of_iterations);
    //     write_water_linear_momentum.writeToFile(number_of_iterations);
    //     write_water_angular_momentum.writeToFile(number_of_iterations);
    //     write_air_linear_momentum.writeToFile(number_of_iterations);
    //     write_air_angular_momentum.writeToFile(number_of_iterations);
    //     write_position_center.writeToFile(number_of_iterations);
    //     write_position_center_v2.writeToFile(number_of_iterations);
    //     write_velocity_center.writeToFile(number_of_iterations);
    //     set_position_center_zero.exec();
    //     calculate_position_center.exec();
    //     //output_position_center.exec();
    //     write_gyration_tensor_S1.writeToFile(number_of_iterations);
    //     write_gyration_tensor_S2.writeToFile(number_of_iterations);
    //     write_position_lower_bond.writeToFile(number_of_iterations);
    //     write_position_upper_bond.writeToFile(number_of_iterations);
    //     write_water_delta_function.writeToFile(number_of_iterations);
    //     write_air_delta_function.writeToFile(number_of_iterations);
    //     TickCount t2 = TickCount::now();
    //     body_states_recording.writeToFile();
    //     TickCount t3 = TickCount::now();
    //     interval += t3 - t2;
    // }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds()
              << " seconds." << std::endl;
    std::cout << std::fixed << std::setprecision(9) << "interval_computing_time_step ="
              << interval_computing_time_step.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_computing_pressure_relaxation = "
              << interval_computing_pressure_relaxation.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_updating_configuration = "
              << interval_updating_configuration.seconds() << "\n";

    return 0;
}