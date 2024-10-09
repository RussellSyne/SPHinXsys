/**
 * @file 	droplet.cpp
 * @brief 	A square droplet deforms to circle due to surface tension.
 * @details A momentum-conservative formulation for surface tension is used here
 *          to reach a long-term stable simulation.
 * @author Shuaihao Zhang and Xiangyu Hu
 */
#include "sphinxsys.h" //SPHinXsys Library.
using namespace SPH;   // Namespace cite here.
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
// Real DL = 1.0;                         /**< Tank length. */
// Real DH = 1.0;                         /**< Tank height. */
// Real LL = 0.6;                         /**< Liquid column length. */
// Real LH = 0.6;                         /**< Liquid column height. */
Real DL = 2.0;                         /**< Tank length. */
Real DH = 2.0;                         /**< Tank height. */
Real LL = 1.0;                         /**< Liquid column length. */
Real LH = 1.0;                         /**< Liquid column height. */
Real particle_spacing_ref = LL / 40.0; /**< Initial reference particle spacing. */
Real BW = particle_spacing_ref * 4;    /**< Extending width for BCs. */
//----------------------------------------------------------------------
//	Material parameters.
//----------------------------------------------------------------------
Real rho0_f = 1.0;       /**< Reference density of water. */
Real rho0_a = 0.001;     /**< Reference density of air. */
Real U_ref = 2.0;        /**< Characteristic velocity. */
Real c_f = 10.0 * U_ref; /**< Reference sound speed. */
Real mu_f = 0.2;         /**< Water viscosity. */
Real mu_a = 0.002;       /**< Air viscosity. */
Real surface_tension_coef = 1.0;

// Real rho0_f = 1.0; /**< Reference density of water. */
// Real rho0_a = 1.0; /**< Reference density of air. */
// // Pressure increases when increasing the U_ref
// Real U_ref = 1.0;        /**< Characteristic velocity. */
// Real c_f = 10.0 * U_ref; /**< Reference sound speed. */
// Real mu_f = 0.2;         /**< Water viscosity. */
// Real mu_a = 0.2;         /**< Air viscosity. */
// Real surface_tension_coef = 1.0;
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
Real droplet_radius = LL / 2;
Vecd droplet_halfsize = Vec2d(droplet_radius, droplet_radius); // local center at origin
Vecd droplet_translation = Vec2d(0, 0);                        // translation to global coordinates
//----------------------------------------------------------------------
// Water body shape definition.
//----------------------------------------------------------------------
class WaterBlock : public ComplexShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TransformShape<GeometricShapeBox>>(Transform(droplet_translation), droplet_halfsize);
    }
};
//----------------------------------------------------------------------
// Air body shape definition.
//----------------------------------------------------------------------/**
class AirBlock : public ComplexShape
{
  public:
    explicit AirBlock(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TransformShape<GeometricShapeBox>>(Transform(air_translation), air_halfsize);
        subtract<TransformShape<GeometricShapeBox>>(Transform(droplet_translation), droplet_halfsize);
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
//
class AveragePressure : public LocalDynamicsReduce<ReduceSum<Real>>
{

  public:
    explicit AveragePressure(SPHBody &sph_body) : LocalDynamicsReduce<ReduceSum<Real>>(sph_body),
                                                  p_(particles_->getVariableDataByName<Real>("Pressure"))
    {
        quantity_name_ = "AveragePressure";
    }
    virtual ~AveragePressure() {};

    Real reduce(size_t index_i, Real dt = 0.0)
    {
        return p_[index_i] / particles_->TotalRealParticles();
    }

  protected:
    Real *p_;
};

class PositionCenter : public LocalDynamicsReduce<ReduceSum<Vecd>>
{

  public:
    explicit PositionCenter(SPHBody &sph_body) : LocalDynamicsReduce<ReduceSum<Vecd>>(sph_body),
                                                 pos_(particles_->getVariableDataByName<Vecd>("Position"))
    {
        quantity_name_ = "PositionCenter";
    }
    virtual ~PositionCenter() {};

    Vecd reduce(size_t index_i, Real dt = 0.0)
    {
        return pos_[index_i] / particles_->TotalRealParticles();
    }

  protected:
    Vecd *pos_;
};
//
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

class ForceByHourglassContact : public LocalDynamicsReduce<ReduceSum<Real>>
{

  public:
    explicit ForceByHourglassContact(SPHBody &sph_body) : LocalDynamicsReduce<ReduceSum<Real>>(sph_body), force_by_hourglass_contact_(particles_->getVariableDataByName<Vecd>("ForceByHourglassContact"))
    {
        quantity_name_ = "ForceByHourglassContact";
    }
    virtual ~ForceByHourglassContact() {};

    Real reduce(size_t index_i, Real dt = 0.0)
    {
        return force_by_hourglass_contact_[index_i].norm();
    }

  protected:
    Vecd *force_by_hourglass_contact_;
};
//-
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
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
    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> corrected_configuration_water(water_inner);
    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> corrected_configuration_air(air_inner);
    Dynamics1Level<fluid_dynamics::MultiPhaseIntegration1stHalfWithWallRiemann> water_pressure_relaxation(water_inner, water_air_contact, water_wall_contact);
    Dynamics1Level<fluid_dynamics::MultiPhaseIntegration2ndHalfWithWallRiemann> water_density_relaxation(water_inner, water_air_contact, water_wall_contact);
    Dynamics1Level<fluid_dynamics::MultiPhaseIntegration1stHalfWithWallRiemann> air_pressure_relaxation(air_inner, air_water_contact, air_wall_contact);
    Dynamics1Level<fluid_dynamics::MultiPhaseIntegration2ndHalfWithWallRiemann> air_density_relaxation(air_inner, air_water_contact, air_wall_contact);

    InteractionWithUpdate<fluid_dynamics::BaseDensitySummationComplex<Inner<>, Contact<>, Contact<>>>
        update_air_density_by_summation(air_inner, air_water_contact, air_wall_contact);
    InteractionWithUpdate<fluid_dynamics::BaseDensitySummationComplex<Inner<>, Contact<>, Contact<>>>
        update_water_density_by_summation(water_inner, water_air_contact, water_wall_contact);
    // The pressure increases when decreasing the coefficient in the transport velocity
    InteractionWithUpdate<fluid_dynamics::MultiPhaseTransportVelocityCorrectionComplex<AllParticles>>
        air_transport_correction(ConstructorArgs(air_inner, 0.02), air_water_contact, air_wall_contact);
    InteractionWithUpdate<fluid_dynamics::MultiPhaseTransportVelocityCorrectionComplex<AllParticles>>
        water_transport_correction(ConstructorArgs(water_inner, 0.02), water_air_contact, water_wall_contact);

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
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> get_water_time_step_size(water_block);
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> get_air_time_step_size(air_block);
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
    ReducedQuantityRecording<AveragePressure> write_water_average_pressure(water_block);
    ReducedQuantityRecording<AveragePressure> write_air_average_pressure(air_block);
    ReducedQuantityRecording<PositionCenter> write_position_center_pressure(water_block);
    ReducedQuantityRecording<DeltaFunction> write_water_delta_function(water_block);
    ReducedQuantityRecording<DeltaFunction> write_air_delta_function(air_block);
    ReducedQuantityRecording<ForceByHourglass> write_water_hourglass_force(water_block);
    ReducedQuantityRecording<ForceByHourglass> write_air_hourglass_force(air_block);
    // ReducedQuantityRecording<ForceByHourglassContact> write_water_contact_hourglass_force(water_block);
    // ReducedQuantityRecording<ForceByHourglassContact> write_air_contact_hourglass_force(air_block);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    wall_boundary_normal_direction.exec();
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    size_t number_of_iterations = 0;
    int screen_output_interval = 100;
    Real end_time = 5.0;
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
    write_water_average_pressure.writeToFile(number_of_iterations);
    write_air_average_pressure.writeToFile(number_of_iterations);
    write_position_center_pressure.writeToFile(number_of_iterations);
    write_water_hourglass_force.writeToFile(number_of_iterations);
    write_air_hourglass_force.writeToFile(number_of_iterations);
    // write_water_contact_hourglass_force.writeToFile(number_of_iterations);
    // write_air_contact_hourglass_force.writeToFile(number_of_iterations);
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
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
            Real Dt = SMIN(Dt_f, Dt_a);

            // The pressure is better when without density summation
            update_air_density_by_summation.exec();
            update_water_density_by_summation.exec();
            air_transport_correction.exec();
            water_transport_correction.exec();

            air_viscous_force.exec();
            water_viscous_force.exec();

            water_surface_tension_stress.exec();
            air_surface_tension_stress.exec();
            water_surface_tension_force.exec();
            air_surface_tension_force.exec();

            interval_computing_time_step += TickCount::now() - time_instance;

            /** Dynamics including pressure relaxation. */
            time_instance = TickCount::now();
            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                // air_viscous_force.exec();
                // water_viscous_force.exec();

                // water_surface_tension_stress.exec();
                // air_surface_tension_stress.exec();
                // water_surface_tension_force.exec();
                // air_surface_tension_force.exec();

                Real dt_f = get_water_time_step_size.exec();
                Real dt_a = get_air_time_step_size.exec();
                dt = SMIN(SMIN(dt_f, dt_a), Dt);

                water_pressure_relaxation.exec(dt);
                air_pressure_relaxation.exec(dt);

                water_density_relaxation.exec(dt);
                air_density_relaxation.exec(dt);

                relaxation_time += dt;
                integration_time += dt;
                physical_time += dt;

                // if (number_of_iterations % screen_output_interval == 0)
                // {
                //     std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                //               << physical_time
                //               << "	Dt = " << Dt << "	dt = " << dt
                //               << "	dt_f = " << dt_f << "	dt_a = " << dt_a << "\n";
                // }
            }
            interval_computing_pressure_relaxation += TickCount::now() - time_instance;

            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << physical_time
                          << "	Dt = " << Dt << "	dt = " << dt << "\n";
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

            interval_updating_configuration += TickCount::now() - time_instance;
        }
        write_water_kinetic_energy.writeToFile(number_of_iterations);
        write_water_linear_momentum.writeToFile(number_of_iterations);
        write_water_angular_momentum.writeToFile(number_of_iterations);
        write_air_linear_momentum.writeToFile(number_of_iterations);
        write_air_angular_momentum.writeToFile(number_of_iterations);
        write_water_average_pressure.writeToFile(number_of_iterations);
        write_air_average_pressure.writeToFile(number_of_iterations);
        write_position_center_pressure.writeToFile(number_of_iterations);
        write_water_delta_function.writeToFile(number_of_iterations);
        write_air_delta_function.writeToFile(number_of_iterations);
        write_water_hourglass_force.writeToFile(number_of_iterations);
        write_air_hourglass_force.writeToFile(number_of_iterations);
        // write_water_contact_hourglass_force.writeToFile(number_of_iterations);
        // write_air_contact_hourglass_force.writeToFile(number_of_iterations);
        TickCount t2 = TickCount::now();
        body_states_recording.writeToFile();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }

    // single time step
    // while (physical_time < end_time)
    // {
    //     Real integration_time = 0.0;
    //     /** Integrate time (loop) until the next output time. */
    //     while (integration_time < output_interval)
    //     {
    //         /** Force Prior due to viscous force and gravity. */
    //         time_instance = TickCount::now();

    //         //update_air_density_by_summation.exec();
    //         //update_water_density_by_summation.exec();
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
    //                       << "	dt = " << dt
    //                       << "	dt_f = " << dt_f
    //                       << "	dt_a = " << dt_a << "\n";
    //             // std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
    //             //           << physical_time
    //             //           << "	dt = " << dt << "\n";
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
    //     write_water_average_pressure.writeToFile(number_of_iterations);
    //     write_air_average_pressure.writeToFile(number_of_iterations);
    //     write_position_center_pressure.writeToFile(number_of_iterations);
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