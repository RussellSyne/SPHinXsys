/**
 * @file    tethered_dead_fish_in_flow.cpp
 * @brief   fish flapping passively in flow
 * @author  Xiangyu Hu and Chi Zhang
 */
#include "sphinxsys.h"
/**
 * Create the shapes for fish and bones.
 */
#include "fish_and_bones.h"
/**
 * @brief Namespace cite here.
 */
using namespace SPH;
/**
 * @brief Basic geometry parameters and numerical setup.
 */
Real DL = 11.0;							/**< Channel length. */
Real DH = 8.0;							/**< Channel height. */
Real resolution_ref = 0.1;				/** Initial particle spacing. */
Real DL_sponge = resolution_ref * 20.0; /**< Sponge region to impose inflow condition. */
Real BW = resolution_ref * 4.0;			/**< Extending width for BCs. */
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(-DL_sponge - BW, -BW), Vec2d(DL + BW, DH + BW));
Real cx = 2.0;			  /**< Center of fish in x direction. */
Real cy = 4.0;			  /**< Center of fish in y direction. */
Real fish_length = 3.738; /**< Length of fish. */
Real fish_shape_resolution = resolution_ref * 0.5;
Vec3d tethering_point(-1.0, cy, 0.0); /**< The tethering point. */
/**
 * Material properties of the fluid.
 */
Real rho0_f = 1.0;
Real U_f = 1.0;
Real c_f = 10.0 * U_f;
Real Re = 5.0e3;
Real mu_f = rho0_f * U_f * (fish_length) / Re;
/**
 * Material properties of the fish body.
 */
Real rho0_s = 1.0;
Real poisson = 0.49;
Real Ae = 2.0e2;
Real Youngs_modulus = Ae * rho0_f * U_f * U_f;
/**
 * Basic geometries for construct SPH bodies.
 */
/**
 * @brief create a water block shape
 */
std::vector<Vecd> createWaterBlockShape()
{
	std::vector<Vecd> pnts_shaping_water_block;
	pnts_shaping_water_block.push_back(Vecd(-DL_sponge, 0.0));
	pnts_shaping_water_block.push_back(Vecd(-DL_sponge, DH));
	pnts_shaping_water_block.push_back(Vecd(DL, DH));
	pnts_shaping_water_block.push_back(Vecd(DL, 0.0));
	pnts_shaping_water_block.push_back(Vecd(-DL_sponge, 0.0));

	return pnts_shaping_water_block;
}
Vec2d buffer_halfsize = Vec2d(0.5 * DL_sponge, 0.5 * DH);
Vec2d buffer_translation = Vec2d(-DL_sponge, 0.0) + buffer_halfsize;
/**
 * @brief create outer wall shape
 */
std::vector<Vecd> createOuterWallShape()
{
	std::vector<Vecd> pnts_shaping_outer_wall;
	pnts_shaping_outer_wall.push_back(Vecd(-DL_sponge - BW, -BW));
	pnts_shaping_outer_wall.push_back(Vecd(-DL_sponge - BW, DH + BW));
	pnts_shaping_outer_wall.push_back(Vecd(DL + BW, DH + BW));
	pnts_shaping_outer_wall.push_back(Vecd(DL + BW, -BW));
	pnts_shaping_outer_wall.push_back(Vecd(-DL_sponge - BW, -BW));

	return pnts_shaping_outer_wall;
}
/**
 * @brief create inner wall shape
 */
std::vector<Vecd> createInnerWallShape()
{
	std::vector<Vecd> pnts_shaping_inner_wall;
	pnts_shaping_inner_wall.push_back(Vecd(-DL_sponge - 2.0 * BW, 0.0));
	pnts_shaping_inner_wall.push_back(Vecd(-DL_sponge - 2.0 * BW, DH));
	pnts_shaping_inner_wall.push_back(Vecd(DL + 2.0 * BW, DH));
	pnts_shaping_inner_wall.push_back(Vecd(DL + 2.0 * BW, 0.0));
	pnts_shaping_inner_wall.push_back(Vecd(-DL_sponge - 2.0 * BW, 0.0));

	return pnts_shaping_inner_wall;
}
/**
 * @brief create blocking shape to separate fish head out
 */
Real head_size = 1.0;
std::vector<Vecd> createFishBlockingShape()
{
	std::vector<Vecd> pnts_blocking_shape;
	pnts_blocking_shape.push_back(Vecd(cx + head_size, cy - 0.4));
	pnts_blocking_shape.push_back(Vecd(cx + head_size, cy + 0.4));
	pnts_blocking_shape.push_back(Vecd(cx + 5.0, cy + 0.4));
	pnts_blocking_shape.push_back(Vecd(cx + 5.0, cy - 0.4));
	pnts_blocking_shape.push_back(Vecd(cx + head_size, cy - 0.4));

	return pnts_blocking_shape;
}
/**
 * Water body shape defintion.
 */
class WaterBlock : public MultiPolygonShape
{
public:
	explicit WaterBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
	{
		multi_polygon_.addAPolygon(createWaterBlockShape(), ShapeBooleanOps::add);
		/** Exclude the fish body. */
		std::vector<Vecd> fish_shape = CreatFishShape(cx, cy, fish_length, fish_shape_resolution);
		multi_polygon_.addAPolygon(fish_shape, ShapeBooleanOps::sub);
	}
};
/**
 * Solid wall shape.
 */
class WallBoundary : public MultiPolygonShape
{
public:
	explicit WallBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
	{
		std::vector<Vecd> outer_shape = createOuterWallShape();
		std::vector<Vecd> inner_shape = createInnerWallShape();
		multi_polygon_.addAPolygon(outer_shape, ShapeBooleanOps::add);
		multi_polygon_.addAPolygon(inner_shape, ShapeBooleanOps::sub);
	}
};
/**
 * Fish body shape 
 */
class FishBody : public MultiPolygonShape
{
public:
	explicit FishBody(const std::string &shape_name) : MultiPolygonShape(shape_name)
	{
		std::vector<Vecd> fish_shape = CreatFishShape(cx, cy, fish_length, fish_shape_resolution);
		multi_polygon_.addAPolygon(fish_shape, ShapeBooleanOps::add);
	}
};
/**
 * @brief create fish head for constraint
 */
MultiPolygon createFishHeadShape(SPHBody &sph_body)
{
	std::vector<Vecd> fish_shape = CreatFishShape(cx, cy, fish_length, sph_body.sph_adaptation_->ReferenceSpacing());
	MultiPolygon multi_polygon;
	multi_polygon.addAPolygon(fish_shape, ShapeBooleanOps::add);
	multi_polygon.addAPolygon(createFishBlockingShape(), ShapeBooleanOps::sub);
	return multi_polygon;
};
/**
 * Observer particle generator.
 */
class FishObserverParticleGenerator : public ObserverParticleGenerator
{
public:
	explicit FishObserverParticleGenerator(SPHBody &sph_body) : ObserverParticleGenerator(sph_body)
	{
		positions_.push_back(Vecd(cx + resolution_ref, cy));
		positions_.push_back(Vecd(cx + fish_length - resolution_ref, cy));
	}
};
/**
 * Inflow boundary condition.
 */
class ParabolicInflow : public fluid_dynamics::InflowVelocityCondition
{
	Real u_ave_, u_ref_, t_ref;

public:
	ParabolicInflow(BodyAlignedBoxByCell &aligned_box_part)
		: InflowVelocityCondition(aligned_box_part),
		  u_ave_(0), u_ref_(1.0), t_ref(4.0) {}

	Vecd getPrescribedVelocity(Vecd &position, Vecd &velocity) override
	{
		Real u = velocity[0];
		Real v = velocity[1];
		if (position[0] < 0.0)
		{
			u = 6.0 * u_ave_ * position[1] * (DH - position[1]) / DH / DH;
			v = 0.0;
		}
		return Vecd(u, v);
	}

	void setupDynamics(Real dt = 0.0) override
	{
		Real run_time = GlobalStaticVariables::physical_time_;
		u_ave_ = run_time < t_ref ? 0.5 * u_ref_ * (1.0 - cos(Pi * run_time / t_ref)) : u_ref_;
	}
};
/**
 * Main program starts here.
 */
int main(int ac, char *av[])
{
	/**
	 * Build up context -- a SPHSystem.
	 */
	SPHSystem system(system_domain_bounds, resolution_ref);
	/** Tag for run particle relaxation for the initial body fitted distribution. */
	system.run_particle_relaxation_ = false;
	/** Tag for computation start with relaxed body fitted particles distribution. */
	system.reload_particles_ = true;
	/** Tag for computation from restart files. 0: start with initial condition. */
	system.restart_step_ = 0;
	system.handleCommandlineOptions(ac, av);
	IOEnvironment io_environment(system);

	/**
	 * @brief   Particles and body creation for water.
	 */
	FluidBody water_block(system, makeShared<WaterBlock>("WaterBody"));
	water_block.defineParticlesAndMaterial<FluidParticles, WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
	water_block.generateParticles<ParticleGeneratorLattice>();
	/**
	 * @brief   Particles and body creation for wall boundary.
	 */
	SolidBody wall_boundary(system, makeShared<WallBoundary>("Wall"));
	wall_boundary.defineParticlesAndMaterial<SolidParticles, Solid>();
	wall_boundary.generateParticles<ParticleGeneratorLattice>();
	/**
	 * @brief   Particles and body creation for fish.
	 */
	SolidBody fish_body(system, makeShared<FishBody>("FishBody"));
	fish_body.defineAdaptationRatios(1.15, 2.0);
	fish_body.defineBodyLevelSetShape();
	fish_body.defineParticlesAndMaterial<ElasticSolidParticles, NeoHookeanSolid>(rho0_s, Youngs_modulus, poisson);
	// Using relaxed particle distribution if needed
	(!system.run_particle_relaxation_ && system.reload_particles_)
		? fish_body.generateParticles<ParticleGeneratorReload>(io_environment, fish_body.getName())
		: fish_body.generateParticles<ParticleGeneratorLattice>();
	/**
	 * @brief   Particle and body creation of fish observer.
	 */
	ObserverBody fish_observer(system, "Observer");
	fish_observer.generateParticles<FishObserverParticleGenerator>();
	/** topology */
	BodyRelationInner water_block_inner(water_block);
	BodyRelationInner fish_body_inner(fish_body);
	ComplexBodyRelation water_block_complex(water_block_inner, {&wall_boundary, &fish_body});
	BodyRelationContact fish_body_contact(fish_body, {&water_block});
	BodyRelationContact fish_observer_contact(fish_observer, {&fish_body});

	BodyStatesRecordingToVtp write_real_body_states(io_environment, system.real_bodies_);
	ReducedQuantityRecording<ReduceDynamics<solid_dynamics::TotalForceOnSolid>> write_total_force_on_fish(io_environment, fish_body);
	ObservedQuantityRecording<Vecd> write_fish_displacement("Position", io_environment, fish_observer_contact);

	/** check whether run particle relaxation for body fitted particle distribution. */
	if (system.run_particle_relaxation_)
	{
		/**
		 * @brief 	Methods used for particle relaxation.
		 */
		/** Random reset the insert body particle position. */
		SimpleDynamics<RandomizeParticlePosition> random_fish_body_particles(fish_body);
		/** Write the body state to Vtp file. */
		BodyStatesRecordingToVtp write_fish_body(io_environment, fish_body);
		/** Write the particle reload files. */
		ReloadParticleIO write_particle_reload_files(io_environment, {&fish_body});

		/** A  Physics relaxation step. */
		relax_dynamics::RelaxationStepInner relaxation_step_inner(fish_body_inner);
		/**
		 * @brief 	Particle relaxation starts here.
		 */
		random_fish_body_particles.parallel_exec(0.25);
		relaxation_step_inner.surface_bounding_.parallel_exec();
		write_fish_body.writeToFile();

		/** relax particles of the insert body. */
		int ite_p = 0;
		while (ite_p < 1000)
		{
			relaxation_step_inner.parallel_exec();
			ite_p += 1;
			if (ite_p % 200 == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the inserted body N = " << ite_p << "\n";
				write_fish_body.writeToFile(ite_p);
			}
		}
		std::cout << "The physics relaxation process of inserted body finish !" << std::endl;

		/** Output results. */
		write_particle_reload_files.writeToFile();
		return 0;
	}

	/**
	 * This section define all numerical methods will be used in this case.
	 */
	/**
	 * @brief   Methods used for updating data structure.
	 */
	/** Periodic BCs in x direction. */
	PeriodicConditionUsingCellLinkedList periodic_condition(water_block, water_block.getBodyShapeBounds(), xAxis);
	SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
	SimpleDynamics<NormalDirectionFromBodyShape> fish_body_normal_direction(fish_body);
	/** Corrected configuration.*/
	InteractionDynamics<solid_dynamics::CorrectConfiguration>
		fish_body_corrected_configuration(fish_body_inner);
	/**
	 * Common particle dynamics.
	 */
	SimpleDynamics<TimeStepInitialization> initialize_a_fluid_step(water_block);
	/** Evaluation of density by summation approach. */
	InteractionWithUpdate<fluid_dynamics::DensitySummationComplex> update_density_by_summation(water_block_complex);
	/** Time step size without considering sound wave speed. */
	ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(water_block, U_f);
	/** Time step size with considering sound wave speed. */
	ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(water_block);
	/** Pressure relaxation using verlet time stepping. */
	Dynamics1Level<fluid_dynamics::PressureRelaxationWithWall> pressure_relaxation(water_block_complex);
	Dynamics1Level<fluid_dynamics::DensityRelaxationRiemannWithWall> density_relaxation(water_block_complex);
	/** Computing viscous acceleration. */
	InteractionDynamics<fluid_dynamics::ViscousAccelerationWithWall> viscous_acceleration(water_block_complex);
	/** Impose transport velocity formulation. */
	InteractionDynamics<fluid_dynamics::TransportVelocityCorrectionComplex> transport_velocity_correction(water_block_complex);
	/** Computing vorticity in the flow. */
	InteractionDynamics<fluid_dynamics::VorticityInner> compute_vorticity(water_block_inner);
	/** Inflow boundary condition. */
	BodyAlignedBoxByCell inflow_buffer(
		water_block, makeShared<AlignedBoxShape>(Transform2d(Vec2d(buffer_translation)), buffer_halfsize));
	SimpleDynamics<ParabolicInflow, BodyAlignedBoxByCell> parabolic_inflow(inflow_buffer);

	/**
	 * Fluid structure interaction model.
	 */
	InteractionDynamics<solid_dynamics::FluidViscousForceOnSolid> viscous_force_on_fish_body(fish_body_contact);
	InteractionDynamics<solid_dynamics::FluidForceOnSolidUpdate> fluid_force_on_fish_body(fish_body_contact, viscous_force_on_fish_body);
	/**
	 * Solid dynamics.
	 */
	/** Time step size calculation. */
	ReduceDynamics<solid_dynamics::AcousticTimeStepSize> fish_body_computing_time_step_size(fish_body);
	/** Process of stress relaxation. */
	Dynamics1Level<solid_dynamics::StressRelaxationFirstHalf>
		fish_body_stress_relaxation_first_half(fish_body_inner);
	Dynamics1Level<solid_dynamics::StressRelaxationSecondHalf>
		fish_body_stress_relaxation_second_half(fish_body_inner);
	/** Update normal direction on fish body.*/
	SimpleDynamics<solid_dynamics::UpdateElasticNormalDirection>
		fish_body_update_normal(fish_body);
	/** Compute the average velocity on fish body. */
	solid_dynamics::AverageVelocityAndAcceleration fish_body_average_velocity(fish_body);
	/**
	 * The multi body system from simbody.
	 */
	SimTK::MultibodySystem MBsystem;
	/** The bodies or matter of the MBsystem. */
	SimTK::SimbodyMatterSubsystem matter(MBsystem);
	/** The forces of the MBsystem.*/
	SimTK::GeneralForceSubsystem forces(MBsystem);
	SimTK::CableTrackerSubsystem cables(MBsystem);
	/** Mass properties of the fixed spot. */
	SimTK::Body::Rigid fixed_spot_info(SimTK::MassProperties(1.0, Vec3d(0), SimTK::UnitInertia(1)));
	SolidBodyPartForSimbody fish_head(fish_body, makeShared<MultiPolygonShape>(createFishHeadShape(fish_body), "FishHead"));
	/** Mass properties of the constrained spot. */
	SimTK::Body::Rigid tethered_spot_info(*fish_head.body_part_mass_properties_);
	/** Mobility of the fixed spot. */
	SimTK::MobilizedBody::Weld fixed_spot(matter.Ground(), SimTK::Transform(tethering_point),
										  fixed_spot_info, SimTK::Transform(Vec3d(0)));
	/** Mobility of the tethered spot.
	 * Set the mass center as the origin location of the planar mobilizer
	 */
	Vec3d displacement0 = fish_head.initial_mass_center_ - tethering_point;
	SimTK::MobilizedBody::Planar tethered_spot(fixed_spot,
											   SimTK::Transform(displacement0), tethered_spot_info, SimTK::Transform(Vec3d(0)));
	/** The tethering line give cable force.
	 * the start point of the cable path is at the origin location of the first mobilizer body,
	 * the end point is the tip of the fish head which has a distance to the origin
	 * location of the second mobilizer body origin location, here, the mass center
	 * of the fish head.
	 */
	Vec3d displacement_cable_end = Vec3d(cx, cy, 0.0) - fish_head.initial_mass_center_;
	SimTK::CablePath tethering_line(cables, fixed_spot, Vec3d(0), tethered_spot, displacement_cable_end);
	SimTK::CableSpring tethering_spring(forces, tethering_line, 100.0, 3.0, 10.0);

	// discreted forces acting on the bodies
	SimTK::Force::DiscreteForces force_on_bodies(forces, matter);
	fixed_spot_info.addDecoration(SimTK::Transform(), SimTK::DecorativeSphere(0.02));
	tethered_spot_info.addDecoration(SimTK::Transform(), SimTK::DecorativeSphere(0.4));
	/** Visualizer from simbody. */
	SimTK::Visualizer viz(MBsystem);
	SimTK::Visualizer::Reporter visualizer_reporter(viz, 0.01);
	MBsystem.addEventReporter(&visualizer_reporter);
	/** Initialize the system and state. */
	SimTK::State state = MBsystem.realizeTopology();
	viz.report(state);
	std::cout << "Hit ENTER to run a short simulation ...";
	getchar();
	/** Time stepping method for multibody system.*/
	SimTK::RungeKuttaMersonIntegrator integ(MBsystem);
	integ.setAccuracy(1e-3);
	integ.setAllowInterpolation(false);
	integ.initialize(state);
	/**
	 * Coupling between SimBody and SPH.
	 */
	ReduceDynamics<solid_dynamics::TotalForceForSimBody, SolidBodyPartForSimbody>
		force_on_tethered_spot(fish_head, MBsystem, tethered_spot, force_on_bodies, integ);
	SimpleDynamics<solid_dynamics::ConstraintBySimBody, SolidBodyPartForSimbody>
		constraint_tethered_spot(fish_head, MBsystem, tethered_spot, force_on_bodies, integ);

	/**
	 * Time steeping starts here.
	 */
	GlobalStaticVariables::physical_time_ = 0.0;
	/**
	 * Initial periodic boundary condition which copies the particle identifies
	 * as extra cell linked list form periodic regions to the corresponding boundaries
	 * for building up of extra configuration.
	 */
	system.initializeSystemCellLinkedLists();
	periodic_condition.update_cell_linked_list_.parallel_exec();
	system.initializeSystemConfigurations();
	/** Prepare quantities, e.g. wall normal, fish body norm,
	 * fluid initial number density and configuration of fish particles, will be used once only.
	 */
	wall_boundary_normal_direction.parallel_exec();
	fish_body_normal_direction.parallel_exec();
	fish_body_corrected_configuration.parallel_exec();
	/** Output for initial condition. */
	write_real_body_states.writeToFile(0);
	write_fish_displacement.writeToFile(0);
	/**
	 * Time parameters
	 */
	int number_of_iterations = 0;
	int screen_output_interval = 100;
	Real end_time = 200.0;
	Real output_interval = end_time / 200.0;
	Real dt = 0.0;	 /**< Default acoustic time step sizes. */
	Real dt_s = 0.0; /**< Default acoustic time step sizes for solid. */
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;

	/**
	 * Main loop starts here.
	 */
	while (GlobalStaticVariables::physical_time_ < end_time)
	{
		Real integration_time = 0.0;
		while (integration_time < output_interval)
		{
			initialize_a_fluid_step.parallel_exec();
			Real Dt = get_fluid_advection_time_step_size.parallel_exec();
			update_density_by_summation.parallel_exec();
			viscous_acceleration.parallel_exec();
			transport_velocity_correction.parallel_exec(Dt);
			/** Viscous force exerting on fish body. */
			viscous_force_on_fish_body.parallel_exec();
			/** Update normal direction on fish body. */
			fish_body_update_normal.parallel_exec();
			Real relaxation_time = 0.0;
			while (relaxation_time < Dt)
			{
				// note that dt needs to sufficiently large to avoid divide zero
				// when computing solid average velocity for FSI
				dt = SMIN(get_fluid_time_step_size.parallel_exec(), Dt);
				/** Fluid dynamics process, first half. */
				pressure_relaxation.parallel_exec(dt);
				/** Fluid pressure force exerting on fish. */
				fluid_force_on_fish_body.parallel_exec();
				/** Fluid dynamics process, second half. */
				density_relaxation.parallel_exec(dt);
				/** Relax fish body by solid dynamics. */
				Real dt_s_sum = 0.0;
				fish_body_average_velocity.initialize_displacement_.parallel_exec();
				while (dt_s_sum < dt)
				{
					dt_s = SMIN(fish_body_computing_time_step_size.parallel_exec(), dt - dt_s_sum);
					fish_body_stress_relaxation_first_half.parallel_exec(dt_s);
					SimTK::State &state_for_update = integ.updAdvancedState();
					force_on_bodies.clearAllBodyForces(state_for_update);
					force_on_bodies.setOneBodyForce(state_for_update, tethered_spot,
													force_on_tethered_spot.parallel_exec());
					integ.stepBy(dt_s);
					constraint_tethered_spot.parallel_exec();
					fish_body_stress_relaxation_second_half.parallel_exec(dt_s);
					dt_s_sum += dt_s;
				}
				// note that dt needs to sufficiently large to avoid divide zero
				fish_body_average_velocity.update_averages_.parallel_exec(dt);
				write_total_force_on_fish.writeToFile(number_of_iterations);

				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
				parabolic_inflow.parallel_exec();
			}
			if (number_of_iterations % screen_output_interval == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
						  << GlobalStaticVariables::physical_time_
						  << "	Dt = " << Dt << "	dt = " << dt << "	dt_s = " << dt_s << "\n";
			}
			number_of_iterations++;

			// visualize the motion of rigid body
			viz.report(integ.getState());
			/** Water block configuration and periodic condition. */
			periodic_condition.bounding_.parallel_exec();
			water_block.updateCellLinkedListWithParticleSort(100);
			fish_body.updateCellLinkedList();
			periodic_condition.update_cell_linked_list_.parallel_exec();
			water_block_complex.updateConfiguration();
			/** Fish body contact configuration. */
			fish_body_contact.updateConfiguration();
			write_fish_displacement.writeToFile(number_of_iterations);
		}
		tick_count t2 = tick_count::now();
		compute_vorticity.parallel_exec();
		write_real_body_states.writeToFile();
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

	return 0;
}