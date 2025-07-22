#include "setup.hpp"
#include "surface_export.hpp"
#include "mesh_streaming.hpp"
#include <memory>

#include "shapes.hpp"

// Global mesh streaming configuration
MeshStreamingConfig mesh_streaming_config;


struct Torus
{
	float x;
	float y;
	float distance_from_sealevel;
	float inner_radius;
	float outer_radius;

	// Constructor
	Torus(float x, float y, float distance_from_sealevel, float inner_radius, float outer_radius)
		: x(x), y(y), distance_from_sealevel(distance_from_sealevel), inner_radius(inner_radius), outer_radius(outer_radius) {}
};


void main_setup_right_hander()
{
	const float f = 0.001f;			 // make smaller
	const float u = 0.17f;			 // peak velocity of speaker membrane
	const float frequency = 0.0007f; // amplitude = u/(2.0f*pif*frequency);
	const float width = 50.0f, depth = 50.0f, height = 10.0f;
	const uint simulation_steps = 50u;
	const uint mem = 4000u;
	const uint3 dims = resolution(float3(width, depth, height), mem);
	LBM lbm(dims, 0.01f, 0.0f, 0.0f, -f);

	const uint nx = lbm.get_Nx(), ny = lbm.get_Ny(), nz = lbm.get_Nz();
	const float outer_radius = nx / 4;
	const float inner_radius = nx / 8;
	const uint D = 400, R = outer_radius;


	const uint num_torii = 4;

	float leading_x = (float)nx / 3, leading_y = (float)ny / 2;
	float z = (float)nz / 2 / 6;

	// straight barrel right hander
	// Torus torii[num_torii] = {
	// 	// leading point
	// 	{leading_x, leading_y, z, inner_radius, outer_radius},
	// 	{leading_x + inner_radius, leading_y + inner_radius, z, inner_radius, outer_radius},
	// 	{leading_x + inner_radius * 1.5, leading_y + inner_radius * 1.5, z, inner_radius, outer_radius},
	// 	{leading_x + inner_radius * 2.0, leading_y + inner_radius * 2.0, z, inner_radius, outer_radius}};

	Torus torii[num_torii] = {
		// leading point
		{leading_x - inner_radius, leading_y, z, inner_radius, outer_radius},
		{leading_x + inner_radius, leading_y + inner_radius * 2, z, inner_radius, outer_radius},
		{leading_x + inner_radius , leading_y + inner_radius * 3.0, z, inner_radius, outer_radius},
		{leading_x + inner_radius , leading_y + inner_radius * 4.0, z, inner_radius, outer_radius}};


	parallel_for(lbm.get_N(), [&](ulong n)
	 {
		 uint x = 0u, y = 0u, z = 0u;
		 lbm.coordinates(n, x, y, z);
		 const uint H = nz / 2u;
		 if (z < H)
		 {
			 lbm.flags[n] = TYPE_F;
			 lbm.rho[n] = units.rho_hydrostatic(f, (float)z, (float)H);
		 }

		 for (int i = 0; i < num_torii; ++i)
		 {
			 Torus torus = torii[i];
			 if (torus_z(x, y, z, float3(torus.x, torus.y, (nz / 2) - ((torus.outer_radius - torus.inner_radius) + torus.distance_from_sealevel)), torus.inner_radius, torus.outer_radius))
			 {
				 lbm.flags[n] = TYPE_S;
			 }
		 }

		 // boundaries
		 if (x == 0u || x == nx - 1u || y == 0u || y == ny - 1u || z == 0u || z == nz - 1u)
			 lbm.flags[n] = TYPE_S;

		 // inflow
		 if (y == 0u && x > 0u && x < nx - 1u && z > 0u && z < nz - 1u)
			 lbm.flags[n] = TYPE_E;

		 // outflow, doesnt work
		 //  if (y == Ny - 1u && x > 0u && x < Nx - 1u && z > 0u && z < Nz - 1u)
		 //  {
		 // 	 lbm.flags[n] = TYPE_E;
		 // 	 lbm.rho[n] = 0.9;
		 //  }
	 });

	#ifdef SURFACE_EXPORT
		// Parse surface export configuration from command line arguments
		SurfaceExportConfig surface_export_config = SurfaceExportConfig::parse_from_arguments(main_arguments);

		// Force enable surface export for testing (comment out for production)
		if(!surface_export_config.enabled) {
			surface_export_config.enabled = true;
			surface_export_config.directory = get_exe_path() + "export/";
			surface_export_config.export_interval = 100u;
			surface_export_config.ascii_format = false;
			println("DEBUG: Force-enabling surface export for testing");
			// Create directory if needed
			create_directory_if_not_exists(surface_export_config.directory);
		}

		println("DEBUG: Beach setup - Surface export config after parsing:");
		println("  Enabled: " + string(surface_export_config.enabled ? "true" : "false"));
		println("  Directory: " + surface_export_config.directory);
		println("  Interval: " + to_string(surface_export_config.export_interval));
		println("  ASCII format: " + string(surface_export_config.ascii_format ? "true" : "false"));
#endif // SURFACE_EXPORT

		// Parse mesh streaming configuration from command line arguments
		MeshStreamingConfig mesh_streaming_config = MeshStreamingConfig::parse_from_arguments(main_arguments);
		
		// Force enable mesh streaming for testing (comment out for production)
		if(!mesh_streaming_config.enabled) {
			mesh_streaming_config.enabled = true;
			mesh_streaming_config.host = "localhost";
			mesh_streaming_config.port = 15703;
			mesh_streaming_config.simulation_id = "fluidx3d_right_hander";
			mesh_streaming_config.stream_interval = 10u; // Stream every 10 timesteps
			mesh_streaming_config.verbose = true;
			println("DEBUG: Force-enabling mesh streaming for testing");
		}
		
		println("DEBUG: Mesh streaming config after parsing:");
		println("  Enabled: " + string(mesh_streaming_config.enabled ? "true" : "false"));
		println("  Target: " + mesh_streaming_config.host + ":" + std::to_string(mesh_streaming_config.port));
		println("  Simulation ID: " + mesh_streaming_config.simulation_id);
		println("  Interval: " + std::to_string(mesh_streaming_config.stream_interval));
		
		// Initialize mesh streamer
		std::unique_ptr<MeshStreamer> mesh_streamer = nullptr;

		while(true) { // main simulation loop
			lbm.u.read_from_device();
			const float uy = u*sinf(2.0f*pif*frequency*(float)lbm.get_t());
			const float uz = 0.5f*u*cosf(2.0f*pif*frequency*(float)lbm.get_t());
			for(uint z=1u; z<nz-1u; z++) {
				for(uint y=0u; y<1u; y++) {
					for(uint x=1u; x<nx-1u; x++) {
						const uint n = x+(y+z*ny)*nx;
						lbm.u.y[n] = uy;
						lbm.u.z[n] = uz;
					}
				}
			}
			lbm.u.write_to_device();
			lbm.run(10u);
  println("HERE");

#ifdef SURFACE_EXPORT
println("SURFACE_EXPORT");

			// Export surface at regular intervals
			if(surface_export_config.enabled && lbm.get_t() > 0u && surface_export_config.should_export(lbm.get_t())) {
			    println("DEBUG: Exporting surface frame at timestep " + to_string(lbm.get_t()));
				export_surface_frame(&lbm, lbm.get_t(), surface_export_config);
			}
#endif // SURFACE_EXPORT

			// Stream mesh to seaview at regular intervals
			if(mesh_streaming_config.enabled && lbm.get_t() > 0u && mesh_streaming_config.should_stream(lbm.get_t())) {
				println("DEBUG: Streaming mesh frame at timestep " + to_string(lbm.get_t()));
				stream_mesh_frame(&lbm, lbm.get_t(), mesh_streamer);
			}
		}


}


void main_setup_beach() { // breaking waves on beach; required extensions in defines.hpp: FP16S, VOLUME_FORCE, EQUILIBRIUM_BOUNDARIES, SURFACE, INTERACTIVE_GRAPHICS, SURFACE_EXPORT
	// ################################################################## define simulation box size, viscosity and volume force ###################################################################
	const float f = 0.001f; // make smaller
	const float u = 0.12f; // peak velocity of speaker membrane
	const float frequency = 0.0007f; // amplitude = u/(2.0f*pif*frequency);
	LBM lbm(128u, 640u, 96u, 0.01f, 0.0f, 0.0f, -f);
	// ###################################################################################### define geometry ######################################################################################
	const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); parallel_for(lbm.get_N(), [&](ulong n) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
		const uint H = Nz/2u;
		if(z<H) {
			lbm.flags[n] = TYPE_F;
			lbm.rho[n] = units.rho_hydrostatic(f, (float)z, (float)H);
		}
		if(plane(x, y, z, float3(lbm.center().x, 128.0f, 0.0f), float3(0.0f, -1.0f, 8.0f))) lbm.flags[n] = TYPE_S;
		if(x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==0u||z==Nz-1u) lbm.flags[n] = TYPE_S; // all non periodic
		if(y==0u && x>0u&&x<Nx-1u&&z>0u&&z<Nz-1u) lbm.flags[n] = TYPE_E;
	}); // ####################################################################### run simulation, export images and data ##########################################################################
#ifdef GRAPHICS
	lbm.graphics.visualization_modes = VIS_FLAG_LATTICE | (lbm.get_D()==1u ? VIS_PHI_RAYTRACE : VIS_PHI_RASTERIZE);
#endif // GRAPHICS
	lbm.run(0u); // initialize simulation

#ifdef SURFACE_EXPORT
	// Parse surface export configuration from command line arguments
	SurfaceExportConfig surface_export_config = SurfaceExportConfig::parse_from_arguments(main_arguments);

	// Force enable surface export for testing (comment out for production)
	if(!surface_export_config.enabled) {
		surface_export_config.enabled = true;
		surface_export_config.directory = get_exe_path() + "export/";
		surface_export_config.export_interval = 100u;
		surface_export_config.ascii_format = false;
		println("DEBUG: Force-enabling surface export for testing");
		// Create directory if needed
		create_directory_if_not_exists(surface_export_config.directory);
	}

	println("DEBUG: Beach setup - Surface export config after parsing:");
	println("  Enabled: " + string(surface_export_config.enabled ? "true" : "false"));
	println("  Directory: " + surface_export_config.directory);
	println("  Interval: " + to_string(surface_export_config.export_interval));
	println("  ASCII format: " + string(surface_export_config.ascii_format ? "true" : "false"));
	if(surface_export_config.enabled && surface_export_config.should_export(0u)) {
		println("DEBUG: Exporting initial surface frame");
		export_surface_frame(&lbm, 0u, surface_export_config);
	}
#endif // SURFACE_EXPORT

	while(true) { // main simulation loop
		lbm.u.read_from_device();
		const float uy = u*sinf(2.0f*pif*frequency*(float)lbm.get_t());
		const float uz = 0.5f*u*cosf(2.0f*pif*frequency*(float)lbm.get_t());
		for(uint z=1u; z<Nz-1u; z++) {
			for(uint y=0u; y<1u; y++) {
				for(uint x=1u; x<Nx-1u; x++) {
					const uint n = x+(y+z*Ny)*Nx;
					lbm.u.y[n] = uy;
					lbm.u.z[n] = uz;
				}
			}
		}
		lbm.u.write_to_device();
		lbm.run(10u);
  println("HERE");

#ifdef SURFACE_EXPORT
println("SURFACE_EXPORT");

		// Export surface at regular intervals
		if(surface_export_config.enabled && lbm.get_t() > 0u && surface_export_config.should_export(lbm.get_t())) {
		    println("DEBUG: Exporting surface frame at timestep " + to_string(lbm.get_t()));
			export_surface_frame(&lbm, lbm.get_t(), surface_export_config);
		}
#endif // SURFACE_EXPORT
	}
}

void main_setup_sphere() {
	const float f = 0.001f;			 // make smaller
	const float u = 0.17f;			 // peak velocity of speaker membrane
	const float frequency = 0.0007f; // amplitude = u/(2.0f*pif*frequency);
	const float width = 50.0f, depth = 50.0f, height = 10.0f;
	const uint simulation_steps = 50u;
	const uint mem = 1000u;
	const uint3 dims = resolution(float3(width, depth, height), mem);
	LBM lbm(dims, 0.01f, 0.0f, 0.0f, -f);

	const uint nx = lbm.get_Nx(), ny = lbm.get_Ny(), nz = lbm.get_Nz();

	// Sphere parameters
	const float sphere_radius = nx / 8.0f;  // Radius is 1/8 of domain width
	const float sphere_x = nx * 0.75f;      // Position sphere to the side (3/4 along x)
	const float sphere_y = ny / 2.0f;       // Center in y direction
	const float sphere_z = sphere_radius;   // Center on floor (z = radius so bottom touches floor)

	parallel_for(lbm.get_N(), [&](ulong n)
	 {
		 uint x = 0u, y = 0u, z = 0u;
		 lbm.coordinates(n, x, y, z);
		 const uint H = nz / 2u;
		 if (z < H)
		 {
			 lbm.flags[n] = TYPE_F;
			 lbm.rho[n] = units.rho_hydrostatic(f, (float)z, (float)H);
		 }

		 // Add sphere
		 if (sphere(x, y, z, float3(sphere_x, sphere_y, sphere_z), sphere_radius))
		 {
			 lbm.flags[n] = TYPE_S;
		 }

		 // boundaries
		 if (x == 0u || x == nx - 1u || y == 0u || y == ny - 1u || z == 0u || z == nz - 1u)
			 lbm.flags[n] = TYPE_S;

		 // inflow
		 if (y == 0u && x > 0u && x < nx - 1u && z > 0u && z < nz - 1u)
			 lbm.flags[n] = TYPE_E;
	 });

	#ifdef SURFACE_EXPORT
		// Parse surface export configuration from command line arguments
		SurfaceExportConfig surface_export_config = SurfaceExportConfig::parse_from_arguments(main_arguments);

		// Force enable surface export for testing (comment out for production)
		if(!surface_export_config.enabled) {
			surface_export_config.enabled = true;
			surface_export_config.directory = get_exe_path() + "export/";
			surface_export_config.export_interval = 100u;
			surface_export_config.ascii_format = false;
			println("DEBUG: Force-enabling surface export for testing");
			// Create directory if needed
			create_directory_if_not_exists(surface_export_config.directory);
		}

		println("DEBUG: Sphere setup - Surface export config after parsing:");
		println("  Enabled: " + string(surface_export_config.enabled ? "true" : "false"));
		println("  Directory: " + surface_export_config.directory);
		println("  Interval: " + to_string(surface_export_config.export_interval));
		println("  ASCII format: " + string(surface_export_config.ascii_format ? "true" : "false"));
		// if(surface_export_config.enabled && surface_export_config.should_export(0u)) {
		// 	println("DEBUG: Exporting initial surface frame");
		// 	export_surface_frame(&lbm, 0u, surface_export_config);
		// }
#endif // SURFACE_EXPORT

		while(true) { // main simulation loop
			lbm.u.read_from_device();
			const float uy = u*sinf(2.0f*pif*frequency*(float)lbm.get_t());
			const float uz = 0.5f*u*cosf(2.0f*pif*frequency*(float)lbm.get_t());
			for(uint z=1u; z<nz-1u; z++) {
				for(uint y=0u; y<1u; y++) {
					for(uint x=1u; x<nx-1u; x++) {
						const uint n = x+(y+z*ny)*nx;
						lbm.u.y[n] = uy;
						lbm.u.z[n] = uz;
					}
				}
			}
			lbm.u.write_to_device();
			lbm.run(10u);
			println("HERE");

#ifdef SURFACE_EXPORT
			println("SURFACE_EXPORT");

			// Export surface at regular intervals
			if(surface_export_config.enabled && lbm.get_t() > 0u && surface_export_config.should_export(lbm.get_t())) {
			    println("DEBUG: Exporting surface frame at timestep " + to_string(lbm.get_t()));
				export_surface_frame(&lbm, lbm.get_t(), surface_export_config);
			}
#endif // SURFACE_EXPORT
		}
}

void main_setup() {
    main_setup_right_hander();
    // main_setup_sphere();
    // main_setup_beach();
    // main_setup_loading_stl_broken();
}



void main_setup_loading_stl_broken() { // breaking waves on beach; required extensions in defines.hpp: FP16S, VOLUME_FORCE, EQUILIBRIUM_BOUNDARIES, SURFACE, INTERACTIVE_GRAPHICS, SURFACE_EXPORT
	// Usage: ./FluidX3D --load-bathymetry path/to/bathymetry.stl
	// The domain size will be automatically set based on the STL bounding box
	// Use --autostart command line argument or uncomment the line below to start simulation automatically
#ifdef INTERACTIVE_GRAPHICS
	key_P = true; // autostart simulation in interactive graphics mode
#endif // INTERACTIVE_GRAPHICS

	/* WAVE GENERATION DEBUGGING SUMMARY:
	 *
	 * Comparing with working main_setup_beach() function (line 1285):
	 *
	 * WORKING BEACH EXAMPLE:
	 * - Uses TYPE_E cells at y=0 boundary
	 * - Direct parameters: u=0.12f, frequency=0.0007f, nu=0.01f, f=0.001f
	 * - Sets velocity on ALL cells at y=0 without checking types
	 * - Uses lbm.run(100u) in main loop
	 * - Runs lbm.run(0u) to initialize before main loop
	 *
	 * KEY FINDINGS:
	 * 1. TYPE_E cells ARE compatible with SURFACE extension
	 * 2. Equilibrium boundaries work by setting fi directly to feq in collision
	 * 3. The velocity we set is used to calculate feq each timestep
	 *
	 * DEBUGGING CHECKS ADDED:
	 * 1. Cell type counts before/after initialization
	 * 2. Velocity values being set and their persistence
	 * 3. Sample cell states showing type, rho, u, and phi
	 * 4. Velocity magnitude warnings for stability
	 * 5. TYPE_E cell creation tracking
	 *
	 * POTENTIAL ISSUES TO CHECK:
	 * 1. Are TYPE_E cells being converted during SURFACE initialization?
	 * 2. Is the velocity magnitude within stable limits (< 0.3)?
	 * 3. Are the unit conversions producing reasonable values?
	 * 4. Is the wave frequency too high/low for the domain?
	 *
	 * CURRENT IMPLEMENTATION:
	 * - USE_BEACH_PARAMS = true: Uses exact beach example values
	 * - USE_BEACH_PARAMS = false: Uses SI unit conversions
	 * - Following beach example structure exactly
	 * - TYPE_E cells replace TYPE_F at y=0 boundary
	 * - Velocity applied to ALL cells at y=0 each timestep
	 */

	// ################################################################## define simulation box size, viscosity and volume force ###################################################################
	// Physical parameters (SI units)
	const float si_g = 9.81f; // gravitational acceleration [m/s^2]
	const float si_nu = 1.0E-6f; // kinematic viscosity of water at 20°C [m^2/s]
	const float si_rho = 1000.0f; // density of water [kg/m^3]
	const float si_wave_height = 0.2f; // wave height [m]
	const float si_wave_period = 4.0f; // wave period [s]
	const float si_cell_size = 0.01f; // physical size of one lattice cell [m]

	// LBM parameters (dimensionless)
	// TEST: Use exact parameters from working beach example
	const bool USE_BEACH_PARAMS = false; // Set to false to use SI unit conversions

	const float lbm_u = USE_BEACH_PARAMS ? 0.12f : 0.15f; // peak velocity
	const float lbm_rho = 1.0f; // reference density in LBM units

	// Check for --load-bathymetry command line argument
	string bathymetry_file = "";
	for(size_t i = 0; i < main_arguments.size(); i++) {
		if(main_arguments[i] == "--load-bathymetry" && i + 1 < main_arguments.size()) {
			bathymetry_file = main_arguments[i + 1];
			break;
		}
	}

	// Default domain parameters
	uint Nx = 16u, Ny = 16u, Nz = 16u;

	if(!bathymetry_file.empty()) {
		// Load bathymetry STL file
		println("Loading bathymetry from: " + bathymetry_file);
		Mesh* bathymetry_mesh = read_stl(bathymetry_file, 1.0f); // Load without auto-scaling

		// Get bounding box of the bathymetry
		const float3 mesh_size = bathymetry_mesh->get_bounding_box_size();
		const float3 mesh_center = bathymetry_mesh->get_bounding_box_center();

		println("Bathymetry bounding box:");
		println("  Min: " + to_string(bathymetry_mesh->pmin.x) + ", " + to_string(bathymetry_mesh->pmin.y) + ", " + to_string(bathymetry_mesh->pmin.z));
		println("  Max: " + to_string(bathymetry_mesh->pmax.x) + ", " + to_string(bathymetry_mesh->pmax.y) + ", " + to_string(bathymetry_mesh->pmax.z));
		println("  Size: " + to_string(mesh_size.x) + " x " + to_string(mesh_size.y) + " x " + to_string(mesh_size.z));

		// Use reasonable domain sizing based on bathymetry scale
		const float padding = 0.5f; // Add 50% padding around the mesh
		const float3 domain_size = mesh_size * (1.0f + padding);

		// Set resolution based on a target of ~5 cells per unit in the largest dimension
		const float max_dim = fmax(fmax(mesh_size.x, mesh_size.y), mesh_size.z);
		const float cells_per_unit = 5.0f; // Fixed 5 cells per unit for consistent resolution

		// Calculate resolution maintaining aspect ratio
		// Nx = max(32u, (uint)(domain_size.x * cells_per_unit));
		// Ny = max(32u, (uint)(domain_size.y * cells_per_unit));
		// Nz = max(32u, (uint)(domain_size.z * cells_per_unit * 1.5f)); // Extra resolution in Z for waves

		// Cap maximum resolution to avoid excessive memory usage
		const uint max_total = 10000000u; // ~10M cells max
		const uint total_cells = Nx * Ny * Nz;
		if(total_cells > max_total) {
			const float scale = cbrt((float)max_total / (float)total_cells);
			Nx = max(64u, (uint)(Nx * scale));
			Ny = max(64u, (uint)(Ny * scale));
			Nz = max(32u, (uint)(Nz * scale));
		}

		println("Domain resolution: " + to_string(Nx) + " x " + to_string(Ny) + " x " + to_string(Nz));

		delete bathymetry_mesh; // Clean up temporary mesh
	}

	// Set up unit conversion based on physicacal cell size
	units.set_m_kg_s(1.0f, lbm_u, lbm_rho, si_cell_size, si_wave_height/si_wave_period, si_rho);

	// Convert physical parameters to LBM units
	const float lbm_nu = USE_BEACH_PARAMS ? 0.01f : units.nu(si_nu); // kinematic viscosity in LBM units
	const float lbm_g = USE_BEACH_PARAMS ? si_g : units.g(si_g); // gravitational acceleration in LBM units
	const float lbm_f = USE_BEACH_PARAMS ? 0.001f : units.f(si_rho, si_g); // force per volume in LBM units
	const float lbm_frequency = USE_BEACH_PARAMS ? 0.0007f : units.frequency(1.0f/si_wave_period); // wave frequency in LBM units

	// DEBUG: Print wave generation parameters
	println("Wave generation parameters:");
	println("  Mode: " + string(USE_BEACH_PARAMS ? "BEACH EXAMPLE PARAMS" : "SI UNIT CONVERSION"));
	println("  lbm_u (peak velocity): " + to_string(lbm_u));
	println("  lbm_frequency: " + to_string(lbm_frequency));
	println("  lbm_nu (viscosity): " + to_string(lbm_nu));
	println("  lbm_f (force): " + to_string(lbm_f));
	if(!USE_BEACH_PARAMS) {
		println("  si_wave_period: " + to_string(si_wave_period) + " s");
		println("  si_wave_height: " + to_string(si_wave_height) + " m");
	}
	println("Working beach example uses: u=0.12f, frequency=0.0007f, nu=0.01f, f=0.001f");

	// Create LBM with determined size and converted parameters
	LBM lbm(Nx, Ny, Nz, lbm_nu, 0.0f, 0.0f, -lbm_f);

	// If bathymetry file was provided, load and voxelize it
	if(!bathymetry_file.empty()) {
		Mesh* bathymetry_mesh = read_stl(bathymetry_file, 1.0f); // Reload the mesh

		// Center bathymetry in the domain
		// Center the bathymetry in the domain
		const float3 mesh_center = bathymetry_mesh->get_bounding_box_center();
		const float3 domain_center = lbm.center();

		// Translate bathymetry so its center aligns with domain center
		// But offset it down so it sits on the bottom part of the domain for waves
		const float3 target_center = float3(domain_center.x, domain_center.y, domain_center.z - 0.2f * lbm.size().z);
		const float3 mesh_offset = target_center - mesh_center;
		bathymetry_mesh->translate(mesh_offset);

		println("Translated bathymetry by: " + to_string(mesh_offset.x) + ", " + to_string(mesh_offset.y) + ", " + to_string(mesh_offset.z));

		// Voxelize the bathymetry as solid geometry
		lbm.voxelize_mesh_on_device(bathymetry_mesh, TYPE_S);
		delete bathymetry_mesh;
	}

	// ###################################################################################### define geometry ######################################################################################
	const uint lbm_Nx=lbm.get_Nx(), lbm_Ny=lbm.get_Ny(), lbm_Nz=lbm.get_Nz(); parallel_for(lbm.get_N(), [&](ulong n) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
		const uint H = lbm_Nz/2u;
		// CRITICAL ISSUE: This sets cells below H as TYPE_F, but later code tries to set
		// some of these same cells (at y=0) as TYPE_E for wave generation.
		// This creates a conflict where TYPE_F cells are overwritten as TYPE_E.
		if(lbm.flags[n]!=TYPE_S && z<H) { // Only set fluid if not already solid from bathymetry
			lbm.flags[n] = TYPE_F;
			lbm.rho[n] = units.rho_hydrostatic(lbm_f, (float)z, (float)H);
		}
		if(!bathymetry_file.empty()) {
			// When using bathymetry, only set boundaries on domain edges
			if(x==0u||x==lbm_Nx-1u||y==lbm_Ny-1u||z==0u||z==lbm_Nz-1u) {
				if(lbm.flags[n]!=TYPE_S) lbm.flags[n] = TYPE_S; // all non periodic except y=0
			}
			// Set equilibrium boundary at y=0 for wave generation (like in working beach example)
			if(y==0u && x>0u&&x<lbm_Nx-1u&&z>0u&&z<lbm_Nz-1u) {
				if(lbm.flags[n]!=TYPE_S) {
					lbm.flags[n] = TYPE_E; // equilibrium boundary for wave generation
					// Debug: Track TYPE_E cell creation
					if(x == lbm_Nx/2u && z == lbm_Nz/2u) {
						println("DEBUG: Setting TYPE_E at middle cell (x=" + to_string(x) + ",y=" + to_string(y) + ",z=" + to_string(z) + ")");
					}
				}
			}
		} else {
			// Default beach slope when no bathymetry is loaded
			if(plane(x, y, z, float3(lbm.center().x, 128.0f, 0.0f), float3(0.0f, -1.0f, 8.0f))) lbm.flags[n] = TYPE_S;
			if(x==0u||x==lbm_Nx-1u||y==0u||y==lbm_Ny-1u||z==0u||z==lbm_Nz-1u) lbm.flags[n] = TYPE_S; // all non periodic
			if(y==0u && x>0u&&x<lbm_Nx-1u&&z>0u&&z<lbm_Nz-1u) lbm.flags[n] = TYPE_E;
		}
	}); // ####################################################################### run simulation, export images and data ##########################################################################
#ifdef GRAPHICS
	lbm.graphics.visualization_modes = VIS_FLAG_LATTICE|VIS_PHI_RAYTRACE; // Show both surface mesh and raytracing
#endif // GRAPHICS

	#ifdef SURFACE_EXPORT
		// Parse surface export configuration from command line arguments
		SurfaceExportConfig surface_export_config = SurfaceExportConfig::parse_from_arguments(main_arguments);
		// if(surface_export_config.should_export(0u)) {
		// 	export_surface_frame(&lbm, 0u, surface_export_config);
		// }
#endif // SURFACE_EXPORT

	const ulong max_timesteps = 2000u; // run for 2000 timesteps to avoid crashes
	lbm.run(0u, max_timesteps); // initialize simulation



	while(true) { // main simulation loop (like beach example)
		// Export frames for visualization
#ifdef GRAPHICS
		if(lbm.graphics.next_frame(max_timesteps, 60.0f)) { // render frames at 60fps
			// Set camera to show the bathymetry and waves
			const float Nx_f = (float)lbm_Nx, Ny_f = (float)lbm_Ny, Nz_f = (float)lbm_Nz;
			// lbm.graphics.set_camera_free(float3(0.7f*Nx_f, -0.7f*Ny_f, 0.5f*Nz_f), -45.0f, 30.0f, 100.0f);
			lbm.graphics.write_frame("tmp/bathymetry_frames/");
		}
#endif // GRAPHICS

		// WAVE GENERATION CODE - Same approach as working beach example
		lbm.u.read_from_device();
		const float uy = lbm_u*sinf(2.0f*pif*lbm_frequency*(float)lbm.get_t());
		const float uz = 0.5f*lbm_u*cosf(2.0f*pif*lbm_frequency*(float)lbm.get_t());

		// Ensure we're not in the first timestep where SURFACE might still be initializing
		if(lbm.get_t() == 0u) {
			println("WARNING: Skipping wave velocity at timestep 0 to allow SURFACE initialization");
			lbm.run(100u); // run simulation step even at timestep 0
			continue;
		}

		// Debug output every 100 timesteps
		if(lbm.get_t() % 100u == 0u || lbm.get_t() < 5u) { // Also debug first few steps
			println("\n=== Wave generation at timestep " + to_string(lbm.get_t()) + " ===");
			println("  lbm_u = " + to_string(lbm_u) + ", lbm_frequency = " + to_string(lbm_frequency));
			println("  uy = " + to_string(uy) + ", uz = " + to_string(uz));

			// Check if TYPE_E cells exist at y=0
			lbm.flags.read_from_device();
			lbm.rho.read_from_device();
			lbm.phi.read_from_device();
			uint type_e_count = 0u;
			uint type_f_count = 0u;
			uint type_s_count = 0u;
			uint type_i_count = 0u;
			uint type_g_count = 0u;

			// Sample a few cells at y=0 to see detailed state
			println("  Sample cells at y=0:");
			for(uint sample_x = lbm_Nx/4u; sample_x < lbm_Nx; sample_x += lbm_Nx/4u) {
				for(uint sample_z = lbm_Nz/4u; sample_z < lbm_Nz; sample_z += lbm_Nz/4u) {
					if(sample_x > 0u && sample_x < lbm_Nx-1u && sample_z > 0u && sample_z < lbm_Nz-1u) {
						const uint n = sample_x + (0u + sample_z*lbm_Ny)*lbm_Nx;
						string type_str = "UNKNOWN";
						if(lbm.flags[n] == TYPE_E) type_str = "TYPE_E";
						else if(lbm.flags[n] == TYPE_F) type_str = "TYPE_F";
						else if(lbm.flags[n] == TYPE_S) type_str = "TYPE_S";
						else if(lbm.flags[n] == TYPE_I) type_str = "TYPE_I";
						else if(lbm.flags[n] == TYPE_G) type_str = "TYPE_G";
						println("    Cell[x=" + to_string(sample_x) + ",y=0,z=" + to_string(sample_z) + "]: "
							+ type_str + ", rho=" + to_string(lbm.rho[n])
							+ ", u=(" + to_string(lbm.u.x[n]) + "," + to_string(lbm.u.y[n]) + "," + to_string(lbm.u.z[n]) + ")"
							+ ", phi=" + to_string(lbm.phi[n]));
					}
				}
			}

			// Count all cell types at y=0
			for(uint z=1u; z<lbm_Nz-1u; z++) {
				for(uint x=1u; x<lbm_Nx-1u; x++) {
					const uint n = x+(0u+z*lbm_Ny)*lbm_Nx;
					if(lbm.flags[n] == TYPE_E) type_e_count++;
					else if(lbm.flags[n] == TYPE_F) type_f_count++;
					else if(lbm.flags[n] == TYPE_S) type_s_count++;
					else if(lbm.flags[n] == TYPE_I) type_i_count++;
					else if(lbm.flags[n] == TYPE_G) type_g_count++;
				}
			}
			println("  Cell counts at y=0: TYPE_E=" + to_string(type_e_count)
				+ ", TYPE_F=" + to_string(type_f_count)
				+ ", TYPE_S=" + to_string(type_s_count)
				+ ", TYPE_I=" + to_string(type_i_count)
				+ ", TYPE_G=" + to_string(type_g_count));

			// Check velocity after setting
			println("  Checking velocity BEFORE setting:");
			uint sample_n = lbm_Nx/2u + (0u + (lbm_Nz/2u)*lbm_Ny)*lbm_Nx;
			println("    Middle cell velocity: u.y=" + to_string(lbm.u.y[sample_n]) + ", u.z=" + to_string(lbm.u.z[sample_n]));
		}
		// Set velocity for ALL cells at y=0, just like beach example
		uint cells_updated = 0u;

		// Velocity magnitude check
		const float velocity_magnitude = sqrt(uy*uy + uz*uz);
		if(lbm.get_t() % 100u == 0u || lbm.get_t() < 5u) {
			println("  Wave velocity magnitude: " + to_string(velocity_magnitude));
			if(velocity_magnitude > 0.3f) {
				println("  WARNING: Velocity magnitude exceeds 0.3 (LBM stability limit)!");
			}
		}

		for(uint z=1u; z<lbm_Nz-1u; z++) {
			for(uint y=0u; y<1u; y++) {
				for(uint x=1u; x<lbm_Nx-1u; x++) {
					const uint n = x+(y+z*lbm_Ny)*lbm_Nx;
					lbm.u.y[n] = uy;
					lbm.u.z[n] = uz;
					cells_updated++;
				}
			}
		}

		// Debug: Check velocity after setting
		if(lbm.get_t() % 100u == 0u || lbm.get_t() < 5u) {
			println("  Updated " + to_string(cells_updated) + " cells with wave velocity");
			println("  Checking velocity AFTER setting:");
			uint sample_n = lbm_Nx/2u + (0u + (lbm_Nz/2u)*lbm_Ny)*lbm_Nx;
			println("    Middle cell velocity: u.y=" + to_string(lbm.u.y[sample_n]) + ", u.z=" + to_string(lbm.u.z[sample_n]));

			// Also check if velocity persists after write_to_device
			lbm.u.write_to_device();
			lbm.u.read_from_device();
			println("  After write/read cycle:");
			println("    Middle cell velocity: u.y=" + to_string(lbm.u.y[sample_n]) + ", u.z=" + to_string(lbm.u.z[sample_n]));
			println("=================================\n");
		} else {
			lbm.u.write_to_device();
		}
		// Use same run pattern as working beach example
		lbm.run(100u); // run 100 steps like the beach example

#ifdef SURFACE_EXPORT
		// Export surface at configured intervals
		if(surface_export_config.should_export(lbm.get_t())) {
			export_surface_frame(&lbm, lbm.get_t(), surface_export_config);
			println("Exported surface at timestep "+to_string(lbm.get_t()));
		}
#endif // SURFACE_EXPORT
	}
}

#ifdef BENCHMARK
#include "info.hpp"
void main_setup() { // benchmark; required extensions in defines.hpp: BENCHMARK, optionally FP16S or FP16C
	// ################################################################## define simulation box size, viscosity and volume force ###################################################################
	uint mlups = 0u; {

		//LBM lbm( 32u,  32u,  32u, 1.0f);
		//LBM lbm( 64u,  64u,  64u, 1.0f);
		//LBM lbm(128u, 128u, 128u, 1.0f);
		LBM lbm(256u, 256u, 256u, 1.0f); // default
		//LBM lbm(384u, 384u, 384u, 1.0f);
		//LBM lbm(512u, 512u, 512u, 1.0f);

		//const uint memory = 1488u; // memory occupation in MB (for multi-GPU benchmarks: make this close to as large as the GPU's VRAM capacity)
		//const uint3 lbm_N = (resolution(float3(1.0f, 1.0f, 1.0f), memory)/4u)*4u; // input: simulation box aspect ratio and VRAM occupation in MB, output: grid resolution
		//LBM lbm(1u*lbm_N.x, 1u*lbm_N.y, 1u*lbm_N.z, 1u, 1u, 1u, 1.0f); // 1 GPU
		//LBM lbm(2u*lbm_N.x, 1u*lbm_N.y, 1u*lbm_N.z, 2u, 1u, 1u, 1.0f); // 2 GPUs
		//LBM lbm(2u*lbm_N.x, 2u*lbm_N.y, 1u*lbm_N.z, 2u, 2u, 1u, 1.0f); // 4 GPUs
		//LBM lbm(2u*lbm_N.x, 2u*lbm_N.y, 2u*lbm_N.z, 2u, 2u, 2u, 1.0f); // 8 GPUs

		// #########################################################################################################################################################################################
		for(uint i=0u; i<1000u; i++) {
			lbm.run(10u, 1000u*10u);
			mlups = max(mlups, to_uint((double)lbm.get_N()*1E-6/info.runtime_lbm_timestep_smooth));
		}
	} // make lbm object go out of scope to free its memory
	print_info("Peak MLUPs/s = "+to_string(mlups));
#if defined(_WIN32)
	wait();
#endif // Windows
} /**/
#endif // BENCHMARK
