#pragma once
#include "utilities.hpp"
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <sys/stat.h>
#ifdef _WIN32
#include <direct.h>
#else
#include <sys/types.h>
#endif

// Function to write surface mesh data to binary STL file
inline void write_stl_binary(const string& filename, const float* vertices, const ulong triangle_count) {
	if(triangle_count == 0u) {
		print_warning("No triangles to export to STL file.");
		return;
	}

	std::ofstream file(filename, std::ios::binary);
	if(!file.is_open()) {
		print_error("Could not open file \"" + filename + "\" for writing.");
		return;
	}

	// STL binary header (80 bytes)
	char header[80];
	memset(header, 0, 80);
	snprintf(header, 80, "FluidX3D surface export - %lu triangles", triangle_count);
	file.write(header, 80);

	// Number of triangles (4 bytes)
	uint32_t num_triangles = (uint32_t)triangle_count;
	file.write(reinterpret_cast<const char*>(&num_triangles), sizeof(uint32_t));

	// Write each triangle
	for(ulong i = 0u; i < triangle_count; i++) {
		const ulong base = i * 9u; // 9 floats per triangle (3 vertices * 3 coordinates)

		// Get triangle vertices
		float3 v0 = float3(vertices[base + 0u], vertices[base + 1u], vertices[base + 2u]);
		float3 v1 = float3(vertices[base + 3u], vertices[base + 4u], vertices[base + 5u]);
		float3 v2 = float3(vertices[base + 6u], vertices[base + 7u], vertices[base + 8u]);

		// Calculate normal (right-hand rule)
		float3 edge1 = v1 - v0;
		float3 edge2 = v2 - v0;
		float3 normal = normalize(cross(edge1, edge2));

		// Write normal (12 bytes)
		file.write(reinterpret_cast<const char*>(&normal.x), sizeof(float));
		file.write(reinterpret_cast<const char*>(&normal.y), sizeof(float));
		file.write(reinterpret_cast<const char*>(&normal.z), sizeof(float));

		// Write vertices (36 bytes)
		file.write(reinterpret_cast<const char*>(&v0.x), sizeof(float));
		file.write(reinterpret_cast<const char*>(&v0.y), sizeof(float));
		file.write(reinterpret_cast<const char*>(&v0.z), sizeof(float));

		file.write(reinterpret_cast<const char*>(&v1.x), sizeof(float));
		file.write(reinterpret_cast<const char*>(&v1.y), sizeof(float));
		file.write(reinterpret_cast<const char*>(&v1.z), sizeof(float));

		file.write(reinterpret_cast<const char*>(&v2.x), sizeof(float));
		file.write(reinterpret_cast<const char*>(&v2.y), sizeof(float));
		file.write(reinterpret_cast<const char*>(&v2.z), sizeof(float));

		// Attribute byte count (2 bytes) - always 0 for standard STL
		uint16_t attribute_count = 0;
		file.write(reinterpret_cast<const char*>(&attribute_count), sizeof(uint16_t));
	}

	file.close();
	println("Exported " + to_string(triangle_count) + " triangles to \"" + filename + "\"");
}

// Function to write surface mesh data to ASCII STL file (for debugging/inspection)
inline void write_stl_ascii(const string& filename, const float* vertices, const ulong triangle_count) {
	if(triangle_count == 0u) {
		print_warning("No triangles to export to STL file.");
		return;
	}

	std::ofstream file(filename);
	if(!file.is_open()) {
		print_error("Could not open file \"" + filename + "\" for writing.");
		return;
	}

	// Set precision for float output
	file << std::fixed << std::setprecision(6);

	// Write header
	file << "solid FluidX3D_surface" << std::endl;

	// Write each triangle
	for(ulong i = 0u; i < triangle_count; i++) {
		const ulong base = i * 9u; // 9 floats per triangle

		// Get triangle vertices
		float3 v0 = float3(vertices[base + 0u], vertices[base + 1u], vertices[base + 2u]);
		float3 v1 = float3(vertices[base + 3u], vertices[base + 4u], vertices[base + 5u]);
		float3 v2 = float3(vertices[base + 6u], vertices[base + 7u], vertices[base + 8u]);

		// Calculate normal
		float3 edge1 = v1 - v0;
		float3 edge2 = v2 - v0;
		float3 normal = normalize(cross(edge1, edge2));

		// Write facet
		file << "  facet normal " << normal.x << " " << normal.y << " " << normal.z << std::endl;
		file << "    outer loop" << std::endl;
		file << "      vertex " << v0.x << " " << v0.y << " " << v0.z << std::endl;
		file << "      vertex " << v1.x << " " << v1.y << " " << v1.z << std::endl;
		file << "      vertex " << v2.x << " " << v2.y << " " << v2.z << std::endl;
		file << "    endloop" << std::endl;
		file << "  endfacet" << std::endl;
	}

	// Write footer
	file << "endsolid FluidX3D_surface" << std::endl;

	file.close();
	println("Exported " + to_string(triangle_count) + " triangles to \"" + filename + "\" (ASCII format)");
}

// Wrapper function that writes binary STL by default
inline void write_stl(const string& filename, const float* vertices, const ulong triangle_count, bool ascii = false) {
	if(ascii) {
		write_stl_ascii(filename, vertices, triangle_count);
	} else {
		write_stl_binary(filename, vertices, triangle_count);
	}
}

// Function to create directory if it doesn't exist
inline bool create_directory_if_not_exists(const string& path) {
	if(path.empty()) return true;

	// Check if directory exists
	struct stat info;
	if(stat(path.c_str(), &info) == 0) {
		if(info.st_mode & S_IFDIR) {
			return true; // Directory already exists
		} else {
			print_error("Path \"" + path + "\" exists but is not a directory.");
			return false;
		}
	}

	// Create directory
#ifdef _WIN32
	if(_mkdir(path.c_str()) != 0) {
#else
	if(mkdir(path.c_str(), 0755) != 0) {
#endif
		print_error("Failed to create directory \"" + path + "\".");
		return false;
	}

	println("Created directory \"" + path + "\"");
	return true;
}

// Structure to hold surface export configuration
struct SurfaceExportConfig {
	string directory = "";
	bool enabled = false;
	uint export_interval = 100u; // Export every N timesteps
	bool ascii_format = false; // Use ASCII STL format (larger files, but human-readable)

	// Parse command line arguments for surface export
	void parse_arguments(const vector<string>& args) {
		for(size_t i = 0; i < args.size(); i++) {
			if(args[i] == "--export-surface-to" && i + 1 < args.size()) {
				directory = args[i + 1];
				enabled = true;

				// Ensure directory ends with separator
				if(!directory.empty() && directory.back() != '/' && directory.back() != '\\') {
					directory += "/";
				}

				// Create directory if needed
				if(!create_directory_if_not_exists(directory)) {
					enabled = false;
					print_error("Surface export disabled due to directory creation failure.");
				}
			}
			else if(args[i] == "--export-surface-interval" && i + 1 < args.size()) {
				export_interval = to_uint(args[i + 1]);
			}
			else if(args[i] == "--export-surface-ascii") {
				ascii_format = true;
			}
		}
	}

	// Check if export should happen at this timestep
	bool should_export(ulong timestep) const {
		return enabled && (timestep % export_interval == 0u);
	}

	// Generate filename for a given timestep
	string get_filename(ulong timestep) const {
		std::stringstream ss;
		ss << directory << "surface_" << std::setfill('0') << std::setw(6) << timestep << ".stl";
		return ss.str();
	}
};
