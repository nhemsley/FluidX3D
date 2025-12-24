#pragma once
#include "utilities.hpp"
#include "lbm.hpp"
#include <string>
#include <memory>
#include <atomic>
#include <chrono>

// Include the generated C header from seaview-network
extern "C" {
    #include "../vendor/seaview/crates/seaview-network/include/seaview_network.h"
}

// Error codes for compatibility
const int SEAVIEW_SUCCESS = 0;
const int SEAVIEW_ERROR_CONNECTION = -1;
const int SEAVIEW_ERROR_SERIALIZATION = -2;
const int SEAVIEW_ERROR_INVALID_DATA = -3;

// Colored output helpers for graceful failure messages
inline void print_warning_yellow(const string& message) {
    print("\033[33m⚠ WARNING: " + message + "\033[0m");
}

inline void print_error_red(const string& message) {
    print("\033[31m✗ ERROR: " + message + "\033[0m");
}

inline void print_success_green(const string& message) {
    print("\033[32m✓ " + message + "\033[0m");
}

/**
 * MeshStreamer - Real-time mesh streaming to seaview visualization
 *
 * This class provides real-time streaming of FluidX3D surface meshes to the
 * seaview visualization application via TCP network protocol.
 *
 * Features:
 * - Zero-copy mesh data transfer where possible
 * - Automatic frame numbering and timestamping
 * - Connection management with retry logic
 * - Performance monitoring and statistics
 * - Thread-safe operation
 */
class MeshStreamer {
private:
    // Network connection
    NetworkSender* sender = nullptr;
    std::string host;
    uint16_t port;
    std::atomic<bool> connected{false};

    // Simulation metadata
    std::string simulation_id;
    std::atomic<uint32_t> frame_counter{0};

    // Domain bounds (updated from LBM)
    float domain_min[3] = {0.0f, 0.0f, 0.0f};
    float domain_max[3] = {1.0f, 1.0f, 1.0f};

    // Performance statistics
    std::atomic<uint64_t> total_frames_sent{0};
    std::atomic<uint64_t> total_bytes_sent{0};
    std::atomic<uint64_t> total_send_time_ms{0};
    std::chrono::steady_clock::time_point last_stats_print;

    // Configuration
    bool enabled = false;
    bool compute_normals = true;
    bool verbose_logging = false;
    uint32_t stats_interval = 100; // Print stats every N frames

    // Error handling
    std::atomic<int> last_error{SEAVIEW_SUCCESS};
    std::atomic<uint32_t> error_count{0};
    std::atomic<bool> connection_failed_notified{false};
    std::chrono::steady_clock::time_point last_connection_attempt;
    uint32_t connection_retry_interval = 30; // seconds

public:
    /**
     * Constructor
     * @param sim_id Unique simulation identifier (human-readable)
     * @param hostname Target seaview host (default: "localhost")
     * @param port_num Target seaview port (default: 15703)
     */
    MeshStreamer(const std::string& sim_id = "fluidx3d",
                 const std::string& hostname = "localhost",
                 uint16_t port_num = 15703)
        : host(hostname), port(port_num), simulation_id(sim_id) {
        last_stats_print = std::chrono::steady_clock::now();

        // Generate unique simulation ID with timestamp if default
        if (simulation_id == "fluidx3d") {
            auto now = std::chrono::system_clock::now();
            auto time_t = std::chrono::system_clock::to_time_t(now);
            simulation_id = "fluidx3d_" + std::to_string(time_t);
        }

        println("MeshStreamer initialized:");
        println("  Simulation ID: " + simulation_id);
        println("  Target: " + host + ":" + std::to_string(port));
    }

    /**
     * Destructor - ensures clean disconnection
     */
    ~MeshStreamer() {
        disconnect();
        print_final_statistics();
    }

    /**
     * Connect to seaview server with graceful failure handling
     * @return true if connection successful
     */
    bool connect() {
        if (connected.load()) {
            if (verbose_logging) println("MeshStreamer: Already connected");
            return true;
        }

        // Check if we should retry connection (rate limiting)
        auto now = std::chrono::steady_clock::now();
        if (connection_failed_notified.load()) {
            auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(now - last_connection_attempt);
            if (elapsed.count() < connection_retry_interval) {
                return false; // Too soon to retry
            }
        }

        if (verbose_logging) println("MeshStreamer: Attempting connection to " + host + ":" + std::to_string(port));

        sender = seaview_network_create_sender(host.c_str(), port);
        last_connection_attempt = now;

        if (sender != nullptr) {
            connected.store(true);
            connection_failed_notified.store(false);
            print_success_green("MeshStreamer connected to seaview at " + host + ":" + std::to_string(port));
            return true;
        } else {
            last_error.store(SEAVIEW_ERROR_CONNECTION);
            error_count.fetch_add(1);

            // Only print warning once, then be silent until retry interval
            if (!connection_failed_notified.load()) {
                print_warning_yellow("Seaview server not running at " + host + ":" + std::to_string(port) +
                                     ". Simulation will continue without streaming. Will retry every " +
                                     std::to_string(connection_retry_interval) + " seconds.");
                connection_failed_notified.store(true);
            }
            return false;
        }
    }

    /**
     * Disconnect from seaview server
     */
    void disconnect() {
        if (sender != nullptr) {
            seaview_network_destroy_sender(sender);
            sender = nullptr;
        }
        connected.store(false);
        if (verbose_logging) println("MeshStreamer: Disconnected");
    }

    /**
     * Set domain bounds for coordinate system
     * @param min_bounds Minimum domain coordinates [x,y,z]
     * @param max_bounds Maximum domain coordinates [x,y,z]
     */
    void set_domain_bounds(const float min_bounds[3], const float max_bounds[3]) {
        for (int i = 0; i < 3; i++) {
            domain_min[i] = min_bounds[i];
            domain_max[i] = max_bounds[i];
        }

        if (verbose_logging) {
            println("MeshStreamer: Domain bounds updated:");
            println("  Min: [" + std::to_string(domain_min[0]) + ", " +
                    std::to_string(domain_min[1]) + ", " + std::to_string(domain_min[2]) + "]");
            println("  Max: [" + std::to_string(domain_max[0]) + ", " +
                    std::to_string(domain_max[1]) + ", " + std::to_string(domain_max[2]) + "]");
        }
    }

    /**
     * Set domain bounds from LBM domain
     * @param lbm LBM instance to extract bounds from
     */
    void set_domain_bounds_from_lbm(const LBM* lbm) {
        if (lbm == nullptr) return;

        // Extract domain size from LBM
        // Note: This assumes uniform domain size - may need adjustment for multi-GPU
        const uint Nx = lbm->get_Nx();
        const uint Ny = lbm->get_Ny();
        const uint Nz = lbm->get_Nz();

        float min_bounds[3] = {0.0f, 0.0f, 0.0f};
        float max_bounds[3] = {static_cast<float>(Nx), static_cast<float>(Ny), static_cast<float>(Nz)};

        set_domain_bounds(min_bounds, max_bounds);
    }

    /**
     * Send mesh frame to seaview
     * @param vertices Vertex data (x,y,z interleaved)
     * @param triangle_count Number of triangles
     * @param normals Optional normal data (can be nullptr)
     * @return true if sent successfully
     */
    bool send_mesh_frame(const float* vertices, ulong triangle_count, const float* normals = nullptr) {
        if (!enabled) return true; // Silent success when disabled

        if (vertices == nullptr || triangle_count == 0) {
            if (verbose_logging) println("MeshStreamer: No mesh data to send");
            return true;
        }

        // Auto-connect if not connected (graceful failure)
        if (!connected.load()) {
            if (!connect()) {
                // Connection failed, but continue simulation gracefully
                if (verbose_logging && !connection_failed_notified.load()) {
                    println("MeshStreamer: Skipping frame " + std::to_string(frame_counter.load()) +
                           " due to connection failure");
                }
                return true; // Return success to continue simulation
            }
        }

        auto start_time = std::chrono::high_resolution_clock::now();

        // Prepare mesh frame
        CMeshFrame mesh_frame = {
            .simulation_id = simulation_id.c_str(),
            .frame_number = frame_counter.fetch_add(1),
            .timestamp = static_cast<uint64_t>(std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::system_clock::now().time_since_epoch()).count()),
            .domain_min = {domain_min[0], domain_min[1], domain_min[2]},
            .domain_max = {domain_max[0], domain_max[1], domain_max[2]},
            .vertex_count = triangle_count * 3, // 3 vertices per triangle
            .vertices = vertices,
            .normals = normals,
            .index_count = 0,
            .indices = nullptr
        };

        // Send mesh
        int result = seaview_network_send_mesh(sender, &mesh_frame);

        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);

        // Update statistics
        total_frames_sent.fetch_add(1);
        size_t frame_bytes = triangle_count * 9 * sizeof(float); // 9 floats per triangle
        if (normals != nullptr) {
            frame_bytes += triangle_count * 9 * sizeof(float); // Add normal data
        }
        total_bytes_sent.fetch_add(frame_bytes);
        total_send_time_ms.fetch_add(duration.count() / 1000);

        if (result == SEAVIEW_SUCCESS) {
            if (verbose_logging) {
                println("MeshStreamer: Sent frame " + std::to_string(mesh_frame.frame_number) +
                        " (" + std::to_string(triangle_count) + " triangles, " +
                        std::to_string(duration.count()) + "μs)");
            }

            // Print periodic statistics
            if (total_frames_sent.load() % stats_interval == 0) {
                print_statistics();
            }

            return true;
        } else {
            last_error.store(result);
            error_count.fetch_add(1);

            print_warning_yellow("MeshStreamer: Failed to send frame " +
                                 std::to_string(mesh_frame.frame_number) +
                                 " (error: " + std::to_string(result) +
                                 "). Simulation continues without streaming.");

            // Disconnect on error to trigger reconnect attempt
            disconnect();
            connection_failed_notified.store(false); // Allow immediate retry attempt
            return true; // Return success to continue simulation gracefully
        }
    }

    /**
     * Send mesh from LBM surface export
     * @param lbm LBM instance to export from
     * @return true if sent successfully
     */
    bool send_mesh_from_lbm(LBM* lbm) {
        if (lbm == nullptr) {
            print_warning("MeshStreamer: LBM instance is null");
            return false;
        }

#if !defined(SURFACE) || !defined(SURFACE_EXPORT)
        print_warning_yellow("MeshStreamer: SURFACE and SURFACE_EXPORT extensions required. Skipping mesh streaming.");
        return true; // Continue simulation gracefully
#else
        // Export surface from LBM
        lbm->enqueue_export_surface();
        float* vertices = lbm->get_surface_vertices();
        ulong triangle_count = lbm->get_triangle_count();

        if (triangle_count == 0 || vertices == nullptr) {
            if (verbose_logging) println("MeshStreamer: No surface data from LBM (empty frame)");
            return true; // Normal case - continue gracefully
        }

        // Update domain bounds from LBM if not set
        if (domain_min[0] == 0.0f && domain_max[0] == 1.0f) {
            set_domain_bounds_from_lbm(lbm);
        }

        // Note: LBM surface export doesn't provide normals separately
        // Normals are computed in the write_stl functions, but not stored
        return send_mesh_frame(vertices, triangle_count, nullptr);
#endif
    }

    /**
     * Configuration methods
     */
    void enable() {
        enabled = true;
        println("MeshStreamer: Enabled");
    }

    void disable() {
        enabled = false;
        println("MeshStreamer: Disabled");
    }

    bool is_enabled() const { return enabled; }
    bool is_connected() const { return connected.load(); }

    void set_verbose(bool verbose) {
        verbose_logging = verbose;
        if (verbose) println("MeshStreamer: Verbose logging enabled");
    }

    void set_connection_retry_interval(uint32_t seconds) {
        connection_retry_interval = seconds;
        if (verbose_logging) println("MeshStreamer: Connection retry interval set to " + std::to_string(seconds) + "s");
    }

    void set_stats_interval(uint32_t interval) {
        stats_interval = interval;
        if (verbose_logging) println("MeshStreamer: Stats interval set to " + std::to_string(interval));
    }

    /**
     * Statistics and monitoring
     */
    void print_statistics() {
        auto now = std::chrono::steady_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(now - last_stats_print);

        uint64_t frames = total_frames_sent.load();
        uint64_t bytes = total_bytes_sent.load();
        uint64_t send_time = total_send_time_ms.load();
        uint32_t errors = error_count.load();

        if (frames > 0 || errors > 0) {
            println("MeshStreamer Statistics:");
            println("  Frames sent: " + std::to_string(frames));
            if (frames > 0) {
                println("  Data sent: " + std::to_string(bytes / (1024*1024)) + " MB");
                double avg_frame_time = static_cast<double>(send_time) / frames;
                double throughput_mbps = (bytes * 8.0) / (elapsed.count() * 1024.0 * 1024.0);
                println("  Avg frame time: " + std::to_string(avg_frame_time) + " ms");
                println("  Throughput: " + std::to_string(throughput_mbps) + " Mbps");
            }

            if (errors > 0) {
                print_warning_yellow("  Connection errors: " + std::to_string(errors));
            }

            string status = connected.load() ? "Connected" :
                           (connection_failed_notified.load() ? "Disconnected (will retry)" : "Connecting...");
            println("  Status: " + status);
        }

        last_stats_print = now;
    }

    void print_final_statistics() {
        println("MeshStreamer Final Statistics:");
        print_statistics();
    }

    /**
     * Get error information
     */
    int get_last_error() const { return last_error.load(); }
    uint32_t get_error_count() const { return error_count.load(); }
    uint64_t get_frames_sent() const { return total_frames_sent.load(); }
    uint64_t get_bytes_sent() const { return total_bytes_sent.load(); }
};

/**
 * Streaming configuration for command-line integration
 */
struct MeshStreamingConfig {
    std::string host = "localhost";
    uint16_t port = 15703;
    std::string simulation_id = "";
    bool enabled = false;
    bool verbose = false;
    uint32_t stream_interval = 1; // Stream every N timesteps
    uint32_t stats_interval = 100;

    /**
     * Parse configuration from command line arguments
     */
    static MeshStreamingConfig parse_from_arguments(const vector<string>& args) {
        MeshStreamingConfig config;
        config.parse_arguments(args);
        return config;
    }

    void parse_arguments(const vector<string>& args) {
        for (size_t i = 0; i < args.size(); i++) {
            if (args[i] == "--stream-to-seaview" && i + 1 < args.size()) {
                // Parse host:port format
                string target = args[i + 1];
                size_t colon_pos = target.find(':');
                if (colon_pos != string::npos) {
                    host = target.substr(0, colon_pos);
                    port = static_cast<uint16_t>(std::stoi(target.substr(colon_pos + 1)));
                } else {
                    host = target;
                }
                enabled = true;
                println("Mesh streaming enabled to " + host + ":" + std::to_string(port));
            }
            else if (args[i] == "--stream-simulation-id" && i + 1 < args.size()) {
                simulation_id = args[i + 1];
                println("Simulation ID set to: " + simulation_id);
            }
            else if (args[i] == "--stream-interval" && i + 1 < args.size()) {
                stream_interval = to_uint(args[i + 1]);
                println("Stream interval set to: " + std::to_string(stream_interval));
            }
            else if (args[i] == "--stream-verbose") {
                verbose = true;
                println("Verbose streaming enabled");
            }
            else if (args[i] == "--stream-stats-interval" && i + 1 < args.size()) {
                stats_interval = to_uint(args[i + 1]);
                println("Stats interval set to: " + std::to_string(stats_interval));
            }
        }
    }

    /**
     * Check if streaming should happen at this timestep
     */
    bool should_stream(ulong timestep) const {
        return enabled && (timestep % stream_interval == 0);
    }
};

// Global streaming configuration
extern MeshStreamingConfig mesh_streaming_config;

/**
 * Stream mesh frame using global configuration
 * @param lbm LBM instance to stream from
 * @param timestep Current simulation timestep
 * @param streamer MeshStreamer instance (will be created if needed)
 */
inline void stream_mesh_frame(LBM* lbm, ulong timestep, std::unique_ptr<MeshStreamer>& streamer) {
    if (!mesh_streaming_config.enabled || !lbm) return;
    if (!mesh_streaming_config.should_stream(timestep)) return;

    // println("DEBUG: stream_mesh_frame called for timestep " + to_string(timestep));

    // Create streamer if needed
    if (!streamer) {
        // println("DEBUG: Creating new MeshStreamer for " + mesh_streaming_config.host + ":" + to_string(mesh_streaming_config.port));
        string sim_id = mesh_streaming_config.simulation_id;
        if (sim_id.empty()) {
            sim_id = "fluidx3d";
        }

        streamer = std::make_unique<MeshStreamer>(
            sim_id,
            mesh_streaming_config.host,
            mesh_streaming_config.port
        );

        streamer->set_verbose(mesh_streaming_config.verbose);
        streamer->set_stats_interval(mesh_streaming_config.stats_interval);
        streamer->enable();
    }

    // Stream the mesh (graceful failure - simulation continues regardless)
    // println("DEBUG: About to call send_mesh_from_lbm");
    streamer->send_mesh_from_lbm(lbm);
    // println("DEBUG: send_mesh_from_lbm completed");
}
