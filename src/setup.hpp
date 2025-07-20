#pragma once

#include "defines.hpp"
#include "lbm.hpp"
#include "shapes.hpp"
#include "surface_export.hpp"
#include "graphics.hpp"

void main_setup(); // main setup script

// #ifdef SURFACE_EXPORT
// void run_with_surface_export(LBM* lbm, const ulong steps=max_ulong); // run simulation with periodic surface export
// void export_surface_frame(LBM* lbm, ulong timestep); // export a single surface frame to STL file

// #endif // SURFACE_EXPORT
