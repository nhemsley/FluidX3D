# Auto-Starting Simulations in FluidX3D

When using FluidX3D with `INTERACTIVE_GRAPHICS` or `INTERACTIVE_GRAPHICS_ASCII` mode, the simulation starts in a paused state by default. You need to press the <kbd>P</kbd> key to start the simulation. This guide shows several methods to automatically start the simulation without manual intervention.

## Method 1: Command Line Argument (Recommended)

The cleanest way is to use the `--autostart` command line argument:

```bash
./FluidX3D --autostart
```

This can be combined with other arguments:

```bash
# With bathymetry loading
./FluidX3D --autostart --load-bathymetry seafloor.stl

# With surface export
./FluidX3D --autostart --export-surface-to ./output/ --export-surface-interval 100

# All combined
./FluidX3D --autostart --load-bathymetry seafloor.stl --export-surface-to ./output/ --export-surface-interval 100
```

## Method 2: Modify Source Code

### Option A: Change Default Initialization

Edit `src/graphics.cpp` and change the `key_P` initialization from `false` to `true`:

```cpp
// Line 12 in src/graphics.cpp
bool key_E=false, key_G=false, key_H=false, key_O=false, key_P=true, key_Q=false, key_T=false, key_Z=false;
//                                                            ^^^^^^^^ changed from false to true
```

This will make ALL simulations start automatically.

### Option B: Set in main() Functions

Uncomment the `key_P = true;` lines in the main() functions in `src/graphics.cpp`:

```cpp
// Around lines 521, 701, 804, and 842
// Alternative: Always autostart by uncommenting the line below
key_P = true;  // <-- uncomment this line
```

## Method 3: Per-Simulation Setup

Add `key_P = true;` at the beginning of your `main_setup()` function:

```cpp
void main_setup() {
    #ifdef INTERACTIVE_GRAPHICS
    key_P = true; // autostart this specific simulation
    #endif
    
    // rest of your setup code...
    LBM lbm(256u, 256u, 256u, 0.02f);
    // ...
}
```

This method is useful when you want some simulations to autostart but not others.

## Method 4: Using a Define

Add a custom define in `src/defines.hpp`:

```cpp
#define AUTOSTART_SIMULATION // uncomment to autostart simulations
```

Then in your setup:

```cpp
void main_setup() {
    #if defined(INTERACTIVE_GRAPHICS) && defined(AUTOSTART_SIMULATION)
    key_P = true;
    #endif
    
    // rest of setup...
}
```

## Comparison of Methods

| Method | Scope | Persistence | Flexibility |
|--------|-------|-------------|-------------|
| Command line argument | Per run | No | High - can choose per execution |
| Source code default | All simulations | Yes | Low - affects all runs |
| Per-simulation setup | Specific simulation | Yes | Medium - per setup function |
| Define option | Compile-time | Yes | Medium - requires recompilation |

## Best Practices

1. **For Development**: Use the command line argument `--autostart` for maximum flexibility
2. **For Production**: Set `key_P = true` in your specific `main_setup()` function
3. **For Batch Processing**: Use the command line argument in your scripts
4. **For Demos**: Modify the source code default to always autostart

## Example Use Cases

### Automated Testing Script
```bash
#!/bin/bash
for stl in *.stl; do
    ./FluidX3D --autostart --load-bathymetry "$stl" --export-surface-to "./results/$stl/"
done
```

### Headless Rendering
When using `GRAPHICS` mode (not interactive) for video rendering, the simulation always runs automatically without needing any modifications.

### Interactive Development
During development, you might want to start paused to set up the view:
```bash
./FluidX3D  # starts paused, press P when ready
```

## Notes

- The `--autostart` argument only affects `INTERACTIVE_GRAPHICS` and `INTERACTIVE_GRAPHICS_ASCII` modes
- In non-interactive `GRAPHICS` mode, simulations always run automatically
- You can still pause/unpause with <kbd>P</kbd> even when autostarted
- The autostart state is just the initial state - all interactive controls remain functional