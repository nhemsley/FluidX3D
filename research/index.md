# FluidX3D Research Documentation Index

This directory contains research documentation for various aspects of FluidX3D development and implementation challenges.

## Available Research Documents

### [Surface Export Implementation Research](surface-export.md)
Analysis of the current state of the surface export feature in FluidX3D, including recommendations for completing the marching cubes implementation. Covers both single-GPU and multi-GPU considerations for exporting fluid surface meshes.

**Key Topics:**
- Current implementation status
- Marching cubes algorithm integration
- Phi field and surface detection
- Implementation recommendations
- Testing strategies

### [Multi-GPU Surface Export Research](multi-gpu-surface-export.md)
Detailed analysis of challenges and solutions for implementing surface export across multiple GPUs in FluidX3D's domain-decomposed architecture.

**Key Topics:**
- FluidX3D multi-GPU architecture overview
- Domain decomposition and halo regions
- Boundary handling challenges
- Proposed solutions including per-GPU export
- Implementation recommendations and testing strategies

## Quick Links

- **Single-GPU Surface Export**: [surface-export.md](surface-export.md)
- **Multi-GPU Surface Export**: [multi-gpu-surface-export.md](multi-gpu-surface-export.md)

## Research Categories

### Performance & Scalability
- [Multi-GPU Surface Export](multi-gpu-surface-export.md) - Handling domain decomposition for surface extraction

### Feature Implementation
- [Surface Export Implementation](surface-export.md) - Completing the marching cubes surface export feature

## Contributing

When adding new research documents:
1. Place the document in this directory with a descriptive filename
2. Add an entry to this index with a brief description
3. Include key topics covered in the document
4. Update the relevant category section