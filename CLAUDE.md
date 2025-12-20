# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a Rust-based computational mathematics project focused on discrete homology theory. The codebase implements graph structures optimized for topological computations, specifically using CSR (Compressed Sparse Row) representation for efficient memory usage and traversal.

## Build and Development Commands

```bash
# Build the project
cargo build

# Run the project
cargo run

# Run tests
cargo test

# Run tests with output
cargo test -- --nocapture

# Build in release mode (optimized)
cargo build --release

# Check code without building
cargo check

# Format code
cargo fmt

# Run clippy linter
cargo clippy
```

## Architecture

### Module Structure

- **src/main.rs**: Entry point; contains top-level API like `n_cube()` for generating n-dimensional hypercubes
- **src/graph.rs**: Core graph data structures and algorithms

### Graph Abstraction (src/graph.rs)

The project uses a trait-based design for undirected graphs:

- **`UGraph` trait**: Defines the interface for undirected graphs with `neighbors(&self, v: u32)` and `n(&self)` methods
- **`CSRGraph` struct**: Implementation using Compressed Sparse Row format
  - `offsets: Vec<u32>`: Prefix sum array where `offsets[i]` marks where vertex i's adjacency list starts
  - `neighbors: Vec<u32>`: Flattened adjacency lists
  - Neighbors are stored sorted in ascending order for deterministic iteration

### Graph Constraints

The `CSRGraph::try_from(AdjMatrix)` constructor enforces:
1. **Reflexivity**: Self-loops required (diagonal must be true)
2. **Symmetry**: Adjacency matrix must be symmetric (undirected graph)
3. **Squareness**: Matrix must be n×n

These constraints suggest the graph represents a relation that must be reflexive and symmetric, typical for discrete topological structures.

### Type Aliases

- `AdjMatrix = Vec<Vec<bool>>`: Adjacency matrix representation (row-major)

## Implementation Notes

### CSR Construction Details

- The `TryFrom<AdjMatrix>` implementation:
  1. Validates matrix is square, symmetric, and reflexive
  2. Computes degree of each vertex
  3. Builds offset array via prefix sum
  4. Fills neighbors array by iterating upper triangle only (since symmetric)
  5. Sorts each adjacency list for deterministic ordering

### Incomplete Features

- `From<CSRGraph> for AdjMatrix` is currently unimplemented (panics)
- The `n_cube()` function in main.rs is defined but empty

### Future Development

The project is in early stages - main only prints "Hello, world!" and key functions like `n_cube()` are stubs. The focus appears to be building graph infrastructure for discrete homology computations on simplicial complexes or cell complexes.
