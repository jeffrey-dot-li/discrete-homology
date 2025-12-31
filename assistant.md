# Assistant Guide for Discrete Homology Project

## Project Context
This is a Rust project focused on discrete homology theory computations. It implements graph structures (CSR, Hypercube) and algorithms for topological analysis.

## Development Environment
- **OS:** MacOS (darwin)
- **Language:** Rust (Edition 2021)
- **Build System:** Cargo

## Key Commands
- `cargo build`: Build the project.
- `cargo test`: Run all tests.
- `cargo test -- --nocapture`: Run tests showing stdout.
- `cargo bench`: Run benchmarks (using criterion).
- `cargo fmt`: Format code.
- `cargo clippy`: Run linter.

## Architecture Overview

### Core Graph Traits (`src/graphs/mod.rs`)
- `AdjMatrix`: Type alias for `Vec<Vec<bool>>`. Represents a reflexive, symmetric relation.
- `UGraph`: Main trait for undirected graphs.
  - `neighbors(v)`: Iterator over neighbors of `v`.
  - `n()`: Number of vertices.
  - `degree(v)`: Degree of vertex `v`.
  - `is_edge(a, b)`: Check existence of edge.

### Implementations
- **CSRGraph** (`src/graphs/mod.rs`): Compressed Sparse Row representation.
  - Constructed via `TryFrom<AdjMatrix>`.
  - Enforces reflexivity (self-loops) and symmetry.
  - Neighbors are sorted.
- **CubeGraph** (`src/graphs/cube.rs`): Represents an N-dimensional hypercube.
  - Generic over dimension `D` (using `Const<N>` or runtime `u32`).
  - Vertices are integers `0..2^N`.
  - Edges exist between vertices differing by exactly one bit (Hamming distance 1) or self-loops.

### Utilities
- `src/shape.rs`: Handles dimension types (const generic vs runtime).
- `src/graph_maps/`: Logic for graph mappings (cube isomorphism, permutations).

## Coding Conventions
- **Style:** Standard Rust fmt.
- **Testing:** Unit tests co-located in modules (`#[cfg(test)] mod tests`).
- **Graph Invariants:**
  - Graphs are **reflexive** (every node has a self-loop).
  - Graphs are **symmetric** (undirected).
