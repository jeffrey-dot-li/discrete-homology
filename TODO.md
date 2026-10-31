
- [x] Generate all boundaries of non degen n+1 maps ($\partial$)
  - `computations::boundary::{partial, boundary_maps}`; `cube_maps::d` is the face map.
- [x] Compute the rank / dimension of all non-degen boundaries.
  - Partial: `computations::homology::{report, rank_mod_p}` — sparse Gaussian
    elimination over `F_p`, degrees 1..3.
  - Validated against Ebel-Kapulkin-Kershaw, arXiv:2608.04262:
    `beta(Gsph) = [1, 0, 1]` (dim H_2 = 1) and `beta(C5) = [1, 1, 0]`.
    `cargo test --lib dump_report -- --nocapture` prints the table.
  - Remaining: the `Poly` capacity dispatch in `homology.rs` is hardcoded to
    CAP = 2, 4, 6. Degree >= 4 needs RedZeD (below), not just more dispatch.
- [ ] Bench amount of time it takes to generate

## Aligning with arXiv:2608.04262 (RedZeD DH)

- [ ] Hyperoctahedral quotient `Aut(Q_n) = S_n |> (Z_2)^n` plus `hyper_orbit` /
      `rep_and_sign` ([GWW25] arXiv:2508.14639); stub in
      `graph_maps/cube_isomorphism.rs`.
- [ ] Semi-degenerate discarding (Def 1.4) and work over `F_p` with
      `p ∤ (n+1)!` (required by Theorem 1.3).
- [ ] RedZeD `R_p` / `R_p^-1` dictionaries, `process_redzed`, and active
      enumeration (Algorithm 2). This is what makes H_3 / H_4 reachable: the
      Greene sphere has ~6.4e12 5-cubes, so the top-dimensional generators must
      never be enumerated.
- [ ] Dominated-vertex preprocessing from [KK24] arXiv:2410.09939.

## Conventions

- `partial` computes `sum_i (-1)^i (A o d_i^+ - A o d_i^-)` with `i`
  zero-indexed; the paper writes `sum_{i=1}^n (-1)^i (A o d_i^0 - A o d_i^1)`.
  The two differ by an overall sign (repo `d = -d_paper`). Ranks, kernels and
  images are unaffected, so the Betti numbers are unaffected.

## Environment

- Relinking was broken on this machine: nix's rustup patch had replaced
  `~/.rustup/toolchains/stable-*/lib/rustlib/x86_64-unknown-linux-gnu/bin/gcc-ld/ld.lld`
  with a shim pointing at a GC'd store path (`rustup-1.28.2`) while 1.29.2 is
  installed, so every relink failed and only pre-built artifacts ran.
  Fixed: the real linker (`gcc-ld-unwrapped/ld.lld`) was copied back over
  `gcc-ld/ld.lld`. Broken shim kept at
  `~/.rustup/backup/ld.lld.broken-nix-shim.bak`. A `rustup` toolchain
  reinstall/update may reintroduce the shim.
- If linking breaks again, `RUSTFLAGS="-C link-arg=-fuse-ld=bfd" cargo test`
  bypasses rust-lld for the system `ld`.
