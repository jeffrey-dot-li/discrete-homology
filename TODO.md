
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
- [x] RedZeD reduction + active enumeration (Algorithm 2), unquotiented.
      `computations::redzed`: `ColumnReducer` is the `R_p` / `R_p^-1` core (rows
      are the *births* of the degree below, i.e. the map `r_{p-1}`), and
      `redzed_dh` runs Phase 1 + Phase 2. Validated: it reproduces
      `homology::betti_numbers` on C3..C6 and Gsph, in both full and active
      modes. Active enumeration skips ~97% of top cubes at degree 2 (Gsph:
      21424 -> 583) and 99.7% at degree 3 (C5: 978350 -> 3207).
  - `beta(C5) = [1, 1, 0, 0]`, so `H_3(C5) = 0`, agreed on by
    `homology::report_to` and by `redzed_dh` in both modes.
  - Active enumeration has to glue along *every* axis, not just the appended
    one (`CubeMap::try_combine_along`). The paper gets this for free by ranging
    the positive face over the whole hyperoctahedral orbit of `B`; without the
    quotient the equivalent is to insert the new coordinate at each of the
    `dim + 1` positions. Pairing on the appended axis alone found a top cube
    only via its two last-coordinate faces, so a cube whose only live face lay
    along an earlier axis was lost -- and its pivot with it, giving
    `beta_3(C5) = 3` instead of `0`.
  - Degrees 1..=3 are wired up. n = 3 builds the degree-3 level by pairing
    degree-2 cubes and then pairs *that* for degree 4, so it is the first
    expensive case; its cross-check against `homology::report_to` is ignored by
    default and wants a release build:
    `cargo test --release --lib -- --ignored --nocapture degree_three`.
    Pick the test graph by degree-3 level size, not by graph size: C4 is worse
    than C5 here (256 degree-2 maps, 64 non-degenerate, so an order of magnitude
    more degree-3 cubes to pair).
  - Remaining: `redzed_dh` and `report_to` hardcode the `Poly` capacity
    dispatch (CAP = 2, 4, 6, 8). n >= 4 needs more dispatch, or a
    non-const-capacity `Poly`. Past that the quotient below is what makes the
    higher degrees affordable.
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
