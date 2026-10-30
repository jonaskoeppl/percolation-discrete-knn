# Planar 2-XnG (XOR two-neighbour graph on Z^2): findings, 30 Oct 2026

The request was to show that the 2-XnG on Z^2 does not percolate, and that its one-arm probability decays at most polynomially.

## Outcome

The model is **not self-dual**, so the natural proof of both statements breaks down. The evidence says the model is slightly **subcritical**. Non-percolation is therefore very likely true but is **not proved**. The polynomial lower bound is very likely **false**: the one-arm probability should decay exponentially, with a correlation length of order 10^3–10^4.

## Rigorous facts, with proofs in the note

- **Characters.** E prod_{e in A} (-1)^{omega(e)} = 3^{-v2(A)} if A is even, and 0 otherwise. Here v2(A) is the number of vertices of degree 2 in A.
  - On any forest the configuration is i.i.d. Bernoulli(1/2).
  - A plaquette has law (1 ± 1/81)/16, with + for an even number of open edges and − for an odd number.
- **Complement symmetry.** 1 - omega has the same law as omega: complement the choices on the even sublattice.
- **No Harris–FKG.** On one plaquette, Cov(1{e1,e3 open}, 1{e2 or e4 open}) = -1/1296.
- **The dual is not a copy of the primal.** Transport the closed dual configuration to Z^2 by any isomorphism. Its 4-point character on a vertex star is 1/81, while for the primal it is 0.
- **Exact crossing probabilities.** Let H_n be the left–right crossing of [0,n+1]x[0,n].
  - P(H_1) = 1/2.
  - P(H_2) = 124345/248832 < 1/2.
  - P(H_3) = 69570879977/139314069504 ≈ 0.4993816.
  - Each value was checked by at least two independent methods: transfer matrix, character formula, and brute force over 6^12 choices.
- **Vertex-XOR family.** Bond percolation at 1/2, AB percolation at 1/2 and the 2-XnG are the XOR graphs with Bin(4,1/2), {0,4} and constant-2 numbers of chosen neighbours, respectively. AB percolation is subcritical.
- **Expansion formula.** E_X f - E_Ber f = sum over nonempty even A of 3^{-v2(A)} fhat(A).
- **Plaquette coefficient of crossing events.** fhat(∂p) = -P_Ber(E_p)/16 ≤ 0.

## Numerics (Monte Carlo)

- **Crossing probability P(H_n):**

  | n | P(H_n) |
  |---|---|
  | 16 | 0.4971 |
  | 256 | 0.4722 |
  | 1024 | 0.421 ± 0.011 |
  | 4096 | 0.30 ± 0.04 |

- **Deficit scaling.** On 32 ≤ n ≤ 1024 the deficit 1/2 - P(H_n) grows like n^0.74.
- **One-arm comparison with Bernoulli(1/2).** At n = 1024 the 2-XnG gives 0.416 and Bernoulli gives 0.459.
- **Adding independent edges.** Adding i.i.d. Bernoulli(delta) edges restores scale invariance only at delta* ≈ 1e-3. That corresponds to an effective edge density of about 1/2 - 5e-4.

## Follow-up: the interpolation ingredient (note version 2, Section 7)

The route to theta = 0 was edge-by-edge interpolation from the 2-XnG to Bernoulli(1/2). Its ingredient: an edge in an intermediate model is not more likely to be open when it is pivotal.

### Exact identity
The interpolation increment equals the sum, over even subgraphs A that contain the edge and avoid the coin set U, of 3^{-v2(A)} fhat(A).

### Strong form (all intermediate models): false in general
- Counterexample: the crossing event H_3 of [0,4]x[0,3], with A a C-shaped simple 12-cycle.
- There fhat(A) = 36·2^-25 > 0.
- The intermediate model with all other edges replaced by coins has crossing probability 1/2 + 36·3^-12·2^-25, strictly above 1/2.
- Mechanism: an open chord inside a cycle, with both arms attached from outside.

### Chain form: open, no counterexample
- The telescoping argument only needs one chain of intermediate models.
- That chain form holds in every exact check, including H_3, using the lexicographic edge order.
- For the one-arm event even the strong form holds in all exactly computable boxes:
  - {0..4}^2;
  - {0..3}x{0..6};
  - {0..4}x{0..5}, with o = (2,2) and with o = (1,2).
- In all of these, P_X < P_Ber, with differences of order 1e-4.

### Two more robust reformulations
Neither is proved; both hold in all computed cases.

**Flux route.** P_nu(f) = E_{mu_nu}[F_f(Phi)], where Phi are the plaquette parities and F_f(Phi) = P_Ber(f | Phi). Two sufficient conditions together give the result:
- (M) F_f is monotone in the set of odd plaquettes.
  - True for crossings and for point-to-boundary events in all tests.
  - False for point-to-point events.
- (D) The 2-XnG parity field is dominated by i.i.d. Bernoulli(1/2).
  - Verified on windows up to 5x4.
  - The margin shrinks with window size: the worst conditional probability rises from 0.4938 to 0.4983.

**kappa-family.** Let nu_lambda = lambda·Unif(2-subsets) + (1-lambda)·Unif{none, all}. Its law depends only on kappa = |1 - 4·lambda/3|:
- kappa = 0 is Bernoulli(1/2). In particular, Bernoulli(1/2) bond percolation is the XOR graph where each vertex picks 2 neighbours w.p. 3/4 and none-or-all w.p. 1/4.
- kappa = 1/3 is the 2-XnG. It equals the 1/2–1/2 mixture.
- kappa = 1 is AB percolation.

The ingredient (K) is that kappa -> P_kappa(o <-> ∂B_n) is non-increasing on [0,1].

## Files

- **Note:** `xor_planar.tex` / `.pdf` (version 2), with code, data (including `data/addendum/`), `prompt.md` and `meta.yml`.
- **Location now:** the session outputs, in the folder `2026-10-30_claude-opus-5.5_2XnG-planar`.
- **Intended repository location:** `ai-outputs/2026-10-30_claude-opus-5.5_2XnG-planar/`.
- **Not yet copied to the repository:** the user's computer was offline in both sessions.

## Open problems

1. Prove theta = 0, for example via the chain form of the ingredient for one-arm events, (M)+(D), or (K).
2. Prove exponential decay (strict subcriticality).