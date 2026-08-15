# Phylogenetic UCA inference in Swig

## Status and purpose

This module estimates the unmutated common ancestor (UCA) of a selected BCR or TCR lineage while representing three kinds of uncertainty that are commonly collapsed into one point estimate:

1. where the UCA attaches to the observed lineage tree;
2. the length of the branch between the UCA and that attachment point; and
3. the UCA nucleotide sequence and V(D)J recombination path.

The method is a **fixed-tree empirical-Bayes approximation**. It is not the MCMC phylo-HMM of Dhar et al. and it does not integrate over observed-tree topology or branch-length uncertainty. Swig first infers one tree from the observed sequences, holds it fixed, and then integrates or optimizes the quantities listed above.

The implementation is isolated under `src/phylo-uca/`. It has no React dependency except for `panel.tsx`; the likelihood, HMM, placement search, worker, and public data contracts are separate modules.

## Inputs and invariants

The analysis uses:

- the exact current user-curated nucleotide alignment from the lineage workbench, including manual Alivibe edits and deleted rows or columns;
- the AIRR rows corresponding to the retained lineage members;
- the active composed V, D, and J reference FASTA sets;
- fixed, user-visible model and search parameters.

The N-masked germline guide row used by the ordinary lineage viewer is **not an observed taxon**. It is removed before tree inference. Columns that are gaps in every remaining observed row are also removed from the phylogenetic calculation and mapped back to their original alignment coordinates in the result.

At least three observed sequences are required. The observed-only tree is inferred with the same double-precision FastTree 2.1.11 WebAssembly executable used by the ordinary lineage-tree action, using nucleotide GTR and the exact retained alignment.

## Character model: four states normally, five only for gapped alignments

This module does **not** use a 5-mer or any other context-dependent likelihood.

In automatic mode:

- an observed alignment with no `-` characters uses an ordinary four-state `A/C/G/T` reversible GTR model;
- an observed alignment containing at least one `-` uses a five-character `A/C/G/T/-` reversible continuous-time Markov chain.

`N` and IUPAC ambiguity symbols are partial observations over nucleotide states. They are not gap observations. A `?` is completely missing information. In GTR5, `-` is an exact fifth character. In an explicitly forced GTR4 analysis, a gap is treated as missing and a warning is recorded.

The GTR5 treatment is a fixed-alignment approximation. It assigns likelihood to a gap character without modeling insertion/deletion history, fragment length, or alignment uncertainty. It is therefore not a TKF-style indel process. The approximation is useful when the user has already committed to an alignment and wants those gap columns to carry information, but the result should not be interpreted as an indel-rate estimate.

### Nucleotide exchangeabilities

There is no established published four-state “antibody GTR” parameter set directly equivalent to the codon-level HLP17/HLP19 models. Swig therefore uses a documented reversible projection of the human HS5F substitution tables as a fixed default:

1. For every HS5F context with central nucleotide `x`, multiply its published mutability by its conditional substitution probability to `y`.
2. Average those directional rates uniformly over the flanking contexts for each central base. Context is discarded after this averaging; it is not part of the Swig likelihood.
3. Solve for the stationary nucleotide frequencies of the resulting directional rate matrix.
4. Symmetrize stationary flux for each pair:

   `F_xy = (pi_x q_xy + pi_y q_yx) / 2`.

5. Form reversible exchangeabilities `r_xy = F_xy / (pi_x pi_y)` and scale them so `r_GT = 1`.

The resulting nucleotide stationary frequencies are:

| State | Frequency |
|---|---:|
| A | 0.2014228 |
| C | 0.3062753 |
| G | 0.2405870 |
| T | 0.2517149 |

The relative exchangeabilities are:

| Pair | Exchangeability |
|---|---:|
| A–C | 1.433127 |
| A–G | 3.121185 |
| A–T | 1.285023 |
| C–G | 1.278638 |
| C–T | 1.938895 |
| G–T | 1.000000 |

For GTR5, the default gap frequency is 0.02 and each nucleotide-gap exchangeability is 0.05 relative to G–T. These gap values are Swig defaults, not estimates from HS5F. All frequencies are normalized and the complete rate matrix is scaled to one expected character change per unit branch length.

This averaged model is a pragmatic site-independent approximation. HLP17/HLP19 and true context-dependent SHM models remain scientifically relevant alternatives, but they do not reduce to a unique four-state GTR matrix. A later model plugin can replace this module without changing the HMM or placement interface.

## Directed Felsenstein messages on the fixed tree

Let `T` be the unrooted observed-sequence tree, `Q` the selected GTR4 or GTR5 generator, and `P(t) = exp(Qt)` its transition matrix. For every directed half-edge `u -> e`, Swig computes the likelihood vector contributed by the component on the `u` side when edge `e` is cut:

`m_(u->e)(x) = l_u(x) product over f != e [ sum_y P_xy(t_f) m_(v_f->f)(y) ]`.

For a leaf, `l_u` is its exact or ambiguous character partial. For an internal node, it is one. Sitewise scaling factors are accumulated in log space, so the calculation does not underflow on large lineages.

Computing both directions of every edge is the all-to-all “Felsenstein up / Felsenstein down” operation. In the terminology used by `MolecularEvolution.jl`, a tip-to-root likelihood transport is `backward`, while a root-to-tip distribution transport is `forward`. Once all directed messages exist, evaluating another point along an edge requires no new traversal of the tree.

## Conditional likelihood at a proposed UCA

Consider an observed-tree edge with endpoints `a` and `b`, total length `L`, and a proposed attachment point `R` at distance `x` from `a`. Transport the two cavity messages to `R` and multiply them:

`L_R(z) = [sum_i P_zi(x) m_(a->e)(i)] [sum_j P_zj(L-x) m_(b->e)(j)]`.

For a proposed UCA pendant length `tau`, transport that likelihood backward from `R` to the sequence at the top of the UCA branch:

`L_UCA(s) = sum_z P_sz(tau) L_R(z)`.

This produces, for every alignment column, a likelihood vector over the four or five possible UCA characters. It is the interface between the phylogeny and the recombination HMM.

## Broad V, D, and J candidate construction

The candidate stage is deliberately broader than the ordinary top-call display.

1. Every selected, co-optimal, and retained near-tied V/J call in every lineage member is collected.
2. Each active reference V and J allele is projected once onto the fixed lineage-guide coordinates with a semi-global affine-like fixed-alignment screen.
3. A broad fixed-alignment V and J score is rerun for **every lineage member**. Candidates within the configured extra-difference window of that member’s best candidate are retained.
4. The union across all members is taken. Observed call hypotheses are never removed by the identity threshold or candidate cap.
5. Candidate caps are applied only after that union, with observed hypotheses first.
6. Every active D record for the selected locus remains a candidate. D genes are not pre-pruned by a top call.

This is designed to exclude only V/J alleles that are extremely implausible while preventing one sequence, one initially reported call, or one inferred guide from determining the complete hypothesis set.

## Recombination HMM

The HMM is a left-to-right automaton over fixed alignment columns. It supplies a structured prior over UCA characters rather than treating sites independently.

### V states

There is one V profile state for each candidate allele. Each candidate is pre-aligned to the lineage columns once. A V state normally advances by staying in the same profile while the alignment column moves right. Its exit probability increases smoothly near the AIRR-derived V endpoint; the width of that logistic transition is `vTrimScale`. This represents uncertainty in V trimming without repeatedly realigning every V at every placement.

### Non-templated states

There is one N state for each number of D segments already used. N emits nucleotides from the configured nucleotide frequencies and has a geometric run-length distribution controlled by `meanNLength`. In GTR5 only, it can emit an alignment gap with `junctionGapProbability`.

Zero N bases are possible because a V or D exit may enter D or J directly. Positive N runs arise by entering and remaining in an N state.

### D states

Every D gene can be entered at any starting nucleotide that leaves at least `minimumDMatch` templated characters. Entry weights use a geometric 5′-trim prior. A D path advances one reference nucleotide for each alignment column. It cannot exit before the minimum match is reached; after that, `dExitProbability` controls its geometric continuation.

After a D exits, the automaton can enter N, enter J, or—at the low configured `additionalDProbability`—enter any D again at any valid start. This loop admits VDDJ, VDDDJ, and higher orders up to `maximumDSegments`. The default maximum is three, but the repeat prior is deliberately small. A direct V-to-J path and a zero-D V-N-J path remain available.

The implementation does not build a dense all-D-to-all-D matrix. It aggregates probability at a D-exit hub and distributes it through normalized D-entry priors. Runtime therefore scales approximately with the number of D profile states, not its square.

### J states

There is one pre-aligned profile state for every candidate J. Entry at different alignment columns represents uncertain 5′ J trimming. The probability of entering J rises smoothly around the AIRR-derived J start with width `jTrimScale`. Once entered, the J profile proceeds rightward to the end of the alignment. Reference positions skipped by the fixed projection are represented directly in that profile.

### HMM emissions and gaps

At alignment column `i`, HMM state `h` defines a prior `p(c | h,i)` over nucleotide characters, and over gap when GTR5 is active. The phylogenetic surface supplies `L_i(c)`. The HMM emission likelihood is:

`e_i(h) = sum_c p(c | h,i) L_i(c)`.

Templated V/D/J states put nearly all mass on their reference character, with a small configurable leakage probability. N states use the non-templated base frequencies. A projected templated gap in GTR5 places most mass on the gap character. In GTR4, a projected reference gap is unknown because the observed alignment contains no explicit gap evidence.

This distinction is important:

- `N` means an unknown or non-templated **nucleotide**;
- `-` means an explicit fixed-alignment **gap character**;
- the HMM does not silently convert one into the other;
- gap runs are not treated as ordinary N-addition runs;
- no context/5-mer state is present.

## Forward likelihood, backward posterior, and Viterbi sequence

For a proposed attachment and pendant length, the HMM forward algorithm marginalizes over:

- V candidate;
- V trimming boundary;
- zero or more N bases;
- D identity, start, end, and number of D segments;
- J candidate and entry boundary;
- UCA character at every site.

The resulting value is the marginal likelihood of the observed tree plus the recombination prior for that placement.

For the selected placement, Swig runs:

- forward/backward HMM inference to obtain marginal `P(UCA_i = c | data)` at every column;
- a max-product Viterbi pass to obtain one joint MAP recombination path and sequence.

These are different estimands. The exported **joint MAP aligned UCA** follows one globally consistent V(D)J path. The **marginal consensus** chooses the highest posterior character independently at each site and need not itself be the highest-probability complete path.

## Placement and branch-length search

Search is coarse to fine.

1. Every observed-tree edge is screened at its midpoint using the existing N-masked lineage guide and several short pendant lengths. This screen is only a computational ranking device; it is not the reported UCA score.
2. The top `fullHmmEdges` edges receive the complete recombination-HMM marginal likelihood.
3. A grid spans attachment fraction along each selected edge and a quadratic grid spans UCA pendant length, giving more resolution near zero.
4. The best placements are refined locally for the configured number of rounds.
5. The reported point estimate maximizes marginal likelihood plus the fixed placement priors.

The default edge prior is proportional to edge length, approximating a uniform prior over location on the continuous tree. A uniform-per-edge alternative is available. Pendant length has an exponential prior with user-configurable mean and a hard user-configurable search maximum.

## Local empirical-Bayes marginalization

When enabled, Swig normalizes posterior weights over the best nearby evaluated `(edge, attachment fraction, pendant length)` points and averages their sitewise UCA nucleotide posteriors. The result reports the effective number of placements:

`N_eff = exp[-sum_k w_k log(w_k)]`.

This represents local uncertainty around the best region of a fixed tree. It is not a continuous quadrature guarantee and does not represent alternative observed-tree topologies. The joint MAP path and placed-tree export remain tied to the single best placement.

## Rooted-tree export

The placed Newick is rooted at the inferred UCA. The named `phylo_UCA` sequence carrier has branch length zero. The **entire** inferred UCA pendant length is assigned to the branch from that root to the observed-tree attachment point; it is not split in half. The original attachment edge is split only according to the optimized attachment distance.

## Default parameters

| Group | Parameter | Default |
|---|---|---:|
| Character model | mode | automatic GTR4/GTR5 |
| Character model | gap equilibrium frequency | 0.02 |
| Candidate screen | V extra differences | 6 |
| Candidate screen | J extra differences | 4 |
| Candidate screen | maximum V / J candidates | 48 / 24 |
| HMM | maximum D segments | 3 |
| HMM | minimum D match | 5 nt |
| HMM | additional-D weight | 0.015 |
| HMM | mean N length | 5 nt |
| HMM | template leakage | 0.003 |
| HMM | junction gap probability, GTR5 only | 0.015 |
| Search | full-HMM edges | 6 |
| Search | maximum UCA branch | 0.30 substitutions/character |
| Search | exponential branch-prior mean | 0.06 |
| Search | local posterior points | 12 |

These are fixed regularizing values, not fitted biological recombination rates. They are exposed in the attached advanced panel so sensitivity analyses are possible. A result JSON stores the complete option set and model provenance.

## Complexity and implementation details

For `n` observed tips and `L` alignment columns, directed tree messages cost `O(n L K^2)` once, where `K` is four or five. Each proposed placement then costs `O(L K^2)` to construct the phylogenetic surface plus one sparse/factorized HMM pass.

If `S` is the number of HMM profile states, a marginal placement evaluation is approximately `O(L S)`. D-repeat transitions do not add `O(S^2)` cost. Full posterior inference stores forward values for the selected local placements; coarse search uses rolling rows. The whole inference runs in a dedicated browser worker. Progress distinguishes observed-tree inference, message construction, edge screening, full-HMM search, local posterior integration, and finalization.

## Exports and session behavior

The panel exports:

- aligned joint-MAP and marginal-consensus UCA FASTA;
- per-column posterior TSV (`A/C/G/T/gap`, entropy, HMM segment, and candidate call);
- UCA-rooted placed Newick;
- complete JSON containing the observed tree, placement set, weights, candidate report, model, parameters, path, sequences, warnings, and alignment fingerprint;
- publication-oriented SVG through the standard lineage tree viewer.

Save session retains the inferred UCA result and all settings because reconstructing it requires a searched placement posterior. The state is bound to the exact lineage ID set and alignment fingerprint. Editing, replacing, or deleting alignment rows invalidates the result rather than silently reusing it.

## Limitations

1. The observed tree is fixed. Tree-topology and ordinary branch-length uncertainty are excluded.
2. FastTree is a practical approximate tree estimator, not a full antibody-specific Bayesian tree model.
3. The HS5F-derived GTR is a context-averaged nucleotide approximation. It is not a 5-mer model and does not reproduce full SHM targeting.
4. The five-character model conditions on a user-curated alignment and is not a generative indel model.
5. Candidate reference annotations and fixed projections must be biologically coherent. Results with incomplete or wrongly aligned germlines require inspection.
6. Local placement marginalization is discrete and local, not full MCMC.
7. Fixed transition parameters influence weakly identified junctions. Important UCAs should be rerun across sensible prior settings.
8. Heavy/light or paired-chain joint UCA inference is not currently performed; each selected locus is modeled separately.

## References and implementation precedents

- MurrellGroup. [MolecularEvolution.jl](https://github.com/MurrellGroup/MolecularEvolution.jl). Directional `forward!`/`backward!` branch-message conventions and generic likelihood partitions.
- MurrellGroup. [EvoOnline](https://github.com/MurrellGroup/EvoOnline). Browser-native reversible phylogenetic model and WebAssembly implementation precedents.
- Matsen FA, Kodner RB, Armbrust EV. [pplacer: linear time maximum-likelihood and Bayesian phylogenetic placement of sequences onto a fixed reference tree](https://doi.org/10.1186/1471-2105-11-538). *BMC Bioinformatics* (2010).
- Yaari G et al. [Models of somatic hypermutation targeting and substitution based on synonymous mutations from high-throughput immunoglobulin sequencing data](https://pmc.ncbi.nlm.nih.gov/articles/PMC3828525/). *Frontiers in Immunology* (2013).
- Hoehn KB, Lunter G, Pybus OG. [A phylogenetic codon substitution model for antibody lineages](https://doi.org/10.1534/genetics.116.196303). *Genetics* (2017). HLP17 is relevant antibody-evolution literature but is a codon model, not the four-state GTR used here.
- Hoehn KB et al. [Repertoire-wide phylogenetic models of B cell molecular evolution reveal evolutionary signatures of aging and vaccination](https://pmc.ncbi.nlm.nih.gov/articles/PMC6842591/). *PNAS* (2019). HLP19 is likewise a codon model.
- Dhar A, Ralph DK, Minin VN, Matsen FA. [A Bayesian phylogenetic hidden Markov model for B cell receptor sequence analysis](https://doi.org/10.1371/journal.pcbi.1008030). *PLOS Computational Biology* (2020). This is the fuller MCMC formulation that motivates, but is not implemented by, Swig’s current fixed-tree empirical-Bayes approximation.
