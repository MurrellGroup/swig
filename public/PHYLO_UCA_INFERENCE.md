# Phylogenetic UCA inference in Swig

## Status and purpose

This module estimates the unmutated common ancestor (UCA) of a selected BCR or TCR lineage while representing three kinds of uncertainty that are commonly collapsed into one point estimate:

1. where the UCA attaches to the observed lineage tree;
2. the length of the branch between the UCA and that attachment point; and
3. the UCA nucleotide sequence and V(D)J recombination path.

The method is a **fixed-observed-tree approximation**. It is not the joint tree/recombination phylo-HMM of Dhar et al. and it does not integrate over observed-tree topology or the fitted branch lengths inside that tree. Swig first infers one observed tree, holds it fixed, and then offers conditional ML, explicit grid integration, or a Metropolis-within-Gibbs sampler for the UCA attachment, UCA pendant length, recombination path, and sequence.

The implementation is isolated under `src/phylo-uca/`. React is confined to `panel.tsx` and the HMM-annotation renderer; the likelihood, HMM, placement search, worker, and public data contracts are separate modules.

## Inputs and invariants

The analysis uses:

- the exact current user-curated nucleotide alignment from the lineage workbench, including manual Alivibe edits and deleted rows or columns;
- the AIRR rows corresponding to the retained lineage members;
- the active composed V, D, and J reference FASTA sets;
- fixed, user-visible model and search parameters.

The N-masked germline guide row used by the ordinary lineage viewer is **not an observed taxon**. It is removed before tree inference. For FastTree fitting, a column is removed only when every remaining tip is missing there (`?` or a leading/trailing gap). A column containing an internal gap is retained even when every tip has that internal gap.

The posterior pass then uses the **full original alignment width**. Columns omitted only from FastTree fitting are restored as missing tip data, so they do not change the fitted observed tree but the V(D)J prior can infer a UCA state there. Keeping these columns also preserves the selected codon phase exactly; deleting a three-column-external padding site must not shift every downstream codon. Result coordinates always refer to the original curated alignment.

At least three observed sequences are required. The observed-only tree is inferred with the same double-precision FastTree 2.1.11 WebAssembly executable used by the ordinary lineage-tree action, using nucleotide GTR and the exact retained alignment.

## Character model: four states normally, five only for gapped alignments

This module does **not** use a 5-mer or any other context-dependent likelihood.

Gap semantics are determined **per observed tip**, not per alignment column. For each sequence, Swig finds the first and last observed nucleotide or IUPAC nucleotide character. Gap runs before the first or after the last such character are fragment boundaries and contribute an all-ones likelihood partial (missing data). Only a `-` between those boundaries is an observed internal gap.

In automatic mode:

- an observed alignment with no **internal** gap characters uses an ordinary four-state `A/C/G/T` reversible GTR model, even when partial reads have leading or trailing gap padding;
- an observed alignment containing at least one internal gap uses a five-character `A/C/G/T/-` reversible continuous-time Markov chain for internal gaps only.

`N` and IUPAC ambiguity symbols are partial observations over nucleotide states. They are not gap observations. A `?` is completely missing information. In GTR5, an internal `-` is an exact fifth character; a leading or trailing `-` remains missing. In an explicitly forced GTR4 analysis, internal gaps are also treated as missing and a warning is recorded.

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

There is one V profile state for each candidate allele. Each candidate is pre-aligned to the lineage columns once. A V state advances in that fixed projection and can exit only from a concrete projected nucleotide. The possible exit sites are assigned the finite geometric 3′-deletion distribution controlled by `vThreePrimeTrimContinuation`; the implementation converts those endpoint masses into the corresponding conditional exit hazards. There is no hidden logistic boundary or unresolved junction-emission state.

### Non-templated states

N emits nucleotides from the configured nucleotide frequencies. Its positive duration distribution has an explicit mass `singleNProbability` at length one and a user-selected one-to-four-phase geometric tail whose continuation is solved so the complete positive-run mean is exactly `meanNLength`. The default two-phase tail avoids the excessive one-nucleotide mass of a memoryless geometric run while adding only a few states. There is one copy of these duration states for each number of D segments already used. In GTR5 only, N can emit an alignment gap with `junctionGapProbability`.

Zero N bases are possible because a V or D exit may enter D or J directly. Positive N runs arise by entering and remaining in an N state.

### D states

Every D allele has equal prior mass before trimming and can be entered at any starting nucleotide that leaves at least `minimumDMatch` templated characters. Within an allele, entry weights use the finite geometric 5′-deletion distribution controlled by `dFivePrimeTrimContinuation`. A D path advances one reference nucleotide for each alignment column. It cannot exit before the minimum match is reached. Thereafter the endpoint hazards implement the finite geometric 3′-deletion distribution controlled by `dThreePrimeTrimContinuation`; this is a prior on deleted terminal bases, not a per-retained-base penalty.

After a D exits, the automaton can enter N, enter J, or—at the low configured `additionalDProbability`—enter any D again at any valid start. This loop admits VDDJ, VDDDJ, and higher orders up to `maximumDSegments`. The default maximum is three, but the repeat prior is deliberately small. A direct V-to-J path and a zero-D V-N-J path remain available.

The implementation does not build a dense all-D-to-all-D matrix. It aggregates probability at a D-exit hub and distributes it through normalized D-entry priors. Runtime therefore scales approximately with the number of D profile states, not its square. D states are evaluated only inside the V/J anchor span plus `junctionSearchFlankColumns`; this exposed window is a computational support bound and can be enlarged for unusually deep trimming.

### J states

There is one pre-aligned profile state for every candidate J. Entry at different concrete projected J nucleotides represents uncertain 5′ J trimming and follows the finite geometric distribution controlled by `jFivePrimeTrimContinuation`. A J path cannot enter on an unknown projection and cannot emit a uniform pseudo-NT column before the real J template. Once entered, the J profile proceeds rightward through its fixed nucleotide projection; only alignment padding after the last concrete J nucleotide may use `terminalPaddingGapProbability`.

### HMM emissions and gaps

At alignment column `i`, HMM state `h` defines a prior `p(c | h,i)` over nucleotide characters, and over gap when GTR5 is active. The phylogenetic surface supplies `L_i(c)`. The HMM emission likelihood is:

`e_i(h) = sum_c p(c | h,i) L_i(c)`.

Templated V/D/J states put mass on their reference character, with optional `templateMismatchProbability` leakage. The default is exactly zero; it is not replaced by an arithmetic floor. Zero leakage does **not** forbid a difference between the UCA template base and the sequence favored at the observed tree attachment, because the exact GTR transition along the separately estimated UCA-to-tree branch still models that mutation. N states use the non-templated base frequencies. A projected templated gap in GTR5 places most mass on the gap character. In GTR4, a projected reference gap is unknown because the observed alignment contains no explicit gap evidence.

This distinction is important:

- `N` means an unknown or non-templated **nucleotide**;
- an internal `-` means an explicit fixed-alignment **gap character**, whereas a terminal `-` is missing coverage;
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

For each retained local placement, Swig runs:

- forward/backward HMM inference to obtain marginal `P(UCA_i = c | data)` at every column;
- an exact three-column forward/backward contraction for every complete codon in the selected alignment frame;
- a max-product Viterbi pass to obtain one joint MAP recombination path and sequence.

These are different estimands. The exported **joint MAP aligned UCA** follows one globally consistent V(D)J path. The **marginal consensus** chooses the highest posterior character independently at each site and need not itself be the highest-probability complete path.

Every reported five-character site vector is normalized after local-placement mixing, with the floating-point residual returned to the largest component, so its entries sum to one to machine precision.

### Exact codon posterior

The codon calculation does **not** multiply three nucleotide marginals. Let `alpha_i(h)` be the HMM forward value after the integrated emission at column `i`, let `beta_i(h)` be the backward value for all later columns, and define

`q_i(c | h) = p(c | h,i) L_i(c) / e_i(h)`.

For a codon beginning at `i`, Swig fixes one character at each of the three columns while retaining the HMM state transitions:

`P(c1,c2,c3,data) = sum_(h1,h2,h3) alpha_i(h1) q_i(c1|h1) T_i(h1,h2) e_(i+1)(h2) q_(i+1)(c2|h2) T_(i+1)(h2,h3) e_(i+2)(h3) q_(i+2)(c3|h3) beta_(i+2)(h3)`.

Normalizing the 125 `A/C/G/T/gap` triples gives the joint codon posterior for that placement. Because the hidden state is not discarded between columns, this sum retains correlation induced by a shared V/J candidate, a particular D identity and position, trimming/recombination path, and N/D/J transitions. Marginalizing this 125-state vector at any of its three positions recovers the corresponding nucleotide posterior (up to floating-point precision).

This is posterior aggregation over three nucleotide columns, **not** a codon substitution process. The expensive tree messages, attachment search, and pendant-length search still use nucleotide GTR4/GTR5. No codon-rate matrix or context/5-mer likelihood has been introduced.

## Posterior frequency logo

The results panel includes a reusable probability-logo component. It is deliberately a **frequency logo**, not an information-content logo:

- every column has total height one;
- character height is its marginal posterior probability;
- entropy does not rescale the stack;
- glyphs use a bold monospace face so narrow characters such as `I` remain visible;
- the browser measures each glyph's actual painted SVG bounds and maps those bounds exactly onto its probability rectangle, removing horizontal side bearings and vertical font padding while retaining descenders such as the tail of `Q`;
- adjacent column rectangles share an edge: the renderer adds no horizontal inter-letter spacing;
- the displayed and SVG-exported stacks use the same normalized vectors as the posterior table.

The viewer has three representations of the same posterior:

- **Nucleotide:** the five single-column states `A/C/G/T/gap`.
- **Codon:** the exact 125-state joint distribution described above. Codon strings are colored by their translated amino-acid class.
- **Amino acid:** an exact deterministic marginal obtained by summing probabilities of synonymous codon states. No nucleotide probabilities are multiplied.

A fully gapped codon maps to `-`; a partially gapped codon maps to `X`. The selected lineage-alignment reading frame determines codon boundaries and is stored with the result. Each view can be exported as SVG.

## HMM-derived V(D)J source tracks

The annotation immediately above the frequency logo comes from the phylogenetic recombination HMM itself. It does not reuse the input AIRR calls. The track surface and logo share one horizontal scroll container and one exact column grid. In codon or amino-acid mode, the three contributing nucleotide columns occupy three equal subcolumns above each codon/AA stack.

Codon display never constrains recombination states to multiples of three. V/D/J/N occupancy and transitions remain nucleotide-column quantities, so a trimming or segment boundary may occur after the first or second nucleotide of a displayed codon. The codon posterior explicitly retains the two intervening HMM transitions rather than replacing three nucleotide states with one codon-level source state.

The user can switch between two estimands.

### Best path

**Best path** conditions on the single highest-posterior evaluated tree attachment and pendant length and displays the max-product Viterbi HMM state path at that placement.

- A V, D, or J row displays the fixed reference character belonging to the chosen template state. It is deliberately not replaced with the UCA MAP character; the track describes source annotation, while the logo below describes UCA character uncertainty.
- An N or uncertain trimming-boundary row has no fixed template character. It displays `q_i(c | h)` as a nucleotide stack at the chosen state.
- Before display aggregation, exactly one source row has total height one at each alignment column, apart from floating-point tolerance.

### Marginalized

For retained placement `k`, let `gamma_(k,i)(h)` be the forward-backward posterior occupancy of HMM state `h` at alignment column `i`, and let `w_k` be the normalized local placement weight.

For a pure template track `r`, the displayed occupancy is

`W_(r,i) = sum_k w_k sum_(h in r) gamma_(k,i)(h)`.

That mass is drawn as the one reference character fixed by track `r` at column `i`. For a non-templated track `n`, character-specific height is

`W_(n,i,c) = sum_k w_k sum_(h in n) gamma_(k,i)(h) q_(k,i)(c | h)`.

Therefore the total height of the N stack is exactly its source occupancy, `sum_c W_(n,i,c)`, rather than a separately normalized logo.

The inference sidecar deliberately keeps route-specific rows. In particular, D start/trim uncertainty is keyed by allele, D-use ordinal, and alignment register (`alignment column - D-reference position`). Those were the repeated same-allele rows shown by earlier versions: they represented alternative registers or D-use routes, not distinct copies of the allele and not repeated calls in the input.

The interactive display now performs a visualization-only aggregation over that sidecar:

- all V, D, or J rows with the same segment kind and allele call are summed into one display row;
- if different routes/registers place different reference nucleotides at one column, the combined row shows their nucleotide mixture rather than choosing one registration;
- in **Best path** only, N states and unresolved V-trim/J-entry boundary mass are summed into spatial `NT1`, `NT2`, and later junction rows. The intended order is V, then the NT preceding each D-use block, that block's D alternatives, the next NT, and finally J. If one collapsed D allele carries mass in multiple D-use blocks, it is placed with the block carrying most of its integrated weight, so perfect ordering is not always possible;
- in **Marginalized** mode, every N and unresolved-boundary route is instead summed into one nucleotide-mixture row at the very top. The remaining rows are grouped V, then D, then J, with alleles inside each segment ordered by their posterior-mass center from left to right (mass breaks near ties).

V, NT, D, and J rows have distinct background colors. The label rail counter-scrolls so it remains visible while the posterior is scrolled horizontally. Hovering a label reports the complete allele/track text, maximum occupancy, weighted center, underlying route/register count, D-use ordinals, and alignment-register offsets. Hovering an occupied nucleotide cell reports its alignment column, total source occupancy, each `A/C/G/T/gap` mass as a percentage of the complete posterior, and the corresponding composition conditional on that display row.

The track toolbar exports either the complete HMM-track canvas or the exact currently visible horizontal and vertical track crop. Both exports include the horizontally aligned UCA probability logo, numbering, and CDR bands directly below the tracks. The visible SVG keeps the floating label rail at the current horizontal scroll position and uses the same horizontal crop for the logo. The standalone nucleotide/codon/amino-acid frequency-logo SVG remains available separately.

An allele group enters the serialized sidecar when its combined occupancy reaches 1% in at least one column. Within a retained D group, raw registers reaching 0.1% are retained; if none does, the strongest register is retained. Raw N and uncertain-boundary rows reaching 0.1% are retained. Display aggregation happens only after these existing sidecar thresholds. The result reports how many subthreshold raw rows were omitted. Neither thresholding nor display aggregation changes the HMM likelihood, UCA posterior, Viterbi path, or exported sequence.

## Shared starting-position screen

All three inference routes evaluate branch interiors rather than restricting attachment candidates to existing tree nodes.

1. Every observed-tree edge is screened at several attachment fractions and short pendant lengths. By default, the screen uses only the fixed-alignment V and J regions. At each column it forms one independent nucleotide mixture across retained V or J allele projections and contracts that mixture with the tree likelihood surface. It does **not** evaluate every complete V/J pairing and does not use D or junction states.
2. The leading V/J-screen edges are refined by alternating continuous attachment-fraction and pendant-length optimizations on that cheap surface. This supplies an initializer only.
3. The legacy single N-masked germline-guide screen remains selectable. Neither cheap screen can become a reported full-HMM likelihood.
4. `fullHmmEdges` controls the breadth admitted to conditional ML and grid integration; zero means every observed-tree edge. In Gibbs/MH mode it controls initializer refinement, while global proposals retain nonzero support on every edge.

## Three inference routes

### Conditional maximum likelihood

For every admitted edge, Swig alternates continuous one-dimensional optimization of attachment fraction and UCA pendant length. Every objective call is the **complete recombination-HMM marginal likelihood**. A fourth-power transformed pendant coordinate concentrates optimizer evaluations near zero without discretizing the parameter. Placement and branch priors do not affect this conditional-ML optimum, and there is no marginalization over placement or length.

### Explicit grid marginalization

The attachment grid is linear in within-edge fraction. The pendant grid contains exact zero followed by user-configurable logarithmically spaced positive values between `minimumPositiveUcaBranchLength` and the maximum. The complete list is shown in the settings before a run and stored in the result. Every grid point receives the full-HMM marginal likelihood. Trapezoid/Voronoi cell widths, the selected edge prior, and the exponential pendant prior provide quadrature weights. The leading user-configured number of quadrature points receive exact nucleotide, codon, and HMM-track posterior calculations; the result warns with their cumulative quadrature mass when lower-mass points are omitted.

### Continuous Metropolis-within-Gibbs (default)

One iteration consists of:

1. an exact forward-filter/backward-sample draw of a coherent HMM path and every UCA character conditional on the current attachment and pendant length; and
2. several cheap MH updates of attachment and length conditional on that sampled UCA; and
3. at the configured interval, an independence proposal that changes edge, continuous within-edge position, and continuous pendant length while marginalizing over the complete HMM.

Within-edge position and pendant length are continuous floating-point states. Reflected random-walk proposals evaluate `P(t)=exp(Qt)` at the proposed pendant length; they do not snap to, index, or interpolate a branch-length grid. Global proposals select an observed-tree edge from a full-support V/J-screen proposal distribution and draw a continuous fraction on it, with the independence-proposal Hastings correction.

Conditional on a sampled UCA, the recombination-path terms cancel from an ordinary MH ratio. Such a proposal therefore needs only the likelihood of that fixed UCA against the two cached directed half-edge messages, plus the attachment/length priors. It does not run another HMM and does not recompute the observed-tree pruning likelihood.

The occasional **collapsed refresh** is deliberately more expensive. It computes `p(data | attachment)` by summing over all HMM paths and UCA characters at the proposed placement and uses that marginal target in an exact Hastings ratio. The edge proposal is a full-support mixture of the already-computed full-HMM initializer weights and the broad V/J screen; position and branch components likewise retain support across their complete allowed ranges. The proposed HMM path/UCA is sampled only after the placement is accepted. This delayed draw is valid because the acceptance probability contains the marginalized likelihood and no proposed latent state; after acceptance Swig immediately draws the latent state from its conditional posterior at the new placement. It would be invalid to accept a collapsed move and retain the old UCA, which Swig does not do.

After burn-in and thinning, retained joint draws estimate nucleotide, exact codon-state, source-track, and number-of-D-segments posteriors. The D-count summary is a direct count from retained paths and adds no HMM evaluations. The UI shows branch-length and full-HMM marginal-likelihood traces, burn-in, acceptance rates, edge switches, simple autocorrelation ESS estimates, ESS per sampling second, and mean timings for a full-HMM draw, collapsed marginal proposal, and fixed-UCA proposal. The seed and all proposal controls are user-visible. Runs with branch or log-target ESS below 20 carry an explicit warning.

For grid or Gibbs/MH inference, the default edge prior is proportional to edge length, corresponding to a uniform attachment location over total tree length when the within-edge coordinate is a fraction. A uniform-per-edge alternative is available. Pendant length has an exponential prior with user-configurable mean and a hard maximum.

The results panel draws every retained point or draw at its exact branch fraction and pendant length. Marker color is `exp(LL_k - LL_best)` from the raw full-HMM marginal likelihood: the best likelihood is red and a relative likelihood approaching zero is blue. The adjacent list reports exact coordinates, likelihood difference, relative likelihood, and marginal/sample weight.

## Rooted-tree export

The placed Newick is rooted at the inferred UCA. The named `phylo_UCA` sequence carrier has branch length zero. The **entire** inferred UCA pendant length is assigned to the branch from that root to the observed-tree attachment point; it is not split in half. The original attachment edge is split only according to the optimized attachment distance.

## Default parameters

| Group | Parameter | Default |
|---|---|---:|
| Character model | mode | automatic GTR4 / internal-gap GTR5 |
| Character model | gap equilibrium frequency | 0.02 |
| Candidate screen | V extra differences | 6 |
| Candidate screen | J extra differences | 4 |
| Candidate screen | maximum V / J candidates | 48 / 24 |
| HMM | maximum D segments | 3 |
| HMM | minimum D match | 5 nt |
| HMM | identifiable first-D probability | 0.934 |
| HMM | additional-D probability | 0.00125 |
| HMM | non-empty N probability | 0.973 |
| HMM | V3 / D5 / D3 / J5 deletion-tail ratios | 0.7527 / 0.8574 / 0.8471 / 0.8708 |
| HMM | positive N mean / one-nt mass / tail phases | 8.8 nt / 0.027 / 2 |
| HMM | N A/C/G/T weights | 0.203 / 0.288 / 0.304 / 0.205 |
| HMM | template leakage | exactly 0 |
| HMM | junction gap probability, GTR5 only | 0.015 |
| HMM | terminal-padding gap probability | 0.01 |
| HMM | D-search flank | 16 columns |
| Search | starting-position screen | independent V/J nucleotide mixture |
| Search | screen points per edge | 5 |
| Search | full-HMM edges | 6 |
| Search | inference route | continuous Gibbs/MH |
| Conditional ML | coordinate rounds / unit tolerance | 2 / 0.002 |
| Grid | attachment / pendant points | 3 / 13 |
| Grid | smallest positive pendant length | 0.00001 |
| Search | maximum UCA branch | 0.30 substitutions/character |
| Grid / Gibbs-MH | exponential branch-prior mean | 0.06 |
| Grid | retained posterior points | 12 |
| Gibbs/MH | iterations / burn-in / thin | 320 / 80 / 2 |
| Gibbs/MH | MH steps per Gibbs draw | 4 |
| Gibbs/MH | pendant / position proposal scale | 0.012 / 0.40 |
| Gibbs/MH | global-jump probability | 0.12 |
| Gibbs/MH | focused global-position mixture / half-width | 0.85 / 0.18 |
| Gibbs/MH | focused collapsed-branch mixture / maximum | 0.90 / 0.03 |
| Gibbs/MH | full-HMM initializer edge mixture | 0.95 |
| Gibbs/MH | collapsed-refresh interval | every 3 iterations |
| Gibbs/MH | reproducible seed | 1729 |

The trim means, positive-N mean, one-base N mass, non-empty N mass, and N composition are compact moment matches to the public human-IGH IGoR/OLGA parameter files. The first-D default is the corresponding probability that at least the identifiable minimum survives. The additional-D value is a conservative starting value informed by the reported rarity of V(DD)J rearrangements. The gap probabilities, maximum-D bound, minimum identifiable D match, finite D-search support, and MCMC proposal controls are Swig regularizers or computational defaults rather than literature estimates.

Partis/ham provides a useful warning against overinterpreting the compact trim tails: its empirical work finds reproducible, often non-parametric allele-specific deletion distributions, and uses tiered aggregation when an allele has too few observations. Swig does **not** import those thousands of fitted parameters or change its topology to match partis; the shared tails above remain a fast, visible starting model. Linearham/partis does support the broader design choice of sampling naive sequence under a recombination HMM with concrete germline states rather than allowing an unresolved early-J pseudo-state.

Every HMM value in the table, including all four N-base weights, is exposed in the advanced panel so sensitivity analyses are possible. The Additional-D control accepts the complete probability interval from 0 through 1 and remains editable after a completed run. A dedicated reset button restores every UCA option to these defaults. The panel also runs a deterministic prior-predictive recombination audit against the active D lengths and reports trim, retained-D, N, junction-span, and D-count summaries. A result JSON stores the complete option set and model provenance.

The supplied lineage-40 zero-leakage regression, all three route timings, seed replication, and collapsed-refresh tuning by ESS per second are reported in [`BENCHMARK_PHYLO_UCA_0.27.0.md`](../BENCHMARK_PHYLO_UCA_0.27.0.md).

## Complexity and implementation details

For `n` observed tips and `L` alignment columns, directed tree messages cost `O(n L K^2)` once, where `K` is four or five. They are never recomputed during placement inference. A full marginal-likelihood placement costs `O(L K^2)` to construct the local phylogenetic surface plus one sparse/factorized HMM pass.

If `S` is the number of HMM profile states, a marginal placement evaluation is approximately `O(L S)`. D-repeat transitions do not add `O(S^2)` cost. The immutable V/D/J state catalog, entry priors, endpoint hazards, and fixed D emission categories are constructed once and reused by ML, grid, and Gibbs/MH. D states are inactive outside the exposed junction support, so those framework columns neither scan D-entry registers nor store D backward values.

Grid/ML search uses rolling rows. An exact Gibbs draw stores backward values for all non-D states across the alignment and D values only across the bounded junction window; this typed storage is reused between iterations rather than reallocated. Each subsequent conditional placement MH proposal costs only `O(L K^2)`. A collapsed refresh costs another rolling full-HMM marginal pass, which is why the default frequency was chosen by ESS per wall-clock second rather than steps per second. Exact analytic codon marginalization in ML/grid mode remains linear in `L` and `S`, with a larger constant than a single-site posterior, and never constructs a dense `S x S` matrix. Gibbs/MH codon output instead counts the exact joint three-character states in retained joint samples.

HMM-source sidecar data is accumulated during the same backward pass. V/J tracks use their fixed projections; D tracks use sparse allele/register keys; N tracks retain five character masses. The display then sums those sparse rows by allele and spatial NT slot. Neither stage reruns placement search or allocates a dense track-by-column-by-state tensor.

The whole inference runs in a dedicated browser worker. Progress distinguishes observed-tree inference, message construction, edge screening, full-HMM search or Gibbs/MH sampling, joint nucleotide/codon posterior integration, and finalization.

## Exports and session behavior

The panel exports:

- aligned joint-MAP and marginal-consensus UCA FASTA;
- per-column posterior TSV (`A/C/G/T/gap`, entropy, HMM segment, and candidate call);
- long-form exact codon-posterior TSV (alignment columns, codon, translated amino acid, probability, and MAP indicator);
- frequency-logo SVG in nucleotide, codon, or amino-acid view;
- complete HMM-source-track-plus-UCA-logo SVG and an SVG cropped to the current horizontal/vertical track viewport with the corresponding logo window below, for either Best path or Marginalized mode;
- UCA-rooted placed Newick;
- complete JSON containing the observed tree, placement set, weights, candidate report, model, parameters, Viterbi and marginalized HMM-source tracks, path, sequences, warnings, and alignment fingerprint;
- publication-oriented SVG through the standard lineage tree viewer.

Save session retains the inferred UCA result and all settings because reconstructing it requires a searched placement posterior. The state is bound to the exact lineage ID set, alignment fingerprint, and selected reading frame. Editing, replacing, deleting alignment rows, or changing the frame invalidates the result rather than silently reusing an incompatible codon posterior.

## Limitations

1. The observed tree is fixed. Tree-topology and ordinary branch-length uncertainty are excluded.
2. FastTree is a practical approximate tree estimator, not a full antibody-specific Bayesian tree model.
3. The HS5F-derived GTR is a context-averaged nucleotide approximation. It is not a 5-mer model and does not reproduce full SHM targeting.
4. The five-character model conditions on user-curated internal alignment gaps and is not a generative indel model; terminal gap padding is missing data.
5. Candidate reference annotations and fixed projections must be biologically coherent. Results with incomplete or wrongly aligned germlines require inspection.
6. Grid integration remains a finite quadrature approximation. Gibbs/MH samples attachment and pendant uncertainty continuously but still conditions on one fixed observed tree and fixed model parameters.
7. Fixed transition parameters influence weakly identified junctions. Important UCAs should be rerun across sensible prior settings.
8. Heavy/light or paired-chain joint UCA inference is not currently performed; each selected locus is modeled separately.

## References and implementation precedents

- MurrellGroup. [MolecularEvolution.jl](https://github.com/MurrellGroup/MolecularEvolution.jl). Directional `forward!`/`backward!` branch-message conventions and generic likelihood partitions.
- MurrellGroup. [EvoOnline](https://github.com/MurrellGroup/EvoOnline). Browser-native reversible phylogenetic model and WebAssembly implementation precedents.
- Matsen FA, Kodner RB, Armbrust EV. [pplacer: linear time maximum-likelihood and Bayesian phylogenetic placement of sequences onto a fixed reference tree](https://doi.org/10.1186/1471-2105-11-538). *BMC Bioinformatics* (2010).
- Marcou Q, Mora T, Walczak AM. [High-throughput immune repertoire analysis with IGoR](https://doi.org/10.1038/s41467-018-02832-w). *Nature Communications* (2018). Source of the public human-IGH rearrangement model used for compact default moment matching; Swig does not run IGoR inference.
- Ralph DK, Matsen FA. [Consistency of VDJ rearrangement and substitution parameters enables accurate B cell receptor sequence annotation](https://doi.org/10.1371/journal.pcbi.1004409). *PLOS Computational Biology* (2016). Partis/ham precedent for allele-specific categorical deletion distributions, N-region modeling, multi-sequence HMM inference, and tiered aggregation in small samples.
- Psathyrella. [ham general-purpose HMM compiler](https://github.com/psathyrella/ham). Software architecture precedent only; Swig contains an independent TypeScript factorized automaton and does not incorporate GPL source.
- Briney BS et al. [Frequency and genetic characterization of V(DD)J recombinants in the human peripheral blood antibody repertoire](https://doi.org/10.1111/j.1365-2567.2012.03605.x). *Immunology* (2012). Empirical rarity reference for the conservative additional-D starting probability.
- Yaari G et al. [Models of somatic hypermutation targeting and substitution based on synonymous mutations from high-throughput immunoglobulin sequencing data](https://pmc.ncbi.nlm.nih.gov/articles/PMC3828525/). *Frontiers in Immunology* (2013).
- Hoehn KB, Lunter G, Pybus OG. [A phylogenetic codon substitution model for antibody lineages](https://doi.org/10.1534/genetics.116.196303). *Genetics* (2017). HLP17 is relevant antibody-evolution literature but is a codon model, not the four-state GTR used here.
- Hoehn KB et al. [Repertoire-wide phylogenetic models of B cell molecular evolution reveal evolutionary signatures of aging and vaccination](https://pmc.ncbi.nlm.nih.gov/articles/PMC6842591/). *PNAS* (2019). HLP19 is likewise a codon model.
- Kepler TB. [Reconstructing a B-cell clonal lineage. I. Statistical inference of unobserved ancestors](https://doi.org/10.12688/f1000research.2-103.v1). *F1000Research* (2013). Related precedent for inference of unobserved B-cell ancestors; Swig is not an implementation of this method.
- Dhar A, Ralph DK, Minin VN, Matsen FA. [A Bayesian phylogenetic hidden Markov model for B cell receptor sequence analysis](https://doi.org/10.1371/journal.pcbi.1008030). *PLOS Computational Biology* (2020), implemented in linearham. This is a fuller joint phylogenetic formulation using partis rearrangement parameters; Swig's sampler is limited to UCA placement/length/path/sequence on one fixed observed tree.
