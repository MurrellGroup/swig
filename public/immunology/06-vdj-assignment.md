# 6. Assigning V, D, J, and C

A rearranged receptor sequence is a mixture of inherited sequence, processed junctions, and any later changes. **V(D)J assignment** asks which reference segments best explain the observed sequence and where their contributions lie.

The analysis begins with a [reference database](https://en.wikipedia.org/wiki/Reference_genome) containing candidate receptor segments. For immunoreceptors, a database is a catalogue of many V, D, J, and sometimes C sequences, rather than one complete rearranged receptor.

## Genes and alleles

A [gene](https://en.wikipedia.org/wiki/Gene) can occur in different inherited versions called [alleles](https://en.wikipedia.org/wiki/Allele). In a name such as `IGHV1-2*02`, `IGH` identifies the heavy-chain locus, `V` identifies a variable segment, `1-2` identifies the gene, and `*02` identifies an allele.

A **gene-level call** such as `IGHV1-2` leaves its allele unresolved. An **allele-level call** such as `IGHV1-2*02` makes a finer distinction. Two reference records can differ at only one position, or even have identical sequences over the part covered by the read. In either case, the available sequence can support several labels. ([R13](references.md#r13), [R16](references.md#r16))

The choice of species and locus determines the relevant reference catalogue. An incomplete catalogue can force a receptor onto the closest available relative. A very broad catalogue increases competition among similar candidates. These issues motivate the reference-refinement analyses in Chapter 8.

## What an alignment establishes

A [sequence alignment](https://en.wikipedia.org/wiki/Sequence_alignment) places an observed sequence beside a reference so corresponding positions can be compared. Matching bases contribute evidence for that reference. Mismatches and gaps require an explanation: mutation, inherited variation, sequencing error, an insertion or deletion, or a different assignment.

An alignment also estimates endpoints. At a V–D junction, a few bases might be attributed to the end of V, to an added tract, or to the beginning of D. The gene labels and the boundaries are therefore related parts of the inference.

Swig uses the SwiftIG caller. Its selectable assignment strategies alter candidate search or refinement; calling profiles alter scoring and some decision rules. These are computational choices applied to the same biological problem. A sensitive setting can recover additional weak matches while also admitting more chance matches. The relevant comparison is performance on material resembling the study's species, read length, reference set, and mutation burden. ([S03](references.md#s03))

## Why V, D, and J differ in certainty

A long, well-covered V segment often contains many informative positions. A D segment may contribute only a short tract after trimming. Added junction bases can match such a tract by chance, and SHM can erase some of its original distinguishing bases.

Suppose a retained D contribution consists of six nucleotides. Multiple D alleles may share those six bases, and an independently added tract may also match them. More confident assignment needs additional evidence from the sequence, segment boundaries, reference alternatives, and the assumptions of the caller. An empty D field can reflect unresolved evidence even for a chain normally assembled through D-containing recombination.

J assignments benefit from conserved structural landmarks but can still be ambiguous when reads are short or references similar. A high-quality alignment should be inspected together with its coverage: 100% identity over a short span answers a narrower question than 100% identity over a complete segment. ([R03](references.md#r03), [R15](references.md#r15))

## Interpreting common result fields

The most useful fields divide into several groups:

| Fields or labels | Biological meaning |
|---|---|
| `locus`, `v_call`, `d_call`, `j_call` | Chain identity and candidate inherited segment labels |
| `c_call` | Constant-region sequence supported by the observed downstream tract |
| `junction`, `junction_aa` | Junction nucleotide sequence and its translation under the declared convention |
| CDR and framework annotations | Subregions of the receptor variable domain |
| `productive`, frame, stop, completeness | Sequence-level coding status and coverage |
| Sequence and germline alignments | The observations and reference comparisons underlying the calls |
| Segment identity and support | Match fraction and evidence measures defined by the caller |

In Swig, segment `*_support` fields are calibrated **expectation values**, or E-values, for chance local alignments under the scoring model. Smaller values indicate stronger evidence against chance matching. They address a different uncertainty from the probability that one allele, among several close relatives, is the true inherited source. Co-optimal call labels preserve a scoring tie. ([S03](references.md#s03))

**Productivity** and **completeness** also answer different questions. A partial sequence can contain an apparently intact coding region while lacking the coverage needed to establish the entire rearrangement. A productive annotation supports sequence-level plausibility; binding and cell function require further evidence. ([R13](references.md#r13))

## Coordinates and numbering

A sequence has several coordinate systems. **Query coordinates** count along the observed read. **Germline coordinates** count along a reference segment. **Alignment coordinates** count columns, including inserted gaps. Amino-acid numbering schemes assign positions in a receptor domain so that related structures can be compared.

[IMGT](https://en.wikipedia.org/wiki/IMGT) conventions define standardized receptor regions. [Kabat numbering](https://en.wikipedia.org/wiki/Kabat_numbering) is another widely used antibody numbering scheme and includes insertion labels at some positions. A mutation label should identify its coordinate convention; a change at alignment column 100 and a change at Kabat position 100 may concern different residues.

For custom references, region metadata supplies the landmarks needed to annotate CDRs, frameworks, and frames. Swig can transfer such metadata from annotated relatives. The reference sequence and the reliability of its region annotation remain separate issues. ([S02](references.md#s02))

## Constant-region evidence

A read extending beyond J may identify a heavy-chain constant region and hence an isotype or subclass. A read ending within the variable region leaves that information unobserved. Closely related constant genes can require particular discriminating bases to be covered.

Swig searches for C evidence downstream of J. Interpreting an isotype call therefore includes checking the aligned span and competing constant references. Inferring membrane-bound versus secreted transcript forms additionally requires coverage of the relevant alternative exons. ([S03](references.md#s03))

## Evidence for two D segments

The optional **Double-D** analysis considers ordered, non-overlapping D matches and compares their support with a single-D explanation. It also checks whether a single D could mimic the proposed pair. Supported calls can be inspected as V–D1–D2–J alignments and used when reviewing unusual junctions or lineages. ([S04](references.md#s04))

For a candidate of particular biological interest, useful follow-up includes examining other members of the same family, checking sequence quality, and reviewing the donor's D references. Agreement across related sequences strengthens the reconstruction while retaining the shared-ancestry dependence of those observations.

---

[Previous: Samples, sequences, and counts](05-samples-sequences-and-counts.md) · [Contents](index.md) · [Next: Errors, duplicates, and chimeras](07-errors-and-chimeras.md)
