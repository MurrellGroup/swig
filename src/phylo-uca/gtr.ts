import {
  PHYLO_UCA_CHARACTERS,
  type PhyloUcaCharacter,
  type PhyloUcaGtrModel,
} from "./types.ts";
import { alignmentGapSemantics, isTerminalAlignmentGap } from "./gaps.ts";

export const PHYLO_UCA_MAX_STATE_COUNT = 5;

const PAIRS: ReadonlyArray<readonly [number, number]> = [
  [0, 1], [0, 2], [0, 3], [0, 4], [1, 2],
  [1, 3], [1, 4], [2, 3], [2, 4], [3, 4],
];

/**
 * Reversible projection of the human HS5F model (Yaari et al. 2013).
 *
 * The original 1,024 context rows were first multiplied by their published
 * mutabilities, averaged uniformly over contexts sharing a central base, and
 * then projected onto a reversible nucleotide flux. No 5-mer/context state is
 * used by the phylogenetic likelihood. Gap exchangeabilities are an explicit
 * fixed-alignment approximation and were not estimated by Yaari et al.
 */
export const HS5F_REVERSIBLE_GTR5: PhyloUcaGtrModel = {
  id: "hs5f-reversible",
  label: "HS5F-averaged reversible GTR5",
  frequencies: [0.197394, 0.300150, 0.235775, 0.246681, 0.02],
  exchangeabilities: [
    1.433127, 3.121185, 1.285023, 0.05,
    1.278638, 1.938895, 0.05,
    1.0, 0.05,
    0.05,
  ],
  provenance: "Human HS5F mutability-weighted substitution flux, uniformly averaged across 5-mer contexts and reversibilized; shared nucleotide-gap exchangeability 0.05 is a Swig fixed-alignment approximation.",
};

export const JC5_CONTROL_MODEL: PhyloUcaGtrModel = {
  id: "jc5",
  label: "Equal-rates five-character control",
  frequencies: [0.2, 0.2, 0.2, 0.2, 0.2],
  exchangeabilities: [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
  provenance: "Five-state equal-frequency/equal-exchangeability control model.",
};

function normalize(values: readonly number[], minimum = 1e-12): number[] {
  const safe = values.map((value) => Number.isFinite(value) && value > minimum ? value : minimum);
  const total = safe.reduce((sum, value) => sum + value, 0);
  return safe.map((value) => value / total);
}

export function empiricalGtr5Model(alignmentSequences: readonly string[], useHs5fExchangeabilities = true): PhyloUcaGtrModel {
  const counts = [1, 1, 1, 1, 1];
  for (const raw of alignmentSequences) {
    const sequence = raw.toUpperCase().replaceAll(".", "-").replaceAll("U", "T");
    const semantics = alignmentGapSemantics(sequence);
    for (let site = 0; site < sequence.length; site += 1) {
      const character = sequence[site];
      if (character === "-" && isTerminalAlignmentGap(site, semantics)) continue;
      const index = PHYLO_UCA_CHARACTERS.indexOf(character as PhyloUcaCharacter);
      if (index >= 0) counts[index] += 1;
    }
  }
  const frequencies = normalize(counts) as [number, number, number, number, number];
  const base = useHs5fExchangeabilities ? HS5F_REVERSIBLE_GTR5.exchangeabilities : JC5_CONTROL_MODEL.exchangeabilities;
  return {
    id: "empirical-gtr",
    label: useHs5fExchangeabilities ? "Empirical frequencies + HS5F-averaged rates" : "Empirical-frequency GTR5",
    frequencies,
    exchangeabilities: [...base] as PhyloUcaGtrModel["exchangeabilities"],
    provenance: `${useHs5fExchangeabilities ? "HS5F-averaged nucleotide" : "Equal"} exchangeabilities with Laplace-smoothed frequencies estimated from the observed-only curated alignment.`,
  };
}

function jacobiSymmetric(matrix: Float64Array, dimension: number): { values: Float64Array; vectors: Float64Array } {
  const values = matrix.slice();
  const vectors = new Float64Array(dimension * dimension);
  for (let index = 0; index < dimension; index += 1) vectors[index * dimension + index] = 1;
  const maximumIterations = 100 * dimension * dimension;
  for (let iteration = 0; iteration < maximumIterations; iteration += 1) {
    let p = 0;
    let q = 1;
    let maximum = 0;
    for (let row = 0; row < dimension; row += 1) {
      for (let column = row + 1; column < dimension; column += 1) {
        const magnitude = Math.abs(values[row * dimension + column]);
        if (magnitude > maximum) {
          maximum = magnitude;
          p = row;
          q = column;
        }
      }
    }
    if (maximum < 1e-14) break;
    const pp = values[p * dimension + p];
    const qq = values[q * dimension + q];
    const pq = values[p * dimension + q];
    const angle = 0.5 * Math.atan2(2 * pq, qq - pp);
    const cosine = Math.cos(angle);
    const sine = Math.sin(angle);
    for (let index = 0; index < dimension; index += 1) {
      if (index === p || index === q) continue;
      const ip = values[index * dimension + p];
      const iq = values[index * dimension + q];
      const nextIp = cosine * ip - sine * iq;
      const nextIq = sine * ip + cosine * iq;
      values[index * dimension + p] = nextIp;
      values[p * dimension + index] = nextIp;
      values[index * dimension + q] = nextIq;
      values[q * dimension + index] = nextIq;
    }
    values[p * dimension + p] = cosine * cosine * pp - 2 * sine * cosine * pq + sine * sine * qq;
    values[q * dimension + q] = sine * sine * pp + 2 * sine * cosine * pq + cosine * cosine * qq;
    values[p * dimension + q] = 0;
    values[q * dimension + p] = 0;
    for (let row = 0; row < dimension; row += 1) {
      const rp = vectors[row * dimension + p];
      const rq = vectors[row * dimension + q];
      vectors[row * dimension + p] = cosine * rp - sine * rq;
      vectors[row * dimension + q] = sine * rp + cosine * rq;
    }
  }
  const eigenvalues = new Float64Array(dimension);
  for (let index = 0; index < dimension; index += 1) eigenvalues[index] = values[index * dimension + index];
  return { values: eigenvalues, vectors };
}

export interface ReversibleCharacterModel {
  dimension: 4 | 5;
  frequencies: Float64Array;
  rateMatrix: Float64Array;
  transition: (length: number) => Float64Array;
}

export function compileGtr(model: PhyloUcaGtrModel, includeGap: boolean): ReversibleCharacterModel {
  const dimension: 4 | 5 = includeGap ? 5 : 4;
  const frequencies = Float64Array.from(normalize(model.frequencies.slice(0, dimension)));
  const rateMatrix = new Float64Array(dimension * dimension);
  for (let pair = 0; pair < PAIRS.length; pair += 1) {
    const [left, right] = PAIRS[pair];
    if (right >= dimension) continue;
    const exchangeability = Number.isFinite(model.exchangeabilities[pair]) && model.exchangeabilities[pair] > 0
      ? model.exchangeabilities[pair]
      : 1e-12;
    rateMatrix[left * dimension + right] = exchangeability * frequencies[right];
    rateMatrix[right * dimension + left] = exchangeability * frequencies[left];
  }
  for (let row = 0; row < dimension; row += 1) {
    let total = 0;
    for (let column = 0; column < dimension; column += 1) if (column !== row) total += rateMatrix[row * dimension + column];
    rateMatrix[row * dimension + row] = -total;
  }
  let expectedRate = 0;
  for (let index = 0; index < dimension; index += 1) expectedRate -= frequencies[index] * rateMatrix[index * dimension + index];
  const scale = expectedRate > 0 ? 1 / expectedRate : 1;
  for (let index = 0; index < rateMatrix.length; index += 1) rateMatrix[index] *= scale;

  const symmetric = new Float64Array(dimension * dimension);
  for (let row = 0; row < dimension; row += 1) {
    for (let column = 0; column < dimension; column += 1) {
      symmetric[row * dimension + column] = Math.sqrt(frequencies[row]) * rateMatrix[row * dimension + column] / Math.sqrt(frequencies[column]);
    }
  }
  const eigensystem = jacobiSymmetric(symmetric, dimension);
  const cache = new Map<string, Float64Array>();
  const transition = (rawLength: number): Float64Array => {
    const length = Number.isFinite(rawLength) && rawLength > 0 ? rawLength : 0;
    const key = length.toPrecision(13);
    const cached = cache.get(key);
    if (cached) return cached;
    const matrix = new Float64Array(dimension * dimension);
    for (let row = 0; row < dimension; row += 1) {
      for (let column = 0; column < dimension; column += 1) {
        let value = 0;
        for (let component = 0; component < dimension; component += 1) {
          value += eigensystem.vectors[row * dimension + component]
            * Math.exp(eigensystem.values[component] * length)
            * eigensystem.vectors[column * dimension + component];
        }
        matrix[row * dimension + column] = Math.max(0, value * Math.sqrt(frequencies[column] / frequencies[row]));
      }
      let rowTotal = 0;
      for (let column = 0; column < dimension; column += 1) rowTotal += matrix[row * dimension + column];
      if (rowTotal <= 0) matrix[row * dimension + row] = 1;
      else for (let column = 0; column < dimension; column += 1) matrix[row * dimension + column] /= rowTotal;
    }
    cache.set(key, matrix);
    return matrix;
  };
  return { dimension, frequencies, rateMatrix, transition };
}

/** Backwards-compatible explicit five-character compiler. */
export function compileGtr5(model: PhyloUcaGtrModel): ReversibleCharacterModel {
  return compileGtr(model, true);
}

/** Likelihood partial for one observed fixed-alignment character. */
export function observedCharacterPartial(character: string, dimension: 4 | 5): Float64Array {
  const normalized = character.toUpperCase().replace("U", "T").replace(".", "-");
  const exact = PHYLO_UCA_CHARACTERS.indexOf(normalized as PhyloUcaCharacter);
  const partial = new Float64Array(dimension);
  if (exact >= 0 && exact < dimension) {
    partial[exact] = 1;
    return partial;
  }
  // IUPAC ambiguity refers to unknown nucleotides, not an unknown gap state.
  const sets: Record<string, string> = {
    N: "ACGT", R: "AG", Y: "CT", K: "GT", M: "AC", S: "CG", W: "AT",
    B: "CGT", D: "AGT", H: "ACT", V: "ACG",
  };
  const allowed = sets[normalized];
  if (allowed) {
    for (const base of allowed) partial[PHYLO_UCA_CHARACTERS.indexOf(base as PhyloUcaCharacter)] = 1;
    return partial;
  }
  // Under an explicit GTR4 override, a gap is missing rather than a fifth
  // character. Auto mode takes this path only when an internal gap is present.
  partial.fill(1);
  return partial;
}

export function transportLikelihood(
  transition: Float64Array,
  source: ArrayLike<number>,
  output = new Float64Array(Math.round(Math.sqrt(transition.length))),
): Float64Array {
  const dimension = Math.round(Math.sqrt(transition.length));
  for (let receiving = 0; receiving < dimension; receiving += 1) {
    let value = 0;
    for (let sourceState = 0; sourceState < dimension; sourceState += 1) {
      value += transition[receiving * dimension + sourceState] * source[sourceState];
    }
    output[receiving] = value;
  }
  return output;
}
