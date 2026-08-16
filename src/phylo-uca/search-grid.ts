import type { PhyloUcaSearchOptions } from "./types.ts";

/**
 * Explicit zero-plus-logarithmic pendant-length grid used by the full HMM.
 *
 * A logarithmic positive part is essential here: a single junction mutation
 * in a few-hundred-column alignment commonly puts the useful scale around
 * 1e-3--1e-2 substitutions/site, even when the allowed maximum is much larger.
 */
export function phyloUcaBranchLengthGrid(options: Pick<PhyloUcaSearchOptions,
  "branchGridPoints" | "minimumPositiveUcaBranchLength" | "maximumUcaBranchLength"
>): number[] {
  const points = Math.max(2, Math.floor(options.branchGridPoints || 2));
  const maximum = Math.max(0, Number.isFinite(options.maximumUcaBranchLength) ? options.maximumUcaBranchLength : 0.3);
  if (!(maximum > 0)) return [0];
  const minimum = Math.min(maximum, Math.max(1e-9, Number.isFinite(options.minimumPositiveUcaBranchLength)
    ? options.minimumPositiveUcaBranchLength
    : Math.min(1e-5, maximum)));
  if (points === 2 || minimum === maximum) return [0, maximum];
  const positive = Array.from({ length: points - 1 }, (_, index) => {
    if (index === 0) return minimum;
    if (index === points - 2) return maximum;
    const fraction = index / (points - 2);
    return Math.exp(Math.log(minimum) + fraction * (Math.log(maximum) - Math.log(minimum)));
  });
  return [0, ...positive].filter((value, index, values) => index === 0 || Math.abs(value - values[index - 1]) > 1e-14);
}

/** Trapezoid/Voronoi cell widths for a sorted one-dimensional quadrature grid. */
export function gridCellWidths(values: readonly number[], lower: number, upper: number): number[] {
  if (!values.length) return [];
  if (values.length === 1) return [Math.max(1e-300, upper - lower)];
  return values.map((value, index) => {
    const left = index === 0 ? lower : (values[index - 1] + value) / 2;
    const right = index === values.length - 1 ? upper : (value + values[index + 1]) / 2;
    return Math.max(1e-300, right - left);
  });
}
