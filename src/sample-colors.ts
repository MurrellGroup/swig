import type { DatasetManifestEntry } from "./study-design.ts";

/** Color-blind-conscious categorical defaults (Okabe–Ito plus extensions). */
export const SAMPLE_COLOR_SEQUENCE = [
  "#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#6F4E7C", "#8C8C2B",
  "#2E8B57", "#B54A4A", "#5B6FB5", "#A66A2C", "#008C95", "#9A5AA5", "#6E7B3D", "#C45D91",
] as const;

export type SampleColorMap = Record<string, string>;

export function sampleIds(datasets: DatasetManifestEntry[]): string[] {
  const seen = new Set<string>();
  const values: string[] = [];
  for (const dataset of datasets) {
    const sample = dataset.sampleId.trim();
    if (!sample || seen.has(sample)) continue;
    seen.add(sample);
    values.push(sample);
  }
  return values;
}

export function createSampleColorMap(datasets: DatasetManifestEntry[], existing: SampleColorMap = {}): SampleColorMap {
  const result: SampleColorMap = {};
  sampleIds(datasets).forEach((sample, index) => {
    const retained = existing[sample];
    result[sample] = /^#[0-9a-f]{6}$/i.test(retained ?? "") ? retained : SAMPLE_COLOR_SEQUENCE[index % SAMPLE_COLOR_SEQUENCE.length];
  });
  return result;
}

export function sampleColor(sample: string, colors: SampleColorMap): string {
  return colors[sample] || "#70817b";
}

export function categoricalLineageColor(lineageId: number): string {
  if (!(lineageId > 0)) return "#70817b";
  return SAMPLE_COLOR_SEQUENCE[(lineageId * 7 + 3) % SAMPLE_COLOR_SEQUENCE.length];
}

export function contrastingText(background: string): string {
  const red = Number.parseInt(background.slice(1, 3), 16);
  const green = Number.parseInt(background.slice(3, 5), 16);
  const blue = Number.parseInt(background.slice(5, 7), 16);
  return red * 0.299 + green * 0.587 + blue * 0.114 < 142 ? "#ffffff" : "#17231f";
}
