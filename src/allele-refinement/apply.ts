export interface CallOverrideVector {
  labels: string[];
  mapNode: Int32Array;
  probability: Float32Array;
}

export interface MutableGermlineCalls {
  ordinal: number;
  vCall: string;
  jCall: string;
  originalVCall?: string;
  originalJCall?: string;
}

/** Restores immutable AIRR calls, then applies thresholded V/J posterior MAP calls. */
export function applyCallOverrides<T extends MutableGermlineCalls>(
  records: T[],
  v: CallOverrideVector | undefined,
  j: CallOverrideVector | undefined,
  minimumPosterior: number,
  intern: (value: string) => string = (value) => value,
): { changedV: number; changedJ: number; threshold: number } {
  const threshold = Math.max(0, Math.min(1, minimumPosterior));
  let changedV = 0;
  let changedJ = 0;
  for (const record of records) {
    const ordinal = record.ordinal;
    record.vCall = record.originalVCall ?? record.vCall;
    record.jCall = record.originalJCall ?? record.jCall;
    const vNode = v?.mapNode[ordinal] ?? -1;
    if (vNode >= 0 && (v?.probability[ordinal] ?? 0) >= threshold && v?.labels[vNode]) {
      const refined = intern(v.labels[vNode]);
      if (refined !== record.vCall) changedV += 1;
      record.vCall = refined;
    }
    const jNode = j?.mapNode[ordinal] ?? -1;
    if (jNode >= 0 && (j?.probability[ordinal] ?? 0) >= threshold && j?.labels[jNode]) {
      const refined = intern(j.labels[jNode]);
      if (refined !== record.jCall) changedJ += 1;
      record.jCall = refined;
    }
  }
  return { changedV, changedJ, threshold };
}
