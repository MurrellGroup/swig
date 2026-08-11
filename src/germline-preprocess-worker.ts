/// <reference lib="webworker" />

import {
  preprocessGermlineFasta,
  type GermlineLocus,
  type GermlinePreprocessReport,
  type GermlineSegment,
  type MetadataAllele,
} from "./germline-preprocess";

interface PreprocessRequest {
  text: string;
  segment: GermlineSegment;
  templateTiers: MetadataAllele[][];
  allowedLoci: GermlineLocus[];
}

self.onmessage = (event: MessageEvent<PreprocessRequest>) => {
  try {
    const { text, segment, templateTiers, allowedLoci } = event.data;
    let report: GermlinePreprocessReport | undefined;
    for (const templates of templateTiers.length ? templateTiers : [[]]) {
      report = preprocessGermlineFasta(report?.fasta ?? text, segment, templates, allowedLoci);
      if (segment !== "V" && segment !== "J") break;
      if (report.annotated === report.count) break;
    }
    self.postMessage({ report });
  } catch (error) {
    self.postMessage({ error: error instanceof Error ? error.message : String(error) });
  }
};

export {};
