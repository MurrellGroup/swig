/// <reference lib="webworker" />

import {
  preprocessGermlineFastaAcrossTiers,
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
    const report: GermlinePreprocessReport = preprocessGermlineFastaAcrossTiers(text, segment, templateTiers, allowedLoci);
    self.postMessage({ report });
  } catch (error) {
    self.postMessage({ error: error instanceof Error ? error.message : String(error) });
  }
};

export {};
