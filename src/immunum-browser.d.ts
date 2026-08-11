declare module "virtual:immunum-browser" {
  export interface NumberingResult {
    chain: "H" | "K" | "L" | "A" | "B" | "D" | "G" | null;
    scheme: "imgt" | "kabat" | null;
    confidence: number | null;
    numbering: Map<string, string> | null;
    query_start: number | null;
    query_end: number | null;
    error: string | null;
  }

  export class Annotator {
    constructor(chains: string[], scheme: string, minimumConfidence?: number | null);
    number(sequence: string): NumberingResult;
    free(): void;
  }

  export function initializeImmunum(): Promise<void>;
}
