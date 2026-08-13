import assert from "node:assert/strict";
import test from "node:test";

import { readAlivibeNucleotideFasta } from "../src/alivibe-roundtrip.ts";

test("Alivibe return captures its frame and exports NT even when the visible view is AA", () => {
  let mode: "nt" | "aa" = "aa";
  const editor = {
    document: {
      getElementById(id: string) {
        if (id === "sel-frame") return { value: "2" };
        if (id === "btn-nt") return { click: () => { mode = "nt"; } };
        return null;
      },
    },
    getClipboardContent() {
      return mode === "nt" ? ">member\nA---ATGGCC\n" : ">member\n-MA\n";
    },
  };
  const returned = readAlivibeNucleotideFasta(editor as unknown as Parameters<typeof readAlivibeNucleotideFasta>[0], 0);
  assert.equal(returned.frameOffset, 1);
  assert.equal(returned.fasta, ">member\nA---ATGGCC\n");
});
