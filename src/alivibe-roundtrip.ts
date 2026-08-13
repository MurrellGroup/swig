import type { AlignmentFrameOffset } from "./lineage-phylogeny.ts";

export interface AlivibeEditorBridge {
  document: Pick<Document, "getElementById">;
  getClipboardContent?: (preferSelection?: boolean) => string;
}

export interface AlivibeNucleotideTransfer {
  fasta: string;
  frameOffset: AlignmentFrameOffset;
}

function validOffset(value: number): AlignmentFrameOffset | undefined {
  return value === 0 || value === 1 || value === 2 ? value : undefined;
}

/** Capture Alivibe's AA phase, then force NT before asking for current-view FASTA. */
export function readAlivibeNucleotideFasta(
  editor: AlivibeEditorBridge,
  fallbackFrameOffset: AlignmentFrameOffset,
): AlivibeNucleotideTransfer {
  const frameValue = Number((editor.document.getElementById("sel-frame") as HTMLSelectElement | null)?.value) - 1;
  const frameOffset = validOffset(frameValue) ?? fallbackFrameOffset;
  (editor.document.getElementById("btn-nt") as HTMLButtonElement | null)?.click();
  return { fasta: editor.getClipboardContent?.(false) ?? "", frameOffset };
}
