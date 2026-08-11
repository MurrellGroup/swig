import type { CSSProperties } from "react";

export type SequenceAlphabet = "nt" | "aa";

const NT_COLORS: Record<string, string> = {
  A: "#81d4fa", G: "#fff176", C: "#a5d6a7", T: "#ff8a80", U: "#ff8a80",
  R: "#e6ee9c", Y: "#80cbc4", M: "#a5d6a7", K: "#ef9a9a", S: "#c5e1a5",
  W: "#ce93d8", H: "#81d4fa", B: "#ef9a9a", V: "#fff59d", D: "#ffcc80",
  N: "#eeeeee", "-": "#ffffff",
};

const AA_COLORS: Record<string, string> = {
  A: "#80a0f0", R: "#f01505", N: "#00ff00", D: "#c048c0", C: "#f08080",
  Q: "#00ff00", E: "#c048c0", G: "#f09048", H: "#15a4a4", I: "#80a0f0",
  L: "#80a0f0", K: "#f01505", M: "#80a0f0", F: "#80a0f0", P: "#ffff00",
  S: "#00ff00", T: "#00ff00", W: "#80a0f0", Y: "#15a4a4", V: "#80a0f0",
  X: "#eeeeee", "*": "#999999", "-": "#ffffff",
};

export function sequenceColor(value: string, alphabet: SequenceAlphabet): string {
  return (alphabet === "nt" ? NT_COLORS : AA_COLORS)[value.toUpperCase()] ?? "#eeeeee";
}

function foreground(background: string): string {
  const red = Number.parseInt(background.slice(1, 3), 16);
  const green = Number.parseInt(background.slice(3, 5), 16);
  const blue = Number.parseInt(background.slice(5, 7), 16);
  return red * 0.299 + green * 0.587 + blue * 0.114 < 125 ? "#ffffff" : "#10201c";
}

export function ColoredSequence({ sequence, alphabet, className = "", style }: {
  sequence: string;
  alphabet: SequenceAlphabet;
  className?: string;
  style?: CSSProperties;
}) {
  return <code className={`colored-sequence ${className}`.trim()} style={style} aria-label={sequence}>{[...sequence].map((value, index) => {
    const background = sequenceColor(value, alphabet);
    return <span key={index} style={{ backgroundColor: background, color: foreground(background) }}>{value}</span>;
  })}</code>;
}
