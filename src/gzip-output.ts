export type OutputCompression = "none" | "gzip";

export type OutputPart = string | Blob | Uint8Array;

export interface OutputSink {
  write: (value: OutputPart) => Promise<void>;
  close?: () => Promise<void>;
  abort?: (reason?: unknown) => Promise<void>;
}

export interface OutputWriter extends OutputSink {
  close: () => Promise<void>;
  abort: (reason?: unknown) => Promise<void>;
}

const encoder = new TextEncoder();

async function writeBlob(
  writer: WritableStreamDefaultWriter<BufferSource>,
  blob: Blob,
): Promise<void> {
  const reader = blob.stream().getReader();
  try {
    while (true) {
      const { done, value } = await reader.read();
      if (done) return;
      await writer.write(value);
    }
  } finally {
    reader.releaseLock();
  }
}

/**
 * Wrap a file or memory sink in a bounded native gzip stream. Blob inputs are
 * read incrementally, which is important when an AIRR result is already backed
 * by a large direct-to-disk file.
 */
export function outputWriter(
  sink: OutputSink,
  compression: OutputCompression,
): OutputWriter {
  if (compression === "none") {
    return {
      write: sink.write,
      close: async () => { await sink.close?.(); },
      abort: async (reason) => { await sink.abort?.(reason); },
    };
  }

  if (typeof CompressionStream === "undefined") {
    throw new Error("Gzip export is unavailable in this browser.");
  }
  const compressor = new CompressionStream("gzip");
  const input = compressor.writable.getWriter();
  const output = compressor.readable.getReader();
  const pump = (async () => {
    while (true) {
      const { done, value } = await output.read();
      if (done) return;
      await sink.write(value);
    }
  })();
  let finished = false;

  return {
    write: async (value) => {
      if (finished) throw new Error("Cannot write after the gzip output has closed.");
      if (typeof value === "string") await input.write(encoder.encode(value));
      else if (value instanceof Blob) await writeBlob(input, value);
      else await input.write(Uint8Array.from(value));
    },
    close: async () => {
      if (finished) return;
      finished = true;
      await input.close();
      await pump;
      await sink.close?.();
    },
    abort: async (reason) => {
      if (finished) return;
      finished = true;
      try { await input.abort(reason); } catch { /* Preserve the original failure. */ }
      try { await output.cancel(reason); } catch { /* Preserve the original failure. */ }
      try { await pump; } catch { /* Preserve the original failure. */ }
      await sink.abort?.(reason);
    },
  };
}

export function compressedName(name: string, compression: OutputCompression): string {
  return compression === "gzip" && !name.toLowerCase().endsWith(".gz") ? `${name}.gz` : name;
}

export function compressedExtension(extension: string, compression: OutputCompression): string {
  return compression === "gzip" ? ".gz" : extension;
}

export function compressedMime(mime: string, compression: OutputCompression): string {
  return compression === "gzip" ? "application/gzip" : mime;
}
