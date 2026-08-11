/* SWIG streams IndexedDB-backed AIRR tables into the browser download manager. */
self.addEventListener("install", () => self.skipWaiting());
self.addEventListener("activate", (event) => event.waitUntil(self.clients.claim()));

function requestResult(request) {
  return new Promise((resolve, reject) => {
    request.onsuccess = () => resolve(request.result);
    request.onerror = () => reject(request.error || new Error("The local AIRR index could not be read."));
  });
}

function openDatabase(name) {
  return new Promise((resolve, reject) => {
    const request = indexedDB.open(name);
    request.onsuccess = () => resolve(request.result);
    request.onerror = () => reject(request.error || new Error("The local AIRR index could not be opened."));
    request.onupgradeneeded = () => {
      request.transaction.abort();
      reject(new Error("The requested AIRR index does not exist."));
    };
  });
}

async function readRecord(database, storeName, key) {
  const transaction = database.transaction(storeName, "readonly");
  return requestResult(transaction.objectStore(storeName).get(key));
}

async function downloadResponse(url) {
  const databaseName = url.searchParams.get("database") || "";
  if (!databaseName.startsWith("swig-results-")) return new Response("Invalid result identifier.", { status: 400 });
  const requestedName = url.searchParams.get("name") || "swig.airr.tsv";
  const safeName = requestedName.replace(/[^a-zA-Z0-9._-]+/g, "_").slice(0, 180) || "swig.airr.tsv";
  let database;
  try {
    database = await openDatabase(databaseName);
    const manifest = await readRecord(database, "meta", "manifest");
    if (!manifest) throw new Error("The AIRR index is not finalized.");
    const encoder = new TextEncoder();
    let chunkIndex = -1;
    let chunkReader = null;
    const body = new ReadableStream({
      async pull(controller) {
        try {
          while (true) {
            if (chunkReader) {
              const { done, value } = await chunkReader.read();
              if (!done) {
                controller.enqueue(value);
                return;
              }
              chunkReader.releaseLock();
              chunkReader = null;
              chunkIndex += 1;
            }
            if (chunkIndex < 0) {
              controller.enqueue(encoder.encode(`${manifest.headerLine}\n`));
              chunkIndex = 0;
              return;
            }
            if (chunkIndex >= manifest.chunks) {
              controller.close();
              database.close();
              return;
            }
            const chunk = await readRecord(database, "chunks", chunkIndex);
            if (!chunk || chunk.storage !== "indexed" || !chunk.data) {
              throw new Error(`AIRR batch ${chunkIndex + 1} is unavailable.`);
            }
            let stream = chunk.data.stream();
            if (chunk.compressed) stream = stream.pipeThrough(new DecompressionStream("gzip"));
            chunkReader = stream.getReader();
          }
        } catch (error) {
          database.close();
          controller.error(error);
        }
      },
      async cancel() {
        await chunkReader?.cancel();
        database.close();
      },
    });
    return new Response(body, {
      headers: {
        "Content-Type": "text/tab-separated-values; charset=utf-8",
        "Content-Disposition": `attachment; filename="${safeName}"`,
        "Cache-Control": "no-store",
      },
    });
  } catch (error) {
    database?.close();
    return new Response(error instanceof Error ? error.message : String(error), { status: 404 });
  }
}

self.addEventListener("fetch", (event) => {
  const url = new URL(event.request.url);
  if (event.request.method === "GET" && url.origin === self.location.origin && url.pathname.endsWith("/__swig_download__")) {
    event.respondWith(downloadResponse(url));
  }
});
