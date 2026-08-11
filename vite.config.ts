import react from "@vitejs/plugin-react";
import { readFileSync } from "node:fs";
import { defineConfig } from "vite";

const IMMUNUM_BROWSER_ID = "\0swig-immunum-browser";

function immunumBrowser() {
  return {
    name: "swig-immunum-browser",
    resolveId(id: string) {
      return id === "virtual:immunum-browser" ? IMMUNUM_BROWSER_ID : null;
    },
    load(id: string) {
      if (id !== IMMUNUM_BROWSER_ID) return null;
      const source = readFileSync(new URL("./node_modules/immunum/immunum.js", import.meta.url), "utf8")
        .replace("exports.Annotator = Annotator;", "")
        .replace(
          /const wasmPath = `[\s\S]*?wasm\.__wbindgen_start\(\);\s*$/,
          `let wasm;\nlet initialization;\nexport async function initializeImmunum() {\n  if (wasm) return;\n  if (!initialization) initialization = (async () => {\n    const response = await fetch(wasmUrl);\n    const imports = __wbg_get_imports();\n    let instance;\n    try {\n      ({ instance } = await WebAssembly.instantiateStreaming(response.clone(), imports));\n    } catch {\n      ({ instance } = await WebAssembly.instantiate(await response.arrayBuffer(), imports));\n    }\n    wasm = instance.exports;\n    wasm.__wbindgen_start();\n  })();\n  await initialization;\n}\nexport { Annotator };\n`,
        );
      return `import wasmUrl from "immunum/immunum_bg.wasm?url";\n${source}`;
    },
  };
}

function githubPagesBase(): string {
  if (process.env.GITHUB_ACTIONS !== "true") return "/";
  const repository = process.env.GITHUB_REPOSITORY?.split("/")[1];
  if (!repository || repository.endsWith(".github.io")) return "/";
  return `/${repository}/`;
}

export default defineConfig({
  base: githubPagesBase(),
  server: {
    host: "0.0.0.0",
    allowedHosts: ["terminal.local"],
  },
  plugins: [immunumBrowser(), react()],
  build: {
    target: "es2022",
    emptyOutDir: true,
  },
});
