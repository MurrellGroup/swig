import react from "@vitejs/plugin-react";
import { defineConfig } from "vite";

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
  plugins: [react()],
  build: {
    target: "es2022",
    emptyOutDir: true,
  },
});
