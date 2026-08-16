#!/usr/bin/env node

import { runCli } from "./swig-cli.mjs";

runCli().catch((error)=>{
  process.stderr.write(`swig-cli: ${error instanceof Error?error.message:String(error)}\n`);
  process.exitCode=1;
});
