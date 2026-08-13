import type { LinkedAirrBinding, SwigSession } from "./session-state";
import type { DatasetManifestEntry } from "./study-design";

export const SWIG_PROJECT_SCHEMA = 1 as const;

export interface ProjectRunManifest {
  id: string;
  name: string;
  createdAt: string;
  updatedAt: string;
  airrPath: string;
  latestSessionPath: string;
  checkpointDirectory: string;
  logPath: string;
  inputManifestPath?: string;
  checkpointCount: number;
  linkedAirr?: LinkedAirrBinding;
}

export interface SwigProjectManifest {
  schema: typeof SWIG_PROJECT_SCHEMA;
  application: "Swig";
  applicationVersion: string;
  projectName: string;
  createdAt: string;
  updatedAt: string;
  activeRunId: string | null;
  runs: ProjectRunManifest[];
}

export interface ProjectWorkspace {
  root: FileSystemDirectoryHandle;
  manifest: SwigProjectManifest;
}

export interface PreparedProjectRun {
  workspace: ProjectWorkspace;
  run: ProjectRunManifest;
  airrHandle: FileSystemFileHandle;
  writable: FileSystemWritableFileStream;
}

type DirectoryPicker = (options?: { id?: string; mode?: "read" | "readwrite" }) => Promise<FileSystemDirectoryHandle>;

const MANIFEST_NAME = "swig-project.json";
const PROJECT_README = `Swig project directory
======================

swig-project.json records the active run and linked state files.
runs/NNN-name/results contains the incrementally written AIRR table.
runs/NNN-name/inputs/datasets.tsv records the editable study metadata used at annotation time.
runs/NNN-name/state/latest.swig-session.json.gz is the current restorable state.
runs/NNN-name/state/checkpoints contains numbered, stage-named state checkpoints.
runs/NNN-name/logs/events.jsonl is an append-only structured event log.

Load this directory through Swig's “Load project directory” control. Do not rename files inside an active run unless swig-project.json is updated accordingly.
`;

export function projectDirectoryPicker(): DirectoryPicker | undefined {
  return (window as Window & { showDirectoryPicker?: DirectoryPicker }).showDirectoryPicker;
}

export function projectDirectoriesSupported(): boolean {
  return typeof window !== "undefined" && window.isSecureContext && Boolean(projectDirectoryPicker());
}

function safeStem(value: string, fallback = "swig-run"): string {
  return value.trim().replace(/[^A-Za-z0-9_.-]+/g, "-").replace(/^-+|-+$/g, "").slice(0, 80) || fallback;
}

async function directoryAt(root: FileSystemDirectoryHandle, path: string, create = false): Promise<FileSystemDirectoryHandle> {
  let current = root;
  for (const part of path.split("/").filter(Boolean)) current = await current.getDirectoryHandle(part, { create });
  return current;
}

async function fileHandleAt(root: FileSystemDirectoryHandle, path: string, create = false): Promise<FileSystemFileHandle> {
  const parts = path.split("/").filter(Boolean);
  const name = parts.pop();
  if (!name) throw new Error("Project file path is empty.");
  return (await directoryAt(root, parts.join("/"), create)).getFileHandle(name, { create });
}

async function writeFile(root: FileSystemDirectoryHandle, path: string, value: FileSystemWriteChunkType): Promise<void> {
  const handle = await fileHandleAt(root, path, true);
  const writable = await handle.createWritable();
  try { await writable.write(value); await writable.close(); }
  catch (error) { await writable.abort(error).catch(() => undefined); throw error; }
}

async function readText(root: FileSystemDirectoryHandle, path: string): Promise<string> {
  return (await (await fileHandleAt(root, path)).getFile()).text();
}

async function maybeReadManifest(root: FileSystemDirectoryHandle): Promise<SwigProjectManifest | null> {
  try {
    const parsed = JSON.parse(await readText(root, MANIFEST_NAME)) as Partial<SwigProjectManifest>;
    if (parsed.schema !== SWIG_PROJECT_SCHEMA || parsed.application !== "Swig" || !Array.isArray(parsed.runs)) throw new Error("This directory contains an unsupported Swig project manifest.");
    return parsed as SwigProjectManifest;
  } catch (error) {
    if (error instanceof DOMException && error.name === "NotFoundError") return null;
    throw error;
  }
}

async function writeManifest(workspace: ProjectWorkspace): Promise<void> {
  workspace.manifest.updatedAt = new Date().toISOString();
  await writeFile(workspace.root, MANIFEST_NAME, `${JSON.stringify(workspace.manifest, null, 2)}\n`);
}

export async function selectProjectDirectory(): Promise<FileSystemDirectoryHandle> {
  const picker = projectDirectoryPicker();
  if (!picker) throw new Error("Project directories are unavailable in this browser. Use a Chromium-based browser on HTTPS, or continue with session files.");
  return picker.call(window, { id: "swig-project", mode: "readwrite" });
}

export async function attachProjectDirectory(root: FileSystemDirectoryHandle, applicationVersion: string): Promise<{ workspace: ProjectWorkspace; existing: boolean }> {
  const current = await maybeReadManifest(root);
  if (current) return { workspace: { root, manifest: current }, existing: true };
  const now = new Date().toISOString();
  const workspace: ProjectWorkspace = { root, manifest: { schema: SWIG_PROJECT_SCHEMA, application: "Swig", applicationVersion, projectName: root.name, createdAt: now, updatedAt: now, activeRunId: null, runs: [] } };
  await writeFile(root, "README.txt", PROJECT_README);
  await writeManifest(workspace);
  return { workspace, existing: false };
}

export function activeProjectRun(workspace: ProjectWorkspace): ProjectRunManifest | null {
  return workspace.manifest.runs.find((run) => run.id === workspace.manifest.activeRunId) ?? null;
}

export async function prepareProjectRun(workspace: ProjectWorkspace, name: string, applicationVersion: string): Promise<PreparedProjectRun> {
  const ordinal = workspace.manifest.runs.reduce((maximum, run) => Math.max(maximum, Number.parseInt(run.id, 10) || 0), 0) + 1;
  const stem = safeStem(name);
  const id = `${String(ordinal).padStart(3, "0")}-${stem}`;
  const base = `runs/${id}`;
  const now = new Date().toISOString();
  const run: ProjectRunManifest = { id, name, createdAt: now, updatedAt: now, airrPath: `${base}/results/${stem}.airr.tsv`, latestSessionPath: `${base}/state/latest.swig-session.json.gz`, checkpointDirectory: `${base}/state/checkpoints`, logPath: `${base}/logs/events.jsonl`, inputManifestPath: `${base}/inputs/datasets.tsv`, checkpointCount: 0 };
  await directoryAt(workspace.root, `${base}/results`, true);
  await directoryAt(workspace.root, `${base}/inputs`, true);
  await directoryAt(workspace.root, run.checkpointDirectory, true);
  await directoryAt(workspace.root, `${base}/logs`, true);
  workspace.manifest.applicationVersion = applicationVersion;
  workspace.manifest.runs.push(run);
  workspace.manifest.activeRunId = id;
  await writeManifest(workspace);
  await appendProjectLog(workspace, run, "run_created", { name });
  const airrHandle = await fileHandleAt(workspace.root, run.airrPath, true);
  return { workspace, run, airrHandle, writable: await airrHandle.createWritable() };
}

function tsv(value:unknown):string{return String(value??"").replace(/[\t\r\n]+/g," ");}

export async function writeProjectDatasetManifest(workspace:ProjectWorkspace,run:ProjectRunManifest,datasets:DatasetManifestEntry[]):Promise<void>{
  const path=run.inputManifestPath??`runs/${run.id}/inputs/datasets.tsv`;
  const lines=["dataset_id\tinput_path\tsample_id\tsubject_id\tcohort\ttimepoint\tcompartment\trecords",...datasets.map((dataset)=>[dataset.datasetId,dataset.inputName,dataset.sampleId,dataset.subjectId,dataset.cohort,dataset.timepoint,dataset.compartment??"",dataset.records??""].map(tsv).join("\t"))];
  await writeFile(workspace.root,path,`${lines.join("\n")}\n`);
  run.inputManifestPath=path;run.updatedAt=new Date().toISOString();await writeManifest(workspace);
  await appendProjectLog(workspace,run,"dataset_manifest_written",{datasets:datasets.length,path});
}

export async function appendProjectLog(workspace: ProjectWorkspace, run: ProjectRunManifest, event: string, detail: Record<string, unknown> = {}): Promise<void> {
  const handle = await fileHandleAt(workspace.root, run.logPath, true);
  const previous = await handle.getFile();
  const writable = await handle.createWritable({ keepExistingData: true });
  try {
    await writable.seek(previous.size);
    await writable.write(`${JSON.stringify({ time: new Date().toISOString(), event, ...detail })}\n`);
    await writable.close();
  } catch (error) { await writable.abort(error).catch(() => undefined); throw error; }
}

function checkpointStage(session: SwigSession): string {
  const post = session.postAnalysis;
  if (post.tree) return "phylogeny";
  if (post.alignment || post.editedAlignments?.length) return "lineage-alignment";
  if (post.missingAlleles) return "missing-alleles";
  if (post.shm) return "shm";
  if (post.lineage) return "lineages";
  if (post.workingStages.some((stage) => stage.id === "selection")) return "repertoire-selection";
  if (post.chimera) return "chimera-filter";
  if (post.collapse) return "collapse-denoise";
  return "annotation";
}

export async function saveProjectCheckpoint(workspace: ProjectWorkspace, session: SwigSession, encoded: Blob, reason = "state_changed"): Promise<ProjectRunManifest> {
  const run = activeProjectRun(workspace);
  if (!run) throw new Error("No active run exists in this project directory.");
  const stage = checkpointStage(session);
  const revision = run.checkpointCount + 1;
  const checkpointName = `${String(revision).padStart(4, "0")}-${stage}.swig-session.json.gz`;
  await writeFile(workspace.root, `${run.checkpointDirectory}/${checkpointName}`, encoded);
  await writeFile(workspace.root, run.latestSessionPath, encoded);
  run.checkpointCount = revision;
  run.updatedAt = new Date().toISOString();
  run.linkedAirr = { ...session.linkedAirr };
  workspace.manifest.applicationVersion = session.applicationVersion;
  await writeManifest(workspace);
  await appendProjectLog(workspace, run, "state_checkpoint", { revision, stage, reason, records: session.linkedAirr.records, checkpoint: checkpointName });
  return run;
}

export async function loadActiveProjectFiles(workspace: ProjectWorkspace): Promise<{ sessionFile: File; airrFile: File; run: ProjectRunManifest }> {
  const run = activeProjectRun(workspace);
  if (!run) throw new Error("This project has no completed or active Swig run.");
  try {
    const [sessionFile, airrFile] = await Promise.all([
      (await fileHandleAt(workspace.root, run.latestSessionPath)).getFile(),
      (await fileHandleAt(workspace.root, run.airrPath)).getFile(),
    ]);
    return { sessionFile, airrFile, run };
  } catch (error) {
    if (error instanceof DOMException && error.name === "NotFoundError") throw new Error("The project manifest points to a missing latest session or AIRR result file.");
    throw error;
  }
}
