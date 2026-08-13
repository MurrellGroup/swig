export interface InputFileCandidate {
  file: File;
  /** Path relative to the selected or dropped ancestor, including its name. */
  relativePath: string;
  rootName: string;
  fromDirectory: boolean;
}

export interface DirectoryDonorSuggestion {
  relativePath: string;
  donor: string | null;
}

export interface DirectoryDonorPlan {
  /** Flat roots require an explicit same/separate donor decision. */
  flatRoots: string[];
  suggestions: DirectoryDonorSuggestion[];
}

interface EntryFile {
  isFile: true;
  isDirectory: false;
  name: string;
  file: (success: (file: File) => void, failure?: (error: DOMException) => void) => void;
}

interface EntryDirectory {
  isFile: false;
  isDirectory: true;
  name: string;
  createReader: () => { readEntries: (success: (entries: LegacyEntry[]) => void, failure?: (error: DOMException) => void) => void };
}

type LegacyEntry = EntryFile | EntryDirectory;

type TransferItemWithHandles = DataTransferItem & {
  getAsFileSystemHandle?: () => Promise<FileSystemHandle | null>;
  webkitGetAsEntry?: () => LegacyEntry | null;
};

function cleanPath(value: string): string {
  return value.replaceAll("\\", "/").split("/").filter(Boolean).join("/");
}

function safeIdentifier(value: string, fallback = "subject"): string {
  return value.trim().replace(/[^A-Za-z0-9_.-]+/g, "_").replace(/^_+|_+$/g, "") || fallback;
}

async function collectModernHandle(handle: FileSystemHandle, parent: string, rootName: string, result: InputFileCandidate[]): Promise<void> {
  const path = cleanPath(`${parent}/${handle.name}`);
  if (handle.kind === "file") {
    result.push({ file: await (handle as FileSystemFileHandle).getFile(), relativePath: path, rootName, fromDirectory: parent.length > 0 });
    return;
  }
  const directory = handle as FileSystemDirectoryHandle & { values: () => AsyncIterableIterator<FileSystemHandle> };
  for await (const child of directory.values()) await collectModernHandle(child, path, rootName, result);
}

function legacyFile(entry: EntryFile): Promise<File> {
  return new Promise((resolve, reject) => entry.file(resolve, reject));
}

function legacyChildren(entry: EntryDirectory): Promise<LegacyEntry[]> {
  const reader = entry.createReader();
  return new Promise((resolve, reject) => {
    const all: LegacyEntry[] = [];
    const next = () => reader.readEntries((entries) => {
      if (!entries.length) resolve(all);
      else { all.push(...entries); next(); }
    }, reject);
    next();
  });
}

async function collectLegacyEntry(entry: LegacyEntry, parent: string, rootName: string, result: InputFileCandidate[]): Promise<void> {
  const path = cleanPath(`${parent}/${entry.name}`);
  if (entry.isFile) {
    result.push({ file: await legacyFile(entry), relativePath: path, rootName, fromDirectory: parent.length > 0 });
    return;
  }
  for (const child of await legacyChildren(entry)) await collectLegacyEntry(child, path, rootName, result);
}

/**
 * Expands files and directories while the drop data store is still available.
 * Modern handles are preferred; the Entries API and FileList are bounded fallbacks.
 */
export async function collectDroppedInput(dataTransfer: DataTransfer): Promise<InputFileCandidate[]> {
  const items = [...dataTransfer.items].filter((item) => item.kind === "file") as TransferItemWithHandles[];
  const result: InputFileCandidate[] = [];
  if (items.some((item) => typeof item.getAsFileSystemHandle === "function")) {
    // Handle acquisition must begin during the drop event's protected data-store window.
    const handles = await Promise.all(items.map((item) => item.getAsFileSystemHandle?.() ?? Promise.resolve(null)));
    for (const handle of handles) {
      if (handle) await collectModernHandle(handle, "", handle.name, result);
    }
    if (result.length) return result;
  }
  if (items.some((item) => typeof item.webkitGetAsEntry === "function")) {
    for (const item of items) {
      const entry = (item as unknown as { webkitGetAsEntry?: () => LegacyEntry | null }).webkitGetAsEntry?.();
      if (entry) await collectLegacyEntry(entry, "", entry.name, result);
    }
    if (result.length) return result;
  }
  return [...dataTransfer.files].map((file) => {
    const relativePath = cleanPath(file.webkitRelativePath || file.name);
    const rootName = relativePath.split("/")[0] || file.name;
    return { file, relativePath, rootName, fromDirectory: relativePath.includes("/") };
  });
}

export function candidatesFromDirectoryPicker(files: File[]): InputFileCandidate[] {
  return files.map((file) => {
    const relativePath = cleanPath(file.webkitRelativePath || file.name);
    const rootName = relativePath.split("/")[0] || file.name;
    return { file, relativePath, rootName, fromDirectory: true };
  });
}

/**
 * The first directory beneath each selected/dropped root is the donor group.
 * Direct children of a flat root remain unresolved for the confirmation dialog.
 */
export function inferDirectoryDonors(candidates: Array<Pick<InputFileCandidate, "relativePath" | "rootName" | "fromDirectory">>): DirectoryDonorPlan {
  const flatRoots = new Set<string>();
  const suggestions = candidates.map((candidate) => {
    if (!candidate.fromDirectory) return { relativePath: candidate.relativePath, donor: null };
    const parts = cleanPath(candidate.relativePath).split("/");
    const rootIndex = parts[0] === candidate.rootName ? 1 : 0;
    const nested = parts.length - rootIndex > 1;
    if (!nested) {
      flatRoots.add(candidate.rootName);
      return { relativePath: candidate.relativePath, donor: null };
    }
    return { relativePath: candidate.relativePath, donor: safeIdentifier(parts[rootIndex], "subject") };
  });
  return { flatRoots: [...flatRoots].sort(), suggestions };
}

export function donorForFlatRoot(rootName: string): string {
  return safeIdentifier(rootName, "subject");
}
