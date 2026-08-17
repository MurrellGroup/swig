import {
  canonicalizeTree,
  collapseShortInternalBranches,
  parseNewick,
  serializeNewick,
  type TreeNode,
} from "./phylogeny.ts";

export interface UploadedLineageTree {
  /** Parsed tree with uploaded child order retained and tip names reconciled. */
  rawNewick: string;
  /** Deterministically ordered form used by the default viewer. */
  canonicalNewick: string;
  /** Optional numerical-floor-collapsed display form. */
  stableNewick: string;
  collapsedEdges: number;
  collapseThreshold: number;
  /** False when the ordinary display tree contains Swig's germline guide tip. */
  observedOnly: boolean;
}

function safeName(value: string): string {
  return value.replace(/[^A-Za-z0-9_.|*+\-]/g, "_");
}

function visitLeaves(node: TreeNode, action: (leaf: TreeNode, index: number) => void, counter = { value: 0 }): void {
  if (!node.children.length) {
    action(node, counter.value++);
    return;
  }
  node.children.forEach((child) => visitLeaves(child, action, counter));
}

function visitBranches(node: TreeNode, action: (node: TreeNode) => void): void {
  node.children.forEach((child) => {
    action(child);
    visitBranches(child, action);
  });
}

/**
 * Validate and reconcile a user Newick tree against Swig's current biological
 * alignment rows. Trees may use the exact alignment identifiers, the original
 * sequence identifiers before Swig's unique `__ordinal` suffix, or a complete
 * zero-based numeric tip set in alignment order.
 */
export function prepareUploadedLineageTree(text: string, observedNames: readonly string[], allowGermlineGuide = false): UploadedLineageTree {
  if (!text.trim()) throw new Error("The uploaded tree file is empty.");
  if (observedNames.length < 2) throw new Error("A lineage tree needs at least two observed alignment rows.");
  if (observedNames.includes("__germline_N_masked__")) throw new Error("Pass only biological alignment identifiers when validating an uploaded tree.");
  if (new Set(observedNames).size !== observedNames.length) throw new Error("The current lineage alignment contains duplicate biological identifiers.");

  let root: TreeNode;
  try {
    root = parseNewick(text);
  } catch (error) {
    throw new Error(`Could not read the uploaded Newick tree: ${error instanceof Error ? error.message : String(error)}`);
  }

  const exact = new Map(observedNames.map((name) => [name, name]));
  const safeExact = new Map<string, string>();
  const short = new Map<string, string[]>();
  observedNames.forEach((name) => {
    safeExact.set(safeName(name), name);
    const base = safeName(name.replace(/__\d+$/, ""));
    short.set(base, [...(short.get(base) ?? []), name]);
  });

  const uploadedLeaves: TreeNode[] = [];
  visitLeaves(root, (leaf) => uploadedLeaves.push(leaf));
  const includesGermlineGuide = uploadedLeaves.some((leaf) => leaf.name === "__germline_N_masked__");
  const expectedTipCount = observedNames.length + (includesGermlineGuide && allowGermlineGuide ? 1 : 0);
  if (includesGermlineGuide && !allowGermlineGuide) {
    throw new Error("The uploaded tree contains __germline_N_masked__. UCA inference accepts an observed-only tree with that synthetic guide removed.");
  }
  if (uploadedLeaves.length !== expectedTipCount) {
    throw new Error(`The uploaded tree has ${uploadedLeaves.length.toLocaleString()} tips, but the current alignment has ${observedNames.length.toLocaleString()} biological rows${allowGermlineGuide ? " (plus at most the optional germline guide)" : ""}.`);
  }
  const numericTips = uploadedLeaves.every((leaf) => /^\d+$/.test(leaf.name))
    && new Set(uploadedLeaves.map((leaf) => Number(leaf.name))).size === uploadedLeaves.length
    && uploadedLeaves.every((leaf) => Number(leaf.name) >= 0 && Number(leaf.name) < observedNames.length);
  const used = new Set<string>();
  uploadedLeaves.forEach((leaf) => {
    if (allowGermlineGuide && leaf.name === "__germline_N_masked__") {
      used.add(leaf.name);
      return;
    }
    const normalized = safeName(leaf.name);
    let resolved = exact.get(leaf.name) ?? safeExact.get(normalized);
    if (!resolved) {
      const candidates = short.get(normalized) ?? [];
      if (candidates.length === 1) resolved = candidates[0];
      else if (candidates.length > 1) throw new Error(`Tree tip ${leaf.name} is ambiguous because that sequence identifier occurs more than once. Use Swig's complete alignment identifier, including its __ordinal suffix.`);
    }
    if (!resolved && numericTips) resolved = observedNames[Number(leaf.name)];
    if (!resolved) throw new Error(`Tree tip ${leaf.name || "(unnamed)"} does not match a biological row in the current alignment.`);
    if (used.has(resolved)) throw new Error(`The uploaded tree contains the alignment row ${resolved} more than once.`);
    used.add(resolved);
    leaf.name = resolved;
  });
  const missing = observedNames.filter((name) => !used.has(name));
  if (missing.length) throw new Error(`The uploaded tree is missing ${missing.slice(0, 4).join(", ")}${missing.length > 4 ? "…" : ""}.`);

  let positiveBranches = 0;
  visitBranches(root, (node) => {
    if (!Number.isFinite(node.length) || node.length < 0) throw new Error("Every uploaded branch length must be a finite, non-negative number.");
    if (node.length > 0) positiveBranches += 1;
  });
  if (!positiveBranches) throw new Error("The uploaded tree needs branch lengths, with at least one positive length. UCA likelihoods use those lengths directly.");

  const rawNewick = `${serializeNewick(root)};`;
  const canonical = canonicalizeTree(root);
  const collapsed = collapseShortInternalBranches(canonical);
  return {
    rawNewick,
    canonicalNewick: `${serializeNewick(canonical)};`,
    stableNewick: `${serializeNewick(collapsed.root)};`,
    collapsedEdges: collapsed.collapsedEdges,
    collapseThreshold: collapsed.threshold,
    observedOnly: !includesGermlineGuide,
  };
}
