export interface TreeNode {
  name: string;
  length: number;
  children: TreeNode[];
}

export interface TreeLayoutNode extends TreeNode {
  x: number;
  y: number;
}

export interface TreeLayout {
  nodes: TreeLayoutNode[];
  edges: Array<{ parent: TreeLayoutNode; child: TreeLayoutNode }>;
  width: number;
  height: number;
  leaves: number;
  mode: TreeLayoutMode;
}

export type TreeLayoutMode = "phylogram" | "cladogram";

export const FASTTREE_DOUBLE_MINIMUM_BRANCH = 5e-9;
export const FASTTREE_AMBIGUOUS_BRANCH_THRESHOLD = 1e-8;

export interface CollapsedTree {
  root: TreeNode;
  collapsedEdges: number;
  threshold: number;
}

export function extractNewick(output: string): string {
  const text = String(output).replace(/\x1b\[[0-9;]*m/g, "");
  let best = "";
  for (let start = 0; start < text.length; start += 1) {
    if (text[start] !== "(") continue;
    let depth = 0;
    let quote = "";
    let commentDepth = 0;
    for (let index = start; index < text.length; index += 1) {
      const value = text[index];
      if (quote) {
        if (value === quote && text[index - 1] !== "\\") quote = "";
        continue;
      }
      if (value === "'" || value === '"') {
        quote = value;
        continue;
      }
      if (value === "[") {
        commentDepth += 1;
        continue;
      }
      if (value === "]" && commentDepth) {
        commentDepth -= 1;
        continue;
      }
      if (commentDepth) continue;
      if (value === "(") depth += 1;
      else if (value === ")") {
        depth -= 1;
        if (depth !== 0) continue;
        let end = index + 1;
        while (end < text.length && text[end] !== ";" && text[end] !== "\n" && text[end] !== "\r") end += 1;
        if (text[end] !== ";") break;
        end += 1;
        const candidate = text.slice(start, end).trim();
        if (candidate.includes(",") && candidate.length > best.length) best = candidate;
        break;
      }
      if (depth < 0) break;
    }
  }
  if (!best) throw new Error("FastTree did not return a complete Newick tree.");
  return `${best.replace(/;+$/, "")};`;
}

export function parseNewick(text: string): TreeNode {
  const source = extractNewick(text).replace(/\[[^\]]*\]/g, "");
  let position = 0;
  const whitespace = () => {
    while (/\s/.test(source[position] ?? "")) position += 1;
  };
  const token = () => {
    whitespace();
    if (source[position] === "'" || source[position] === '"') {
      const quote = source[position++];
      let value = "";
      while (position < source.length && source[position] !== quote) value += source[position++];
      if (source[position] === quote) position += 1;
      whitespace();
      return value;
    }
    const start = position;
    while (position < source.length && !"(),:;".includes(source[position]) && !/\s/.test(source[position])) position += 1;
    const value = source.slice(start, position).trim();
    whitespace();
    return value;
  };
  const length = () => {
    whitespace();
    if (source[position] !== ":") return 0;
    position += 1;
    const value = token();
    const parsed = Number(value);
    return Number.isFinite(parsed) ? parsed : 0;
  };
  const parse = (): TreeNode => {
    whitespace();
    const children: TreeNode[] = [];
    if (source[position] === "(") {
      position += 1;
      while (position < source.length) {
        children.push(parse());
        if (source[position] === ",") {
          position += 1;
          continue;
        }
        if (source[position] === ")") {
          position += 1;
          break;
        }
        throw new Error("Could not parse the complete Newick tree.");
      }
    }
    const name = token();
    const branchLength = length();
    return { name, length: branchLength, children };
  };
  const root = parse();
  whitespace();
  if (source[position] === ";") position += 1;
  whitespace();
  if (position < source.length) throw new Error("Could not parse the complete Newick tree.");
  return root;
}

export function serializeNewick(node: TreeNode): string {
  const children = node.children.length ? `(${node.children.map(serializeNewick).join(",")})` : "";
  const name = node.name.replace(/[\s():;,]/g, "_");
  const length = node.length > 0 ? `:${node.length.toPrecision(8).replace(/0+$/, "").replace(/\.$/, "")}` : "";
  return `${children}${name}${length}`;
}

export function rootOnOutgroup(root: TreeNode, outgroupName: string): TreeNode {
  const adjacency = new Map<TreeNode, Array<{ node: TreeNode; length: number }>>();
  let outgroup: TreeNode | null = null;
  const visit = (node: TreeNode, parent?: TreeNode) => {
    if (!node.children.length && node.name === outgroupName) outgroup = node;
    if (!adjacency.has(node)) adjacency.set(node, []);
    for (const child of node.children) {
      adjacency.get(node)!.push({ node: child, length: child.length });
      adjacency.set(child, [{ node, length: child.length }]);
      visit(child, node);
    }
    void parent;
  };
  visit(root);
  if (!outgroup) throw new Error(`The tree does not contain the outgroup ${outgroupName}.`);
  const neighbor = adjacency.get(outgroup)?.[0];
  if (!neighbor) return root;
  const orient = (node: TreeNode, parent: TreeNode, branchLength: number): TreeNode => ({
    name: node.name,
    length: branchLength,
    children: (adjacency.get(node) ?? []).filter((edge) => edge.node !== parent).map((edge) => orient(edge.node, node, edge.length)),
  });
  const half = neighbor.length / 2;
  return {
    name: "",
    length: 0,
    children: [orient(outgroup, neighbor.node, half), orient(neighbor.node, outgroup, half)],
  };
}

/**
 * Return a deterministic child order without changing any splits or lengths.
 * Newick child order is arbitrary, but a stable order prevents equivalent
 * trees from jumping around in the lineage viewer between runtimes.
 */
export function canonicalizeTree(root: TreeNode): TreeNode {
  const visit = (node: TreeNode): { node: TreeNode; key: string } => {
    if (!node.children.length) return { node: { ...node, children: [] }, key: node.name };
    const children = node.children.map(visit).sort((left, right) => left.key.localeCompare(right.key));
    return {
      node: { ...node, children: children.map((child) => child.node) },
      key: children.map((child) => child.key).sort().join("\0"),
    };
  };
  return visit(root).node;
}

/**
 * Collapse internal branches at FastTreeDbl's numerical floor. FastTree can
 * resolve these zero-information edges differently across native and WASM
 * floating-point targets. Leaves are never collapsed, and the untouched raw
 * FastTree Newick remains available separately for audit.
 */
export function collapseShortInternalBranches(
  root: TreeNode,
  threshold = FASTTREE_AMBIGUOUS_BRANCH_THRESHOLD,
): CollapsedTree {
  if (!Number.isFinite(threshold) || threshold < 0) throw new Error("The branch-collapse threshold must be a non-negative number.");
  let collapsedEdges = 0;
  const visit = (node: TreeNode): TreeNode => {
    const children = node.children.flatMap((child) => {
      const next = visit(child);
      if (next.children.length && next.length <= threshold) {
        collapsedEdges += 1;
        return next.children;
      }
      return [next];
    });
    return { ...node, children };
  };
  return { root: canonicalizeTree(visit(root)), collapsedEdges, threshold };
}

export function layoutTree(root: TreeNode, width = 900, rowHeight = 24, requestedMode: TreeLayoutMode = "phylogram"): TreeLayout {
  let hasPositiveLength = false;
  const detectLengths = (node: TreeNode) => {
    for (const child of node.children) {
      if (child.length > 0) hasPositiveLength = true;
      detectLengths(child);
    }
  };
  detectLengths(root);
  const mode: TreeLayoutMode = requestedMode === "phylogram" && hasPositiveLength ? "phylogram" : "cladogram";
  const distances = new Map<TreeNode, number>();
  let maximum = 0;
  const measure = (node: TreeNode, distance: number) => {
    distances.set(node, distance);
    maximum = Math.max(maximum, distance);
    // A zero-length FastTree edge is biologically meaningful. Never inflate it
    // to one unit in a phylogram; use unit depth only in explicit cladogram mode.
    node.children.forEach((child) => measure(child, distance + (mode === "phylogram" ? Math.max(0, child.length) : 1)));
  };
  measure(root, 0);
  const leaves = [...distances.keys()].filter((node) => !node.children.length);
  const leafY = new Map(leaves.map((leaf, index) => [leaf, 28 + index * rowHeight]));
  const nodes: TreeLayoutNode[] = [];
  const edges: Array<{ parent: TreeLayoutNode; child: TreeLayoutNode }> = [];
  const build = (node: TreeNode): TreeLayoutNode => {
    const children = node.children.map(build);
    const y = children.length ? children.reduce((sum, child) => sum + child.y, 0) / children.length : leafY.get(node) ?? 0;
    const layout: TreeLayoutNode = {
      ...node,
      children,
      x: 24 + (distances.get(node) ?? 0) / Math.max(maximum, 1) * (width - 220),
      y,
    };
    nodes.push(layout);
    children.forEach((child) => edges.push({ parent: layout, child }));
    return layout;
  };
  build(root);
  return { nodes, edges, width, height: Math.max(90, leaves.length * rowHeight + 50), leaves: leaves.length, mode };
}
