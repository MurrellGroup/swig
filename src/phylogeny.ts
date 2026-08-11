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
}

export function parseNewick(text: string): TreeNode {
  const source = text.trim().replace(/;$/, "");
  let position = 0;
  const token = () => {
    const start = position;
    while (position < source.length && !"(),:;".includes(source[position])) position += 1;
    return source.slice(start, position).trim();
  };
  const length = () => {
    if (source[position] !== ":") return 0;
    position += 1;
    const value = token();
    const parsed = Number(value);
    return Number.isFinite(parsed) ? parsed : 0;
  };
  const parse = (): TreeNode => {
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
      }
    }
    const name = token();
    const branchLength = length();
    return { name, length: branchLength, children };
  };
  const root = parse();
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

export function layoutTree(root: TreeNode, width = 900, rowHeight = 24): TreeLayout {
  const distances = new Map<TreeNode, number>();
  let maximum = 0;
  const measure = (node: TreeNode, distance: number) => {
    distances.set(node, distance);
    maximum = Math.max(maximum, distance);
    node.children.forEach((child) => measure(child, distance + (child.length > 0 ? child.length : 1)));
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
  return { nodes, edges, width, height: Math.max(90, leaves.length * rowHeight + 50), leaves: leaves.length };
}
