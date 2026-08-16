import { useMemo } from "react";

import { layoutTree, parseNewick, type TreeNode } from "../phylogeny.ts";
import type { PhyloUcaPlacement } from "./types.ts";

interface PlacementPoint extends PhyloUcaPlacement {
  index: number;
  deltaLogLikelihood: number;
  relativeLikelihood: number;
  color: string;
  x: number;
  y: number;
  pendantPixels: number;
}

function interpolateColor(probability: number): string {
  const amount = Math.max(0, Math.min(1, probability));
  const blue = [33, 102, 172];
  const red = [178, 24, 43];
  const rgb = blue.map((value, index) => Math.round(value + amount * (red[index] - value)));
  return `rgb(${rgb.join(",")})`;
}

/** Match the deterministic pre-order internal labels used by PhyloUcaTreeMessages. */
function labelInternalNodes(root: TreeNode): TreeNode {
  let nextId = 0;
  const visit = (node: TreeNode): TreeNode => {
    const id = nextId++;
    return {
      ...node,
      name: node.children.length ? `internal_${id + 1}` : node.name,
      children: node.children.map(visit),
    };
  };
  return visit(root);
}

export function placementRelativeLikelihoods(placements: readonly PhyloUcaPlacement[]): Array<Pick<PlacementPoint, "deltaLogLikelihood" | "relativeLikelihood" | "color">> {
  const best = placements.reduce((maximum, placement) => Math.max(maximum, placement.logMarginalLikelihood), Number.NEGATIVE_INFINITY);
  return placements.map((placement) => {
    const deltaLogLikelihood = placement.logMarginalLikelihood - best;
    const relativeLikelihood = Math.exp(Math.min(0, deltaLogLikelihood));
    return { deltaLogLikelihood, relativeLikelihood, color: interpolateColor(relativeLikelihood) };
  });
}

function number(value: number, digits = 4): string {
  return Number.isFinite(value) ? value.toFixed(digits) : "—";
}

export function PhyloUcaPlacementMap({ newick, placements, inferenceMode, title }: {
  newick: string;
  placements: readonly PhyloUcaPlacement[];
  inferenceMode: "maximum-likelihood" | "grid-marginalization" | "gibbs-mh";
  title: string;
}) {
  const model = useMemo(() => {
    const retained = placements.filter((placement) => placement.localPosteriorWeight > 0);
    if (!retained.length) return null;
    const root = labelInternalNodes(parseNewick(newick));
    const layout = layoutTree(root, 920, 27, "phylogram", 34, 68);
    const edges = new Map(layout.edges.map((edge) => [`${edge.parent.name}\u0000${edge.child.name}`, edge]));
    const relative = placementRelativeLikelihoods(retained);
    const maximumPendant = Math.max(1e-12, ...retained.map((placement) => placement.ucaBranchLength));
    const points = retained.flatMap((placement, index): PlacementPoint[] => {
      const edge = edges.get(`${placement.endpointA}\u0000${placement.endpointB}`)
        ?? edges.get(`${placement.endpointB}\u0000${placement.endpointA}`);
      if (!edge) return [];
      const forward = edge.parent.name === placement.endpointA;
      const fraction = forward ? placement.edgeFraction : 1 - placement.edgeFraction;
      return [{
        ...placement,
        index: index + 1,
        ...relative[index],
        x: edge.parent.x + Math.max(0, Math.min(1, fraction)) * (edge.child.x - edge.parent.x),
        y: edge.child.y,
        pendantPixels: 10 + 38 * placement.ucaBranchLength / maximumPendant,
      }];
    });
    return { layout, points };
  }, [newick, placements]);

  if (!model) return null;
  const height = model.layout.height + 38;
  const heading = inferenceMode === "maximum-likelihood" ? "Conditional-ML tree attachment"
    : inferenceMode === "grid-marginalization" ? "Full-HMM grid marginalization"
      : "Continuous Gibbs/MH placement samples";
  const description = inferenceMode === "maximum-likelihood"
    ? "The marker is the single maximum-likelihood attachment and pendant length after continuous full-HMM optimization."
    : inferenceMode === "grid-marginalization"
      ? "Markers are the explicit attachment/pendant-length grid points retained in quadrature marginalization."
      : "Markers are retained posterior draws. Within-edge attachment fraction and pendant length are continuous; repeated points can occur when an MH proposal is rejected.";
  return <section className="phylo-uca-placement-map">
    <header><div><span className="section-kicker">{heading}</span><h5>Tree attachment points used in the UCA result</h5><p>{description} Markers sit along branches, not merely at nodes. Pendant length is proportional to the inferred branch from the attachment to the UCA. Color uses exp(ΔLL) from the best raw full-HMM marginal likelihood: red = 1; near-zero relative likelihood = blue.</p></div><div className="phylo-placement-color-key"><span>≈0</span><i /><span>1</span></div></header>
    <div className="phylo-uca-placement-scroll"><svg xmlns="http://www.w3.org/2000/svg" width={model.layout.width} height={height} viewBox={`0 0 ${model.layout.width} ${height}`} role="img" aria-label={title}>
      <title>{title}</title>
      <desc>Observed phylogeny with every full-HMM attachment and UCA branch-length point used by the selected inference route.</desc>
      <rect width={model.layout.width} height={height} fill="#fffdf9" />
      {model.layout.edges.map((edge, index) => <path key={index} d={`M ${edge.parent.x} ${edge.parent.y} V ${edge.child.y} H ${edge.child.x}`} fill="none" stroke="#788681" strokeWidth="1.1" />)}
      {model.layout.nodes.filter((node) => !node.children.length).map((node) => <text key={node.name} x={node.x + 5} y={node.y + 3} fill="#52605b" fontFamily="Inter,Arial,sans-serif" fontSize="8">{node.name.length > 24 ? `${node.name.slice(0, 22)}…` : node.name}</text>)}
      {model.points.map((point) => {
        const direction = point.index % 2 ? -1 : 1;
        const endpointY = point.y + direction * point.pendantPixels;
        const tooltip = `Point ${point.index}; ${point.endpointA} ↔ ${point.endpointB}; edge fraction ${(point.edgeFraction * 100).toFixed(3)}%; UCA branch ${number(point.ucaBranchLength, 6)}; log marginal ${number(point.logMarginalLikelihood)}; ΔLL ${number(point.deltaLogLikelihood)}; exp(ΔLL) ${point.relativeLikelihood.toExponential(4)}; local posterior weight ${(point.localPosteriorWeight * 100).toFixed(3)}%`;
        return <g key={`${point.index}-${point.edgeId}-${point.distanceFromA}-${point.ucaBranchLength}`} aria-label={tooltip}><title>{tooltip}</title>
          <circle cx={point.x} cy={point.y} r="2.4" fill="#182722" />
          <line x1={point.x} x2={point.x} y1={point.y} y2={endpointY} stroke={point.color} strokeWidth="1.4" />
          <circle cx={point.x} cy={endpointY} r="5" fill={point.color} stroke="#ffffff" strokeWidth="1.2" />
          <text x={point.x} y={endpointY + 2.4} textAnchor="middle" fill="#ffffff" fontFamily="Inter,Arial,sans-serif" fontSize="6" fontWeight="800">{point.index}</text>
        </g>;
      })}
    </svg></div>
    <div className="phylo-uca-placement-list" role="list" aria-label="Exact local placement likelihoods">{model.points.map((point) => <div role="listitem" key={`row-${point.index}`}><i style={{ background: point.color }} /><b>{point.index}</b><span>{point.endpointA} ↔ {point.endpointB} @ {(point.edgeFraction * 100).toFixed(2)}%</span><span>UCA {number(point.ucaBranchLength, 6)}</span><span>ΔLL {number(point.deltaLogLikelihood)}</span><span>e<sup>ΔLL</sup> {point.relativeLikelihood.toExponential(3)}</span><span>weight {(point.localPosteriorWeight * 100).toFixed(2)}%</span></div>)}</div>
  </section>;
}
