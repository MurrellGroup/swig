#!/usr/bin/env node
import { _ as annotateDoubleDBatch, a as deduplicate, c as parseFasta$2, d as runDenoisePartitionJob, f as runExactDedupJob, g as annotateAirrBatch, h as DEFAULT_PIPELINE_PLAN, i as createExactDedupPlan, l as prepareReferenceMsa, n as assignLineages, o as finishExactDedupPlan, p as sequenceFingerprint, s as normalizeNt, t as DenoiseAccumulator, v as datasetScopeValue, y as stableDatasetSeed } from "./chunks/post-analysis-core-tTuqWlbd.mjs";
import { createReadStream, createWriteStream, readFileSync, writeFileSync } from "node:fs";
import { mkdir, readFile, writeFile } from "node:fs/promises";
import { createHash } from "node:crypto";
import { basename, dirname, extname, isAbsolute, join, resolve } from "node:path";
import { createInterface } from "node:readline";
import { fileURLToPath } from "node:url";
import { Worker } from "node:worker_threads";
import { createGunzip, gunzipSync, gzipSync } from "node:zlib";
import { once } from "node:events";
import { availableParallelism } from "node:os";
import { Readable } from "node:stream";
//#region src/allele-refinement/evidence.ts
function callTokens(value) {
	return [...new Set(value.split(",").map((call) => call.trim()).filter(Boolean))];
}
function parseAlternativeEvidence(value) {
	return value.split(";").map((entry) => {
		const fields = entry.split("|");
		const score = Number(fields[1]);
		return {
			call: fields[0]?.trim() ?? "",
			score: Number.isFinite(score) ? score : null
		};
	}).filter((entry) => entry.call);
}
function modelScopeValue(row, options) {
	return datasetScopeValue({
		datasetId: row.datasetId,
		sampleId: row.sampleId,
		subjectId: row.subjectId
	}, options.scope);
}
/**
* Per-SNP evidence odds assigned to an unreported reference neighbour.
*
* If the true and reported alleles differ at one nucleotide, a read-level
* substitution rate mu gives an equal-substitution likelihood ratio of
* mu / [3(1-mu)] for changing the true diagnostic base into the reported
* base rather than retaining it. The baseline term deliberately represents
* unmodelled assignment uncertainty, not a sequencing-error probability.
*/
function adaptiveNeighbourOdds(readShm, options) {
	const baseline = Math.max(0, Math.min(.999999, options.baselineNeighbourOdds));
	const maximumShm = Math.max(0, Math.min(.95, options.maximumShm));
	const mu = Math.max(0, Math.min(maximumShm, readShm ?? 0));
	const sensitivity = Math.max(0, options.shmLeakageSensitivity);
	const somaticOdds = mu > 0 ? sensitivity * mu / (3 * (1 - mu)) : 0;
	const cap = Math.max(0, Math.min(.999999, options.maximumNeighbourOdds));
	return Math.min(cap, baseline + somaticOdds);
}
/**
* Converts a literal SwiftIG call into a sparse, database-aware evidence row.
* Co-optimal calls start at exactly equal weight. Unreported nucleotide
* neighbours receive geometric leakage; this is an explicit evidence kernel,
* not a claim that the affine alignment score is a calibrated read likelihood.
*/
function buildSparseEvidenceRow(row, graph, options) {
	const direct = /* @__PURE__ */ new Map();
	const selected = callTokens(row.call);
	for (const call of selected) {
		const node = graph.callToNode.get(call);
		if (node !== void 0) direct.set(node, 1);
	}
	const temperature = Math.max(1e-6, options.alternativeScoreTemperature);
	for (const alternative of parseAlternativeEvidence(row.alternatives)) {
		const node = graph.callToNode.get(alternative.call);
		if (node === void 0 || direct.has(node)) continue;
		const relative = alternative.score !== null && row.score !== null ? Math.min(1, Math.exp((alternative.score - row.score) / temperature)) : options.unscoredAlternativeWeight;
		direct.set(node, Math.max(direct.get(node) ?? 0, relative));
	}
	if (!direct.size) return null;
	const expanded = new Map(direct);
	const substitutionOdds = adaptiveNeighbourOdds(row.shm, options);
	const baselineOdds = Math.max(0, Math.min(options.maximumNeighbourOdds, options.baselineNeighbourOdds));
	if ((substitutionOdds > 0 || baselineOdds > 0) && options.neighbourRadius > 0) for (const [source, sourceWeight] of direct) for (const neighbour of graph.neighbours[source]) {
		if (neighbour.distance > options.neighbourRadius) break;
		const leaked = sourceWeight * (neighbour.substitutionOnly ? substitutionOdds : baselineOdds) ** neighbour.distance;
		if (leaked > (expanded.get(neighbour.index) ?? 0)) expanded.set(neighbour.index, leaked);
	}
	let entries = [...expanded.entries()].map(([node, weight]) => ({
		node,
		weight
	})).filter((entry) => entry.weight > 0).sort((left, right) => right.weight - left.weight || left.node - right.node);
	const retainedCap = Math.max(1, Math.floor(options.maxCandidatesPerRow), direct.size);
	const truncated = entries.length > retainedCap;
	if (truncated) entries = entries.slice(0, retainedCap);
	const sum = entries.reduce((total, entry) => total + entry.weight, 0);
	if (!(sum > 0)) return null;
	entries.forEach((entry) => {
		entry.weight /= sum;
	});
	const top = entries[0];
	const groupKey = `${modelScopeValue(row, options)}\u0000${row.locus}\u0000${graph.segment}`;
	return {
		ordinal: row.ordinal,
		groupKey,
		entries,
		weight: options.weighting === "abundance" ? Math.max(1, row.abundance) : 1,
		localTop: top.node,
		localTopProbability: top.weight,
		truncated
	};
}
/** Fixed-row/dynamically-grown-NNZ builder used by the browser worker. */
var SparseEvidenceAccumulator = class {
	graph;
	options;
	rowOffsets;
	ordinals;
	weights;
	rowGroups;
	localTop;
	localTopProbability;
	columns = new Uint32Array(4096);
	logEvidence = new Float32Array(4096);
	rowCount = 0;
	nonZeros = 0;
	groupIndex = /* @__PURE__ */ new Map();
	groupKeys = [];
	skippedRows = 0;
	truncatedRows = 0;
	constructor(maximumRows, graph, options) {
		this.graph = graph;
		this.options = options;
		this.rowOffsets = new Uint32Array(maximumRows + 1);
		this.ordinals = new Uint32Array(maximumRows);
		this.weights = new Float32Array(maximumRows);
		this.rowGroups = new Uint32Array(maximumRows);
		this.localTop = new Int32Array(maximumRows);
		this.localTop.fill(-1);
		this.localTopProbability = new Float32Array(maximumRows);
	}
	reserve(required) {
		if (required <= this.columns.length) return;
		let capacity = this.columns.length;
		while (capacity < required) capacity = Math.max(capacity + 4096, Math.ceil(capacity * 1.6));
		const columns = new Uint32Array(capacity);
		const evidence = new Float32Array(capacity);
		columns.set(this.columns);
		evidence.set(this.logEvidence);
		this.columns = columns;
		this.logEvidence = evidence;
	}
	add(input) {
		const row = buildSparseEvidenceRow(input, this.graph, this.options);
		if (!row) {
			this.skippedRows += 1;
			return false;
		}
		if (this.rowCount >= this.ordinals.length) throw new Error("Allele-refinement input exceeded the declared AIRR record count.");
		this.reserve(this.nonZeros + row.entries.length);
		let group = this.groupIndex.get(row.groupKey);
		if (group === void 0) {
			group = this.groupKeys.length;
			this.groupKeys.push(row.groupKey);
			this.groupIndex.set(row.groupKey, group);
		}
		this.rowOffsets[this.rowCount] = this.nonZeros;
		this.ordinals[this.rowCount] = row.ordinal;
		this.weights[this.rowCount] = row.weight;
		this.rowGroups[this.rowCount] = group;
		this.localTop[this.rowCount] = row.localTop;
		this.localTopProbability[this.rowCount] = row.localTopProbability;
		for (const entry of row.entries) {
			this.columns[this.nonZeros] = entry.node;
			this.logEvidence[this.nonZeros] = Math.log(Math.max(Number.MIN_VALUE, entry.weight));
			this.nonZeros += 1;
		}
		this.rowCount += 1;
		this.rowOffsets[this.rowCount] = this.nonZeros;
		if (row.truncated) this.truncatedRows += 1;
		return true;
	}
	finish() {
		return {
			rowOffsets: this.rowOffsets.subarray(0, this.rowCount + 1),
			columns: this.columns.subarray(0, this.nonZeros),
			logEvidence: this.logEvidence.subarray(0, this.nonZeros),
			ordinals: this.ordinals.subarray(0, this.rowCount),
			weights: this.weights.subarray(0, this.rowCount),
			groupKeys: [...this.groupKeys],
			rowGroups: this.rowGroups.subarray(0, this.rowCount),
			localTop: this.localTop.subarray(0, this.rowCount),
			localTopProbability: this.localTopProbability.subarray(0, this.rowCount),
			skippedRows: this.skippedRows,
			truncatedRows: this.truncatedRows
		};
	}
};
//#endregion
//#region src/allele-refinement/apply.ts
function posteriorMapPassesPolicy(policy, probability, minimumPosterior) {
	return policy === "best" || probability >= Math.max(0, Math.min(1, minimumPosterior));
}
//#endregion
//#region src/allele-refinement/reference-graph.ts
function normalizeReferenceSequence(value) {
	return value.toUpperCase().replaceAll("U", "T").replace(/[^ACGTN]/g, "");
}
function referenceLocus(name) {
	return name.toUpperCase().match(/^(IGH|IGK|IGL|TRA|TRB|TRD|TRG)/)?.[1] ?? "";
}
/** Banded Levenshtein distance with an early cutoff; only tiny radii are used here. */
function boundedReferenceDistance(left, right, maximum) {
	if (Math.abs(left.length - right.length) > maximum) return null;
	if (left === right) return 0;
	if (left.length === right.length) {
		let distance = 0;
		for (let index = 0; index < left.length; index += 1) {
			if (left[index] === right[index]) continue;
			distance += 1;
			if (distance > maximum) return null;
		}
		return distance;
	}
	const columns = right.length + 1;
	let previous = new Uint16Array(columns);
	let current = new Uint16Array(columns);
	for (let column = 0; column < columns; column += 1) previous[column] = column;
	for (let row = 1; row <= left.length; row += 1) {
		current.fill(maximum + 1);
		current[0] = row;
		const begin = Math.max(1, row - maximum);
		const end = Math.min(right.length, row + maximum);
		let rowMinimum = maximum + 1;
		for (let column = begin; column <= end; column += 1) {
			const substitution = previous[column - 1] + (left[row - 1] === right[column - 1] ? 0 : 1);
			const deletion = previous[column] + 1;
			const insertion = current[column - 1] + 1;
			const value = Math.min(substitution, deletion, insertion);
			current[column] = value;
			rowMinimum = Math.min(rowMinimum, value);
		}
		if (rowMinimum > maximum) return null;
		[previous, current] = [current, previous];
	}
	const distance = previous[right.length];
	return distance <= maximum ? distance : null;
}
function buildReferenceAlleleGraph(fasta, segment, maximumDistance) {
	const bySequence = /* @__PURE__ */ new Map();
	const callToNode = /* @__PURE__ */ new Map();
	let exactDuplicateLabels = 0;
	for (const record of parseFasta$2(fasta, true)) {
		const name = record.name.split(/\s+/, 1)[0]?.trim();
		const sequence = normalizeReferenceSequence(record.sequence);
		if (!name || !sequence) continue;
		const locus = referenceLocus(name);
		const key = `${locus}\u0000${sequence}`;
		let node = bySequence.get(key);
		if (!node) {
			node = {
				index: bySequence.size,
				segment,
				locus,
				names: [],
				sequence
			};
			bySequence.set(key, node);
		} else exactDuplicateLabels += 1;
		if (!node.names.includes(name)) node.names.push(name);
	}
	const nodes = [...bySequence.values()];
	nodes.forEach((node, index) => {
		node.index = index;
		node.names.sort();
		node.names.forEach((name) => callToNode.set(name, index));
	});
	const neighbours = nodes.map(() => []);
	const radius = Math.max(0, Math.min(5, Math.floor(maximumDistance)));
	if (radius > 0) {
		const byLocusLength = /* @__PURE__ */ new Map();
		nodes.forEach((node, index) => {
			const key = `${node.locus}\u0000${node.sequence.length}`;
			const values = byLocusLength.get(key);
			if (values) values.push(index);
			else byLocusLength.set(key, [index]);
		});
		const candidates = /* @__PURE__ */ new Set();
		const addCandidate = (left, right) => {
			const low = Math.min(left, right);
			const high = Math.max(left, right);
			candidates.add(low * nodes.length + high);
		};
		for (const indices of byLocusLength.values()) {
			if (indices.length < 2) continue;
			const length = nodes[indices[0]].sequence.length;
			const blocks = Math.min(radius + 1, Math.max(1, length));
			for (let block = 0; block < blocks; block += 1) {
				const begin = Math.floor(block * length / blocks);
				const end = Math.floor((block + 1) * length / blocks);
				const buckets = /* @__PURE__ */ new Map();
				for (const index of indices) {
					const key = nodes[index].sequence.slice(begin, end);
					const previous = buckets.get(key);
					if (previous) {
						for (const other of previous) addCandidate(other, index);
						previous.push(index);
					} else buckets.set(key, [index]);
				}
			}
		}
		const groups = [...byLocusLength.entries()].map(([key, indices]) => {
			const split = key.lastIndexOf("\0");
			return {
				locus: key.slice(0, split),
				length: Number(key.slice(split + 1)),
				indices
			};
		});
		for (let leftGroup = 0; leftGroup < groups.length; leftGroup += 1) for (let rightGroup = leftGroup + 1; rightGroup < groups.length; rightGroup += 1) {
			const left = groups[leftGroup];
			const right = groups[rightGroup];
			if (left.locus !== right.locus || left.length === right.length || Math.abs(left.length - right.length) > radius) continue;
			for (const leftIndex of left.indices) for (const rightIndex of right.indices) addCandidate(leftIndex, rightIndex);
		}
		for (const encoded of candidates) {
			const left = Math.floor(encoded / nodes.length);
			const right = encoded % nodes.length;
			const distance = boundedReferenceDistance(nodes[left].sequence, nodes[right].sequence, radius);
			if (distance === null || distance === 0) continue;
			const substitutionOnly = nodes[left].sequence.length === nodes[right].sequence.length;
			neighbours[left].push({
				index: right,
				distance,
				substitutionOnly
			});
			neighbours[right].push({
				index: left,
				distance,
				substitutionOnly
			});
		}
	}
	neighbours.forEach((values) => values.sort((left, right) => left.distance - right.distance || left.index - right.index));
	return {
		segment,
		nodes,
		callToNode,
		neighbours,
		exactDuplicateLabels
	};
}
//#endregion
//#region src/allele-refinement/input.ts
function numeric$1(value) {
	const parsed = Number(value);
	return value !== "" && Number.isFinite(parsed) ? parsed : null;
}
function normalizedIdentity(value) {
	if (value === null) return null;
	const normalized = value > 1 ? value / 100 : value;
	return normalized >= 0 && normalized <= 1 ? normalized : null;
}
function bestAlternativeIdentity(value) {
	let best = null;
	for (const entry of value.split(";")) {
		const identity = normalizedIdentity(numeric$1(entry.split("|")[2] ?? ""));
		if (identity !== null && (best === null || identity > best)) best = identity;
	}
	return best;
}
function bestIdentity(primary, alternatives) {
	const direct = normalizedIdentity(numeric$1(primary));
	const alternate = bestAlternativeIdentity(alternatives);
	if (direct === null) return alternate;
	if (alternate === null) return direct;
	return Math.max(direct, alternate);
}
function toRefinementInputRow(row, segment) {
	const prefix = segment.toLowerCase();
	const vIdentity = bestIdentity(row.values.v_identity ?? "", row.values.v_alternatives ?? "");
	const segmentIdentity = bestIdentity(row.values[`${prefix}_identity`] ?? "", row.values[`${prefix}_alternatives`] ?? "");
	const closestIdentity = vIdentity ?? segmentIdentity;
	return {
		ordinal: row.ordinal,
		sequenceId: row.values.sequence_id ?? "",
		datasetId: row.values.swig_dataset_id ?? "",
		sampleId: row.values.sample_id ?? "",
		subjectId: row.values.subject_id ?? "",
		locus: row.values.locus ?? "",
		call: row.values[`${prefix}_call`] ?? "",
		score: numeric$1(row.values[`${prefix}_score`] ?? ""),
		identity: numeric$1(row.values[`${prefix}_identity`] ?? ""),
		shm: closestIdentity === null ? null : Math.max(0, 1 - closestIdentity),
		alternatives: row.values[`${prefix}_alternatives`] ?? "",
		abundance: Math.max(1, Math.floor(Number(row.values.duplicate_count) || 1))
	};
}
//#endregion
//#region src/allele-refinement/active-set-model.ts
function parseGroupKey$1(key) {
	const [scopeValue = "", locus = "", segment = "V"] = key.split("\0");
	return {
		scopeValue,
		locus,
		segment
	};
}
function logistic(value) {
	if (value >= 0) return 1 / (1 + Math.exp(-value));
	const exponential = Math.exp(value);
	return exponential / (1 + exponential);
}
function logSumExp(values) {
	const maximum = Math.max(...values);
	if (!Number.isFinite(maximum)) return maximum;
	let sum = 0;
	for (const value of values) sum += Math.exp(value - maximum);
	return maximum + Math.log(sum);
}
/**
* Small deterministic quadrature for the positive-frequency slab. A gamma
* shape below one gives the empirically required long left tail. The prior
* mean is tied to the expected active fraction of the complete locus-matched
* database, while the integration floor is numerical rather than a count or
* retention threshold.
*/
function positiveFrequencyQuadrature(databaseNodes, options) {
	const points = Math.max(4, Math.min(64, Math.floor(options.activeSetQuadraturePoints)));
	const floor = Math.max(1e-12, Math.min(.05, options.activeSetFrequencyFloor));
	const ceiling = .95;
	const shape = Math.max(.05, Math.min(10, options.activeSetTailShape));
	const activeFraction = Math.max(1e-6, Math.min(.999999, options.activeSetPriorActiveFraction));
	const expectedActive = Math.max(2, databaseNodes * activeFraction);
	const scale = Math.min(.5, 1 / expectedActive) / shape;
	const logFloor = Math.log(floor);
	const logCeiling = Math.log(ceiling);
	const raw = [];
	for (let index = 0; index < points; index += 1) {
		const lower = Math.exp(logFloor + (logCeiling - logFloor) * index / points);
		const upper = Math.exp(logFloor + (logCeiling - logFloor) * (index + 1) / points);
		const frequency = Math.sqrt(lower * upper);
		const logDensity = (shape - 1) * Math.log(frequency) - frequency / scale;
		raw.push({
			frequency,
			logWeight: logDensity + Math.log(Math.max(Number.MIN_VALUE, upper - lower))
		});
	}
	const normalization = logSumExp(raw.map((entry) => entry.logWeight));
	return raw.map((entry) => ({
		...entry,
		logWeight: entry.logWeight - normalization
	}));
}
function patternMixtures(patterns, active, theta) {
	const mixtures = new Float64Array(patterns.length);
	patterns.forEach((pattern, patternIndex) => {
		let mixture = 0;
		for (let offset = 0; offset < pattern.columns.length; offset += 1) {
			const column = pattern.columns[offset];
			if (active[column]) mixture += theta[column] * pattern.evidence[offset];
		}
		mixtures[patternIndex] = mixture;
	});
	return mixtures;
}
/** Ordinary sparse EM restricted to the exact current active set. */
function fitActiveFrequencies(patterns, active, initialTheta, options) {
	const theta = Float64Array.from(initialTheta);
	let thetaSum = 0;
	for (let index = 0; index < theta.length; index += 1) {
		if (!active[index]) theta[index] = 0;
		thetaSum += theta[index];
	}
	if (!(thetaSum > 0)) {
		for (let index = 0; index < theta.length; index += 1) if (active[index]) theta[index] = 1;
		thetaSum = theta.reduce((sum, value) => sum + value, 0);
	}
	for (let index = 0; index < theta.length; index += 1) theta[index] /= thetaSum;
	const maximumIterations = Math.max(1, Math.floor(options.maxIterations));
	const tolerance = Math.max(1e-12, options.convergenceTolerance);
	let counts = new Float64Array(theta.length);
	let converged = false;
	let finalMaximumChange = Number.POSITIVE_INFINITY;
	let iterations = 0;
	for (; iterations < maximumIterations; iterations += 1) {
		counts = new Float64Array(theta.length);
		for (const pattern of patterns) {
			let normalizer = 0;
			for (let offset = 0; offset < pattern.columns.length; offset += 1) {
				const column = pattern.columns[offset];
				if (active[column]) normalizer += theta[column] * pattern.evidence[offset];
			}
			if (!(normalizer > 0)) {
				let fallback = -1;
				let fallbackEvidence = -1;
				for (let offset = 0; offset < pattern.columns.length; offset += 1) {
					const column = pattern.columns[offset];
					if (active[column] && pattern.evidence[offset] > fallbackEvidence) {
						fallback = column;
						fallbackEvidence = pattern.evidence[offset];
					}
				}
				if (fallback >= 0) counts[fallback] += pattern.weight;
				continue;
			}
			for (let offset = 0; offset < pattern.columns.length; offset += 1) {
				const column = pattern.columns[offset];
				if (!active[column]) continue;
				counts[column] += pattern.weight * theta[column] * pattern.evidence[offset] / normalizer;
			}
		}
		let total = 0;
		for (let index = 0; index < counts.length; index += 1) if (active[index]) total += counts[index];
		if (!(total > 0)) break;
		finalMaximumChange = 0;
		for (let index = 0; index < theta.length; index += 1) {
			const next = active[index] ? counts[index] / total : 0;
			finalMaximumChange = Math.max(finalMaximumChange, Math.abs(next - theta[index]));
			theta[index] = next;
		}
		if (finalMaximumChange <= tolerance) {
			converged = true;
			iterations += 1;
			break;
		}
	}
	return {
		theta,
		counts,
		iterations,
		converged,
		finalMaximumChange
	};
}
/**
* Approximate posterior inclusion probability from an exact one-dimensional
* conditional likelihood and the heavy-tailed positive-frequency quadrature.
* Only patterns containing the candidate are visited; all other records enter
* through one aggregated log(1-f) term.
*/
function inclusionProbability(candidate, active, theta, mixtures, patterns, incident, totalWeight, quadrature, priorLogOdds) {
	const isActive = Boolean(active[candidate]);
	const candidateWeight = isActive ? theta[candidate] : 0;
	const remaining = 1 - candidateWeight;
	if (isActive && remaining <= 1e-12) return 1;
	let incidentWeight = 0;
	const ratios = [];
	for (const entry of incident) {
		const pattern = patterns[entry.pattern];
		incidentWeight += pattern.weight;
		const baseline = isActive ? (mixtures[entry.pattern] - candidateWeight * entry.evidence) / remaining : mixtures[entry.pattern];
		if (!(baseline > 1e-300)) return 1;
		ratios.push({
			weight: pattern.weight,
			ratio: entry.evidence / baseline
		});
	}
	const outsideWeight = Math.max(0, totalWeight - incidentWeight);
	return logistic(priorLogOdds + logSumExp(quadrature.map(({ frequency, logWeight }) => {
		let delta = outsideWeight * Math.log1p(-frequency);
		for (const entry of ratios) {
			const relative = 1 - frequency + frequency * entry.ratio;
			if (!(relative > 0)) return Number.NEGATIVE_INFINITY;
			delta += entry.weight * Math.log(relative);
		}
		return logWeight + delta;
	})));
}
function fitGroup$1(matrix, graph, options, key, rowIndices) {
	const parsed = parseGroupKey$1(key);
	const nodeSet = /* @__PURE__ */ new Set();
	let nonZeros = 0;
	let effectiveRows = 0;
	for (const row of rowIndices) {
		effectiveRows += matrix.weights[row];
		const begin = matrix.rowOffsets[row];
		const end = matrix.rowOffsets[row + 1];
		nonZeros += end - begin;
		for (let offset = begin; offset < end; offset += 1) nodeSet.add(matrix.columns[offset]);
	}
	const nodes = [...nodeSet].sort((left, right) => left - right);
	const localByNode = new Map(nodes.map((node, index) => [node, index]));
	const databaseNodes = graph.nodes.reduce((count, node) => count + Number(node.locus === parsed.locus || !node.locus && nodeSet.has(node.index)), 0);
	const inactivePriorNodes = Math.max(0, databaseNodes - nodes.length);
	const localEvidenceCounts = new Float64Array(nodes.length);
	const patternsByKey = /* @__PURE__ */ new Map();
	for (const row of rowIndices) {
		const columns = [];
		const evidence = [];
		let patternKey = "";
		for (let offset = matrix.rowOffsets[row]; offset < matrix.rowOffsets[row + 1]; offset += 1) {
			const local = localByNode.get(matrix.columns[offset]);
			if (local === void 0) continue;
			const value = Math.exp(matrix.logEvidence[offset]);
			columns.push(local);
			evidence.push(value);
			localEvidenceCounts[local] += value * matrix.weights[row];
			patternKey += `${local}:${matrix.logEvidence[offset].toFixed(6)},`;
		}
		const previous = patternsByKey.get(patternKey);
		if (previous) previous.weight += matrix.weights[row];
		else patternsByKey.set(patternKey, {
			columns,
			evidence,
			weight: matrix.weights[row]
		});
	}
	const patterns = [...patternsByKey.values()];
	const incident = Array.from({ length: nodes.length }, () => []);
	patterns.forEach((pattern, patternIndex) => pattern.columns.forEach((column, offset) => {
		incident[column].push({
			pattern: patternIndex,
			evidence: pattern.evidence[offset]
		});
	}));
	const active = new Uint8Array(nodes.length);
	active.fill(1);
	let theta = Float64Array.from(localEvidenceCounts);
	let totalInitial = theta.reduce((sum, value) => sum + value, 0);
	if (!(totalInitial > 0)) totalInitial = Math.max(1, theta.length);
	for (let index = 0; index < theta.length; index += 1) theta[index] = theta[index] > 0 ? theta[index] / totalInitial : 1 / Math.max(1, theta.length);
	const activeFraction = Math.max(1e-6, Math.min(.999999, options.activeSetPriorActiveFraction));
	const priorLogOdds = Math.log(activeFraction) - Math.log1p(-activeFraction);
	const threshold = Math.max(0, Math.min(1, options.activeSetInclusionThreshold));
	const quadrature = positiveFrequencyQuadrature(databaseNodes, options);
	const totalWeight = patterns.reduce((sum, pattern) => sum + pattern.weight, 0);
	const inclusion = new Float64Array(nodes.length);
	inclusion.fill(1);
	let totalIterations = 0;
	let selectionConverged = false;
	let lastFit = fitActiveFrequencies(patterns, active, theta, options);
	totalIterations += lastFit.iterations;
	theta = lastFit.theta;
	const maximumSweeps = Math.max(1, Math.min(50, Math.floor(options.activeSetMaxSweeps)));
	for (let sweep = 0; sweep < maximumSweeps; sweep += 1) {
		const mixtures = patternMixtures(patterns, active, theta);
		const candidates = [];
		for (let index = 0; index < nodes.length; index += 1) {
			if (!active[index]) continue;
			const probability = inclusionProbability(index, active, theta, mixtures, patterns, incident[index], totalWeight, quadrature, priorLogOdds);
			inclusion[index] = probability;
			if (probability < threshold) candidates.push({
				index,
				probability
			});
		}
		candidates.sort((left, right) => left.probability - right.probability || theta[left.index] - theta[right.index] || left.index - right.index);
		const activePerPattern = patterns.map((pattern) => pattern.columns.reduce((count, column) => count + Number(active[column]), 0));
		let removed = 0;
		for (const candidate of candidates) {
			if (!active[candidate.index]) continue;
			if (incident[candidate.index].some((entry) => activePerPattern[entry.pattern] <= 1)) continue;
			active[candidate.index] = 0;
			theta[candidate.index] = 0;
			for (const entry of incident[candidate.index]) activePerPattern[entry.pattern] -= 1;
			removed += 1;
		}
		if (!removed) {
			selectionConverged = true;
			break;
		}
		lastFit = fitActiveFrequencies(patterns, active, theta, options);
		totalIterations += lastFit.iterations;
		theta = lastFit.theta;
	}
	lastFit = fitActiveFrequencies(patterns, active, theta, options);
	totalIterations += lastFit.iterations;
	theta = lastFit.theta;
	const finalMixtures = patternMixtures(patterns, active, theta);
	for (let index = 0; index < nodes.length; index += 1) if (active[index]) inclusion[index] = inclusionProbability(index, active, theta, finalMixtures, patterns, incident[index], totalWeight, quadrature, priorLogOdds);
	const alleles = nodes.map((node, index) => {
		const mean = active[index] ? theta[index] : 0;
		return {
			nodeIndex: node,
			names: [...graph.nodes[node].names],
			sequenceLength: graph.nodes[node].sequence.length,
			posteriorMean: mean,
			expectedAssignments: active[index] ? lastFit.counts[index] : 0,
			localEvidenceAssignments: localEvidenceCounts[index],
			posteriorSd: Math.sqrt(Math.max(0, mean * (1 - mean) / Math.max(1, effectiveRows + 1))),
			inclusionProbability: inclusion[index],
			active: Boolean(active[index])
		};
	}).sort((left, right) => Number(right.active) - Number(left.active) || right.posteriorMean - left.posteriorMean || (right.inclusionProbability ?? 0) - (left.inclusionProbability ?? 0) || left.names.join(",").localeCompare(right.names.join(",")));
	return {
		rowIndices,
		localByNode,
		active,
		theta,
		summary: {
			key,
			scopeValue: parsed.scopeValue,
			locus: parsed.locus,
			segment: parsed.segment,
			rows: rowIndices.length,
			effectiveRows,
			nonZeros,
			databaseNodes,
			inactivePriorNodes,
			inferenceModel: "active-set",
			activeAlleles: active.reduce((sum, value) => sum + Number(value), 0),
			alleles,
			iterations: totalIterations,
			converged: lastFit.converged && selectionConverged,
			finalMaximumChange: lastFit.finalMaximumChange
		}
	};
}
/** Fast exact-zero hurdle alternative to the continuous Dirichlet model. */
function fitSparseActiveSetAlleleModel(matrix, graph, options, totalRecords, onProgress) {
	const mapNode = new Int32Array(totalRecords);
	const localTopNode = new Int32Array(totalRecords);
	mapNode.fill(-1);
	localTopNode.fill(-1);
	const mapProbability = new Float32Array(totalRecords);
	const localTopProbability = new Float32Array(totalRecords);
	const posteriorEntropy = new Float32Array(totalRecords);
	const modelIndex = new Int32Array(totalRecords);
	modelIndex.fill(-1);
	const assignmentWeight = new Float32Array(totalRecords);
	const rowsByGroup = matrix.groupKeys.map(() => []);
	for (let row = 0; row < matrix.ordinals.length; row += 1) rowsByGroup[matrix.rowGroups[row]].push(row);
	const models = [];
	let changedMapRows = 0;
	for (let group = 0; group < rowsByGroup.length; group += 1) {
		const fit = fitGroup$1(matrix, graph, options, matrix.groupKeys[group], rowsByGroup[group]);
		models.push(fit.summary);
		for (const row of fit.rowIndices) {
			const ordinal = matrix.ordinals[row];
			modelIndex[ordinal] = group;
			assignmentWeight[ordinal] = matrix.weights[row];
			let maximum = Number.NEGATIVE_INFINITY;
			for (let offset = matrix.rowOffsets[row]; offset < matrix.rowOffsets[row + 1]; offset += 1) {
				const local = fit.localByNode.get(matrix.columns[offset]);
				if (local === void 0 || !fit.active[local] || !(fit.theta[local] > 0)) continue;
				maximum = Math.max(maximum, matrix.logEvidence[offset] + Math.log(fit.theta[local]));
			}
			let normalizer = 0;
			if (Number.isFinite(maximum)) for (let offset = matrix.rowOffsets[row]; offset < matrix.rowOffsets[row + 1]; offset += 1) {
				const local = fit.localByNode.get(matrix.columns[offset]);
				if (local === void 0 || !fit.active[local] || !(fit.theta[local] > 0)) continue;
				normalizer += Math.exp(matrix.logEvidence[offset] + Math.log(fit.theta[local]) - maximum);
			}
			let topNode = matrix.localTop[row];
			let topProbability = 1;
			let entropy = 0;
			if (normalizer > 0) {
				topNode = -1;
				topProbability = -1;
				for (let offset = matrix.rowOffsets[row]; offset < matrix.rowOffsets[row + 1]; offset += 1) {
					const node = matrix.columns[offset];
					const local = fit.localByNode.get(node);
					if (local === void 0 || !fit.active[local] || !(fit.theta[local] > 0)) continue;
					const probability = Math.exp(matrix.logEvidence[offset] + Math.log(fit.theta[local]) - maximum) / normalizer;
					if (probability > topProbability || probability === topProbability && node < topNode) {
						topNode = node;
						topProbability = probability;
					}
					if (probability > 0) entropy -= probability * Math.log(probability);
				}
			}
			mapNode[ordinal] = topNode;
			mapProbability[ordinal] = Math.max(0, topProbability);
			posteriorEntropy[ordinal] = entropy;
			localTopNode[ordinal] = matrix.localTop[row];
			localTopProbability[ordinal] = matrix.localTopProbability[row];
			if (topNode !== matrix.localTop[row]) changedMapRows += 1;
		}
		onProgress?.(group + 1, rowsByGroup.length);
	}
	return {
		segment: graph.segment,
		nodes: graph.nodes,
		mapNode,
		mapProbability,
		posteriorEntropy,
		localTopNode,
		localTopProbability,
		modelIndex,
		assignmentWeight,
		models,
		modeledRows: matrix.ordinals.length,
		changedMapRows,
		skippedRows: matrix.skippedRows,
		matrixNonZeros: matrix.columns.length,
		truncatedRows: matrix.truncatedRows,
		exactDuplicateLabels: graph.exactDuplicateLabels
	};
}
//#endregion
//#region src/allele-refinement/model.ts
/** Accurate enough for variational mixture updates over positive Dirichlet parameters. */
function digamma(value) {
	let x = value;
	let result = 0;
	while (x < 8) {
		result -= 1 / x;
		x += 1;
	}
	const inverse = 1 / x;
	const square = inverse * inverse;
	return result + Math.log(x) - .5 * inverse - square * (1 / 12 - square * (1 / 120 - square / 252));
}
function parseGroupKey(key) {
	const [scopeValue = "", locus = "", segment = "V"] = key.split("\0");
	return {
		scopeValue,
		locus,
		segment
	};
}
function fitGroup(matrix, graph, options, key, rowIndices) {
	const parsed = parseGroupKey(key);
	const nodeSet = /* @__PURE__ */ new Set();
	let nonZeros = 0;
	let effectiveRows = 0;
	for (const row of rowIndices) {
		effectiveRows += matrix.weights[row];
		const begin = matrix.rowOffsets[row];
		const end = matrix.rowOffsets[row + 1];
		nonZeros += end - begin;
		for (let offset = begin; offset < end; offset += 1) nodeSet.add(matrix.columns[offset]);
	}
	const nodes = [...nodeSet].sort((left, right) => left - right);
	const databaseNodes = graph.nodes.reduce((count, node) => count + Number(node.locus === parsed.locus || !node.locus && nodeSet.has(node.index)), 0);
	const inactivePriorNodes = Math.max(0, databaseNodes - nodes.length);
	const local = new Map(nodes.map((node, index) => [node, index]));
	const alpha = Math.max(1e-9, options.alphaPerAllele);
	const gamma = new Float64Array(nodes.length);
	const localEvidenceCounts = new Float64Array(nodes.length);
	gamma.fill(alpha);
	for (const row of rowIndices) {
		const begin = matrix.rowOffsets[row];
		const end = matrix.rowOffsets[row + 1];
		const weight = matrix.weights[row];
		for (let offset = begin; offset < end; offset += 1) {
			const index = local.get(matrix.columns[offset]);
			if (index === void 0) continue;
			const value = Math.exp(matrix.logEvidence[offset]) * weight;
			gamma[index] += value;
			localEvidenceCounts[index] += value;
		}
	}
	const patterns = /* @__PURE__ */ new Map();
	for (const row of rowIndices) {
		const begin = matrix.rowOffsets[row];
		const end = matrix.rowOffsets[row + 1];
		const columns = [];
		const evidence = [];
		let key = "";
		for (let offset = begin; offset < end; offset += 1) {
			const column = matrix.columns[offset];
			const logValue = matrix.logEvidence[offset];
			columns.push(column);
			evidence.push(logValue);
			key += `${column}:${logValue.toFixed(6)},`;
		}
		const previous = patterns.get(key);
		if (previous) previous.weight += matrix.weights[row];
		else patterns.set(key, {
			columns,
			logEvidence: evidence,
			weight: matrix.weights[row]
		});
	}
	let converged = false;
	let finalMaximumChange = Number.POSITIVE_INFINITY;
	let iterations = 0;
	for (; iterations < options.maxIterations; iterations += 1) {
		let gammaSum = inactivePriorNodes * alpha;
		for (const value of gamma) gammaSum += value;
		const expectedLog = Float64Array.from(gamma, (value) => digamma(value) - digamma(gammaSum));
		const counts = new Float64Array(nodes.length);
		for (const pattern of patterns.values()) {
			let maximum = Number.NEGATIVE_INFINITY;
			for (let offset = 0; offset < pattern.columns.length; offset += 1) {
				const index = local.get(pattern.columns[offset]);
				maximum = Math.max(maximum, pattern.logEvidence[offset] + expectedLog[index]);
			}
			let normalizer = 0;
			for (let offset = 0; offset < pattern.columns.length; offset += 1) normalizer += Math.exp(pattern.logEvidence[offset] + expectedLog[local.get(pattern.columns[offset])] - maximum);
			for (let offset = 0; offset < pattern.columns.length; offset += 1) {
				const probability = Math.exp(pattern.logEvidence[offset] + expectedLog[local.get(pattern.columns[offset])] - maximum) / normalizer;
				counts[local.get(pattern.columns[offset])] += pattern.weight * probability;
			}
		}
		finalMaximumChange = 0;
		for (let index = 0; index < gamma.length; index += 1) {
			const next = alpha + counts[index];
			finalMaximumChange = Math.max(finalMaximumChange, Math.abs(next - gamma[index]) / Math.max(1, gamma[index]));
			gamma[index] = next;
		}
		if (finalMaximumChange <= options.convergenceTolerance) {
			converged = true;
			iterations += 1;
			break;
		}
	}
	let gammaSum = inactivePriorNodes * alpha;
	for (const value of gamma) gammaSum += value;
	const expectedLog = Float64Array.from(gamma, (value) => digamma(value) - digamma(gammaSum));
	const alleles = nodes.map((node, index) => {
		const mean = gamma[index] / gammaSum;
		const variance = gamma[index] * (gammaSum - gamma[index]) / (gammaSum * gammaSum * (gammaSum + 1));
		return {
			nodeIndex: node,
			names: [...graph.nodes[node].names],
			sequenceLength: graph.nodes[node].sequence.length,
			posteriorMean: mean,
			expectedAssignments: Math.max(0, gamma[index] - alpha),
			localEvidenceAssignments: localEvidenceCounts[index],
			posteriorSd: Math.sqrt(Math.max(0, variance))
		};
	}).sort((left, right) => right.posteriorMean - left.posteriorMean || left.names.join(",").localeCompare(right.names.join(",")));
	return {
		rowIndices,
		nodes,
		expectedLog,
		summary: {
			key,
			scopeValue: parsed.scopeValue,
			locus: parsed.locus,
			segment: parsed.segment,
			rows: rowIndices.length,
			effectiveRows,
			nonZeros,
			databaseNodes,
			inactivePriorNodes,
			inferenceModel: "dirichlet",
			alleles,
			iterations,
			converged,
			finalMaximumChange
		}
	};
}
function fitSparseAlleleModel(matrix, graph, options, totalRecords, onProgress) {
	if (options.model === "active-set") return fitSparseActiveSetAlleleModel(matrix, graph, options, totalRecords, onProgress);
	const mapNode = new Int32Array(totalRecords);
	const localTopNode = new Int32Array(totalRecords);
	mapNode.fill(-1);
	localTopNode.fill(-1);
	const mapProbability = new Float32Array(totalRecords);
	const localTopProbability = new Float32Array(totalRecords);
	const posteriorEntropy = new Float32Array(totalRecords);
	const modelIndex = new Int32Array(totalRecords);
	modelIndex.fill(-1);
	const assignmentWeight = new Float32Array(totalRecords);
	const rowsByGroup = matrix.groupKeys.map(() => []);
	for (let row = 0; row < matrix.ordinals.length; row += 1) rowsByGroup[matrix.rowGroups[row]].push(row);
	const models = [];
	let changedMapRows = 0;
	for (let group = 0; group < rowsByGroup.length; group += 1) {
		const fit = fitGroup(matrix, graph, options, matrix.groupKeys[group], rowsByGroup[group]);
		models.push(fit.summary);
		const local = new Map(fit.nodes.map((node, index) => [node, index]));
		fit.rowIndices.forEach((row) => {
			const ordinal = matrix.ordinals[row];
			modelIndex[ordinal] = group;
			assignmentWeight[ordinal] = matrix.weights[row];
			const begin = matrix.rowOffsets[row];
			let topNode = -1;
			let topProbability = -1;
			let entropy = 0;
			let maximum = Number.NEGATIVE_INFINITY;
			for (let offset = begin; offset < matrix.rowOffsets[row + 1]; offset += 1) maximum = Math.max(maximum, matrix.logEvidence[offset] + fit.expectedLog[local.get(matrix.columns[offset])]);
			let normalizer = 0;
			for (let offset = begin; offset < matrix.rowOffsets[row + 1]; offset += 1) normalizer += Math.exp(matrix.logEvidence[offset] + fit.expectedLog[local.get(matrix.columns[offset])] - maximum);
			for (let offset = begin; offset < matrix.rowOffsets[row + 1]; offset += 1) {
				const probability = Math.exp(matrix.logEvidence[offset] + fit.expectedLog[local.get(matrix.columns[offset])] - maximum) / normalizer;
				const node = matrix.columns[offset];
				if (probability > topProbability || probability === topProbability && node < topNode) {
					topNode = node;
					topProbability = probability;
				}
				if (probability > 0) entropy -= probability * Math.log(probability);
			}
			mapNode[ordinal] = topNode;
			mapProbability[ordinal] = Math.max(0, topProbability);
			posteriorEntropy[ordinal] = entropy;
			localTopNode[ordinal] = matrix.localTop[row];
			localTopProbability[ordinal] = matrix.localTopProbability[row];
			if (topNode !== matrix.localTop[row]) changedMapRows += 1;
		});
		onProgress?.(group + 1, rowsByGroup.length);
	}
	return {
		segment: graph.segment,
		nodes: graph.nodes,
		mapNode,
		mapProbability,
		posteriorEntropy,
		localTopNode,
		localTopProbability,
		modelIndex,
		assignmentWeight,
		models,
		modeledRows: matrix.ordinals.length,
		changedMapRows,
		skippedRows: matrix.skippedRows,
		matrixNonZeros: matrix.columns.length,
		truncatedRows: matrix.truncatedRows,
		exactDuplicateLabels: graph.exactDuplicateLabels
	};
}
//#endregion
//#region src/allele-refinement/export.ts
function assignment(result, ordinal) {
	if (!result) return null;
	const node = result.mapNode[ordinal] ?? -1;
	if (node < 0 || !result.nodes[node]) return null;
	const local = result.localTopNode[ordinal] ?? -1;
	return {
		call: result.nodes[node].names.join(","),
		probability: result.mapProbability[ordinal] ?? 0,
		entropy: result.posteriorEntropy[ordinal] ?? 0,
		localCall: local >= 0 && result.nodes[local] ? result.nodes[local].names.join(",") : "",
		localProbability: result.localTopProbability[ordinal] ?? 0
	};
}
function refinedCall(result, segment, ordinal, policy, minimumPosterior) {
	const value = assignment(result?.segments[segment], ordinal);
	return value && posteriorMapPassesPolicy(policy, value.probability, minimumPosterior) ? value.call : null;
}
//#endregion
//#region src/allele-refinement/types.ts
const DEFAULT_ALLELE_REFINEMENT_OPTIONS = {
	model: "dirichlet",
	scope: "subject",
	segments: ["V", "J"],
	alphaPerAllele: .1,
	activeSetPriorActiveFraction: .15,
	activeSetInclusionThreshold: .5,
	activeSetTailShape: .35,
	activeSetFrequencyFloor: 1e-6,
	activeSetQuadraturePoints: 12,
	activeSetMaxSweeps: 8,
	baselineNeighbourOdds: .01,
	shmLeakageSensitivity: 1,
	maximumNeighbourOdds: .25,
	maximumShm: .3,
	neighbourRadius: 2,
	alternativeScoreTemperature: 2,
	unscoredAlternativeWeight: .25,
	maxCandidatesPerRow: 32,
	weighting: "unique",
	maxIterations: 100,
	convergenceTolerance: 1e-6
};
//#endregion
//#region src/germline-evidence.ts
const DEFAULT_MISSING_ALLELE_OPTIONS = {
	unit: "lineage",
	minimumIndependentUnits: 6,
	minimumCoveredUnits: 20,
	minimumAlleleFraction: .2,
	maximumShmRate: .08,
	minimumAlignedBases: 180,
	maximumCandidateSnps: 6,
	maximumPValue: 1e-6,
	minimumDistinctJCalls: 3,
	minimumDistinctJunctionLengths: 3,
	minimumDistinctCdr3s: 6,
	minimumLinkedFraction: .9,
	minimumNearGermlineUnits: 3,
	maximumOtherAlternateFraction: .02
};
const FLAG_JOINT_COVERAGE = 1;
const FLAG_JOINT_ALTERNATE = 2;
const FLAG_REFERENCE_PRESENT = 4;
const FLAG_CONFLICT_PRESENT = 8;
const MAX_PROPOSALS_PER_V = 64;
const SUPPORTING_UNIT_EXPORT_LIMIT = 1e3;
const EMPTY_MUTATIONS = [];
function topCall$1(value) {
	return value.split(",")[0]?.trim() ?? "";
}
function cleanAlignment(value) {
	return value.toUpperCase().replace(/\./g, "-").replace(/[^ACGTN-]/g, "N");
}
function cleanCdr3(value) {
	return value.toUpperCase().replace(/[^ACGTN]/g, "");
}
function cleanReference(value) {
	return value.toUpperCase().replace(/[^ACGTN]/g, "");
}
function lineageUnit(lineageId) {
	return Number.isFinite(lineageId) && lineageId > 0 ? Math.floor(lineageId) : 0;
}
function displayUnit(item) {
	return `${item.subject} · lineage:${item.unit}`;
}
function cdr3Fingerprint(value) {
	let first = 2166136261, second = 2246822519;
	for (let index = 0; index < value.length; index += 1) {
		const code = value.charCodeAt(index);
		first = Math.imul(first ^ code, 16777619);
		second = Math.imul(second ^ code, 3266489917);
	}
	return [first >>> 0, second >>> 0];
}
function parseFasta$1(value) {
	const result = /* @__PURE__ */ new Map();
	let name = "", sequence = "";
	const commit = () => {
		if (name && sequence) result.set(name, cleanReference(sequence));
	};
	for (const line of value.split(/\r?\n/)) if (line.startsWith(">")) {
		commit();
		name = line.slice(1).trim().split(/\s+/)[0];
		sequence = "";
	} else sequence += line.trim();
	commit();
	return result;
}
function logGamma(value) {
	const coefficients = [
		76.18009172947146,
		-86.50532032941678,
		24.01409824083091,
		-1.231739572450155,
		.001208650973866179,
		-5395239384953e-18
	];
	let x = value, y = value, tmp = x + 5.5;
	tmp -= (x + .5) * Math.log(tmp);
	let series = 1.000000000190015;
	for (const coefficient of coefficients) {
		y += 1;
		series += coefficient / y;
	}
	return -tmp + Math.log(2.5066282746310007 * series / x);
}
function betaFraction(a, b, x) {
	const maximum = 200, epsilon = 3e-12, tiny = 1e-300;
	const qab = a + b, qap = a + 1, qam = a - 1;
	let c = 1, d = 1 - qab * x / qap;
	if (Math.abs(d) < tiny) d = tiny;
	d = 1 / d;
	let h = d;
	for (let m = 1; m <= maximum; m += 1) {
		const m2 = 2 * m;
		let aa = m * (b - m) * x / ((qam + m2) * (a + m2));
		d = 1 + aa * d;
		if (Math.abs(d) < tiny) d = tiny;
		c = 1 + aa / c;
		if (Math.abs(c) < tiny) c = tiny;
		d = 1 / d;
		h *= d * c;
		aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
		d = 1 + aa * d;
		if (Math.abs(d) < tiny) d = tiny;
		c = 1 + aa / c;
		if (Math.abs(c) < tiny) c = tiny;
		d = 1 / d;
		const delta = d * c;
		h *= delta;
		if (Math.abs(delta - 1) < epsilon) break;
	}
	return h;
}
function regularizedBeta(x, a, b) {
	if (x <= 0) return 0;
	if (x >= 1) return 1;
	const bt = Math.exp(logGamma(a + b) - logGamma(a) - logGamma(b) + a * Math.log(x) + b * Math.log(1 - x));
	return x < (a + 1) / (a + b + 2) ? bt * betaFraction(a, b, x) / a : 1 - bt * betaFraction(b, a, 1 - x) / b;
}
function binomialSurvival(successes, trials, probability) {
	if (successes <= 0) return 1;
	if (successes > trials) return 0;
	return Math.min(1, Math.max(0, regularizedBeta(Math.min(.999999999999, Math.max(1e-12, probability)), successes, trials - successes + 1)));
}
/** AID WRCY/RGYW context at the substituted base. */
function hotspotAt(sequence, zeroBasedPosition) {
	const base = sequence[zeroBasedPosition];
	const isW = (value = "") => value === "A" || value === "T";
	const isR = (value = "") => value === "A" || value === "G";
	const isY = (value = "") => value === "C" || value === "T";
	if (base === "C") return isW(sequence[zeroBasedPosition - 2]) && isR(sequence[zeroBasedPosition - 1]) && isY(sequence[zeroBasedPosition + 1]);
	if (base === "G") return isR(sequence[zeroBasedPosition - 1]) && isY(sequence[zeroBasedPosition + 1]) && isW(sequence[zeroBasedPosition + 2]);
	return false;
}
function parseObservation(row, ordinal, unit, retainBases = false) {
	const vCall = topCall$1(row.v_call);
	const query = cleanAlignment(row.v_sequence_alignment || "");
	const germline = cleanAlignment(row.v_germline_alignment || "");
	if (!unit || !vCall || !query || !germline) return null;
	const length = Math.min(query.length, germline.length);
	let mutations = EMPTY_MUTATIONS;
	const covered = [];
	const baseByPosition = retainBases ? /* @__PURE__ */ new Map() : void 0;
	const reportedStart = Math.floor(Number(row.v_germline_start));
	let position = Number.isFinite(reportedStart) && reportedStart > 0 ? reportedStart - 1 : 0;
	let aligned = 0, intervalStart = 0, lastCovered = 0;
	for (let column = 0; column < length; column += 1) {
		const q = query[column], g = germline[column];
		if (g !== "-") position += 1;
		if (/^[ACGT]$/.test(q) && /^[ACGT]$/.test(g)) {
			aligned += 1;
			if (!intervalStart || position !== lastCovered + 1) {
				if (intervalStart) covered.push(intervalStart, lastCovered);
				intervalStart = position;
			}
			lastCovered = position;
			baseByPosition?.set(position, {
				reference: g,
				query: q
			});
			if (q !== g) {
				if (mutations === EMPTY_MUTATIONS) mutations = [];
				mutations.push({
					position,
					reference: g,
					alternate: q
				});
			}
		}
	}
	if (intervalStart) covered.push(intervalStart, lastCovered);
	const cdr3 = cleanCdr3(row.cdr3 || row.junction || ""), [cdr3HashA, cdr3HashB] = cdr3Fingerprint(cdr3);
	return {
		ordinal,
		unit,
		subject: (row.subject_id || "unassigned-subject").trim() || "unassigned-subject",
		vCall,
		jCall: topCall$1(row.j_call),
		cdr3HashA,
		cdr3HashB,
		junctionLength: cdr3.length,
		aligned,
		rate: aligned ? mutations.length / aligned : 1,
		mutations,
		coverageStart: covered[0] ?? 0,
		coverageEnd: covered[covered.length - 1] ?? 0,
		coverageBreaks: covered.length > 2 ? covered : void 0,
		baseByPosition
	};
}
function eventKey(event) {
	return `${event.position}:${event.reference}>${event.alternate}`;
}
function eventAt(item, event) {
	const mutation = item.mutations.find((value) => value.position === event.position);
	return mutation?.reference === event.reference && mutation.alternate === event.alternate;
}
function coversEvent(item, event) {
	const intervals = item.coverageBreaks;
	if (!intervals) return event.position >= item.coverageStart && event.position <= item.coverageEnd;
	for (let index = 0; index < intervals.length; index += 2) if (event.position >= intervals[index] && event.position <= intervals[index + 1]) return true;
	return false;
}
function setIntersectionSize(left, right) {
	let result = 0;
	for (const value of left) if (right.has(value)) result += 1;
	return result;
}
function isSubset(left, right) {
	const rightEvents = "substitutions" in right ? right.substitutions : right.events;
	const leftEvents = "substitutions" in left ? left.substitutions : left.events;
	const rightKeys = new Set(rightEvents.map(eventKey));
	return leftEvents.every((event) => rightKeys.has(eventKey(event)));
}
function mutationNullProbability(rate, hotspot) {
	return Math.min(.25, Math.max(1e-5, rate / 3) * (hotspot ? 5 : 1));
}
function diversity(items) {
	return {
		jCalls: new Set(items.map((item) => item.jCall).filter(Boolean)),
		junctionLengths: new Set(items.map((item) => item.junctionLength).filter((value) => value > 0)),
		cdr3s: new Set(items.filter((item) => item.junctionLength > 0).map((item) => `${item.cdr3HashA.toString(36)}:${item.cdr3HashB.toString(36)}`)),
		subjects: new Set(items.map((item) => item.subject).filter(Boolean))
	};
}
var MissingAlleleAccumulator = class {
	options;
	best = /* @__PURE__ */ new Map();
	inputRecords = 0;
	eligibleRecords = 0;
	constructor(options = {}) {
		this.options = {
			...DEFAULT_MISSING_ALLELE_OPTIONS,
			...options,
			unit: "lineage"
		};
	}
	add(row, ordinal, lineageId = 0) {
		this.inputRecords += 1;
		const unit = lineageUnit(lineageId);
		const item = parseObservation(row, ordinal, unit);
		if (!item) return;
		if (item.aligned < this.options.minimumAlignedBases || item.rate > this.options.maximumShmRate) return;
		this.eligibleRecords += 1;
		const existing = this.best.get(unit);
		if (!existing || item.rate < existing.rate || item.rate === existing.rate && item.aligned > existing.aligned || item.rate === existing.rate && item.aligned === existing.aligned && item.ordinal < existing.ordinal) this.best.set(unit, item);
	}
	prepareValidation(referenceFasta = "") {
		const references = parseFasta$1(referenceFasta);
		const allKnown = new Set([...references.values()].map(cleanReference));
		const groups = /* @__PURE__ */ new Map();
		for (const item of this.best.values()) {
			const group = groups.get(item.vCall) ?? [];
			group.push(item);
			groups.set(item.vCall, group);
		}
		const proposals = [];
		let proposalTruncations = 0;
		for (const [vCall, items] of groups) {
			if (items.length < this.options.minimumCoveredUnits) continue;
			const parent = references.get(vCall) || "";
			if (!parent) continue;
			const coverageDelta = new Int32Array(parent.length + 2);
			const changes = /* @__PURE__ */ new Map();
			let summedRates = 0;
			for (const item of items) {
				summedRates += item.rate;
				if (item.coverageBreaks) for (let interval = 0; interval < item.coverageBreaks.length; interval += 2) {
					const start = Math.max(1, item.coverageBreaks[interval]), end = Math.min(parent.length, item.coverageBreaks[interval + 1]);
					if (start <= end) {
						coverageDelta[start] += 1;
						coverageDelta[end + 1] -= 1;
					}
				}
				else {
					const start = Math.max(1, item.coverageStart), end = Math.min(parent.length, item.coverageEnd);
					if (start <= end) {
						coverageDelta[start] += 1;
						coverageDelta[end + 1] -= 1;
					}
				}
				for (const mutation of item.mutations) {
					if (mutation.position < 1 || mutation.position > parent.length || parent[mutation.position - 1] !== mutation.reference) continue;
					const key = eventKey(mutation);
					const found = changes.get(key) ?? {
						...mutation,
						supporting: []
					};
					found.supporting.push(item);
					changes.set(key, found);
				}
			}
			const coverageByPosition = new Uint32Array(parent.length + 1);
			let runningCoverage = 0;
			for (let position = 1; position <= parent.length; position += 1) {
				runningCoverage += coverageDelta[position];
				coverageByPosition[position] = runningCoverage;
			}
			const passing = [];
			for (const change of changes.values()) {
				const coverage = coverageByPosition[change.position];
				const support = change.supporting.length;
				const fraction = support / Math.max(1, coverage);
				const background = Math.max(0, (summedRates - change.supporting.reduce((sum, item) => sum + 1 / Math.max(1, item.aligned), 0)) / items.length);
				const hotspot = hotspotAt(parent, change.position - 1);
				const pValue = binomialSurvival(support, coverage, mutationNullProbability(background, hotspot));
				if (support >= this.options.minimumIndependentUnits && coverage >= this.options.minimumCoveredUnits && fraction >= this.options.minimumAlleleFraction && pValue <= this.options.maximumPValue) passing.push({
					position: change.position,
					reference: change.reference,
					alternate: change.alternate,
					support,
					coverage,
					fraction,
					pValue,
					hotspot,
					supportingUnits: new Set(change.supporting.map((item) => item.unit))
				});
			}
			passing.sort((left, right) => left.pValue - right.pValue || right.fraction - left.fraction || right.support - left.support || left.position - right.position);
			const signatures = /* @__PURE__ */ new Set();
			const vProposals = [];
			const itemByUnit = new Map(items.map((item) => [item.unit, item]));
			for (const seed of passing) {
				const events = [seed];
				let currentSupport = new Set(seed.supportingUnits);
				for (const candidate of passing) {
					if (events.length >= this.options.maximumCandidateSnps || events.some((event) => eventKey(event) === eventKey(candidate))) continue;
					let currentEligible = 0, candidateEligible = 0, both = 0;
					for (const unit of currentSupport) if (coversEvent(itemByUnit.get(unit), candidate)) {
						currentEligible += 1;
						if (candidate.supportingUnits.has(unit)) both += 1;
					}
					for (const unit of candidate.supportingUnits) {
						const item = itemByUnit.get(unit);
						if (events.every((event) => coversEvent(item, event))) candidateEligible += 1;
					}
					if (both < this.options.minimumIndependentUnits || both / Math.max(1, currentEligible) < this.options.minimumLinkedFraction || both / Math.max(1, candidateEligible) < this.options.minimumLinkedFraction) continue;
					events.push(candidate);
					currentSupport = new Set([...currentSupport].filter((unit) => candidate.supportingUnits.has(unit)));
				}
				events.sort((left, right) => left.position - right.position);
				const signature = events.map(eventKey).join("|");
				if (signatures.has(signature)) continue;
				signatures.add(signature);
				const covered = items.filter((item) => events.every((event) => coversEvent(item, event)));
				const supporting = [...currentSupport].map((unit) => itemByUnit.get(unit));
				const fraction = supporting.length / Math.max(1, covered.length);
				const groupsSeen = diversity(supporting);
				if (supporting.length < this.options.minimumIndependentUnits || covered.length < this.options.minimumCoveredUnits || fraction < this.options.minimumAlleleFraction || groupsSeen.jCalls.size < this.options.minimumDistinctJCalls || groupsSeen.junctionLengths.size < this.options.minimumDistinctJunctionLengths || groupsSeen.cdr3s.size < this.options.minimumDistinctCdr3s) continue;
				const meanBackgroundShm = items.reduce((sum, item) => sum + Math.max(0, item.mutations.length - events.filter((event) => eventAt(item, event)).length) / Math.max(1, item.aligned), 0) / items.length;
				const nullProbability = events.reduce((probability, event) => probability * mutationNullProbability(meanBackgroundShm, event.hotspot), 1);
				const pValue = binomialSurvival(supporting.length, covered.length, nullProbability);
				if (pValue > this.options.maximumPValue) continue;
				const sequence = [...parent];
				for (const event of events) if (sequence[event.position - 1] === event.reference) sequence[event.position - 1] = event.alternate;
				const candidateSequence = sequence.join("");
				if (allKnown.has(candidateSequence)) continue;
				const storedEvents = events.map((event) => ({
					position: event.position,
					reference: event.reference,
					alternate: event.alternate,
					support: event.support,
					coverage: event.coverage,
					fraction: event.fraction,
					pValue: event.pValue,
					hotspot: event.hotspot
				}));
				vProposals.push({
					vCall,
					sequence: candidateSequence,
					events: storedEvents,
					preliminarySupport: new Set(supporting.map((item) => item.unit)),
					meanBackgroundShm,
					nullProbability,
					pValue
				});
			}
			vProposals.sort((left, right) => left.pValue - right.pValue || right.preliminarySupport.size - left.preliminarySupport.size || right.events.length - left.events.length);
			const maximal = vProposals.filter((candidate, index) => !vProposals.some((other, otherIndex) => otherIndex !== index && other.events.length > candidate.events.length && isSubset(candidate, other) && setIntersectionSize(candidate.preliminarySupport, other.preliminarySupport) / Math.max(1, candidate.preliminarySupport.size) >= this.options.minimumLinkedFraction));
			if (maximal.length > MAX_PROPOSALS_PER_V) proposalTruncations += maximal.length - MAX_PROPOSALS_PER_V;
			for (const proposal of maximal.slice(0, MAX_PROPOSALS_PER_V)) {
				proposal.preliminarySupport.clear();
				proposals.push(proposal);
			}
		}
		return new MissingAlleleValidator(this.options, this.inputRecords, this.eligibleRecords, this.best, groups.size, proposals, proposalTruncations, allKnown);
	}
};
var MissingAlleleValidator = class {
	groups = /* @__PURE__ */ new Map();
	options;
	inputRecords;
	eligibleRecords;
	discoveryUnits;
	vAllelesTested;
	proposals;
	proposalTruncations;
	allKnown;
	constructor(options, inputRecords, eligibleRecords, discoveryUnits, vAllelesTested, proposals, proposalTruncations, allKnown) {
		this.options = options;
		this.inputRecords = inputRecords;
		this.eligibleRecords = eligibleRecords;
		this.discoveryUnits = discoveryUnits;
		this.vAllelesTested = vAllelesTested;
		this.proposals = proposals;
		this.proposalTruncations = proposalTruncations;
		this.allKnown = allKnown;
		for (const proposal of proposals) {
			const group = this.groups.get(proposal.vCall) ?? {
				proposals: [],
				observations: [],
				flags: new Uint8Array(0),
				mixedVCall: new Uint8Array(0)
			};
			group.proposals.push(proposal);
			this.groups.set(proposal.vCall, group);
		}
		for (const item of discoveryUnits.values()) {
			const group = this.groups.get(item.vCall);
			if (!group) continue;
			item.validationIndex = group.observations.length;
			group.observations.push(item);
		}
		for (const group of this.groups.values()) {
			group.flags = new Uint8Array(group.observations.length * group.proposals.length);
			group.mixedVCall = new Uint8Array(group.observations.length);
		}
	}
	add(row, ordinal, lineageId = 0) {
		const unit = lineageUnit(lineageId);
		const discovery = this.discoveryUnits.get(unit);
		if (!discovery) return;
		const group = this.groups.get(discovery.vCall);
		const unitIndex = discovery.validationIndex;
		if (!group || unitIndex === void 0) return;
		const observedVCall = topCall$1(row.v_call);
		if (observedVCall && observedVCall !== discovery.vCall) {
			group.mixedVCall[unitIndex] = 1;
			return;
		}
		const item = parseObservation(row, ordinal, unit, true);
		if (!item) return;
		group.proposals.forEach((proposal, proposalIndex) => {
			let allCovered = true, allAlternate = true, referencePresent = false, conflictPresent = false;
			for (const event of proposal.events) {
				const base = item.baseByPosition?.get(event.position);
				if (!base) {
					allCovered = false;
					allAlternate = false;
					continue;
				}
				if (base.reference !== event.reference) {
					allCovered = false;
					allAlternate = false;
					conflictPresent = true;
					continue;
				}
				if (base.query === event.reference) {
					referencePresent = true;
					allAlternate = false;
				} else if (base.query !== event.alternate) {
					conflictPresent = true;
					allAlternate = false;
				}
			}
			const offset = unitIndex * group.proposals.length + proposalIndex;
			if (allCovered) group.flags[offset] |= FLAG_JOINT_COVERAGE;
			if (allCovered && allAlternate) group.flags[offset] |= FLAG_JOINT_ALTERNATE;
			if (referencePresent) group.flags[offset] |= FLAG_REFERENCE_PRESENT;
			if (conflictPresent) group.flags[offset] |= FLAG_CONFLICT_PRESENT;
		});
	}
	finish() {
		const candidates = [];
		const candidateSupportSets = /* @__PURE__ */ new Map();
		let referenceVetoedLineagePatterns = 0, conflictingLineagePatterns = 0;
		for (const [vCall, group] of this.groups) group.proposals.forEach((proposal, proposalIndex) => {
			const supporting = [];
			let coveredUnits = 0, referencePresentUnits = 0, conflictingUnits = 0;
			for (let unitIndex = 0; unitIndex < group.observations.length; unitIndex += 1) {
				const item = group.observations[unitIndex], flags = group.flags[unitIndex * group.proposals.length + proposalIndex], mixedVCall = Boolean(group.mixedVCall[unitIndex]);
				if (flags & FLAG_JOINT_COVERAGE && !mixedVCall) coveredUnits += 1;
				const hadAlternate = Boolean(flags & FLAG_JOINT_ALTERNATE);
				const hadReference = Boolean(flags & FLAG_REFERENCE_PRESENT);
				const hadConflict = Boolean(flags & FLAG_CONFLICT_PRESENT) || mixedVCall;
				if (hadAlternate && hadReference) referencePresentUnits += 1;
				if (flags & FLAG_CONFLICT_PRESENT && !mixedVCall) conflictingUnits += 1;
				if (hadAlternate && !hadReference && !hadConflict) supporting.push(item);
			}
			referenceVetoedLineagePatterns += referencePresentUnits;
			conflictingLineagePatterns += conflictingUnits;
			const alleleFraction = supporting.length / Math.max(1, coveredUnits);
			const groupsSeen = diversity(supporting);
			const candidateEventKeys = new Set(proposal.events.map(eventKey));
			const nearGermlineUnits = supporting.filter((item) => item.mutations.filter((mutation) => !candidateEventKeys.has(eventKey(mutation))).length <= 2).length;
			const pValue = binomialSurvival(supporting.length, coveredUnits, proposal.nullProbability);
			const otherAlternateFraction = conflictingUnits / Math.max(1, coveredUnits);
			if (supporting.length < this.options.minimumIndependentUnits || coveredUnits < this.options.minimumCoveredUnits || alleleFraction < this.options.minimumAlleleFraction || otherAlternateFraction > this.options.maximumOtherAlternateFraction || pValue > this.options.maximumPValue || groupsSeen.jCalls.size < this.options.minimumDistinctJCalls || groupsSeen.junctionLengths.size < this.options.minimumDistinctJunctionLengths || groupsSeen.cdr3s.size < this.options.minimumDistinctCdr3s || nearGermlineUnits < this.options.minimumNearGermlineUnits || this.allKnown.has(proposal.sequence)) return;
			const units = supporting.map(displayUnit).sort();
			const candidate = {
				id: "",
				vCall,
				parentAllele: vCall,
				substitutions: proposal.events.map((event) => ({
					position: event.position,
					reference: event.reference,
					alternate: event.alternate,
					hotspot: event.hotspot,
					support: supporting.length,
					coverage: coveredUnits,
					fraction: alleleFraction,
					pValue
				})),
				independentUnits: supporting.length,
				coveredUnits,
				alleleFraction,
				pValue,
				distinctJCalls: groupsSeen.jCalls.size,
				distinctJunctionLengths: groupsSeen.junctionLengths.size,
				distinctCdr3s: groupsSeen.cdr3s.size,
				distinctSubjects: groupsSeen.subjects.size,
				nearGermlineUnits,
				referencePresentUnits,
				conflictingUnits,
				otherAlternateFraction,
				meanBackgroundShm: proposal.meanBackgroundShm,
				sequence: proposal.sequence,
				supportingUnits: units.slice(0, SUPPORTING_UNIT_EXPORT_LIMIT),
				supportingUnitsTruncated: units.length > SUPPORTING_UNIT_EXPORT_LIMIT,
				caution: "Diagnostic candidate only. Confirm with genotype-aware germline inference and independent data before adding it to a reference set."
			};
			candidates.push(candidate);
			candidateSupportSets.set(candidate, new Set(supporting.map((item) => item.unit)));
		});
		candidates.sort((left, right) => left.pValue - right.pValue || right.independentUnits - left.independentUnits || right.substitutions.length - left.substitutions.length);
		const maximal = candidates.filter((candidate, index) => !candidates.some((other, otherIndex) => otherIndex !== index && other.vCall === candidate.vCall && other.substitutions.length > candidate.substitutions.length && isSubset(candidate, other) && setIntersectionSize(candidateSupportSets.get(candidate), candidateSupportSets.get(other)) / Math.max(1, candidate.independentUnits) >= this.options.minimumLinkedFraction));
		maximal.forEach((candidate, index) => {
			candidate.id = `${candidate.vCall.replace(/[^A-Za-z0-9_.-]/g, "_")}_candidate_${index + 1}`;
		});
		const warnings = [];
		if (referenceVetoedLineagePatterns) warnings.push(`${referenceVetoedLineagePatterns.toLocaleString()} lineage–candidate pattern${referenceVetoedLineagePatterns === 1 ? " was" : "s were"} vetoed because another covered member of the same lineage contained a parent-reference nucleotide at a proposed site.`);
		if (conflictingLineagePatterns) warnings.push(`${conflictingLineagePatterns.toLocaleString()} lineage–candidate pattern${conflictingLineagePatterns === 1 ? " contained" : "s contained"} a third nucleotide state; candidates above the configured cross-lineage fraction were rejected.`);
		if (this.proposalTruncations) warnings.push(`${this.proposalTruncations.toLocaleString()} lower-ranked preliminary pattern${this.proposalTruncations === 1 ? " was" : "s were"} omitted by the 64-pattern-per-V browser memory guard.`);
		if (maximal.length) warnings.push(`${maximal.length} linked low-SHM haplotype${maximal.length === 1 ? "" : "s"} passed the two-pass diagnostic. These are referral candidates, not inferred genotypes.`);
		return {
			mode: "lineage",
			validationPasses: 2,
			inputRecords: this.inputRecords,
			eligibleRecords: this.eligibleRecords,
			independentUnits: this.discoveryUnits.size,
			vAllelesTested: this.vAllelesTested,
			candidatePatternsTested: this.proposals.length,
			referenceVetoedLineagePatterns,
			conflictingLineagePatterns,
			proposalTruncations: this.proposalTruncations,
			candidates: maximal,
			warnings
		};
	}
};
//#endregion
//#region src/repertoire-selection.ts
const DEFAULT_REPERTOIRE_SELECTION = {
	sequenceId: "",
	datasetId: "",
	sampleId: "",
	subjectId: "",
	cohort: "",
	timepoint: "",
	compartment: "",
	locus: "",
	vCall: "",
	vCallIncludeAmbiguous: false,
	d1Call: "",
	d1CallIncludeAmbiguous: false,
	d2Call: "",
	jCall: "",
	jCallIncludeAmbiguous: false,
	cCall: "",
	cCallIncludeAmbiguous: false,
	isotype: "",
	cdr3Nt: "",
	cdr3Aa: "",
	motif: "",
	motifTarget: "cdr3_aa",
	motifSyntax: "substring",
	motifMode: "any",
	productive: "any",
	completeVdj: "any",
	vjInFrame: "any",
	stopCodon: "any",
	hasD: "any",
	hasCdr3: "any",
	doubleD: "any",
	minCdr3NtLength: 0,
	maxCdr3NtLength: 0,
	minCdr3AaLength: 0,
	maxCdr3AaLength: 0,
	minVIdentity: 0,
	minJIdentity: 0,
	minVMutation: 0,
	maxVMutation: 0
};
const IUPAC = {
	A: "A",
	C: "C",
	G: "G",
	T: "[TU]",
	U: "[TU]",
	R: "[AG]",
	Y: "[CTU]",
	S: "[GC]",
	W: "[ATU]",
	K: "[GTU]",
	M: "[AC]",
	B: "[CGTU]",
	D: "[AGTU]",
	H: "[ACTU]",
	V: "[ACG]",
	N: "[ACGTU]",
	X: "."
};
function tokens(value) {
	return value.split(/[\n,;]+/).map((item) => item.trim()).filter(Boolean);
}
function callMatches(value, query, includeAmbiguous) {
	const wanted = tokens(query).map((item) => item.toUpperCase());
	if (!wanted.length) return true;
	const assignments = value.split(",").map((item) => item.trim().toUpperCase()).filter(Boolean);
	if (!includeAmbiguous && assignments.length > 1) return false;
	return wanted.some((item) => assignments.some((assignment) => assignment === item || assignment.replace(/\*.*$/, "") === item || assignment.includes(item)));
}
function facetMatches(value, query) {
	const wanted = tokens(query).map((item) => item.toUpperCase());
	if (!wanted.length) return true;
	const observed = value.toUpperCase();
	return wanted.some((item) => observed === item || observed.includes(item));
}
function boolValue(value) {
	return /^(T|TRUE|YES|1)$/i.test(value.trim());
}
function triMatches(value, wanted) {
	return wanted === "any" || boolValue(value) === (wanted === "yes");
}
function motifRegexes(options) {
	return tokens(options.motif).map((motif) => {
		if (options.motifSyntax === "regex") return new RegExp(motif, "i");
		if (options.motifSyntax === "iupac") {
			const source = [...motif.toUpperCase()].map((symbol) => IUPAC[symbol] ?? symbol.replace(/[.*+?^${}()|[\]\\]/g, "\\$&")).join("");
			return new RegExp(source, "i");
		}
		return new RegExp(motif.replace(/[.*+?^${}()|[\]\\]/g, "\\$&"), "i");
	});
}
function validateRepertoireSelection(options) {
	const errors = [];
	for (const [minimum, maximum, label] of [
		[
			options.minCdr3NtLength,
			options.maxCdr3NtLength,
			"CDR3 nucleotide length"
		],
		[
			options.minCdr3AaLength,
			options.maxCdr3AaLength,
			"CDR3 amino-acid length"
		],
		[
			options.minVMutation,
			options.maxVMutation,
			"V mutation fraction"
		]
	]) if (maximum > 0 && minimum > maximum) errors.push(`${label}: the minimum exceeds the maximum.`);
	try {
		motifRegexes(options);
	} catch (error) {
		errors.push(`Motif expression: ${error instanceof Error ? error.message : String(error)}`);
	}
	return errors;
}
function number(value) {
	const parsed = Number(value);
	return Number.isFinite(parsed) ? parsed : 0;
}
function repertoireRowMatches(row, options, doubleD) {
	const text = (name) => row[name] ?? "";
	if (options.sequenceId && !text("sequence_id").toLowerCase().includes(options.sequenceId.trim().toLowerCase())) return false;
	if (!facetMatches(text("swig_dataset_id"), options.datasetId) || !facetMatches(text("sample_id"), options.sampleId) || !facetMatches(text("subject_id"), options.subjectId) || !facetMatches(text("swig_cohort"), options.cohort) || !facetMatches(text("swig_timepoint"), options.timepoint) || !facetMatches(text("swig_compartment"), options.compartment)) return false;
	if (options.locus && !facetMatches(text("locus"), options.locus)) return false;
	if (!callMatches(text("v_call"), options.vCall, options.vCallIncludeAmbiguous) || !callMatches(doubleD?.values.d_call || text("d_call"), options.d1Call, options.d1CallIncludeAmbiguous) || !facetMatches(doubleD?.values.d2_call || text("d2_call"), options.d2Call) || !callMatches(text("j_call"), options.jCall, options.jCallIncludeAmbiguous) || !callMatches(text("c_call"), options.cCall, options.cCallIncludeAmbiguous) || !facetMatches(text("isotype"), options.isotype)) return false;
	const cdr3Nt = (text("cdr3") || text("junction")).toUpperCase();
	const cdr3Aa = (text("cdr3_aa") || text("junction_aa")).toUpperCase();
	if (options.cdr3Nt && !cdr3Nt.includes(options.cdr3Nt.replace(/\s/g, "").toUpperCase())) return false;
	if (options.cdr3Aa && !cdr3Aa.includes(options.cdr3Aa.replace(/\s/g, "").toUpperCase())) return false;
	if (!triMatches(text("productive"), options.productive) || !triMatches(text("complete_vdj"), options.completeVdj) || !triMatches(text("vj_in_frame"), options.vjInFrame) || !triMatches(text("stop_codon"), options.stopCodon)) return false;
	const dPresent = Boolean(doubleD?.values.d_call || text("d_call"));
	if (options.hasD !== "any" && dPresent !== (options.hasD === "yes")) return false;
	if (options.hasCdr3 !== "any" && Boolean(cdr3Nt || cdr3Aa) !== (options.hasCdr3 === "yes")) return false;
	const isDoubleD = Boolean(doubleD?.values.d2_call || text("d2_call"));
	if (options.doubleD === "positive" && !isDoubleD) return false;
	if (options.doubleD === "negative" && isDoubleD) return false;
	if (options.minCdr3NtLength && cdr3Nt.length < options.minCdr3NtLength) return false;
	if (options.maxCdr3NtLength && cdr3Nt.length > options.maxCdr3NtLength) return false;
	if (options.minCdr3AaLength && cdr3Aa.length < options.minCdr3AaLength) return false;
	if (options.maxCdr3AaLength && cdr3Aa.length > options.maxCdr3AaLength) return false;
	const vIdentity = number(text("v_identity"));
	const jIdentity = number(text("j_identity"));
	if (options.minVIdentity && vIdentity < options.minVIdentity) return false;
	if (options.minJIdentity && jIdentity < options.minJIdentity) return false;
	const vMutation = vIdentity ? 1 - vIdentity : 0;
	if (options.minVMutation && vMutation < options.minVMutation) return false;
	if (options.maxVMutation && (!vIdentity || vMutation > options.maxVMutation)) return false;
	const expressions = motifRegexes(options);
	if (expressions.length) {
		const target = options.motifTarget === "sequence" ? text("sequence_alignment") || text("sequence") : options.motifTarget === "cdr3_nt" ? cdr3Nt : options.motifTarget === "cdr3_aa" ? cdr3Aa : text("junction_aa");
		const matched = expressions.map((expression) => expression.test(target));
		if (options.motifMode === "all" ? matched.some((value) => !value) : matched.every((value) => !value)) return false;
	}
	return true;
}
//#endregion
//#region src/sequence-stream.ts
const DEFAULT_FASTQ_QUALITY_FILTER = {
	enabled: false,
	maximumExpectedErrors: .01,
	phredOffset: 33,
	trim3Prime: {
		enabled: false,
		windowSize: 4,
		minimumMeanPhred: 20,
		minimumLength: 50
	}
};
function emptyFastqQualityFilterStats(enabled = false, applicable = false) {
	return {
		enabled,
		applicable,
		recordsEvaluated: 0,
		recordsRetained: 0,
		recordsPassedThrough: 0,
		recordsRejectedExpectedErrors: 0,
		recordsRejectedMinimumLength: 0,
		recordsTrimmed: 0,
		basesTrimmed: 0
	};
}
function addFastqQualityFilterStats(target, source) {
	return {
		enabled: target.enabled || source.enabled,
		applicable: target.applicable || source.applicable,
		recordsEvaluated: target.recordsEvaluated + source.recordsEvaluated,
		recordsRetained: target.recordsRetained + source.recordsRetained,
		recordsPassedThrough: target.recordsPassedThrough + source.recordsPassedThrough,
		recordsRejectedExpectedErrors: target.recordsRejectedExpectedErrors + source.recordsRejectedExpectedErrors,
		recordsRejectedMinimumLength: target.recordsRejectedMinimumLength + source.recordsRejectedMinimumLength,
		recordsTrimmed: target.recordsTrimmed + source.recordsTrimmed,
		basesTrimmed: target.basesTrimmed + source.basesTrimmed
	};
}
async function* fastaRecords(source) {
	let header = "";
	let sequence = [];
	for await (const line of source) {
		if (line.startsWith(">")) {
			if (header) {
				const joined = sequence.join("");
				if (!joined) throw new Error(`Empty FASTA record: ${header.slice(1).trim() || "unnamed sequence"}.`);
				yield `${header}\n${joined}\n`;
			}
			header = line;
			sequence = [];
			continue;
		}
		if (!line.trim()) continue;
		if (!header) throw new Error("Expected a FASTA header beginning with '>'.");
		sequence.push(line.replace(/\s/g, ""));
	}
	if (header) {
		const joined = sequence.join("");
		if (!joined) throw new Error(`Empty FASTA record: ${header.slice(1).trim() || "unnamed sequence"}.`);
		yield `${header}\n${joined}\n`;
	}
}
async function* fastqRecords(source) {
	let header = "";
	let plus = "";
	let sequence = [];
	let quality = [];
	let sequenceLength = 0;
	let qualityLength = 0;
	let state = "header";
	for await (const line of source) {
		if (state === "header") {
			if (!line.trim()) continue;
			if (!line.startsWith("@")) throw new Error("Expected a FASTQ header beginning with '@'.");
			header = line;
			plus = "";
			sequence = [];
			quality = [];
			sequenceLength = 0;
			qualityLength = 0;
			state = "sequence";
			continue;
		}
		if (state === "sequence") {
			if (line.startsWith("+")) {
				if (!sequenceLength) throw new Error(`Empty FASTQ sequence: ${header.slice(1).trim() || "unnamed sequence"}.`);
				plus = line;
				state = "quality";
			} else {
				const normalized = line.replace(/\s/g, "");
				sequence.push(normalized);
				sequenceLength += normalized.length;
			}
			continue;
		}
		quality.push(line);
		qualityLength += line.length;
		if (qualityLength > sequenceLength) throw new Error(`FASTQ sequence/quality length mismatch: ${header.slice(1).trim() || "unnamed sequence"}.`);
		if (qualityLength === sequenceLength) {
			yield {
				header,
				sequence: sequence.join(""),
				plus,
				quality: quality.join("")
			};
			state = "header";
		}
	}
	if (state === "sequence") throw new Error(`The FASTQ input ended before a '+' line for ${header.slice(1).trim() || "a sequence"}.`);
	if (state === "quality") throw new Error(`The FASTQ quality string is truncated for ${header.slice(1).trim() || "a sequence"}.`);
}
async function* airrRecords(source) {
	let header = "";
	let delimiter = "	";
	let sequenceColumn = -1;
	for await (const line of source) {
		if (!header) {
			if (!line.trim()) continue;
			header = line;
			delimiter = header.includes("	") ? "	" : ",";
			sequenceColumn = header.split(delimiter).indexOf("sequence");
			if (sequenceColumn < 0) throw new Error("AIRR input requires a 'sequence' column.");
			continue;
		}
		if (!line.trim()) continue;
		if (line.split(delimiter)[sequenceColumn]) yield {
			header,
			row: line
		};
	}
}
function expectedErrorTable(offset) {
	const table = new Float64Array(127);
	for (let code = offset; code < table.length; code += 1) table[code] = 10 ** (-(code - offset) / 10);
	return table;
}
const PHRED33_EXPECTED_ERRORS = expectedErrorTable(33);
const PHRED64_EXPECTED_ERRORS = expectedErrorTable(64);
function canonicalFastq(record, end = record.sequence.length) {
	if (end === record.sequence.length) return `${record.header}\n${record.sequence}\n${record.plus}\n${record.quality}\n`;
	return `${record.header}\n${record.sequence.slice(0, end)}\n${record.plus}\n${record.quality.slice(0, end)}\n`;
}
/**
* Apply the FASTQ filter without constructing an output record for rejected
* reads. The full-read expected-error sum is accumulated once; when 3' trim
* removes a base, its contribution is subtracted before the threshold test.
*/
function filterFastqRecord(record, options, stats) {
	stats.recordsEvaluated += 1;
	const offset = options.phredOffset;
	const errorTable = offset === 64 ? PHRED64_EXPECTED_ERRORS : PHRED33_EXPECTED_ERRORS;
	const quality = record.quality;
	const trim = options.trim3Prime;
	const windowSize = Math.max(1, Math.floor(trim.windowSize));
	let end = quality.length;
	const initialWindowStart = Math.max(0, end - windowSize);
	let windowStart = initialWindowStart;
	let windowPhredSum = 0;
	let expectedErrors = 0;
	for (let index = 0; index < quality.length; index += 1) {
		const code = quality.charCodeAt(index);
		if (code < offset || code >= errorTable.length) {
			const name = record.header.slice(1).trim() || "unnamed sequence";
			throw new Error(`FASTQ quality character outside Phred+${offset} range in ${name} at base ${index + 1}.`);
		}
		expectedErrors += errorTable[code];
		if (trim.enabled && index >= initialWindowStart) windowPhredSum += code - offset;
	}
	if (trim.enabled) {
		while (end > 0 && windowPhredSum < trim.minimumMeanPhred * (end - windowStart)) {
			const removedCode = quality.charCodeAt(end - 1);
			expectedErrors -= errorTable[removedCode];
			windowPhredSum -= removedCode - offset;
			end -= 1;
			const nextWindowStart = Math.max(0, end - windowSize);
			if (nextWindowStart < windowStart) {
				const addedCode = quality.charCodeAt(nextWindowStart);
				windowPhredSum += addedCode - offset;
			}
			windowStart = nextWindowStart;
		}
		if (end < Math.max(1, Math.floor(trim.minimumLength))) {
			stats.recordsRejectedMinimumLength += 1;
			return null;
		}
	}
	if (Math.max(0, expectedErrors) > options.maximumExpectedErrors) {
		stats.recordsRejectedExpectedErrors += 1;
		return null;
	}
	const trimmedBases = quality.length - end;
	if (trimmedBases) {
		stats.recordsTrimmed += 1;
		stats.basesTrimmed += trimmedBases;
	}
	stats.recordsRetained += 1;
	return canonicalFastq(record, end);
}
function seededRandom(seed) {
	let value = Math.trunc(seed) >>> 0;
	return () => {
		value = value + 1831565813 >>> 0;
		let mixed = value;
		mixed = Math.imul(mixed ^ mixed >>> 15, mixed | 1);
		mixed ^= mixed + Math.imul(mixed ^ mixed >>> 7, mixed | 61);
		return ((mixed ^ mixed >>> 14) >>> 0) / 4294967296;
	};
}
function validateFastqQualityFilter(options) {
	if (!options.enabled) return;
	if (!Number.isFinite(options.maximumExpectedErrors) || options.maximumExpectedErrors < 0) throw new Error("Maximum FASTQ expected errors must be a non-negative number.");
	if (options.phredOffset !== 33 && options.phredOffset !== 64) throw new Error("FASTQ quality encoding must be Phred+33 or Phred+64.");
	if (options.trim3Prime.enabled && (!Number.isFinite(options.trim3Prime.windowSize) || options.trim3Prime.windowSize < 1 || !Number.isFinite(options.trim3Prime.minimumMeanPhred) || options.trim3Prime.minimumMeanPhred < 0 || !Number.isFinite(options.trim3Prime.minimumLength) || options.trim3Prime.minimumLength < 1)) throw new Error("FASTQ 3' trimming requires a positive window and retained length, and a non-negative mean Phred threshold.");
}
const DEFAULT_DENOISE = {
	mode: "conservative",
	errorRate: .00473,
	alpha: .01,
	callResolution: "allele",
	ambiguity: "strict",
	minimumParentCount: 2,
	ambiguousPolicy: "exclude",
	unresolvedPolicy: "discard",
	fadNeighborThreshold: 1,
	fadMethod: 2,
	expectedZeroErrorFraction: 1,
	maximumHammingDistance: 1,
	maximumEditDistance: 2,
	minimumIndelParentRatio: 2,
	maxCandidatesPerVariant: 5e4,
	scope: "sample",
	respectConstantCall: true
};
const DEFAULT_CLI_CONFIG = {
	schema: 1,
	application: "swig-cli",
	studyName: "swig-study",
	studyDesign: "independent",
	inputs: [],
	references: {
		species: "Homo sapiens",
		scope: "BCR",
		prepareMetadata: true
	},
	preprocessing: {
		fastqFilter: {
			...DEFAULT_FASTQ_QUALITY_FILTER,
			trim3Prime: { ...DEFAULT_FASTQ_QUALITY_FILTER.trim3Prime }
		},
		subsample: {
			enabled: false,
			size: 1e4,
			seed: 1
		}
	},
	annotation: {
		workers: 0,
		batchRecords: 0,
		callingProfile: "truth_optimized",
		assignerStrategy: "riat_mp",
		minimumIdentity: .6,
		strand: 0,
		airrMode: "preserve",
		doubleD: {
			mode: "off",
			minimumVjSpan: 40,
			seedLength: 11,
			pseudoTrim: 5,
			maximumPseudoMismatches: 3,
			minimumScoreGain: 8
		}
	},
	pipeline: {
		collapse: {
			...DEFAULT_PIPELINE_PLAN.collapse,
			enabled: true,
			denoise: { ...DEFAULT_DENOISE }
		},
		chimera: {
			...DEFAULT_PIPELINE_PLAN.chimera,
			enabled: false,
			priorProbability: .05,
			baseMutationProbability: .05,
			mutationRates: [
				0,
				.0179,
				.0357,
				.0536,
				.0714,
				.0893,
				.1071,
				.125,
				.1429,
				.1607,
				.1786,
				.1964,
				.2143,
				.2321,
				.25
			],
			mutationSwitchProbability: 0,
			minimumDfr: 1,
			detailed: false
		},
		selection: {
			...DEFAULT_REPERTOIRE_SELECTION,
			enabled: false
		},
		alleleRefinement: {
			...DEFAULT_ALLELE_REFINEMENT_OPTIONS,
			segments: [...DEFAULT_ALLELE_REFINEMENT_OPTIONS.segments],
			enabled: false,
			reassignmentPolicy: "confidence",
			applyMinimumPosterior: .8
		},
		lineage: {
			...DEFAULT_PIPELINE_PLAN.lineage,
			enabled: true,
			maxCandidateComparisons: 5e4
		},
		shm: {
			...DEFAULT_PIPELINE_PLAN.shm,
			enabled: true
		},
		missingAlleles: {
			...DEFAULT_MISSING_ALLELE_OPTIONS,
			enabled: false
		}
	},
	output: {
		directory: "swig-output",
		prefix: "swig",
		writeAnnotatedAirr: true,
		writeLineageStudy: true
	}
};
function finite(value, fallback) {
	const parsed = Number(value);
	return Number.isFinite(parsed) ? parsed : fallback;
}
function safeIdentifier(value, fallback) {
	return value.trim().replace(/[^A-Za-z0-9_.-]+/g, "-").replace(/^-+|-+$/g, "") || fallback;
}
function normalizeCliConfig(value) {
	const collapse = value.pipeline?.collapse ?? {};
	const chimera = value.pipeline?.chimera ?? {};
	const selection = value.pipeline?.selection ?? {};
	const allele = value.pipeline?.alleleRefinement ?? {};
	const lineage = value.pipeline?.lineage ?? {};
	const shm = value.pipeline?.shm ?? {};
	const missing = value.pipeline?.missingAlleles ?? {};
	const fastqFilter = value.preprocessing?.fastqFilter ?? {};
	const trim3Prime = fastqFilter.trim3Prime ?? {};
	const subsample = value.preprocessing?.subsample ?? {};
	const inputs = (value.inputs ?? []).map((input, index) => {
		const stem = input.path.split(/[\\/]/).pop()?.replace(/\.(?:fasta?|fna|fas|fastq|fq|airr\.tsv|tsv)(?:\.gz)?$/i, "") || `dataset-${index + 1}`;
		const datasetId = safeIdentifier(input.datasetId ?? stem, `dataset-${index + 1}`);
		const sampleId = safeIdentifier(input.sampleId ?? datasetId, `sample-${index + 1}`);
		return {
			...input,
			format: input.format ?? "auto",
			datasetId,
			sampleId,
			subjectId: safeIdentifier(input.subjectId ?? sampleId, `subject-${index + 1}`),
			cohort: input.cohort ?? "",
			timepoint: input.timepoint ?? "",
			compartment: input.compartment ?? ""
		};
	});
	const annotation = {
		...DEFAULT_CLI_CONFIG.annotation,
		...value.annotation,
		doubleD: {
			...DEFAULT_CLI_CONFIG.annotation.doubleD,
			...value.annotation?.doubleD
		}
	};
	annotation.workers = Math.max(0, Math.floor(finite(annotation.workers, 0)));
	annotation.batchRecords = Math.max(0, Math.floor(finite(annotation.batchRecords, 0)));
	annotation.minimumIdentity = Math.max(0, Math.min(1, finite(annotation.minimumIdentity, .6)));
	const normalizedSubsample = {
		...DEFAULT_CLI_CONFIG.preprocessing.subsample,
		...subsample
	};
	normalizedSubsample.enabled = Boolean(normalizedSubsample.enabled);
	normalizedSubsample.size = Math.max(1, Math.floor(finite(normalizedSubsample.size, 1e4)));
	normalizedSubsample.seed = Math.trunc(finite(normalizedSubsample.seed, 1));
	return {
		...DEFAULT_CLI_CONFIG,
		...value,
		schema: 1,
		application: "swig-cli",
		inputs,
		references: {
			...DEFAULT_CLI_CONFIG.references,
			...value.references,
			prepareMetadata: value.references?.prepareMetadata !== false,
			inline: value.references?.inline ? { ...value.references.inline } : void 0,
			files: value.references?.files ? { ...value.references.files } : void 0
		},
		preprocessing: {
			fastqFilter: {
				...DEFAULT_CLI_CONFIG.preprocessing.fastqFilter,
				...fastqFilter,
				trim3Prime: {
					...DEFAULT_CLI_CONFIG.preprocessing.fastqFilter.trim3Prime,
					...trim3Prime
				}
			},
			subsample: normalizedSubsample
		},
		annotation,
		pipeline: {
			collapse: {
				...DEFAULT_CLI_CONFIG.pipeline.collapse,
				...collapse,
				denoise: {
					...DEFAULT_DENOISE,
					...collapse.denoise,
					mode: collapse.mode ?? collapse.denoise?.mode ?? DEFAULT_DENOISE.mode,
					scope: collapse.scope ?? collapse.denoise?.scope ?? DEFAULT_DENOISE.scope,
					respectConstantCall: collapse.respectConstantCall ?? collapse.denoise?.respectConstantCall ?? true
				}
			},
			chimera: {
				...DEFAULT_CLI_CONFIG.pipeline.chimera,
				...chimera,
				mutationRates: [...chimera.mutationRates ?? DEFAULT_CLI_CONFIG.pipeline.chimera.mutationRates]
			},
			selection: {
				...DEFAULT_CLI_CONFIG.pipeline.selection,
				...selection
			},
			alleleRefinement: {
				...DEFAULT_CLI_CONFIG.pipeline.alleleRefinement,
				...allele,
				segments: [...allele.segments ?? DEFAULT_CLI_CONFIG.pipeline.alleleRefinement.segments]
			},
			lineage: {
				...DEFAULT_CLI_CONFIG.pipeline.lineage,
				...lineage
			},
			shm: {
				...DEFAULT_CLI_CONFIG.pipeline.shm,
				...shm
			},
			missingAlleles: {
				...DEFAULT_CLI_CONFIG.pipeline.missingAlleles,
				...missing
			}
		},
		output: {
			...DEFAULT_CLI_CONFIG.output,
			...value.output
		}
	};
}
//#endregion
//#region src/germline-preprocess.ts
const LOCI = [
	"IGH",
	"IGK",
	"IGL",
	"TRA",
	"TRB",
	"TRD",
	"TRG"
];
const IMGT_V_GAPPED_ENDS = [
	78,
	114,
	165,
	195,
	312
];
const EMPTY_BOUNDS = [
	-1,
	-1,
	-1,
	-1,
	-1,
	-1,
	-1,
	-1,
	-1,
	-1
];
const SOURCE_IMGT_GAPPED = 1;
const SOURCE_AIRR_C = 2;
const SOURCE_TRANSFERRED_IMGT = 3;
const SOURCE_VALIDATED_J_MOTIF = 4;
const SOURCE_PROVIDED = 5;
const SOURCE_TRANSFERRED_J = 6;
const MATCH_PRESETS = {
	strict: {
		vSameGeneMinIdentity: .8,
		vNearestMinIdentity: .72,
		jSameGeneMinIdentity: .75,
		jNearestMinIdentity: .68,
		nearestCandidates: 12
	},
	permissive: {
		vSameGeneMinIdentity: .65,
		vNearestMinIdentity: .55,
		jSameGeneMinIdentity: .6,
		jNearestMinIdentity: .5,
		nearestCandidates: 32
	},
	best_guess: {
		vSameGeneMinIdentity: 0,
		vNearestMinIdentity: 0,
		jSameGeneMinIdentity: 0,
		jNearestMinIdentity: 0,
		nearestCandidates: 64
	}
};
function identityOption(value, fallback, label) {
	const resolved = value ?? fallback;
	if (!Number.isFinite(resolved) || resolved < 0 || resolved > 1) throw new Error(`${label} must be between 0 and 1.`);
	return resolved;
}
function resolveGermlineMatchOptions(options = {}) {
	const mode = options.mode ?? "strict";
	const preset = MATCH_PRESETS[mode];
	if (!preset) throw new Error(`Unsupported germline metadata match mode: ${String(mode)}.`);
	const nearestCandidates = options.nearestCandidates ?? preset.nearestCandidates;
	if (!Number.isSafeInteger(nearestCandidates) || nearestCandidates < 1 || nearestCandidates > 1e4) throw new Error("nearestCandidates must be an integer between 1 and 10000.");
	return {
		mode,
		vSameGeneMinIdentity: identityOption(options.vSameGeneMinIdentity, preset.vSameGeneMinIdentity, "vSameGeneMinIdentity"),
		vNearestMinIdentity: identityOption(options.vNearestMinIdentity, preset.vNearestMinIdentity, "vNearestMinIdentity"),
		jSameGeneMinIdentity: identityOption(options.jSameGeneMinIdentity, preset.jSameGeneMinIdentity, "jSameGeneMinIdentity"),
		jNearestMinIdentity: identityOption(options.jNearestMinIdentity, preset.jNearestMinIdentity, "jNearestMinIdentity"),
		nearestCandidates,
		includeDiagnostics: Boolean(options.includeDiagnostics)
	};
}
const V_CHAIN_LOCI = {
	VH: "IGH",
	VK: "IGK",
	VL: "IGL",
	VA: "TRA",
	VB: "TRB",
	VD: "TRD",
	VG: "TRG"
};
const J_CHAIN_LOCI = {
	JH: "IGH",
	JK: "IGK",
	JL: "IGL",
	JA: "TRA",
	JB: "TRB",
	JD: "TRD",
	JG: "TRG"
};
const METADATA_SOURCE_LABELS = {
	[SOURCE_IMGT_GAPPED]: "IMGT-gapped delineation",
	[SOURCE_AIRR_C]: "AIRR-C annotation",
	[SOURCE_TRANSFERRED_IMGT]: "validated IMGT-boundary transfer",
	[SOURCE_VALIDATED_J_MOTIF]: "frame-validated J motif",
	[SOURCE_PROVIDED]: "provided annotation",
	[SOURCE_TRANSFERRED_J]: "validated J-anchor transfer"
};
function positiveModulo(value, modulus) {
	return (value % modulus + modulus) % modulus;
}
function compactMetadata(frame = -1, cdr3Stop = -1, bounds = EMPTY_BOUNDS, source = 0) {
	return [
		frame,
		cdr3Stop,
		...bounds,
		source
	];
}
function inferLocus(value) {
	const upper = value.toUpperCase();
	return LOCI.find((locus) => upper.includes(locus)) ?? null;
}
function inferSegment(name) {
	const match = /^(?:IGH|IGK|IGL|TRA|TRB|TRD|TRG)([VDJC])/.exec(name.toUpperCase());
	return match ? match[1] : null;
}
function germlineName(header) {
	const identifier = header.trim().split(/\s+/, 1)[0] ?? "";
	const fields = identifier.split("|");
	return fields.find((field) => inferLocus(field) && field.includes("*")) ?? fields.find((field) => inferLocus(field)) ?? identifier;
}
function normalizeSequence(raw) {
	const compact = raw.toUpperCase().replaceAll("U", "T").replace(/\s/g, "");
	if (/[^ACGTNRYKMSWBDHV.\-]/.test(compact)) throw new Error("Germline FASTA contains characters outside the IUPAC nucleotide alphabet.");
	const withoutGaps = compact.replace(/[.\-]/g, "");
	const ambiguous = (withoutGaps.match(/[^ACGT]/g) ?? []).length;
	return {
		sequence: withoutGaps.replace(/[^ACGT]/g, "N"),
		ambiguous
	};
}
function normalizeIndexSequence(raw) {
	const compact = raw.toUpperCase().replaceAll("U", "T").replace(/\s/g, "");
	if (/[^ACGTNRYKMSWBDHV.\-]/.test(compact)) throw new Error("Germline FASTA contains characters outside the IUPAC nucleotide alphabet.");
	const sequence = compact.replace(/[.\-]/g, "");
	return {
		sequence,
		ambiguous: (sequence.match(/[^ACGT]/g) ?? []).length
	};
}
function parseFasta(text) {
	const records = [];
	let header = "";
	let sequence = [];
	for (const line of text.split(/\r?\n/)) if (line.startsWith(">")) {
		if (header) records.push({
			header,
			rawSequence: sequence.join("")
		});
		header = line.slice(1).trim();
		sequence = [];
	} else if (header && line.trim()) sequence.push(line.trim());
	else if (!header && line.trim()) throw new Error("Expected a FASTA header beginning with '>'.");
	if (header) records.push({
		header,
		rawSequence: sequence.join("")
	});
	if (!records.length) throw new Error("No germline FASTA records were found.");
	return records;
}
function parseImgtFields(header) {
	const fields = (header.split(/\s+/, 1)[0] ?? "").split("|");
	return fields.length >= 8 && fields[1] && fields[4] ? fields : null;
}
function imgtFrame(fields) {
	const value = Number.parseInt(fields?.[7] ?? "", 10);
	return value >= 1 && value <= 3 ? value - 1 : -1;
}
function metadataFromHeader(header) {
	const match = /(?:^|\s)SWIGMETA=([\d,\-]+)/i.exec(header);
	if (!match) return void 0;
	const values = match[1].split(",").map(Number);
	if (values.length !== 13 || values.some((value) => !Number.isInteger(value))) return void 0;
	return values;
}
function metadataHeader(metadata) {
	return metadata ? ` SWIGMETA=${metadata.join(",")}` : "";
}
function hasRegionMetadata(metadata) {
	return Boolean(metadata && metadata.slice(2, 12).every((value) => value >= 0));
}
function hasJMetadata(metadata) {
	return Boolean(metadata && metadata[0] >= 0 && metadata[1] >= 0);
}
function validateMetadata(metadata, sequenceLength, segment) {
	if (metadata[0] < -1 || metadata[0] > 2 || metadata[1] < -1 || metadata[1] >= sequenceLength) return false;
	if (segment === "V" && hasRegionMetadata(metadata)) {
		const bounds = metadata.slice(2, 12);
		if (bounds.some((value, index) => value < 0 || value > sequenceLength || index && value < bounds[index - 1])) return false;
		for (let index = 0; index < bounds.length; index += 2) if (bounds[index + 1] <= bounds[index]) return false;
	}
	return true;
}
function nearestFrameCysEnd(sequence, predictedEnd, frame, maximumDistance = 60) {
	const candidates = [];
	const start = Math.max(0, predictedEnd - maximumDistance);
	const stop = Math.min(sequence.length - 2, predictedEnd + maximumDistance);
	for (let index = start; index <= stop; index += 1) {
		const codon = sequence.slice(index, index + 3);
		if (codon !== "TGT" && codon !== "TGC") continue;
		if (frame >= 0 && positiveModulo(index - frame, 3) !== 0) continue;
		candidates.push(index + 3);
	}
	return candidates.sort((left, right) => Math.abs(left - predictedEnd) - Math.abs(right - predictedEnd) || right - left)[0];
}
function imgtVMetadata(rawSequence, fields) {
	if (fields?.[4]?.toUpperCase() !== "V-REGION") return void 0;
	const gapped = rawSequence.toUpperCase().replaceAll("U", "T").replace(/\s/g, "");
	if (gapped.length < IMGT_V_GAPPED_ENDS.at(-1)) return void 0;
	const sequence = gapped.replace(/[.\-]/g, "");
	const ends = IMGT_V_GAPPED_ENDS.map((end) => gapped.slice(0, end).replace(/[.\-]/g, "").length);
	const anchorEnd = nearestFrameCysEnd(sequence, ends[4], imgtFrame(fields));
	if (anchorEnd && anchorEnd > ends[3]) ends[4] = anchorEnd;
	const bounds = [
		0,
		ends[0],
		ends[0],
		ends[1],
		ends[1],
		ends[2],
		ends[2],
		ends[3],
		ends[3],
		ends[4]
	];
	const result = compactMetadata(imgtFrame(fields), -1, bounds, SOURCE_IMGT_GAPPED);
	return validateMetadata(result, sequence.length, "V") ? result : void 0;
}
function isJAnchor(sequence, index) {
	const aromatic = sequence.slice(index, index + 3);
	return (aromatic === "TGG" || aromatic === "TTT" || aromatic === "TTC") && /^GG[ACGT]$/.test(sequence.slice(index + 3, index + 6));
}
function jMetadata(sequence, fields) {
	const providedFrame = imgtFrame(fields);
	const candidates = [];
	const limit = Math.min(sequence.length - 6, 54);
	for (let index = Math.max(0, providedFrame); index <= limit; index += providedFrame >= 0 ? 3 : 1) if (isJAnchor(sequence, index)) candidates.push(index);
	if (!candidates.length || providedFrame < 0 && candidates.length !== 1) return void 0;
	const anchor = candidates[0];
	const frame = providedFrame >= 0 ? providedFrame : anchor % 3;
	if (anchor % 3 !== frame) return void 0;
	return compactMetadata(frame, anchor - 1, EMPTY_BOUNDS, SOURCE_VALIDATED_J_MOTIF);
}
function globalAlignment(query, reference) {
	const columns = reference.length + 1;
	const rows = query.length + 1;
	const scores = new Int16Array(columns);
	const next = new Int16Array(columns);
	const trace = new Uint8Array(rows * columns);
	for (let column = 1; column < columns; column += 1) {
		scores[column] = -3 * column;
		trace[column] = 2;
	}
	for (let row = 1; row < rows; row += 1) {
		next[0] = -3 * row;
		trace[row * columns] = 1;
		for (let column = 1; column < columns; column += 1) {
			const diagonal = scores[column - 1] + (query[row - 1] === reference[column - 1] ? 2 : -2);
			const up = scores[column] - 3;
			const left = next[column - 1] - 3;
			const best = Math.max(diagonal, up, left);
			next[column] = best;
			trace[row * columns + column] = best === diagonal ? 0 : best === up ? 1 : 2;
		}
		scores.set(next);
	}
	let row = query.length;
	let column = reference.length;
	const alignedQuery = [];
	const alignedReference = [];
	let matches = 0;
	while (row || column) {
		const direction = trace[row * columns + column];
		if (row && column && direction === 0) {
			const queryBase = query[--row];
			const referenceBase = reference[--column];
			alignedQuery.push(queryBase);
			alignedReference.push(referenceBase);
			if (queryBase === referenceBase) matches += 1;
		} else if (row && (direction === 1 || !column)) {
			alignedQuery.push(query[--row]);
			alignedReference.push("-");
		} else {
			alignedQuery.push("-");
			alignedReference.push(reference[--column]);
		}
	}
	alignedQuery.reverse();
	alignedReference.reverse();
	return {
		query: alignedQuery.join(""),
		reference: alignedReference.join(""),
		identity: matches / Math.max(1, alignedQuery.length)
	};
}
function canonicalAllele(name) {
	return name.replace(/_S\d+$/i, "").toUpperCase();
}
function canonicalGene(name) {
	return canonicalAllele(name).split("*", 1)[0];
}
function kmerSet(sequence, size = 9) {
	const output = /* @__PURE__ */ new Set();
	for (let index = 0; index + size <= sequence.length; index += 1) {
		const kmer = sequence.slice(index, index + size);
		if (!kmer.includes("N")) output.add(kmer);
	}
	return output;
}
const TEMPLATE_KMER_CACHE = /* @__PURE__ */ new WeakMap();
function templateKmerSet(template) {
	const cached = TEMPLATE_KMER_CACHE.get(template);
	if (cached) return cached;
	const kmers = kmerSet(template[1]);
	TEMPLATE_KMER_CACHE.set(template, kmers);
	return kmers;
}
function nearestTemplates(query, templates, limit = 12) {
	const queryKmers = kmerSet(query);
	return templates.map((template) => {
		const templateKmers = templateKmerSet(template);
		let shared = 0;
		for (const kmer of templateKmers) if (queryKmers.has(kmer)) shared += 1;
		const union = queryKmers.size + templateKmers.size - shared;
		return {
			template,
			similarity: union ? shared / union : 0,
			lengthDelta: Math.abs(query.length - template[1].length)
		};
	}).sort((a, b) => b.similarity - a.similarity || a.lengthDelta - b.lengthDelta).slice(0, limit).map(({ template }) => template);
}
function templateCandidateGroups(queryName, templates, eligible = hasRegionMetadata) {
	const delineated = templates.filter((template) => eligible(template[2]));
	const alleleName = canonicalAllele(queryName);
	const geneName = canonicalGene(queryName);
	const exactAllele = delineated.filter(([name]) => canonicalAllele(name) === alleleName);
	const exactGene = delineated.filter(([name]) => canonicalGene(name) === geneName);
	const preferred = exactAllele.length ? exactAllele : exactGene;
	const preferredKeys = new Set(preferred.map(([name, sequence]) => `${name}\u0000${sequence}`));
	return {
		preferred,
		preferredKind: exactAllele.length ? "same_allele" : "same_gene",
		remaining: delineated.filter(([name, sequence]) => !preferredKeys.has(`${name}\u0000${sequence}`))
	};
}
function mapReferenceCoordinates(alignment, referenceLength) {
	const mapped = new Array(referenceLength + 1).fill(-1);
	let queryPosition = 0;
	let referencePosition = 0;
	mapped[0] = 0;
	for (let column = 0; column < alignment.query.length; column += 1) {
		const queryBase = alignment.query[column];
		if (alignment.reference[column] === "-") {
			if (queryBase !== "-") queryPosition += 1;
			mapped[referencePosition] = queryPosition;
			continue;
		}
		if (mapped[referencePosition] < 0) mapped[referencePosition] = queryPosition;
		if (queryBase !== "-") queryPosition += 1;
		referencePosition += 1;
		mapped[referencePosition] = queryPosition;
	}
	return mapped;
}
function incrementRejection(result, reason) {
	result.rejectionCounts[reason] = (result.rejectionCounts[reason] ?? 0) + 1;
}
function updateBest(result, template, identity) {
	if (result.bestIdentity === void 0 || identity > result.bestIdentity) {
		result.bestIdentity = identity;
		result.bestCandidate = template[0];
	}
}
function matchKind(name, templateName, fallback) {
	if (canonicalAllele(templateName) === canonicalAllele(name)) return "same_allele";
	if (canonicalGene(templateName) === canonicalGene(name)) return "same_gene";
	return fallback;
}
function transferMetadata(name, sequence, templates, options) {
	const result = {
		attemptedCandidates: 0,
		rejectionCounts: {}
	};
	const groups = templateCandidateGroups(name, templates);
	const candidates = [{
		values: groups.preferred,
		fallback: groups.preferredKind
	}, {
		values: null,
		fallback: "nearest"
	}];
	for (const group of candidates) {
		const values = group.values ?? nearestTemplates(sequence, groups.remaining, options.nearestCandidates);
		if (!values.length) continue;
		let selected;
		for (const template of values) {
			const alignment = globalAlignment(sequence, template[1]);
			const named = canonicalGene(template[0]) === canonicalGene(name);
			const relation = matchKind(name, template[0], group.fallback);
			result.attemptedCandidates += 1;
			updateBest(result, template, alignment.identity);
			if (alignment.identity < (named ? options.vSameGeneMinIdentity : options.vNearestMinIdentity)) {
				incrementRejection(result, "below_identity");
				continue;
			}
			const templateMetadata = template[2];
			const templateBounds = templateMetadata.slice(2, 12);
			const mapped = mapReferenceCoordinates(alignment, template[1].length);
			const bounds = templateBounds.map((boundary) => mapped[boundary]);
			if (bounds.some((value, index) => value < 0 || value > sequence.length || index && value < bounds[index - 1])) {
				incrementRejection(result, "unmapped_or_nonmonotonic_boundary");
				continue;
			}
			let valid = true;
			for (let index = 0; index < bounds.length; index += 2) if (bounds[index + 1] <= bounds[index]) valid = false;
			if (!valid) {
				incrementRejection(result, "empty_region");
				continue;
			}
			const templateFrame = templateMetadata[0] >= 0 ? templateMetadata[0] : templateBounds[0] % 3;
			const frame = positiveModulo(bounds[0] + templateFrame - templateBounds[0], 3);
			const anchorEnd = nearestFrameCysEnd(sequence, bounds[9], frame, 24);
			if (!anchorEnd || anchorEnd <= bounds[8]) {
				incrementRejection(result, "missing_frame_consistent_v_anchor");
				continue;
			}
			bounds[9] = anchorEnd;
			const metadata = compactMetadata(frame, -1, bounds, SOURCE_TRANSFERRED_IMGT);
			if (!validateMetadata(metadata, sequence.length, "V")) {
				incrementRejection(result, "invalid_projected_metadata");
				continue;
			}
			if (!selected || alignment.identity > selected.identity) selected = {
				metadata,
				identity: alignment.identity,
				template: template[0],
				matchKind: relation
			};
		}
		if (selected) return {
			...result,
			...selected
		};
	}
	if (!result.attemptedCandidates) incrementRejection(result, "no_annotated_template_candidates");
	return result;
}
function transferJMetadata(name, sequence, templates, options) {
	const result = {
		attemptedCandidates: 0,
		rejectionCounts: {}
	};
	const groups = templateCandidateGroups(name, templates, hasJMetadata);
	const candidates = [{
		values: groups.preferred,
		fallback: groups.preferredKind
	}, {
		values: null,
		fallback: "nearest"
	}];
	for (const group of candidates) {
		const values = group.values ?? nearestTemplates(sequence, groups.remaining, options.nearestCandidates);
		if (!values.length) continue;
		let selected;
		for (const template of values) {
			const alignment = globalAlignment(sequence, template[1]);
			const named = canonicalGene(template[0]) === canonicalGene(name);
			const relation = matchKind(name, template[0], group.fallback);
			result.attemptedCandidates += 1;
			updateBest(result, template, alignment.identity);
			if (alignment.identity < (named ? options.jSameGeneMinIdentity : options.jNearestMinIdentity)) {
				incrementRejection(result, "below_identity");
				continue;
			}
			const referenceAnchor = template[2][1] + 1;
			if (referenceAnchor < 0 || referenceAnchor + 6 > template[1].length) {
				incrementRejection(result, "invalid_template_j_anchor");
				continue;
			}
			const mapped = mapReferenceCoordinates(alignment, template[1].length);
			const anchor = mapped[referenceAnchor];
			const anchorEnd = mapped[referenceAnchor + 6];
			if (anchor < 0 || anchorEnd - anchor !== 6) {
				incrementRejection(result, "incomplete_j_anchor_projection");
				continue;
			}
			if (!isJAnchor(sequence, anchor)) {
				incrementRejection(result, "target_j_motif_mismatch");
				continue;
			}
			const transferred = compactMetadata(anchor % 3, anchor - 1, EMPTY_BOUNDS, SOURCE_TRANSFERRED_J);
			if (!validateMetadata(transferred, sequence.length, "J")) {
				incrementRejection(result, "invalid_projected_metadata");
				continue;
			}
			if (!selected || alignment.identity > selected.identity) selected = {
				metadata: transferred,
				identity: alignment.identity,
				template: template[0],
				matchKind: relation
			};
		}
		if (selected) return {
			...result,
			...selected
		};
	}
	if (!result.attemptedCandidates) incrementRejection(result, "no_annotated_template_candidates");
	return result;
}
function annotationPresent(segment, metadata) {
	if (segment === "V") return hasRegionMetadata(metadata);
	if (segment === "J") return hasJMetadata(metadata);
	return true;
}
function integerField(value, context) {
	if (value === void 0 || !/^-?\d+$/.test(value)) throw new Error(`${context} must be an integer.`);
	return Number.parseInt(value, 10);
}
function annotationDataLines(text) {
	const result = [];
	text.split(/\r?\n/).forEach((raw, index) => {
		const line = raw.trim();
		if (!line || line.startsWith("#")) return;
		result.push({
			fields: line.split(/\s+/),
			line: index + 1
		});
	});
	return result;
}
function assertChainLocus(name, locus, chain, chainLoci, source) {
	const declared = chainLoci[chain.toUpperCase()];
	if (declared && declared !== locus) throw new Error(`${source} declares ${name} as ${chain} (${declared}), but its germline identifier is ${locus}.`);
}
/**
* Perform the FASTA preparation done by IgBLAST's `edit_imgt_file.pl` where it
* is relevant to SwiftIG: canonicalize the germline identifier and remove IMGT
* alignment gaps before indexing. Existing SWIGMETA is deliberately removed.
*/
function prepareIgblastStyleGermlineFasta(text, segment, allowedLoci = LOCI) {
	const inputRecords = parseFasta(text);
	const seen = /* @__PURE__ */ new Set();
	const loci = /* @__PURE__ */ new Set();
	let ambiguousBases = 0;
	return {
		fasta: inputRecords.map((input) => {
			const name = germlineName(input.header);
			if (!name) throw new Error(`A ${segment} germline record has an empty identifier.`);
			if (seen.has(name)) throw new Error(`Duplicate ${segment} germline identifier: ${name}.`);
			seen.add(name);
			const locus = inferLocus(`${name} ${input.header}`);
			if (!locus) throw new Error(`${name} does not identify one of the supported IG/TR loci.`);
			if (!allowedLoci.includes(locus)) throw new Error(`${name} belongs to ${locus}, outside the selected ${allowedLoci.join("/")} search space.`);
			const detectedSegment = inferSegment(name);
			if (detectedSegment && detectedSegment !== segment) throw new Error(`${name} looks like a ${detectedSegment} gene, not a ${segment} gene.`);
			let normalized;
			try {
				normalized = normalizeIndexSequence(input.rawSequence);
			} catch (error) {
				throw new Error(`${name}: ${error instanceof Error ? error.message : String(error)}`);
			}
			if (!normalized.sequence) throw new Error(`${name} has an empty nucleotide sequence.`);
			ambiguousBases += normalized.ambiguous;
			loci.add(locus);
			return `>${name}\n${normalized.sequence}\n`;
		}).join(""),
		count: inputRecords.length,
		ambiguousBases,
		loci: [...loci].sort()
	};
}
/** Apply the exact, 1-based inclusive V-domain coordinates in IgBLAST `.ndm.imgt` data. */
function applyIgblastInternalData(fasta, data) {
	const entries = /* @__PURE__ */ new Map();
	for (const record of annotationDataLines(data)) {
		if (record.fields.length !== 13) throw new Error(`Invalid IgBLAST internal-data record at line ${record.line}; expected 13 fields.`);
		const bounds = [];
		for (let index = 0; index < 10; index += 2) {
			const start = integerField(record.fields[index + 1], `IgBLAST internal-data line ${record.line} start`);
			const stop = integerField(record.fields[index + 2], `IgBLAST internal-data line ${record.line} stop`);
			if (start < 1 || stop < start) throw new Error(`Invalid 1-based FWR/CDR interval at IgBLAST internal-data line ${record.line}.`);
			bounds.push(start - 1, stop);
		}
		const frame = integerField(record.fields[12], `IgBLAST internal-data line ${record.line} coding frame`);
		if (frame < -1 || frame > 2) throw new Error(`Invalid coding frame at IgBLAST internal-data line ${record.line}.`);
		entries.set(record.fields[0], {
			bounds,
			chain: record.fields[11],
			frame
		});
	}
	if (!entries.size) throw new Error("The IgBLAST internal-data file contains no annotation records.");
	const prepared = prepareIgblastStyleGermlineFasta(fasta, "V");
	let matched = 0;
	const unmatched = [];
	return {
		fasta: parseFasta(prepared.fasta).map((record) => {
			const name = germlineName(record.header);
			const normalized = normalizeIndexSequence(record.rawSequence);
			const locus = inferLocus(name);
			const entry = entries.get(name);
			if (!entry) {
				unmatched.push(name);
				return `>${name}\n${normalized.sequence}\n`;
			}
			assertChainLocus(name, locus, entry.chain, V_CHAIN_LOCI, "IgBLAST internal data");
			const metadata = compactMetadata(entry.frame, -1, entry.bounds, SOURCE_PROVIDED);
			if (!validateMetadata(metadata, normalized.sequence.length, "V")) throw new Error(`IgBLAST internal-data coordinates for ${name} exceed or contradict its ungapped V sequence.`);
			matched += 1;
			return `>${name}${metadataHeader(metadata)}\n${normalized.sequence}\n`;
		}).join(""),
		matched,
		annotated: matched,
		total: prepared.count,
		unmatched
	};
}
/**
* Apply IgBLAST J auxiliary data. Frame and CDR3-stop coordinates are 0-based;
* the optional fifth field is retained separately for FWR4 end trimming.
*/
function applyIgblastAuxiliaryData(fasta, data) {
	const entries = /* @__PURE__ */ new Map();
	for (const record of annotationDataLines(data)) {
		if (record.fields.length < 3 || record.fields.length > 5) throw new Error(`Invalid IgBLAST auxiliary record at line ${record.line}; expected 3 to 5 fields.`);
		const frame = integerField(record.fields[1], `IgBLAST auxiliary line ${record.line} coding frame`);
		if (frame < -1 || frame > 2) throw new Error(`Invalid coding frame at IgBLAST auxiliary line ${record.line}.`);
		const cdr3Stop = record.fields.length >= 4 ? integerField(record.fields[3], `IgBLAST auxiliary line ${record.line} CDR3 stop`) : -1;
		if (cdr3Stop < -1) throw new Error(`Invalid CDR3 stop at IgBLAST auxiliary line ${record.line}.`);
		const fwr4EndOffset = record.fields.length === 5 ? integerField(record.fields[4], `IgBLAST auxiliary line ${record.line} FWR4 end offset`) : void 0;
		if (fwr4EndOffset !== void 0 && fwr4EndOffset < 0) throw new Error(`Invalid FWR4 end offset at IgBLAST auxiliary line ${record.line}.`);
		entries.set(record.fields[0], {
			frame,
			chain: record.fields[2],
			cdr3Stop,
			fwr4EndOffset
		});
	}
	if (!entries.size) throw new Error("The IgBLAST auxiliary file contains no annotation records.");
	const prepared = prepareIgblastStyleGermlineFasta(fasta, "J");
	let matched = 0;
	let annotated = 0;
	const unmatched = [];
	const fwr4EndOffsets = {};
	return {
		fasta: parseFasta(prepared.fasta).map((record) => {
			const name = germlineName(record.header);
			const normalized = normalizeIndexSequence(record.rawSequence);
			const locus = inferLocus(name);
			const entry = entries.get(name);
			if (!entry) {
				unmatched.push(name);
				return `>${name}\n${normalized.sequence}\n`;
			}
			assertChainLocus(name, locus, entry.chain, J_CHAIN_LOCI, "IgBLAST auxiliary data");
			const metadata = compactMetadata(entry.frame, entry.cdr3Stop, EMPTY_BOUNDS, SOURCE_PROVIDED);
			if (!validateMetadata(metadata, normalized.sequence.length, "J")) throw new Error(`IgBLAST auxiliary coordinates for ${name} exceed its ungapped J sequence.`);
			matched += 1;
			if (entry.cdr3Stop >= 0) annotated += 1;
			if (entry.fwr4EndOffset !== void 0) fwr4EndOffsets[name] = entry.fwr4EndOffset;
			return `>${name}${metadataHeader(metadata)}\n${normalized.sequence}\n`;
		}).join(""),
		matched,
		annotated,
		total: prepared.count,
		unmatched,
		fwr4EndOffsets
	};
}
/** Apply the exact 0-based coding-frame-one starts in IgBLAST `-d_frame_data`. */
function applyIgblastDFrameData(fasta, data) {
	const entries = /* @__PURE__ */ new Map();
	for (const record of annotationDataLines(data)) {
		if (record.fields.length !== 2) throw new Error(`Invalid IgBLAST D-frame record at line ${record.line}; expected 2 fields.`);
		const frame = integerField(record.fields[1], `IgBLAST D-frame line ${record.line} start`);
		if (frame < -1 || frame > 2) throw new Error(`Invalid D-frame start at line ${record.line}.`);
		entries.set(record.fields[0], frame);
	}
	if (!entries.size) throw new Error("The IgBLAST D-frame file contains no annotation records.");
	const prepared = prepareIgblastStyleGermlineFasta(fasta, "D");
	let matched = 0;
	let annotated = 0;
	const unmatched = [];
	return {
		fasta: parseFasta(prepared.fasta).map((record) => {
			const name = germlineName(record.header);
			const normalized = normalizeIndexSequence(record.rawSequence);
			const frame = entries.get(name);
			if (frame === void 0) {
				unmatched.push(name);
				return `>${name}\n${normalized.sequence}\n`;
			}
			matched += 1;
			if (frame < 0) return `>${name}\n${normalized.sequence}\n`;
			annotated += 1;
			return `>${name}${metadataHeader(compactMetadata(frame, -1, EMPTY_BOUNDS, SOURCE_PROVIDED))}\n${normalized.sequence}\n`;
		}).join(""),
		matched,
		annotated,
		total: prepared.count,
		unmatched
	};
}
function preprocessGermlineFasta(text, segment, templates = [], allowedLoci = LOCI, match = {}) {
	const matchOptions = resolveGermlineMatchOptions(match);
	const inputRecords = parseFasta(text);
	const records = [];
	const diagnostics = [];
	const seen = /* @__PURE__ */ new Set();
	const loci = /* @__PURE__ */ new Set();
	const warnings = [];
	let exactImgt = 0;
	let transferred = 0;
	let motifValidated = 0;
	let ambiguousBases = 0;
	for (const input of inputRecords) {
		const name = germlineName(input.header);
		if (!name) throw new Error(`A ${segment} germline record has an empty identifier.`);
		if (seen.has(name)) throw new Error(`Duplicate ${segment} germline identifier: ${name}.`);
		seen.add(name);
		const locus = inferLocus(`${name} ${input.header}`);
		if (!locus) throw new Error(`${name} does not identify one of the supported IG/TR loci.`);
		if (!allowedLoci.includes(locus)) throw new Error(`${name} belongs to ${locus}, outside the selected ${allowedLoci.join("/")} search space.`);
		const detectedSegment = inferSegment(name);
		if (detectedSegment && detectedSegment !== segment) throw new Error(`${name} looks like a ${detectedSegment} gene, not a ${segment} gene.`);
		let normalized;
		try {
			normalized = normalizeSequence(input.rawSequence);
		} catch (error) {
			throw new Error(`${name}: ${error instanceof Error ? error.message : String(error)}`);
		}
		if (!normalized.sequence) throw new Error(`${name} has an empty nucleotide sequence.`);
		ambiguousBases += normalized.ambiguous;
		const fields = parseImgtFields(input.header);
		let metadata = metadataFromHeader(input.header);
		let status = metadata ? "retained" : "unresolved";
		let transfer;
		if (metadata && !validateMetadata(metadata, normalized.sequence.length, segment)) throw new Error(`${name} contains invalid SWIGMETA coordinates.`);
		if (metadata) metadata = [...metadata.slice(0, 12), metadata[12] || SOURCE_PROVIDED];
		if (!metadata && segment === "V") {
			metadata = imgtVMetadata(input.rawSequence, fields);
			if (metadata) status = "imgt";
		}
		if (!metadata && segment === "V" && templates.length) {
			transfer = transferMetadata(name, normalized.sequence, templates, matchOptions);
			metadata = transfer.metadata;
			if (metadata) status = "transferred";
		}
		if (!metadata && segment === "J" && templates.length) {
			transfer = transferJMetadata(name, normalized.sequence, templates, matchOptions);
			metadata = transfer.metadata;
			if (metadata) status = "transferred";
		}
		if (!metadata && segment === "J") {
			metadata = jMetadata(normalized.sequence, fields);
			if (metadata) status = "motif";
		}
		if (!metadata && segment === "D") {
			const frame = imgtFrame(fields);
			if (frame >= 0) {
				metadata = compactMetadata(frame, -1, EMPTY_BOUNDS, SOURCE_IMGT_GAPPED);
				status = "imgt";
			}
		}
		if (!metadata && (segment === "D" || segment === "C")) status = "normalized";
		if (metadata?.[12] === SOURCE_IMGT_GAPPED) exactImgt += 1;
		if (metadata?.[12] === SOURCE_TRANSFERRED_IMGT || metadata?.[12] === SOURCE_TRANSFERRED_J) transferred += 1;
		if (metadata?.[12] === SOURCE_VALIDATED_J_MOTIF) motifValidated += 1;
		if (!annotationPresent(segment, metadata) && warnings.length < 8) warnings.push(`${name}: ${segment === "V" ? "no validated IMGT FWR/CDR delineation" : segment === "J" ? "no unique frame-consistent F/W-G J motif" : "no segment metadata"}.`);
		loci.add(locus);
		records.push({
			header: input.header,
			name,
			rawSequence: input.rawSequence,
			sequence: normalized.sequence,
			locus,
			metadata
		});
		if (matchOptions.includeDiagnostics) diagnostics.push({
			segment,
			name,
			locus,
			status,
			source: metadataSource(metadata),
			template: transfer?.template,
			identity: transfer?.identity,
			matchKind: transfer?.matchKind,
			taxonomicTier: transfer?.metadata ? match.taxonomicTier : void 0,
			attemptedCandidates: transfer?.attemptedCandidates ?? 0,
			bestCandidate: transfer?.bestCandidate,
			bestIdentity: transfer?.bestIdentity,
			rejectionCounts: transfer && Object.keys(transfer.rejectionCounts).length ? transfer.rejectionCounts : void 0
		});
	}
	const annotated = records.filter((record) => annotationPresent(segment, record.metadata)).length;
	return {
		fasta: records.map((record) => `>${record.name}${metadataHeader(record.metadata)}\n${record.sequence}\n`).join(""),
		count: records.length,
		annotated,
		unannotated: records.length - annotated,
		exactImgt,
		transferred,
		motifValidated,
		ambiguousBases,
		loci: [...loci].sort(),
		warnings,
		diagnostics: matchOptions.includeDiagnostics ? diagnostics : void 0
	};
}
function mergeDiagnosticHistory(previous, current) {
	if (!previous) return current;
	if (current.status === "retained" && previous.status !== "unresolved") return previous;
	const rejectionCounts = { ...previous.rejectionCounts ?? {} };
	for (const [reason, count] of Object.entries(current.rejectionCounts ?? {})) rejectionCounts[reason] = (rejectionCounts[reason] ?? 0) + count;
	const previousBest = previous.bestIdentity ?? -1;
	const currentBest = current.bestIdentity ?? -1;
	return {
		...current.status !== "unresolved" ? current : previous,
		attemptedCandidates: previous.attemptedCandidates + current.attemptedCandidates,
		bestCandidate: currentBest > previousBest ? current.bestCandidate : previous.bestCandidate,
		bestIdentity: Math.max(previousBest, currentBest) >= 0 ? Math.max(previousBest, currentBest) : void 0,
		rejectionCounts: Object.keys(rejectionCounts).length ? rejectionCounts : void 0
	};
}
/**
* Apply the same progressively broadened template search used by the browser
* worker. Existing valid SWIGMETA is retained, while records still lacking V
* or J metadata are offered to each successive template tier.
*/
function preprocessGermlineFastaAcrossTiers(text, segment, templateTiers, allowedLoci = LOCI, match = {}) {
	let report;
	const diagnosticHistory = /* @__PURE__ */ new Map();
	const tiers = templateTiers.length ? templateTiers : [[]];
	for (let tier = 0; tier < tiers.length; tier += 1) {
		report = preprocessGermlineFasta(report?.fasta ?? text, segment, [...tiers[tier]], allowedLoci, {
			...match,
			includeDiagnostics: Boolean(match.includeDiagnostics),
			taxonomicTier: tier
		});
		for (const diagnostic of report.diagnostics ?? []) diagnosticHistory.set(diagnostic.name, mergeDiagnosticHistory(diagnosticHistory.get(diagnostic.name), diagnostic));
		if (segment !== "V" && segment !== "J") break;
		if (report.annotated === report.count) break;
	}
	if (match.includeDiagnostics) report.diagnostics = [...diagnosticHistory.values()];
	return report;
}
function annotationCoverage(fasta, segment) {
	const records = parseFasta(fasta);
	let annotated = 0;
	for (const record of records) {
		const normalized = normalizeSequence(record.rawSequence);
		const metadata = metadataFromHeader(record.header);
		if (metadata && validateMetadata(metadata, normalized.sequence.length, segment) && annotationPresent(segment, metadata)) annotated += 1;
	}
	return {
		annotated,
		total: records.length
	};
}
function alleleMetadataHeader(allele) {
	return `>${allele[0]}${metadataHeader(allele[2])}\n${allele[1]}\n`;
}
function metadataSource(metadata) {
	return metadata ? METADATA_SOURCE_LABELS[metadata[12]] ?? "unclassified annotation" : "none";
}
//#endregion
//#region src/post-analysis-record.ts
/** Shared browser-worker/CLI conversion into the compact scientific record. */
function airrRowToPostAnalysisRecord(row, intern = (value) => value) {
	const value = (name) => row.values[name] ?? "";
	const sequence = normalizeNt(value("sequence"));
	return {
		ordinal: row.ordinal,
		sequenceId: value("sequence_id"),
		datasetId: intern(value("swig_dataset_id")),
		sampleId: intern(value("sample_id")),
		subjectId: intern(value("subject_id")),
		cohort: intern(value("swig_cohort")),
		timepoint: intern(value("swig_timepoint")),
		compartment: intern(value("swig_compartment")),
		locus: intern(value("locus")),
		vCall: intern(value("v_call")),
		jCall: intern(value("j_call")),
		originalVCall: intern(value("v_call")),
		originalJCall: intern(value("j_call")),
		cCall: intern(value("c_call")),
		cdr3Nt: normalizeNt(value("cdr3")),
		cdr3Aa: value("cdr3_aa").toUpperCase().replace(/[^A-Z*]/g, ""),
		productive: /^(?:T|TRUE)$/i.test(value("productive")),
		sequenceFingerprint: sequenceFingerprint(sequence),
		trimmedFingerprint: sequenceFingerprint(normalizeNt(value("sequence_alignment") || value("sequence"))),
		inputCount: Math.max(1, Math.floor(Number(value("duplicate_count")) || 1))
	};
}
/** Exact V-to-J slice used by every denoising front end. AIRR coordinates are one-based. */
function denoiseVdjSequence(row) {
	const raw = row.values.sequence ?? "";
	const start = Math.floor(Number(row.values.v_sequence_start));
	const end = Math.floor(Number(row.values.j_sequence_end));
	if (raw && Number.isFinite(start) && Number.isFinite(end) && start >= 1 && end >= start && end <= raw.length) return raw.slice(start - 1, end);
	return row.values.sequence_alignment || raw;
}
//#endregion
//#region src/reference-pack.ts
const BCR_LOCI = [
	"IGH",
	"IGK",
	"IGL"
];
const TCR_LOCI = [
	"TRA",
	"TRB",
	"TRD",
	"TRG"
];
function lociForScope(species, scope) {
	return (scope === "BCR" ? BCR_LOCI : scope === "TCR" ? TCR_LOCI : [scope]).filter((locus) => species.loci[locus]);
}
function allelesForScope(species, scope, segment) {
	return lociForScope(species, scope).flatMap((locus) => species.loci[locus]?.[segment] ?? []);
}
/** Exact species-first taxonomic template tiers used for custom V/J metadata transfer. */
function germlineTemplateTiers(pack, selected, scope, segment) {
	const baseTaxon = selected.name.split("_", 1)[0];
	const genus = baseTaxon.split(" ", 1)[0];
	const groups = [
		[selected],
		pack.species.filter((entry) => entry.name !== selected.name && entry.name.split("_", 1)[0] === baseTaxon),
		pack.species.filter((entry) => entry.name.split("_", 1)[0] !== baseTaxon && entry.name.startsWith(`${genus} `)),
		pack.species.filter((entry) => !entry.name.startsWith(`${genus} `))
	];
	const seen = /* @__PURE__ */ new Set();
	return groups.map((group) => group.flatMap((entry) => allelesForScope(entry, scope, segment)).filter((allele) => {
		const key = `${allele[0]}\u0000${allele[1]}`;
		if (seen.has(key)) return false;
		seen.add(key);
		return true;
	})).filter((group) => group.length);
}
function allelesToFasta(alleles) {
	return alleles.map(alleleMetadataHeader).join("");
}
function compileReferences(species, scope, overrides = {}) {
	const loci = lociForScope(species, scope);
	const segments = {
		V: allelesForScope(species, scope, "V"),
		D: allelesForScope(species, scope, "D"),
		J: allelesForScope(species, scope, "J"),
		C: allelesForScope(species, scope, "C")
	};
	return {
		V: overrides.V ?? allelesToFasta(segments.V),
		D: overrides.D ?? allelesToFasta(segments.D),
		J: overrides.J ?? allelesToFasta(segments.J),
		C: overrides.C ?? allelesToFasta(segments.C),
		counts: {
			V: overrides.V ? countFastaRecords(overrides.V) : segments.V.length,
			D: overrides.D ? countFastaRecords(overrides.D) : segments.D.length,
			J: overrides.J ? countFastaRecords(overrides.J) : segments.J.length,
			C: overrides.C ? countFastaRecords(overrides.C) : segments.C.length
		},
		annotation: {
			V: overrides.V ? annotationCoverage(overrides.V, "V") : {
				annotated: segments.V.filter((allele) => allele[2]?.slice(2, 12).every((value) => value >= 0)).length,
				total: segments.V.length
			},
			J: overrides.J ? annotationCoverage(overrides.J, "J") : {
				annotated: segments.J.filter((allele) => Boolean(allele[2] && allele[2][0] >= 0 && allele[2][1] >= 0)).length,
				total: segments.J.length
			}
		},
		loci
	};
}
function countFastaRecords(value) {
	return value.split(/\r?\n/).filter((line) => line.startsWith(">")).length;
}
//#endregion
//#region src/shm-analysis.ts
const CODON = {
	TTT: "F",
	TTC: "F",
	TTA: "L",
	TTG: "L",
	TCT: "S",
	TCC: "S",
	TCA: "S",
	TCG: "S",
	TAT: "Y",
	TAC: "Y",
	TAA: "*",
	TAG: "*",
	TGT: "C",
	TGC: "C",
	TGA: "*",
	TGG: "W",
	CTT: "L",
	CTC: "L",
	CTA: "L",
	CTG: "L",
	CCT: "P",
	CCC: "P",
	CCA: "P",
	CCG: "P",
	CAT: "H",
	CAC: "H",
	CAA: "Q",
	CAG: "Q",
	CGT: "R",
	CGC: "R",
	CGA: "R",
	CGG: "R",
	ATT: "I",
	ATC: "I",
	ATA: "I",
	ATG: "M",
	ACT: "T",
	ACC: "T",
	ACA: "T",
	ACG: "T",
	AAT: "N",
	AAC: "N",
	AAA: "K",
	AAG: "K",
	AGT: "S",
	AGC: "S",
	AGA: "R",
	AGG: "R",
	GTT: "V",
	GTC: "V",
	GTA: "V",
	GTG: "V",
	GCT: "A",
	GCC: "A",
	GCA: "A",
	GCG: "A",
	GAT: "D",
	GAC: "D",
	GAA: "E",
	GAG: "E",
	GGT: "G",
	GGC: "G",
	GGA: "G",
	GGG: "G"
};
function topCall(value) {
	return value.split(",")[0]?.trim() ?? "";
}
function numeric(value) {
	const result = Number(value);
	return Number.isFinite(result) ? result : 0;
}
function positiveCount(value) {
	const result = Math.round(numeric(value));
	return result > 0 ? result : 1;
}
function valid(base) {
	return /^[ACGT]$/.test(base);
}
function rangeValue(row, name) {
	return Math.max(0, Math.floor(numeric(row[name])));
}
function inRegion(queryPosition, row, prefix) {
	const start = rangeValue(row, `${prefix}_start`);
	const end = rangeValue(row, `${prefix}_end`);
	return start > 0 && end >= start && queryPosition >= start && queryPosition <= end;
}
function computeShmMetric(row, ordinal = 0, lineageId = 0, stratum = "all") {
	const query = (row.v_sequence_alignment || row.sequence_alignment || "").toUpperCase().replace(/\./g, "-");
	const germline = (row.v_germline_alignment || row.germline_alignment || "").toUpperCase().replace(/\./g, "-");
	if (!query || !germline) return null;
	const length = Math.min(query.length, germline.length);
	const queryStart = Math.max(1, rangeValue(row, "v_sequence_start") || 1);
	const frame = Math.max(1, rangeValue(row, "sequence_frame") || rangeValue(row, "v_frame") || 1);
	let queryPosition = queryStart - 1;
	let comparedNt = 0, vNtMutations = 0, cdrNtCompared = 0, cdrNtMutations = 0, frameworkNtCompared = 0, frameworkNtMutations = 0;
	const codonQuery = [], codonGermline = [];
	let synonymous = 0, vAaReplacements = 0, comparedCodons = 0;
	const flushCodon = () => {
		if (codonQuery.length !== 3 || codonGermline.length !== 3) return;
		const q = codonQuery.join(""), g = codonGermline.join("");
		if (!/^[ACGT]{3}$/.test(q) || !/^[ACGT]{3}$/.test(g)) return;
		comparedCodons += 1;
		if (q !== g) if (CODON[q] === CODON[g]) synonymous += 1;
		else vAaReplacements += 1;
	};
	let codingBases = (frame - 1) % 3;
	for (let column = 0; column < length; column += 1) {
		const q = query[column], g = germline[column];
		if (q !== "-") queryPosition += 1;
		if (valid(q) && valid(g)) {
			comparedNt += 1;
			const mismatch = q !== g;
			if (mismatch) vNtMutations += 1;
			const cdr = inRegion(queryPosition, row, "cdr1") || inRegion(queryPosition, row, "cdr2");
			const framework = inRegion(queryPosition, row, "fwr1") || inRegion(queryPosition, row, "fwr2") || inRegion(queryPosition, row, "fwr3");
			if (cdr) {
				cdrNtCompared += 1;
				if (mismatch) cdrNtMutations += 1;
			}
			if (framework) {
				frameworkNtCompared += 1;
				if (mismatch) frameworkNtMutations += 1;
			}
		}
		if (q !== "-" && g !== "-") {
			if (codingBases < 0) codingBases = 0;
			codonQuery.push(q);
			codonGermline.push(g);
			codingBases += 1;
			if (codingBases % 3 === 0) {
				flushCodon();
				codonQuery.length = 0;
				codonGermline.length = 0;
			}
		} else if (q !== "-" || g !== "-") {
			codonQuery.length = 0;
			codonGermline.length = 0;
			codingBases += q !== "-" ? 1 : 0;
		}
	}
	if (!comparedNt) return null;
	return {
		ordinal,
		lineageId,
		sequenceId: row.sequence_id || `record_${ordinal + 1}`,
		vCall: topCall(row.v_call),
		jCall: topCall(row.j_call),
		sampleId: row.sample_id || "Unassigned sample",
		subjectId: row.subject_id || "",
		timepoint: row.swig_timepoint || "",
		compartment: row.swig_compartment || "",
		cdr3Nt: row.cdr3 || row.junction || "",
		cdr3Aa: row.cdr3_aa || row.junction_aa || "",
		duplicateCount: positiveCount(row.duplicate_count || row.consensus_count || "1"),
		comparedNt,
		vNtMutations,
		vNtRate: vNtMutations / comparedNt,
		comparedCodons,
		vAaReplacements,
		vAaRate: comparedCodons ? vAaReplacements / comparedCodons : 0,
		synonymous,
		cdrNtCompared,
		cdrNtMutations,
		cdrNtRate: cdrNtCompared ? cdrNtMutations / cdrNtCompared : 0,
		frameworkNtCompared,
		frameworkNtMutations,
		frameworkNtRate: frameworkNtCompared ? frameworkNtMutations / frameworkNtCompared : 0,
		stratum
	};
}
function metricValue(record, metric) {
	return record[metric];
}
function weightedQuantile(values, quantile) {
	if (!values.length) return 0;
	const sorted = [...values].sort((a, b) => a.value - b.value);
	const threshold = sorted.reduce((sum, item) => sum + item.weight, 0) * Math.max(0, Math.min(1, quantile));
	let total = 0;
	for (const item of sorted) {
		total += item.weight;
		if (total >= threshold) return item.value;
	}
	return sorted.at(-1)?.value ?? 0;
}
function distributions(records, metric, key) {
	const groups = /* @__PURE__ */ new Map();
	for (const record of records) {
		const label = key(record);
		if (!label) continue;
		const group = groups.get(label) ?? [];
		group.push(record);
		groups.set(label, group);
	}
	return [...groups.entries()].map(([label, group]) => {
		const values = group.map((record) => metricValue(record, metric));
		const weights = group.map((record) => record.duplicateCount);
		const abundance = weights.reduce((a, b) => a + b, 0);
		const weighted = values.map((value, index) => ({
			value,
			weight: weights[index]
		}));
		return {
			label,
			values,
			weights,
			records: group.length,
			abundance,
			median: weightedQuantile(weighted, .5),
			mean: abundance ? weighted.reduce((sum, item) => sum + item.value * item.weight, 0) / abundance : 0,
			maximum: Math.max(...values),
			p95: weightedQuantile(weighted, .95)
		};
	}).sort((a, b) => b.abundance - a.abundance || a.label.localeCompare(b.label));
}
function lineageSampleDistributions(records, metric) {
	const distributionsByPair = distributions(records, metric, (record) => record.lineageId > 0 ? `${record.lineageId}\u0000${record.sampleId}` : "");
	const metadata = /* @__PURE__ */ new Map();
	for (const record of records) if (record.lineageId > 0) metadata.set(`${record.lineageId}\u0000${record.sampleId}`, record);
	return distributionsByPair.map((distribution) => {
		const record = metadata.get(distribution.label);
		return {
			...distribution,
			label: `Lineage ${record.lineageId} · ${record.sampleId}`,
			lineageId: record.lineageId,
			sampleId: record.sampleId,
			subjectId: record.subjectId,
			timepoint: record.timepoint,
			compartment: record.compartment
		};
	});
}
var ShmAccumulator = class {
	options;
	samples = [];
	sampledByLineageSample = /* @__PURE__ */ new Map();
	lowestByLineage = /* @__PURE__ */ new Map();
	skipped = 0;
	analyzed = 0;
	abundance = 0;
	constructor(options = {}) {
		this.options = {
			metric: options.metric ?? "vNtRate",
			maxSamplesPerLineage: Math.max(10, options.maxSamplesPerLineage ?? 2e3),
			maxGlobalSamples: Math.max(1e3, options.maxGlobalSamples ?? 1e5)
		};
	}
	add(row, ordinal, lineageId = 0, stratum = "all") {
		const metric = computeShmMetric(row, ordinal, lineageId, stratum);
		if (!metric) {
			this.skipped += 1;
			return;
		}
		this.analyzed += 1;
		this.abundance += metric.duplicateCount;
		if (lineageId > 0) {
			const previous = this.lowestByLineage.get(lineageId);
			if (!previous || metric.vNtRate < previous.vNtRate || metric.vNtRate === previous.vNtRate && (metric.vNtMutations < previous.vNtMutations || metric.vNtMutations === previous.vNtMutations && metric.ordinal < previous.ordinal)) this.lowestByLineage.set(lineageId, metric);
		}
		const pair = `${lineageId}\u0000${metric.sampleId}`;
		const pairSamples = this.sampledByLineageSample.get(pair) ?? 0;
		if (pairSamples >= this.options.maxSamplesPerLineage) return;
		if (this.samples.length < this.options.maxGlobalSamples) {
			this.samples.push(metric);
			this.sampledByLineageSample.set(pair, pairSamples + 1);
			return;
		}
		const slot = (ordinal * 2654435761 >>> 0) % this.analyzed;
		if (slot < this.options.maxGlobalSamples) {
			const removed = this.samples[slot];
			const removedPair = `${removed.lineageId}\u0000${removed.sampleId}`;
			const removedCount = this.sampledByLineageSample.get(removedPair) ?? 1;
			if (removedCount <= 1) this.sampledByLineageSample.delete(removedPair);
			else this.sampledByLineageSample.set(removedPair, removedCount - 1);
			this.samples[slot] = metric;
			this.sampledByLineageSample.set(pair, (this.sampledByLineageSample.get(pair) ?? 0) + 1);
		}
	}
	finish() {
		const records = this.samples;
		const metric = this.options.metric;
		const bins = Array.from({ length: 20 }, () => ({
			count: 0,
			abundance: 0
		}));
		for (const record of records) {
			const value = metricValue(record, metric);
			const normalized = metric.toLowerCase().includes("rate") ? Math.min(.999999, Math.max(0, value)) : Math.min(.999999, Math.max(0, value / 50));
			const bin = Math.floor(normalized * bins.length);
			bins[bin].count += 1;
			bins[bin].abundance += record.duplicateCount;
		}
		return {
			analyzedRecords: this.analyzed,
			analyzedAbundance: this.abundance,
			skippedRecords: this.skipped,
			sampledRecords: records.length,
			metric,
			records,
			lowestByLineage: [...this.lowestByLineage.values()].sort((left, right) => left.lineageId - right.lineageId).map((record) => ({
				lineageId: record.lineageId,
				ordinal: record.ordinal,
				vNtRate: record.vNtRate,
				vNtMutations: record.vNtMutations,
				cdr3Nt: record.cdr3Nt,
				cdr3Aa: record.cdr3Aa
			})),
			lineages: distributions(records, metric, (record) => record.lineageId ? `Lineage ${record.lineageId}` : "Unassigned"),
			vGenes: distributions(records, metric, (record) => record.vCall.replace(/\*.*$/, "")),
			strata: distributions(records, metric, (record) => record.stratum),
			lineageSamples: lineageSampleDistributions(records, metric),
			histogram: bins.map((bin, index) => ({
				label: `${index * 5}–${(index + 1) * 5}%`,
				...bin
			}))
		};
	}
};
//#endregion
//#region cli-src/swig-cli.mjs
const VERSION = "0.37.6";
const CLI_STREAM_HIGH_WATER_MARK = 8 * 1024 * 1024;
const CLI_GZIP_CHUNK_SIZE = 1024 * 1024;
const CLI_DIRECTORY = dirname(fileURLToPath(import.meta.url));
function defaultCliAssets() {
	const directory = join(CLI_DIRECTORY, "assets");
	return {
		wasmPath: join(directory, "swiftig.wasm"),
		referencePackPath: join(directory, "imgt-reference-pack.json.gz")
	};
}
function usage() {
	return `swig-cli ${VERSION}\n\nRun a complete non-phylogenetic Swig pipeline:\n  swig-cli run reads.fastq.gz --out swig-output\n  swig-cli run --config swig.config.json [--out DIRECTORY] [--workers N]\n\nRun only streaming V(D)J assignment (AIRR outfmt 19):\n  swig-cli --vdj -query reads.fasta -germline_db_V V.fasta -germline_db_D D.fasta \\\n    -germline_db_J J.fasta -out calls.airr.tsv\n\nPrepare custom germlines once and reuse their inferred annotations:\n  swig-cli prepare-reference -germline_db_V V.fasta -germline_db_D D.fasta \\\n    -germline_db_J J.fasta -organism human -ig_seqtype Ig --out-prefix refs/custom\n\nDisplay bundled-data attribution and license:\n  swig-cli notices\n\nCreate an editable config:\n  swig-cli init swig.config.json\n\nSingle-input metadata options:\n  --sample SAMPLE_ID  --donor SUBJECT_ID  --dataset DATASET_ID\n\nSamples with the same subjectId/--donor are treated as the same donor.\nLineage phylogenetics is intentionally not run by swig-cli.`;
}
function prepareReferenceUsage() {
	return `swig-cli ${VERSION} prepare-reference\n\nInfer, validate, and persist reusable SWIGMETA germline annotations.\n\nRequired:\n  -germline_db_V FASTA  -germline_db_J FASTA  --out-prefix PREFIX\n\nOptional references and exact metadata:\n  -germline_db_D FASTA  -c_region_db FASTA\n  -custom_internal_data FILE  -auxiliary_data FILE  -d_frame_data FILE\n  -organism NAME (default human)  -ig_seqtype Ig|TCR (default Ig)\n\nMatching controls:\n  --match-mode strict|permissive|best-guess  (default strict)\n  --best-guess             Alias for --match-mode best-guess; disables identity floors\n  --nearest-candidates N   Non-gene candidates aligned after named candidates fail\n  --v-same-gene-min-identity X  --v-nearest-min-identity X\n  --j-same-gene-min-identity X  --j-nearest-min-identity X\n  --require-complete       Exit nonzero if any V/J record remains unresolved\n\nOutputs are PREFIX.V/D/J/C.fasta, PREFIX.swig-reference.json, and\nPREFIX.annotation-diagnostics.tsv. The manifest can be passed directly to\nswig-cli --vdj with --prepared-reference.`;
}
function vdjUsage() {
	return `swig-cli ${VERSION} --vdj\n\nLow-overhead, streaming SwiftIG V(D)J assignment with IgBLAST-style option names.\n\nRequired:\n  -out AIRR_TSV, plus either --prepared-reference MANIFEST or\n  -germline_db_V FASTA and -germline_db_J FASTA\n\nInput and optional references:\n  -query FASTA            Query FASTA or '-' for stdin (default '-')\n  -germline_db_D FASTA    D germline FASTA\n  -c_region_db FASTA      Constant-region FASTA\n\nAnnotation modes (default: assignments only; CDR/FWR fields remain empty):\n  -custom_internal_data FILE  IgBLAST V .ndm.imgt data (1-based inclusive intervals)\n  -auxiliary_data FILE        IgBLAST J .aux data (0-based frame/CDR3 stop)\n  -d_frame_data FILE          IgBLAST D frame-one starts\n  --swigannots                Infer/validate metadata as in Swig Web\n\n  --prepared-reference FILE   Reuse a prepare-reference manifest and its FASTAs\n  --match-mode MODE           Metadata transfer only: strict, permissive, best-guess\n  --best-guess                Disable metadata-transfer identity floors (not read mapping)\nExecution:\n  -num_threads N          Exact worker count; --workers N overrides it\n  --workers N             Exact workers with no CLI cap; 0 chooses up to 8\n  --batch-records N       Records per bounded WASM batch; 0/omitted selects 2000, 1000,\n                          or 500 according to worker count\n  --assigner NAME         riat_mp (default), aer, aer_robust, or standard\n  -strand both|plus|minus -outfmt 19 -organism NAME -ig_seqtype Ig|TCR\n\nThe germline options take FASTA files, not makeblastdb binary prefixes. Output is SwiftIG AIRR,\nnot IgBLAST pairwise/tabular formatting. The output path is mandatory and is written incrementally.`;
}
function thirdPartyNotices() {
	return "Bundled IMGT/GENE-DB reference data\n\nSource: IMGT/GENE-DB release 202632-7, retrieved 2026-08-08.\nCopyright © 1995-2026 IMGT®, the international ImMunoGeneTics information system®.\nAttribution: IMGT®, the international ImMunoGeneTics information system®, https://www.imgt.org/, Institute of Human Genetics, Université de Montpellier and CNRS.\nLicense: CC BY 4.0, https://creativecommons.org/licenses/by/4.0/\nTerms: https://www.imgt.org/about/termsofuse.php\nCitation: Giudicelli V, Chaume D, Lefranc M-P. Nucleic Acids Research. 2005;33:D593-D597. https://doi.org/10.1093/nar/gki010\n\nSwig modifies the source data by selecting and reorganizing IG/TR V/D/J/C records, normalizing and ungapping nucleotide sequences, deriving compact coordinate metadata, selecting one source sequence per allele identifier, and joining selected coding IGH/TR constant exons. Membrane-only and untranslated constant exons are omitted. IMGT, Université de Montpellier, and CNRS do not endorse Swig or warrant the modified pack or its use.";
}
function argumentValue(args, name) {
	const index = args.indexOf(name);
	if (index >= 0) return args[index + 1];
	return args.find((value) => value.startsWith(`${name}=`))?.slice(name.length + 1);
}
function hasFlag(args, name) {
	return args.includes(name);
}
function positional(args) {
	const value = [];
	for (let i = 0; i < args.length; i += 1) {
		if (args[i].startsWith("--")) {
			if (!args[i].includes("=") && !["--help", "--version"].includes(args[i])) i += 1;
			continue;
		}
		value.push(args[i]);
	}
	return value;
}
function cleanCell(value) {
	return String(value ?? "").replace(/[\t\r\n]+/g, " ");
}
function parseIntegerOption(value, label, { minimum = 0, allowZero = true } = {}) {
	if (value === void 0 || !/^\d+$/.test(value)) throw new Error(`${label} requires an integer value.`);
	const parsed = Number(value);
	if (!Number.isSafeInteger(parsed) || parsed < minimum || !allowZero && parsed === 0) throw new Error(`${label} has an invalid value: ${value}.`);
	return parsed;
}
function parseFiniteOption(value, label) {
	const parsed = Number(value);
	if (value === void 0 || !Number.isFinite(parsed)) throw new Error(`${label} requires a numeric value.`);
	return parsed;
}
function parseVdjArguments(rawArgs, context = "vdj") {
	const aliases = new Map([
		["--query", "-query"],
		["--output", "-out"],
		["--out", "-out"],
		["--germline-db-v", "-germline_db_V"],
		["--germline-db-d", "-germline_db_D"],
		["--germline-db-j", "-germline_db_J"],
		["--c-region-db", "-c_region_db"],
		["--out-prefix", "-out"]
	]);
	const valued = new Set([
		"-query",
		"-out",
		"-germline_db_V",
		"-germline_db_D",
		"-germline_db_J",
		"-c_region_db",
		"-custom_internal_data",
		"-auxiliary_data",
		"-d_frame_data",
		"-organism",
		"-domain_system",
		"-ig_seqtype",
		"-strand",
		"-outfmt",
		"-num_threads",
		"--workers",
		"--batch-records",
		"--minimum-identity",
		"--assigner",
		"--calling-profile",
		"-min_D_match",
		"-min_J_length",
		"-num_alignments_D",
		"-num_alignments_J",
		"-D_penalty",
		"-J_penalty",
		"--prepared-reference",
		"--match-mode",
		"--nearest-candidates",
		"--v-same-gene-min-identity",
		"--v-nearest-min-identity",
		"--j-same-gene-min-identity",
		"--j-nearest-min-identity"
	]);
	const flags = new Set([
		"--swigannots",
		"-show_translation",
		"--best-guess",
		"--require-complete"
	]);
	const options = {};
	for (let index = 0; index < rawArgs.length; index += 1) {
		let token = rawArgs[index];
		if ([
			"--vdj",
			"--precompute_aux",
			"--precompute-aux"
		].includes(token)) continue;
		if ([
			"-h",
			"-help",
			"--help"
		].includes(token)) {
			options.help = true;
			continue;
		}
		if (["-version", "--version"].includes(token)) {
			options.version = true;
			continue;
		}
		const equals = token.indexOf("=");
		let inline;
		if (equals > 0) {
			inline = token.slice(equals + 1);
			token = token.slice(0, equals);
		}
		token = aliases.get(token) ?? token;
		if (flags.has(token)) {
			if (inline !== void 0) throw new Error(`${token} does not take a value.`);
			options[token] = true;
			continue;
		}
		if (!valued.has(token)) throw new Error(`Unsupported ${context === "prepare" ? "prepare-reference" : "--vdj"} option ${token}.\n\n${context === "prepare" ? prepareReferenceUsage() : vdjUsage()}`);
		const value = inline !== void 0 ? inline : rawArgs[++index];
		if (value === void 0) throw new Error(`${token} requires a value.`);
		options[token] = value;
	}
	return options;
}
function germlineMatchOptions(options, { diagnostics = false } = {}) {
	let mode = String(options["--match-mode"] ?? "strict").replaceAll("-", "_");
	if (options["--best-guess"]) {
		if (options["--match-mode"] && mode !== "best_guess") throw new Error("--best-guess conflicts with a different --match-mode value.");
		mode = "best_guess";
	}
	const optionalIdentity = (name) => options[name] === void 0 ? void 0 : parseFiniteOption(options[name], name);
	return resolveGermlineMatchOptions({
		mode,
		nearestCandidates: options["--nearest-candidates"] === void 0 ? void 0 : parseIntegerOption(options["--nearest-candidates"], "--nearest-candidates", {
			minimum: 1,
			allowZero: false
		}),
		vSameGeneMinIdentity: optionalIdentity("--v-same-gene-min-identity"),
		vNearestMinIdentity: optionalIdentity("--v-nearest-min-identity"),
		jSameGeneMinIdentity: optionalIdentity("--j-same-gene-min-identity"),
		jNearestMinIdentity: optionalIdentity("--j-nearest-min-identity"),
		includeDiagnostics: diagnostics
	});
}
function resolveFrom(base, value) {
	return isAbsolute(value) ? value : resolve(base, value);
}
function detectFormat(path, requested = "auto") {
	if (requested !== "auto") return requested;
	const plain = path.replace(/\.gz$/i, "");
	const extension = extname(plain).toLowerCase();
	if ([
		".fa",
		".fasta",
		".fna",
		".fas"
	].includes(extension)) return "fasta";
	if ([".fq", ".fastq"].includes(extension)) return "fastq";
	if ([".tsv", ".airr"].includes(extension) || plain.toLowerCase().endsWith(".airr.tsv")) return "airr";
	throw new Error(`Cannot infer the input format for ${path}; set format to fasta, fastq, or airr in the config.`);
}
function inputLines(input) {
	const range = input.gzipRange;
	if (typeof input.inline === "string") {
		if (range) throw new Error(`${input.path} cannot combine inline input with gzipRange.`);
		return createInterface({
			input: Readable.from([input.inline]),
			crlfDelay: Infinity
		});
	}
	if (input.path === "-" && range) throw new Error("Standard input cannot be combined with gzipRange.");
	if (range && (!Number.isSafeInteger(range.start) || !Number.isSafeInteger(range.end) || range.start < 0 || range.end <= range.start)) throw new Error(`${input.path} has an invalid gzipRange.`);
	const readOptions = range ? {
		start: range.start,
		end: range.end - 1,
		highWaterMark: CLI_GZIP_CHUNK_SIZE
	} : { highWaterMark: CLI_GZIP_CHUNK_SIZE };
	const raw = input.path === "-" ? process.stdin : createReadStream(input.path, readOptions);
	return createInterface({
		input: /\.gz$/i.test(input.path) || range ? raw.pipe(createGunzip({ chunkSize: CLI_GZIP_CHUNK_SIZE })) : raw,
		crlfDelay: Infinity
	});
}
async function* sequenceBatches(input, batchRecords, preprocessing, airrMode, datasetIndex, state) {
	const format = detectFormat(input.path, input.format);
	const lines = inputLines(input);
	const filter = preprocessing.fastqFilter;
	validateFastqQualityFilter(filter);
	state.fastqFilter = emptyFastqQualityFilterStats(filter.enabled, filter.enabled && format === "fastq");
	const parsed = async function* () {
		if (format === "airr") {
			for await (const record of airrRecords(lines)) {
				state.inputRecords += 1;
				state.eligibleRecords += 1;
				if (filter.enabled) {
					state.fastqFilter.recordsRetained += 1;
					state.fastqFilter.recordsPassedThrough += 1;
				}
				const delimiter = record.header.includes("	") ? "	" : ",";
				const header = delimiter === "	" ? record.header : record.header.split(delimiter).join("	");
				const row = delimiter === "	" ? record.row : record.row.split(delimiter).join("	");
				yield {
					ordinal: state.inputRecords - 1,
					text: `${row}\n`,
					header
				};
			}
			return;
		}
		if (format === "fasta") {
			for await (const text of fastaRecords(lines)) {
				state.inputRecords += 1;
				state.eligibleRecords += 1;
				if (filter.enabled) {
					state.fastqFilter.recordsRetained += 1;
					state.fastqFilter.recordsPassedThrough += 1;
				}
				yield {
					ordinal: state.inputRecords - 1,
					text
				};
			}
			return;
		}
		for await (const record of fastqRecords(lines)) {
			state.inputRecords += 1;
			const text = filter.enabled ? filterFastqRecord(record, filter, state.fastqFilter) : canonicalFastq(record);
			if (text !== null) {
				state.eligibleRecords += 1;
				yield {
					ordinal: state.inputRecords - 1,
					text
				};
			}
		}
	};
	let selected = parsed();
	if (preprocessing.subsample.enabled) {
		const size = preprocessing.subsample.size;
		const random = seededRandom(stableDatasetSeed(preprocessing.subsample.seed, datasetIndex));
		const reservoir = [];
		for await (const record of selected) if (reservoir.length < size) reservoir.push(record);
		else {
			const replacement = Math.floor(random() * state.eligibleRecords);
			if (replacement < size) reservoir[replacement] = record;
		}
		reservoir.sort((left, right) => left.ordinal - right.ordinal);
		selected = (async function* () {
			yield* reservoir;
		})();
	}
	let batch = [];
	const emit = () => {
		const header = format === "airr" ? batch[0]?.header : void 0;
		const body = batch.map((record) => record.text).join("");
		const result = format === "airr" && airrMode === "preserve" ? {
			format,
			count: batch.length,
			header,
			body
		} : {
			format,
			count: batch.length,
			text: header ? `${header}\n${body}` : body
		};
		batch = [];
		return result;
	};
	for await (const record of selected) {
		batch.push(record);
		state.selectedRecords += 1;
		if (batch.length >= batchRecords) yield emit();
	}
	if (batch.length) yield emit();
	if (!state.inputRecords) throw new Error(`No sequence records were found in ${input.path}.`);
	if (!state.eligibleRecords && state.fastqFilter.applicable) throw new Error(`The FASTQ quality filter rejected all ${state.inputRecords.toLocaleString()} reads in ${input.path}.`);
}
var WasmPool = class {
	constructor(size, init) {
		this.size = size;
		this.init = init;
		this.workers = [];
		this.available = [];
		this.queued = [];
		this.pending = /* @__PURE__ */ new Map();
		this.nextId = 1;
		this.failed = null;
	}
	async start() {
		for (let index = 0; index < this.size; index += 1) {
			const webWorker = Boolean(process.versions.bun && globalThis.Worker);
			const worker = webWorker ? new globalThis.Worker(new URL("./swig-worker.js", import.meta.url)) : new Worker(new URL("./swig-worker.mjs", import.meta.url));
			const receive = (message) => {
				const pending = this.pending.get(message.id);
				if (!pending) return;
				this.pending.delete(message.id);
				const { id, error, ...result } = message;
				if (error) pending.reject(new Error(error));
				else pending.resolve(result);
				if (pending.release) this.release(worker);
			};
			const fail = (error) => {
				const failure = error instanceof Error ? error : new Error(error?.message ?? String(error));
				this.failed = failure;
				for (const [id, pending] of this.pending) {
					this.pending.delete(id);
					pending.reject(failure);
				}
				while (this.queued.length) this.queued.shift().reject(failure);
			};
			if (webWorker) {
				worker.addEventListener("message", (event) => receive(event.data));
				worker.addEventListener("error", fail);
			} else {
				worker.on("message", receive);
				worker.on("error", fail);
			}
			this.workers.push(worker);
		}
		await Promise.all(this.workers.map((worker) => this.request(worker, {
			type: "init",
			...this.init
		})));
		this.available.push(...this.workers);
	}
	request(worker, message, release = false) {
		const id = this.nextId++;
		return new Promise((resolvePromise, reject) => {
			this.pending.set(id, {
				resolve: resolvePromise,
				reject,
				worker,
				release
			});
			worker.postMessage({
				id,
				...message
			});
		});
	}
	dispatch(worker, job) {
		const id = this.nextId++;
		this.pending.set(id, {
			resolve: job.resolve,
			reject: job.reject,
			worker,
			release: true
		});
		worker.postMessage({
			id,
			type: "annotate",
			...job.message
		});
	}
	release(worker) {
		const job = this.queued.shift();
		if (job) this.dispatch(worker, job);
		else this.available.push(worker);
	}
	run(message) {
		if (this.failed) return Promise.reject(this.failed);
		return new Promise((resolvePromise, reject) => {
			const job = {
				message,
				resolve: resolvePromise,
				reject
			};
			const worker = this.available.shift();
			if (worker) this.dispatch(worker, job);
			else this.queued.push(job);
		});
	}
	async close() {
		const error = /* @__PURE__ */ new Error("SwiftIG worker pool closed before queued work completed.");
		while (this.queued.length) this.queued.shift().reject(error);
		for (const worker of this.workers) await worker.terminate();
		this.workers = [];
		this.available = [];
	}
};
async function runPostAnalysisTasks(tasks, requestedWorkers, onCompleted) {
	if (!tasks.length) return [];
	const workerCount = Math.max(1, Math.min(Math.floor(requestedWorkers) || 1, tasks.length));
	if (workerCount === 1) return tasks.map((task) => {
		const result = task.kind === "denoise" ? runDenoisePartitionJob(task.job) : runExactDedupJob(task.job);
		onCompleted?.(result);
		return result;
	});
	return new Promise((resolvePromise, reject) => {
		const workers = [];
		const results = new Array(tasks.length);
		let next = 0, completed = 0, settled = false;
		const close = () => {
			for (const worker of workers) worker.terminate();
		};
		const fail = (error) => {
			if (settled) return;
			settled = true;
			close();
			reject(error instanceof Error ? error : new Error(error?.message ?? String(error)));
		};
		const dispatch = (worker) => {
			if (next >= tasks.length) return;
			const id = next++;
			worker.postMessage({
				id,
				...tasks[id]
			});
		};
		for (let index = 0; index < workerCount; index += 1) {
			const webWorker = Boolean(process.versions.bun && globalThis.Worker);
			const worker = webWorker ? new globalThis.Worker(new URL("./post-analysis-worker.js", import.meta.url)) : new Worker(new URL("./post-analysis-worker.mjs", import.meta.url));
			const receive = (message) => {
				if (settled) return;
				if (message.error || !message.result) {
					fail(new Error(message.error || "A post-analysis worker returned no result."));
					return;
				}
				results[message.id] = message.result;
				onCompleted?.(message.result);
				completed += 1;
				if (completed === tasks.length) {
					settled = true;
					close();
					resolvePromise(results);
				} else dispatch(worker);
			};
			if (webWorker) {
				worker.addEventListener("message", (event) => receive(event.data));
				worker.addEventListener("error", fail);
			} else {
				worker.on("message", receive);
				worker.on("error", fail);
			}
			workers.push(worker);
			dispatch(worker);
		}
	});
}
var ChmmPool = class {
	constructor(size, init) {
		this.size = size;
		this.init = init;
		this.workers = [];
		this.available = [];
		this.queued = [];
		this.pending = /* @__PURE__ */ new Map();
		this.nextId = 1;
		this.failed = null;
	}
	async start() {
		for (let index = 0; index < this.size; index += 1) {
			const webWorker = Boolean(process.versions.bun && globalThis.Worker);
			const worker = webWorker ? new globalThis.Worker(new URL("./chmmairra-worker.js", import.meta.url)) : new Worker(new URL("./chmmairra-worker.mjs", import.meta.url));
			const receive = (message) => {
				const pending = this.pending.get(message.id);
				if (!pending) return;
				this.pending.delete(message.id);
				if (message.error) pending.reject(new Error(message.error));
				else pending.resolve(message.result);
				if (pending.release) this.release(worker);
			};
			const fail = (error) => {
				const failure = error instanceof Error ? error : new Error(error?.message ?? String(error));
				this.failed = failure;
				for (const [id, pending] of this.pending) {
					this.pending.delete(id);
					pending.reject(failure);
				}
				while (this.queued.length) this.queued.shift().reject(failure);
			};
			if (webWorker) {
				worker.addEventListener("message", (event) => receive(event.data));
				worker.addEventListener("error", fail);
			} else {
				worker.on("message", receive);
				worker.on("error", fail);
			}
			this.workers.push(worker);
		}
		await Promise.all(this.workers.map((worker) => this.request(worker, {
			type: "init",
			...this.init
		})));
		this.available.push(...this.workers);
	}
	request(worker, message, release = false) {
		const id = this.nextId++;
		return new Promise((resolvePromise, reject) => {
			this.pending.set(id, {
				resolve: resolvePromise,
				reject,
				release
			});
			worker.postMessage({
				id,
				...message
			});
		});
	}
	dispatch(worker, job) {
		const id = this.nextId++;
		this.pending.set(id, {
			resolve: job.resolve,
			reject: job.reject,
			release: true
		});
		worker.postMessage({
			id,
			type: "batch",
			rows: job.rows
		});
	}
	release(worker) {
		const job = this.queued.shift();
		if (job) this.dispatch(worker, job);
		else this.available.push(worker);
	}
	run(rows) {
		if (this.failed) return Promise.reject(this.failed);
		return new Promise((resolvePromise, reject) => {
			const job = {
				rows,
				resolve: resolvePromise,
				reject
			};
			const worker = this.available.shift();
			if (worker) this.dispatch(worker, job);
			else this.queued.push(job);
		});
	}
	async close() {
		const error = /* @__PURE__ */ new Error("CHMMAIRRa worker pool closed before queued work completed.");
		while (this.queued.length) this.queued.shift().reject(error);
		await Promise.all(this.workers.map((worker) => Promise.resolve(worker.terminate())));
		this.workers = [];
		this.available = [];
	}
};
function workerInitialization(wasmPath, references, callingProfile, assignerStrategy, tuning) {
	return {
		wasmPath,
		referenceV: references.V,
		referenceD: references.D,
		referenceJ: references.J,
		referenceC: references.C,
		callingProfile,
		assignerStrategy,
		hasTuning: Boolean(tuning),
		tuningDMatch: tuning?.dMatch ?? 0,
		tuningDMismatch: tuning?.dMismatch ?? 0,
		tuningDGapOpen: tuning?.dGapOpen ?? 0,
		tuningDGapExtend: tuning?.dGapExtend ?? 0,
		tuningTopD: tuning?.topD ?? 0,
		tuningMinDMatch: tuning?.minDMatch ?? 0,
		tuningJMatch: tuning?.jMatch ?? 0,
		tuningJMismatch: tuning?.jMismatch ?? 0,
		tuningJGapOpen: tuning?.jGapOpen ?? 0,
		tuningJGapExtend: tuning?.jGapExtend ?? 0,
		tuningTopJ: tuning?.topJ ?? 0,
		tuningMinJLength: tuning?.minJLength ?? 0
	};
}
function workerAnnotation(text, count, format, minimumIdentity, strand, doubleD) {
	return {
		text,
		count,
		format,
		minimumIdentity,
		strand,
		doubleDMode: doubleD.mode,
		doubleDMinimumVjSpan: doubleD.minimumVjSpan ?? 0,
		doubleDSeedLength: doubleD.seedLength ?? 0,
		doubleDPseudoTrim: doubleD.pseudoTrim ?? 0,
		doubleDMaximumPseudoMismatches: doubleD.maximumPseudoMismatches ?? 0,
		doubleDMinimumScoreGain: doubleD.minimumScoreGain ?? 0
	};
}
function automaticBatchRecords(workers) {
	return workers <= 2 ? 2e3 : workers <= 4 ? 1e3 : 500;
}
function parseTable(headerText, bodyText) {
	const headers = headerText.replace(/\r$/, "").split("	");
	const rows = [];
	for (const line of bodyText.split(/\r?\n/)) {
		if (!line) continue;
		const fields = line.split("	");
		const values = {};
		headers.forEach((header, index) => {
			values[header] = fields[index] ?? "";
		});
		rows.push(values);
	}
	return {
		headers,
		rows
	};
}
function trimAuxiliaryFwr4(headerText, bodyText, endOffsets, jLengths) {
	if (!Object.keys(endOffsets).length) return bodyText;
	const headers = headerText.split("	");
	const at = new Map(headers.map((header, index) => [header, index]));
	if ([
		"j_call",
		"j_sequence_start",
		"j_germline_start",
		"j_sequence_alignment",
		"j_germline_alignment",
		"fwr4",
		"fwr4_aa",
		"fwr4_end"
	].some((header) => !at.has(header))) return bodyText;
	const lines = [];
	for (const line of bodyText.split(/\r?\n/)) {
		if (!line) continue;
		const values = line.split("	");
		const call = (values[at.get("j_call")] ?? "").split(",", 1)[0];
		const offset = endOffsets[call], referenceLength = jLengths.get(call);
		if (offset === void 0 || !referenceLength || offset === 0) {
			lines.push(line);
			continue;
		}
		const queryStart = Number(values[at.get("j_sequence_start")]);
		const referenceStart = Number(values[at.get("j_germline_start")]);
		const queryAlignment = values[at.get("j_sequence_alignment")] ?? "";
		const germlineAlignment = values[at.get("j_germline_alignment")] ?? "";
		const currentEnd = Number(values[at.get("fwr4_end")]);
		const targetReferenceEnd = referenceLength - offset;
		if (!Number.isFinite(queryStart) || !Number.isFinite(referenceStart) || !Number.isFinite(currentEnd) || targetReferenceEnd < referenceStart - 1) {
			lines.push(line);
			continue;
		}
		let queryPosition = queryStart - 1, referencePosition = referenceStart - 1;
		const columns = Math.min(queryAlignment.length, germlineAlignment.length);
		for (let column = 0; column < columns; column += 1) {
			if (referencePosition >= targetReferenceEnd) break;
			if (queryAlignment[column] !== "-") queryPosition += 1;
			if (germlineAlignment[column] !== "-") referencePosition += 1;
		}
		if (referencePosition < targetReferenceEnd || queryPosition >= currentEnd) {
			lines.push(line);
			continue;
		}
		const trim = currentEnd - queryPosition;
		const nucleotide = values[at.get("fwr4")] ?? "";
		const kept = Math.max(0, nucleotide.length - trim);
		values[at.get("fwr4")] = nucleotide.slice(0, kept);
		values[at.get("fwr4_aa")] = (values[at.get("fwr4_aa")] ?? "").slice(0, Math.floor(kept / 3));
		values[at.get("fwr4_end")] = String(queryPosition);
		lines.push(values.join("	"));
	}
	return lines.length ? `${lines.join("\n")}\n` : "";
}
function serializeRowBody(headers, rows) {
	return `${rows.map((row) => headers.map((header) => cleanCell(row[header])).join("	")).join("\n")}${rows.length ? "\n" : ""}`;
}
function addHeaders(headers, names) {
	for (const name of names) if (!headers.includes(name)) headers.push(name);
}
async function loadReferences(config, base, referencePackPath) {
	const packed = JSON.parse(gunzipSync(readFileSync(referencePackPath)).toString("utf8"));
	const species = packed.species.find((entry) => entry.name === config.references.species);
	if (!species) throw new Error(`The bundled reference pack has no exact species named ${config.references.species}.`);
	const overrides = { ...config.references.inline };
	const preparation = {};
	const allowedLoci = lociForScope(species, config.references.scope);
	if (!allowedLoci.length) throw new Error(`The bundled reference pack has no ${config.references.scope} loci for ${config.references.species}.`);
	for (const segment of [
		"V",
		"D",
		"J",
		"C"
	]) {
		const path = config.references.files?.[segment];
		if (!path) continue;
		const raw = (await readFile(resolveFrom(base, path))).toString("utf8");
		if (config.references.prepareMetadata) {
			process.stderr.write(`Preparing custom ${segment} reference metadata from ${basename(path)}…\n`);
			const report = preprocessGermlineFastaAcrossTiers(raw, segment, germlineTemplateTiers(packed, species, config.references.scope, segment), allowedLoci);
			overrides[segment] = report.fasta;
			const { fasta, ...summary } = report;
			preparation[segment] = {
				file: path,
				...summary
			};
			const detail = segment === "V" ? `${report.annotated.toLocaleString()} with validated FWR/CDR metadata` : segment === "J" ? `${report.annotated.toLocaleString()} with validated frame/CDR3-anchor metadata` : "validated and normalized";
			process.stderr.write(`Prepared ${report.count.toLocaleString()} ${segment} record${report.count === 1 ? "" : "s"}; ${detail}.\n`);
			for (const warning of report.warnings) process.stderr.write(`Reference warning: ${warning}\n`);
		} else overrides[segment] = raw;
	}
	const references = compileReferences(species, config.references.scope, overrides);
	if (!references.V.trim() || !references.J.trim()) throw new Error("The selected reference composition must contain V and J records.");
	return {
		references,
		preparation
	};
}
async function readMaybeCompressedText(path, label) {
	if (path === "-") throw new Error(`${label} must be a FASTA/data file; standard input is reserved for -query.`);
	let bytes;
	try {
		bytes = await readFile(resolve(path));
	} catch (error) {
		if (error?.code === "ENOENT") throw new Error(`${label} was not found at ${path}. The IgBLAST-style germline options require source FASTA files, not makeblastdb prefixes.`);
		throw error;
	}
	if (bytes[0] === 31 && bytes[1] === 139) bytes = gunzipSync(bytes);
	return bytes.toString("utf8");
}
async function readVdjFasta(path, label) {
	const text = await readMaybeCompressedText(path, label);
	if (!text.trimStart().startsWith(">")) throw new Error(`${label} must point to a nucleotide FASTA file. Binary BLAST database files/prefixes are not decoded by swig-cli.`);
	return text;
}
function vdjScope(options) {
	const value = String(options["-ig_seqtype"] ?? "Ig").toUpperCase();
	if (value === "IG") return "BCR";
	if (value === "TCR") return "TCR";
	throw new Error("-ig_seqtype must be Ig or TCR.");
}
function vdjSpecies(pack, requested) {
	const aliases = {
		human: "Homo sapiens",
		mouse: "Mus musculus_C57BL/6",
		rat: "Rattus norvegicus_BN; Sprague-Dawley",
		rabbit: "Oryctolagus cuniculus",
		rhesus_monkey: "Macaca mulatta_AG07107"
	};
	const value = String(requested ?? "human");
	const target = aliases[value.toLowerCase()] ?? value;
	const species = pack.species.find((entry) => entry.name.toLowerCase() === target.toLowerCase());
	if (!species) throw new Error(`The embedded Swig reference pack has no species matching -organism ${value}. Use one of human, mouse, rat, rabbit, rhesus_monkey, or an exact pack species name.`);
	return species;
}
function fastaLengths(fasta) {
	const result = /* @__PURE__ */ new Map();
	let name = "", sequence = "";
	const finish = () => {
		if (name) result.set(name, sequence.replace(/\s/g, "").length);
	};
	for (const line of fasta.split(/\r?\n/)) if (line.startsWith(">")) {
		finish();
		name = line.slice(1).trim().split(/\s+/, 1)[0];
		sequence = "";
	} else if (name) sequence += line.trim();
	finish();
	return result;
}
function unmatchedMessage(kind, application) {
	const preview = application.unmatched.slice(0, 8).join(", ");
	return `${kind} did not match ${application.unmatched.length.toLocaleString()} selected germline identifier${application.unmatched.length === 1 ? "" : "s"}${preview ? `: ${preview}${application.unmatched.length > 8 ? ", …" : ""}` : ""}.`;
}
function sha256(value) {
	return createHash("sha256").update(value).digest("hex");
}
function diagnosticTsv(diagnostics) {
	const header = [
		"segment",
		"name",
		"locus",
		"status",
		"annotation_source",
		"taxonomic_tier",
		"match_kind",
		"template",
		"identity",
		"attempted_candidates",
		"best_candidate",
		"best_identity",
		"rejections"
	];
	const rows = diagnostics.map((item) => [
		item.segment,
		item.name,
		item.locus,
		item.status,
		item.source,
		item.taxonomicTier ?? "",
		item.matchKind ?? "",
		item.template ?? "",
		item.identity === void 0 ? "" : item.identity.toFixed(6),
		item.attemptedCandidates,
		item.bestCandidate ?? "",
		item.bestIdentity === void 0 ? "" : item.bestIdentity.toFixed(6),
		Object.entries(item.rejectionCounts ?? {}).sort(([left], [right]) => left.localeCompare(right)).map(([reason, count]) => `${reason}:${count}`).join(";")
	].map(cleanCell).join("	"));
	return `${header.join("	")}\n${rows.join("\n")}${rows.length ? "\n" : ""}`;
}
function referencePreparationSummary(report) {
	const statusCounts = {};
	for (const diagnostic of report.diagnostics ?? []) statusCounts[diagnostic.status] = (statusCounts[diagnostic.status] ?? 0) + 1;
	const { fasta, diagnostics, ...summary } = report;
	return {
		...summary,
		statusCounts
	};
}
function logReferenceFailures(segment, report) {
	for (const diagnostic of report.diagnostics ?? []) {
		if (diagnostic.status !== "unresolved") continue;
		const best = diagnostic.bestCandidate ? `; best candidate ${diagnostic.bestCandidate}${diagnostic.bestIdentity === void 0 ? "" : ` at ${(diagnostic.bestIdentity * 100).toFixed(1)}% identity`}` : "; no annotated candidate";
		const failures = Object.entries(diagnostic.rejectionCounts ?? {}).sort((left, right) => right[1] - left[1] || left[0].localeCompare(right[0])).map(([reason, count]) => `${reason}=${count}`).join(", ");
		process.stderr.write(`[prepare:${segment}] unresolved ${diagnostic.name}${best}${failures ? `; rejected: ${failures}` : ""}.\n`);
	}
}
async function runPrepareReference(rawArgs, assets) {
	const options = parseVdjArguments(rawArgs, "prepare");
	if (options.help) {
		process.stdout.write(`${prepareReferenceUsage()}\n`);
		return;
	}
	if (options.version) {
		process.stdout.write(`${VERSION}\n`);
		return;
	}
	const outputValue = options["-out"];
	if (!outputValue || outputValue === "-") throw new Error("prepare-reference requires --out-prefix with a filesystem path.");
	const paths = {
		V: options["-germline_db_V"],
		D: options["-germline_db_D"],
		J: options["-germline_db_J"],
		C: options["-c_region_db"]
	};
	if (!paths.V || !paths.J) throw new Error("prepare-reference requires -germline_db_V and -germline_db_J source FASTA files.");
	const packBytes = readFileSync(assets.referencePackPath);
	const pack = JSON.parse(gunzipSync(packBytes).toString("utf8"));
	const species = vdjSpecies(pack, options["-organism"]);
	const scope = vdjScope(options);
	const allowedLoci = lociForScope(species, scope);
	if (!allowedLoci.length) throw new Error(`The embedded reference pack has no ${scope} loci for ${species.name}.`);
	const match = germlineMatchOptions(options, { diagnostics: true });
	const outputPrefix = resolve(String(outputValue));
	await mkdir(dirname(outputPrefix), { recursive: true });
	process.stderr.write(`[prepare] ${species.name} ${scope}; loci ${allowedLoci.join(",")}; ${match.mode.replaceAll("_", "-")} matching; ${match.nearestCandidates} nearest candidates.\n`);
	process.stderr.write(`[prepare] V identity floors same-gene=${match.vSameGeneMinIdentity.toFixed(3)}, nearest=${match.vNearestMinIdentity.toFixed(3)}; J same-gene=${match.jSameGeneMinIdentity.toFixed(3)}, nearest=${match.jNearestMinIdentity.toFixed(3)}.\n`);
	if (match.mode === "best_guess") process.stderr.write("[prepare] Best-guess mode disables identity rejection; coordinate, frame, and conserved-anchor validation remain mandatory.\n");
	const exactData = {
		V: options["-custom_internal_data"] ? await readMaybeCompressedText(options["-custom_internal_data"], "-custom_internal_data") : null,
		J: options["-auxiliary_data"] ? await readMaybeCompressedText(options["-auxiliary_data"], "-auxiliary_data") : null,
		D: options["-d_frame_data"] ? await readMaybeCompressedText(options["-d_frame_data"], "-d_frame_data") : null
	};
	if (exactData.D && !paths.D) throw new Error("-d_frame_data requires -germline_db_D.");
	const references = {
		V: "",
		D: "",
		J: "",
		C: ""
	};
	const files = {};
	const segments = {};
	const allDiagnostics = [];
	let fwr4EndOffsets = {};
	for (const segment of [
		"V",
		"D",
		"J",
		"C"
	]) {
		const path = paths[segment];
		if (!path) continue;
		const raw = await readVdjFasta(path, segment === "C" ? "-c_region_db" : `-germline_db_${segment}`);
		let input = raw;
		process.stderr.write(`[prepare:${segment}] Reading and validating ${basename(path)}.\n`);
		if (segment === "V" && exactData.V) {
			const application = applyIgblastInternalData(input, exactData.V);
			input = application.fasta;
			process.stderr.write(`[prepare:V] IgBLAST internal data matched ${application.matched.toLocaleString()}/${application.total.toLocaleString()} records; unmatched records continue to homology transfer.\n`);
		}
		if (segment === "J" && exactData.J) {
			const application = applyIgblastAuxiliaryData(input, exactData.J);
			input = application.fasta;
			fwr4EndOffsets = application.fwr4EndOffsets ?? {};
			process.stderr.write(`[prepare:J] IgBLAST auxiliary data matched ${application.matched.toLocaleString()}/${application.total.toLocaleString()} records and annotated ${application.annotated.toLocaleString()} CDR3 stops; unmatched records continue to homology/motif inference.\n`);
		}
		if (segment === "D" && exactData.D) {
			const application = applyIgblastDFrameData(input, exactData.D);
			input = application.fasta;
			process.stderr.write(`[prepare:D] IgBLAST D-frame data matched ${application.matched.toLocaleString()}/${application.total.toLocaleString()} records.\n`);
		}
		const started = performance.now();
		const report = preprocessGermlineFastaAcrossTiers(input, segment, germlineTemplateTiers(pack, species, scope, segment), allowedLoci, match);
		const seconds = (performance.now() - started) / 1e3;
		references[segment] = report.fasta;
		allDiagnostics.push(...report.diagnostics ?? []);
		logReferenceFailures(segment, report);
		files[segment] = {
			path: basename(`${outputPrefix}.${segment}.fasta`),
			sha256: sha256(report.fasta),
			records: report.count,
			annotated: report.annotated
		};
		segments[segment] = {
			source: basename(path),
			sourceSha256: sha256(raw),
			...referencePreparationSummary(report),
			seconds: Number(seconds.toFixed(3))
		};
		process.stderr.write(`[prepare:${segment}] ${report.annotated.toLocaleString()}/${report.count.toLocaleString()} annotated; ${report.unannotated.toLocaleString()} unresolved; ${seconds.toFixed(2)} s.\n`);
	}
	const diagnosticsText = diagnosticTsv(allDiagnostics);
	const diagnosticsPath = `${outputPrefix}.annotation-diagnostics.tsv`;
	const manifestPath = `${outputPrefix}.swig-reference.json`;
	const incomplete = (segments.V?.unannotated ?? 0) + (segments.J?.unannotated ?? 0);
	const manifest = {
		schema: 1,
		application: "Swig prepared germline reference",
		applicationVersion: VERSION,
		createdAt: (/* @__PURE__ */ new Date()).toISOString(),
		complete: incomplete === 0,
		species: species.name,
		scope,
		loci: allowedLoci,
		referencePack: {
			source: pack.source,
			release: pack.release,
			retrieved: pack.retrieved,
			sha256: sha256(packBytes)
		},
		match: {
			...match,
			includeDiagnostics: void 0
		},
		files,
		diagnostics: {
			path: basename(diagnosticsPath),
			sha256: sha256(diagnosticsText),
			records: allDiagnostics.length
		},
		segments,
		fwr4EndOffsets
	};
	for (const segment of Object.keys(files)) await writeFile(join(dirname(outputPrefix), files[segment].path), references[segment]);
	await writeFile(diagnosticsPath, diagnosticsText);
	await writeFile(manifestPath, `${JSON.stringify(manifest, null, 2)}\n`);
	process.stderr.write(`[prepare] Wrote ${Object.keys(files).length} prepared FASTA${Object.keys(files).length === 1 ? "" : "s"}, ${basename(diagnosticsPath)}, and ${basename(manifestPath)}.\n`);
	process.stdout.write(`${JSON.stringify({
		manifest: manifestPath,
		diagnostics: diagnosticsPath,
		complete: manifest.complete,
		segments: Object.fromEntries(Object.entries(segments).map(([segment, value]) => [segment, {
			count: value.count,
			annotated: value.annotated,
			unannotated: value.unannotated,
			seconds: value.seconds
		}]))
	}, null, 2)}\n`);
	if (options["--require-complete"] && incomplete) throw new Error(`${incomplete.toLocaleString()} V/J germline record${incomplete === 1 ? " remains" : "s remain"} unresolved; inspect ${diagnosticsPath}.`);
}
async function loadPreparedReference(path) {
	const manifestPath = resolve(path);
	let manifest;
	try {
		manifest = JSON.parse(await readFile(manifestPath, "utf8"));
	} catch (error) {
		throw new Error(`Could not read prepared-reference manifest ${path}: ${error instanceof Error ? error.message : String(error)}`);
	}
	if (manifest?.schema !== 1 || manifest?.application !== "Swig prepared germline reference") throw new Error(`${path} is not a supported Swig prepared-reference manifest.`);
	const base = dirname(manifestPath);
	const references = {
		V: "",
		D: "",
		J: "",
		C: ""
	};
	for (const segment of [
		"V",
		"D",
		"J",
		"C"
	]) {
		const entry = manifest.files?.[segment];
		if (!entry) continue;
		const relative = typeof entry === "string" ? entry : entry.path;
		if (typeof relative !== "string" || !relative) throw new Error(`Prepared-reference manifest has no valid ${segment} FASTA path.`);
		const fasta = await readVdjFasta(resolve(base, relative), `prepared ${segment} reference`);
		const expected = typeof entry === "object" ? entry.sha256 : void 0;
		if (expected && sha256(fasta) !== expected) throw new Error(`Prepared ${segment} FASTA failed its SHA-256 check: ${resolve(base, relative)}.`);
		references[segment] = fasta;
		process.stderr.write(`[prepared-reference:${segment}] Loaded and verified ${basename(relative)}.\n`);
	}
	if (!references.V || !references.J) throw new Error("Prepared-reference manifest must provide V and J FASTAs.");
	if (manifest.complete === false) process.stderr.write("Reference warning: this prepared-reference manifest contains unresolved V/J metadata; assignments remain available but some region annotations may be blank.\n");
	return {
		references,
		mode: "prepared-reference",
		fwr4EndOffsets: manifest.fwr4EndOffsets ?? {},
		jLengths: fastaLengths(references.J)
	};
}
async function prepareVdjReferences(options, assets) {
	const preparedPath = options["--prepared-reference"];
	if (preparedPath) {
		const conflicts = [
			"-germline_db_V",
			"-germline_db_D",
			"-germline_db_J",
			"-c_region_db",
			"-custom_internal_data",
			"-auxiliary_data",
			"-d_frame_data",
			"--swigannots",
			"--match-mode",
			"--best-guess",
			"--nearest-candidates",
			"--v-same-gene-min-identity",
			"--v-nearest-min-identity",
			"--j-same-gene-min-identity",
			"--j-nearest-min-identity"
		].filter((name) => options[name]);
		if (conflicts.length) throw new Error(`--prepared-reference cannot be combined with ${conflicts.join(", ")}.`);
		return loadPreparedReference(preparedPath);
	}
	const paths = {
		V: options["-germline_db_V"],
		D: options["-germline_db_D"],
		J: options["-germline_db_J"],
		C: options["-c_region_db"]
	};
	if (!paths.V || !paths.J) throw new Error("--vdj requires -germline_db_V and -germline_db_J source FASTA files.");
	const raw = {
		V: await readVdjFasta(paths.V, "-germline_db_V"),
		D: paths.D ? await readVdjFasta(paths.D, "-germline_db_D") : "",
		J: await readVdjFasta(paths.J, "-germline_db_J"),
		C: paths.C ? await readVdjFasta(paths.C, "-c_region_db") : ""
	};
	const useSwig = Boolean(options["--swigannots"]);
	const internalPath = options["-custom_internal_data"];
	const auxiliaryPath = options["-auxiliary_data"];
	if (useSwig && (internalPath || auxiliaryPath)) throw new Error("--swigannots cannot be combined with -custom_internal_data or -auxiliary_data; choose one annotation source.");
	const scope = vdjScope(options);
	const match = germlineMatchOptions(options);
	let pack;
	let species;
	let allowedLoci;
	if (useSwig || Boolean(auxiliaryPath && !internalPath)) {
		pack = JSON.parse(gunzipSync(readFileSync(assets.referencePackPath)).toString("utf8"));
		species = vdjSpecies(pack, options["-organism"]);
		allowedLoci = lociForScope(species, scope);
		if (!allowedLoci.length) throw new Error(`The embedded reference pack has no ${scope} loci for ${species.name}.`);
	}
	const references = {
		V: "",
		D: "",
		J: "",
		C: ""
	};
	let mode = "assignments-only";
	if (useSwig) {
		mode = "swig-metadata";
		for (const segment of [
			"V",
			"D",
			"J",
			"C"
		]) {
			if (!raw[segment]) continue;
			const report = preprocessGermlineFastaAcrossTiers(raw[segment], segment, germlineTemplateTiers(pack, species, scope, segment), allowedLoci, match);
			references[segment] = report.fasta;
			process.stderr.write(`Swig metadata preparation: ${segment} ${report.annotated.toLocaleString()}/${report.count.toLocaleString()} annotated.\n`);
			for (const warning of report.warnings) process.stderr.write(`Reference warning: ${warning}\n`);
		}
	} else for (const segment of [
		"V",
		"D",
		"J",
		"C"
	]) if (raw[segment]) references[segment] = prepareIgblastStyleGermlineFasta(raw[segment], segment).fasta;
	if (internalPath) {
		mode = "igblast-data";
		const application = applyIgblastInternalData(references.V, await readMaybeCompressedText(internalPath, "-custom_internal_data"));
		if (application.matched !== application.total) throw new Error(`${unmatchedMessage("-custom_internal_data", application)} IgBLAST requires custom internal data for every selected V sequence.`);
		references.V = application.fasta;
		process.stderr.write(`IgBLAST internal data: ${application.matched.toLocaleString()}/${application.total.toLocaleString()} V records annotated.\n`);
	} else if (auxiliaryPath) {
		mode = "igblast-data";
		const report = preprocessGermlineFastaAcrossTiers(raw.V, "V", germlineTemplateTiers(pack, species, scope, "V"), allowedLoci, match);
		references.V = report.fasta;
		process.stderr.write(`Embedded ${species.name} V metadata: ${report.annotated.toLocaleString()}/${report.count.toLocaleString()} records annotated.\n`);
		for (const warning of report.warnings) process.stderr.write(`Reference warning: ${warning}\n`);
	}
	let fwr4EndOffsets = {};
	if (auxiliaryPath) {
		mode = "igblast-data";
		const application = applyIgblastAuxiliaryData(references.J, await readMaybeCompressedText(auxiliaryPath, "-auxiliary_data"));
		if (!application.matched) throw new Error("-auxiliary_data has no exact identifiers matching the selected J FASTA.");
		references.J = application.fasta;
		fwr4EndOffsets = application.fwr4EndOffsets ?? {};
		process.stderr.write(`IgBLAST auxiliary data: ${application.annotated.toLocaleString()}/${application.total.toLocaleString()} J records have CDR3-stop annotations.\n`);
		if (application.unmatched.length) process.stderr.write(`Reference warning: ${unmatchedMessage("-auxiliary_data", application)}\n`);
	}
	const dFramePath = options["-d_frame_data"];
	if (dFramePath) {
		if (!references.D) throw new Error("-d_frame_data requires -germline_db_D.");
		const application = applyIgblastDFrameData(references.D, await readMaybeCompressedText(dFramePath, "-d_frame_data"));
		if (!application.matched) throw new Error("-d_frame_data has no exact identifiers matching the selected D FASTA.");
		references.D = application.fasta;
		process.stderr.write(`IgBLAST D-frame data: ${application.annotated.toLocaleString()}/${application.total.toLocaleString()} D records annotated.\n`);
		if (application.unmatched.length) process.stderr.write(`Reference warning: ${unmatchedMessage("-d_frame_data", application)}\n`);
	}
	return {
		references,
		mode,
		fwr4EndOffsets,
		jLengths: fastaLengths(references.J)
	};
}
function compactAlleleResult(result) {
	return {
		version: result.version,
		options: result.options,
		totalRecords: result.totalRecords,
		activeRecords: result.activeRecords,
		runAt: result.runAt,
		warnings: result.warnings,
		segments: Object.fromEntries(Object.entries(result.segments).map(([segment, value]) => [segment, {
			segment: value.segment,
			models: value.models,
			modeledRows: value.modeledRows,
			changedMapRows: value.changedMapRows,
			skippedRows: value.skippedRows,
			matrixNonZeros: value.matrixNonZeros,
			truncatedRows: value.truncatedRows,
			exactDuplicateLabels: value.exactDuplicateLabels
		}]))
	};
}
function fitAlleles(rows, references, options) {
	const segments = {};
	for (const segment of options.segments) {
		const fasta = references[segment];
		if (!fasta?.trim()) continue;
		const graph = buildReferenceAlleleGraph(fasta, segment, options.neighbourRadius);
		const accumulator = new SparseEvidenceAccumulator(rows.length, graph, options);
		rows.forEach((values, ordinal) => accumulator.add(toRefinementInputRow({
			ordinal,
			values
		}, segment)));
		segments[segment] = fitSparseAlleleModel(accumulator.finish(), graph, options, rows.length);
	}
	const warnings = [];
	if (options.model === "active-set") warnings.push("The hurdle model estimates repertoire-active usage, not literal genomic presence.");
	if (options.weighting === "abundance") warnings.push("Abundance weighting allows clonal expansion to influence the fitted mixture.");
	return {
		version: 1,
		options,
		totalRecords: rows.length,
		activeRecords: rows.length,
		segments,
		runAt: (/* @__PURE__ */ new Date()).toISOString(),
		warnings
	};
}
function applyAlleles(rows, result, policy, threshold, headers) {
	for (const segment of [
		"V",
		"D",
		"J"
	]) {
		const prefix = segment.toLowerCase();
		if (!result.segments[segment]) continue;
		addHeaders(headers, [`swig_original_${prefix}_call`, `swig_repertoire_${prefix}_call`]);
		rows.forEach((row, ordinal) => {
			const call = refinedCall(result, segment, ordinal, policy, threshold);
			if (!call) return;
			row[`swig_original_${prefix}_call`] = row[`${prefix}_call`] ?? "";
			row[`swig_repertoire_${prefix}_call`] = call;
			row[`${prefix}_call`] = call;
		});
	}
}
async function runChimera(rows, activeMask, config, scope, headers, references, workers) {
	const selectedSegment = config.segment.toUpperCase();
	const msaText = config.uploadedMsa?.trim() || (config.msaSource === "selected" ? references[selectedSegment]?.trim() : "");
	if (!msaText) throw new Error("CLI chimera filtering requires either pipeline.chimera.uploadedMsa or an aligned selected-segment reference.");
	try {
		prepareReferenceMsa(msaText);
	} catch (error) {
		if (config.msaSource === "selected" && !config.uploadedMsa?.trim()) throw new Error(`The selected ${selectedSegment} references are not already aligned. Export this pipeline from Swig Web to embed its Kalign MSA, or set pipeline.chimera.uploadedMsa explicitly. ${error instanceof Error ? error.message : String(error)}`);
		throw error;
	}
	const segment = config.segment.toLowerCase();
	const options = {
		method: config.model === "auto" ? scope === "TCR" || String(scope).startsWith("TR") ? "DB" : "BW" : config.model,
		priorProbability: config.priorProbability,
		baseMutationProbability: config.baseMutationProbability,
		mutationRates: config.mutationRates,
		mutationSwitchProbability: config.mutationSwitchProbability,
		detailed: config.detailed
	};
	let evaluated = 0, flagged = 0, unevaluated = 0;
	addHeaders(headers, ["swig_chimera_probability", "swig_chimera_status"]);
	if (config.detailed) addHeaders(headers, ["swig_chimera_starting_reference", "swig_chimera_recombinations"]);
	const batches = [];
	let batch = [];
	for (let ordinal = 0; ordinal < rows.length; ordinal += 1) {
		if (!activeMask[ordinal]) continue;
		const row = rows[ordinal], call = row[`${segment}_call`] ?? "", sequence = row[`${segment}_sequence_alignment`] ?? "", germline = row[`${segment}_germline_alignment`] ?? "";
		batch.push({
			ordinal,
			call,
			sequenceAlignment: sequence,
			germlineAlignment: germline
		});
		if (batch.length >= 250) {
			batches.push(batch);
			batch = [];
		}
	}
	if (batch.length) batches.push(batch);
	if (!batches.length) return {
		evaluated,
		flagged,
		unevaluated,
		threshold: config.posteriorThreshold
	};
	const workerCount = Math.max(1, Math.min(Math.floor(workers) || 1, batches.length));
	process.stderr.write(`CHMMAIRRa: ${workerCount} worker${workerCount === 1 ? "" : "s"} across ${batches.length.toLocaleString()} row batches.\n`);
	const pool = new ChmmPool(workerCount, {
		msa: msaText,
		options,
		minDfr: config.minimumDfr
	});
	try {
		await pool.start();
		const outputs = await Promise.all(batches.map((rows) => pool.run(rows)));
		for (const output of outputs) for (const result of output.results) {
			const row = rows[result.ordinal];
			row.swig_chimera_status = result.status;
			if (result.status === "evaluated") {
				row.swig_chimera_probability = String(result.probability);
				evaluated += 1;
				if (config.detailed) {
					row.swig_chimera_starting_reference = result.startingReference;
					row.swig_chimera_recombinations = result.recombinations.map((event) => `${event.left}->${event.right}@${event.position}`).join(";");
				}
				if (result.probability >= config.posteriorThreshold) {
					activeMask[result.ordinal] = 0;
					flagged += 1;
				}
			} else {
				unevaluated += 1;
				if (!config.retainUnevaluated) activeMask[result.ordinal] = 0;
			}
		}
	} finally {
		await pool.close();
	}
	return {
		evaluated,
		flagged,
		unevaluated,
		threshold: config.posteriorThreshold
	};
}
function writeChunk(stream, chunk) {
	if (stream.write(chunk)) return Promise.resolve();
	return once(stream, "drain");
}
function createCliWriteStream(path) {
	return createWriteStream(path, { highWaterMark: CLI_STREAM_HIGH_WATER_MARK });
}
async function finishWritable(stream) {
	stream.end();
	await once(stream, "finish");
}
async function writeRowsFile(path, headers, rows, include = () => true) {
	const stream = createCliWriteStream(path);
	try {
		await once(stream, "open");
		await writeChunk(stream, `${headers.join("	")}\n`);
		let batch = [];
		for (let ordinal = 0; ordinal < rows.length; ordinal += 1) {
			if (!include(rows[ordinal], ordinal)) continue;
			batch.push(rows[ordinal]);
			if (batch.length >= 2e3) {
				await writeChunk(stream, serializeRowBody(headers, batch));
				batch = [];
			}
		}
		if (batch.length) await writeChunk(stream, serializeRowBody(headers, batch));
		await finishWritable(stream);
	} catch (error) {
		stream.destroy();
		throw error;
	}
}
async function writeLineageStudy(path, manifestPath, headers, rows, activeMask, lineages, config, references, shm) {
	const stream = createCliWriteStream(path);
	const hash = createHash("sha256");
	let offset = 0;
	const append = async (text) => {
		const chunk = Buffer.from(text, "utf8");
		hash.update(chunk);
		offset += chunk.byteLength;
		await writeChunk(stream, chunk);
	};
	await append(`${headers.join("	")}\n`);
	const retainedIds = new Set(lineages.summaries.map((summary) => summary.id));
	const indexed = rows.map((row, ordinal) => ({
		row,
		ordinal,
		lineageId: lineages.assignments[ordinal] ?? 0
	})).filter((item) => activeMask[item.ordinal] && retainedIds.has(item.lineageId)).sort((left, right) => left.lineageId - right.lineageId || left.ordinal - right.ordinal);
	const ranges = [];
	let current = 0, start = offset, count = 0, rangeHash = createHash("sha256");
	const finishRange = () => {
		if (current > 0) ranges.push({
			lineageId: current,
			start,
			end: offset,
			records: count,
			sha256: rangeHash.digest("hex")
		});
	};
	for (const item of indexed) {
		if (item.lineageId !== current) {
			finishRange();
			current = item.lineageId;
			start = offset;
			count = 0;
			rangeHash = createHash("sha256");
		}
		const text = `${headers.map((header) => cleanCell(item.row[header])).join("	")}\n`;
		rangeHash.update(text, "utf8");
		await append(text);
		count += 1;
	}
	finishRange();
	stream.end();
	await once(stream, "finish");
	const lowestShmByLineage = new Map((shm?.lowestByLineage ?? []).map((value) => [value.lineageId, value]));
	const shmSummaries = shm?.lineages.flatMap((group) => {
		const match = /^Lineage\s+(\d+)$/.exec(group.label);
		if (!match) return [];
		const lineageId = Number(match[1]), lowest = lowestShmByLineage.get(lineageId);
		return [{
			lineageId,
			mean: group.mean,
			p95: group.p95 ?? 0,
			ordinal: lowest?.ordinal,
			cdr3Nt: lowest?.cdr3Nt,
			cdr3Aa: lowest?.cdr3Aa
		}];
	}) ?? [];
	const manifest = {
		schema: 1,
		application: "Swig lineage study",
		applicationVersion: VERSION,
		createdAt: (/* @__PURE__ */ new Date()).toISOString(),
		linkedAirr: {
			name: basename(path),
			size: offset,
			records: indexed.length,
			headers,
			sha256: hash.digest("hex")
		},
		analysis: {
			inputName: config.studyName,
			species: config.references.species,
			scope: config.references.scope,
			references,
			datasets: config.inputs.map((input) => ({
				datasetId: input.datasetId,
				inputName: input.inputName || basename(input.path),
				sampleId: input.sampleId,
				subjectId: input.subjectId,
				cohort: input.cohort,
				timepoint: input.timepoint,
				compartment: input.compartment
			})),
			callingProfile: config.annotation.callingProfile,
			assignerStrategy: config.annotation.assignerStrategy,
			minimumIdentity: config.annotation.minimumIdentity,
			strand: config.annotation.strand,
			lineage: {
				scope: config.pipeline.lineage.scope,
				identity: config.pipeline.lineage.identity,
				resolution: config.pipeline.lineage.resolution,
				ambiguity: config.pipeline.lineage.ambiguity,
				productiveOnly: config.pipeline.lineage.productiveOnly,
				maxCandidateComparisons: config.pipeline.lineage.maxCandidateComparisons
			}
		},
		summaries: lineages.summaries,
		ranges,
		shm: shm ? {
			metric: shm.metric,
			summaries: shmSummaries
		} : void 0
	};
	await writeFile(manifestPath, gzipSync(JSON.stringify(manifest)));
	return manifest;
}
async function runPipeline(config, base, assets) {
	if (!config.inputs.length) throw new Error("No input datasets were specified.");
	const outputDirectory = resolveFrom(base, config.output.directory);
	await mkdir(outputDirectory, { recursive: true });
	const prefix = config.output.prefix;
	const loadedReferences = await loadReferences(config, base, assets.referencePackPath);
	const references = loadedReferences.references;
	if (config.annotation.airrMode === "preserve" && config.annotation.doubleD.mode !== "off" && config.inputs.some((input) => detectFormat(input.path, input.format) === "airr")) throw new Error("Double-D screening of AIRR input requires annotation.airrMode = \"reannotate\".");
	const needsAnnotation = config.inputs.some((input) => detectFormat(input.path, input.format) !== "airr" || config.annotation.airrMode === "reannotate");
	const annotationBatchRecords = config.annotation.batchRecords || automaticBatchRecords(config.annotation.workers);
	const pool = needsAnnotation ? new WasmPool(config.annotation.workers, workerInitialization(assets.wasmPath, references, config.annotation.callingProfile, config.annotation.assignerStrategy)) : null;
	if (pool) {
		process.stderr.write(`Starting SwiftIG pool (${config.annotation.workers} worker${config.annotation.workers === 1 ? "" : "s"}; ${annotationBatchRecords.toLocaleString()} records/batch).\n`);
		await pool.start();
	}
	const annotatedStream = config.output.writeAnnotatedAirr ? createCliWriteStream(join(outputDirectory, `${prefix}.annotated.airr.tsv`)) : null;
	let annotatedOutputHeaders = null;
	const rows = [];
	const headers = [];
	let annotatedRecords = 0;
	let inputRecords = 0;
	let eligibleRecords = 0;
	let fastqFilterStats = emptyFastqQualityFilterStats(config.preprocessing.fastqFilter.enabled, false);
	try {
		if (annotatedStream) await once(annotatedStream, "open");
		for (let datasetIndex = 0; datasetIndex < config.inputs.length; datasetIndex += 1) {
			const input = config.inputs[datasetIndex];
			const preserveAirr = detectFormat(input.path, input.format) === "airr" && config.annotation.airrMode === "preserve";
			process.stderr.write(`${preserveAirr ? "Loading existing AIRR calls from" : "Annotating"} ${input.inputName || basename(input.path)} (${input.sampleId}; donor ${input.subjectId})…\n`);
			const preprocessingState = {
				inputRecords: 0,
				eligibleRecords: 0,
				selectedRecords: 0,
				fastqFilter: emptyFastqQualityFilterStats(false, false)
			};
			const pending = [];
			const consume = async (item) => {
				const result = await item.promise;
				const manifest = {
					datasetId: input.datasetId,
					inputName: input.inputName || basename(input.path),
					sampleId: input.sampleId,
					subjectId: input.subjectId,
					cohort: input.cohort,
					timepoint: input.timepoint,
					compartment: input.compartment
				};
				const annotated = annotateAirrBatch(result.header, result.body, manifest);
				const table = parseTable(annotated.header, annotated.body);
				const doubleDAnnotated = result.doubleDHeader ? annotateDoubleDBatch(result.doubleDHeader, result.doubleDBody, manifest) : null;
				const doubleD = doubleDAnnotated ? parseTable(doubleDAnnotated.header, doubleDAnnotated.body) : null;
				const ddById = new Map((doubleD?.rows ?? []).map((row) => [row.sequence_id, row]));
				const batchHeaders = [...table.headers];
				if (doubleD) addHeaders(batchHeaders, doubleD.headers);
				for (const row of table.rows) {
					const dd = ddById.get(row.sequence_id);
					if (dd) Object.assign(row, dd);
					rows.push(row);
				}
				addHeaders(headers, batchHeaders);
				if (annotatedStream) {
					if (!annotatedOutputHeaders) {
						annotatedOutputHeaders = batchHeaders;
						await writeChunk(annotatedStream, `${annotatedOutputHeaders.join("	")}\n`);
					} else {
						const newHeaders = batchHeaders.filter((header) => !annotatedOutputHeaders.includes(header));
						if (newHeaders.length) throw new Error(`Streaming annotated AIRR output cannot add columns after its header was written (${newHeaders.join(", ")}). Use matching AIRR schemas for all preserved inputs or reannotate them.`);
					}
					await writeChunk(annotatedStream, serializeRowBody(annotatedOutputHeaders, table.rows));
				}
				annotatedRecords += table.rows.length;
			};
			for await (const batch of sequenceBatches(input, annotationBatchRecords, config.preprocessing, config.annotation.airrMode, datasetIndex, preprocessingState)) {
				const promise = batch.format === "airr" && config.annotation.airrMode === "preserve" ? Promise.resolve({
					direct: true,
					header: batch.header,
					body: batch.body,
					doubleDHeader: "",
					doubleDBody: ""
				}) : pool.run(workerAnnotation(batch.text, batch.count, batch.format === "fasta" ? 1 : batch.format === "fastq" ? 2 : 3, config.annotation.minimumIdentity, config.annotation.strand, config.annotation.doubleD));
				pending.push({
					promise,
					count: batch.count
				});
				if (pending.length >= Math.max(2, config.annotation.workers * 2)) await consume(pending.shift());
			}
			while (pending.length) await consume(pending.shift());
			inputRecords += preprocessingState.inputRecords;
			eligibleRecords += preprocessingState.eligibleRecords;
			fastqFilterStats = addFastqQualityFilterStats(fastqFilterStats, preprocessingState.fastqFilter);
		}
	} catch (error) {
		annotatedStream?.destroy();
		throw error;
	} finally {
		await pool?.close();
	}
	if (annotatedStream) await finishWritable(annotatedStream);
	rows.forEach((row, ordinal) => {
		if (!row.sequence_id) row.sequence_id = `swig_${ordinal + 1}`;
	});
	let alleleResult = null;
	if (config.pipeline.alleleRefinement.enabled) {
		process.stderr.write("Fitting repertoire-level allele model…\n");
		alleleResult = fitAlleles(rows, references, config.pipeline.alleleRefinement);
		applyAlleles(rows, alleleResult, config.pipeline.alleleRefinement.reassignmentPolicy, config.pipeline.alleleRefinement.applyMinimumPosterior, headers);
		await writeFile(join(outputDirectory, `${prefix}.allele-models.json`), JSON.stringify(compactAlleleResult(alleleResult), null, 2));
	}
	const records = rows.map((values, ordinal) => airrRowToPostAnalysisRecord({
		ordinal,
		values
	}));
	let collapseResult = null;
	let activeMask = new Uint8Array(rows.length);
	activeMask.fill(1);
	if (config.pipeline.collapse.enabled) {
		process.stderr.write(`${config.pipeline.collapse.mode === "exact" ? "Collapsing exact duplicates" : "Denoising reads"}…\n`);
		if (config.pipeline.collapse.mode === "exact") {
			const shardCount = records.length >= 1e4 ? config.annotation.workers : 1;
			if (shardCount === 1) collapseResult = deduplicate(records, config.pipeline.collapse.key, config.pipeline.collapse.unresolvedPolicy, config.pipeline.collapse.scope, config.pipeline.collapse.respectConstantCall);
			else {
				const plan = createExactDedupPlan(records, config.pipeline.collapse.key, config.pipeline.collapse.unresolvedPolicy, config.pipeline.collapse.scope, config.pipeline.collapse.respectConstantCall, shardCount);
				const workerCount = Math.max(1, Math.min(config.annotation.workers, plan.jobs.length));
				process.stderr.write(`Exact collapse: ${workerCount} worker${workerCount === 1 ? "" : "s"} across ${plan.jobs.length.toLocaleString()} deterministic key shards.\n`);
				collapseResult = finishExactDedupPlan(plan, await runPostAnalysisTasks(plan.jobs.map((job) => ({
					kind: "exact",
					job
				})), workerCount));
			}
		} else {
			const accumulator = new DenoiseAccumulator(records, {
				...config.pipeline.collapse.denoise,
				mode: config.pipeline.collapse.mode,
				scope: config.pipeline.collapse.scope,
				respectConstantCall: config.pipeline.collapse.respectConstantCall,
				unresolvedPolicy: config.pipeline.collapse.unresolvedPolicy
			});
			rows.forEach((values, ordinal) => accumulator.add(ordinal, denoiseVdjSequence({
				ordinal,
				values
			})));
			const jobs = accumulator.preparePartitionJobs();
			const requested = jobs.reduce((total, job) => total + job.variants.length, 0) >= 500 ? config.annotation.workers : 1;
			const workerCount = Math.max(1, Math.min(requested, jobs.length || 1));
			process.stderr.write(`${config.pipeline.collapse.mode.toUpperCase()} denoising: ${workerCount} worker${workerCount === 1 ? "" : "s"} across ${jobs.length.toLocaleString()} independent V/J partition${jobs.length === 1 ? "" : "s"}.\n`);
			const results = await runPostAnalysisTasks(jobs.map((job) => ({
				kind: "denoise",
				job
			})), workerCount);
			collapseResult = accumulator.finishWithPartitionResults(results);
		}
		activeMask = Uint8Array.from(collapseResult.counts, (count) => count > 0 ? 1 : 0);
		addHeaders(headers, ["duplicate_count"]);
		rows.forEach((row, ordinal) => {
			if (collapseResult.counts[ordinal]) row.duplicate_count = String(collapseResult.counts[ordinal]);
		});
	}
	let chimeraSummary = null;
	if (config.pipeline.chimera.enabled) {
		process.stderr.write("Filtering candidate chimeras…\n");
		chimeraSummary = await runChimera(rows, activeMask, config.pipeline.chimera, config.references.scope, headers, references, config.annotation.workers);
	}
	if (config.pipeline.selection.enabled) {
		const errors = validateRepertoireSelection(config.pipeline.selection);
		if (errors.length) throw new Error(errors.join(" "));
		process.stderr.write("Applying explicit repertoire selection…\n");
		rows.forEach((row, ordinal) => {
			if (activeMask[ordinal] && !repertoireRowMatches(row, config.pipeline.selection)) activeMask[ordinal] = 0;
		});
	}
	let lineages = null;
	if (config.pipeline.lineage.enabled) {
		process.stderr.write("Assigning lineages…\n");
		const doubleDMask = Uint8Array.from(rows, (row) => row.d2_call ? 1 : 0);
		lineages = assignLineages(records, {
			identity: config.pipeline.lineage.identity,
			callResolution: config.pipeline.lineage.resolution,
			ambiguity: config.pipeline.lineage.ambiguity,
			productiveOnly: config.pipeline.lineage.productiveOnly,
			requireSameLocus: true,
			maxCandidateComparisons: config.pipeline.lineage.maxCandidateComparisons,
			scope: config.pipeline.lineage.scope
		}, collapseResult ?? void 0, activeMask, doubleDMask);
		addHeaders(headers, ["clone_id"]);
		rows.forEach((row, ordinal) => {
			row.clone_id = lineages.assignments[ordinal] > 0 ? String(lineages.assignments[ordinal]) : "";
		});
	}
	let shm = null;
	if (config.pipeline.shm.enabled) {
		process.stderr.write("Calculating SHM summaries…\n");
		const accumulator = new ShmAccumulator({
			metric: config.pipeline.shm.metric,
			maxSamplesPerLineage: 2e3
		});
		rows.forEach((row, ordinal) => {
			if (activeMask[ordinal]) accumulator.add(row, ordinal, lineages?.assignments[ordinal] ?? 0, "All selected");
		});
		shm = accumulator.finish();
		await writeFile(join(outputDirectory, `${prefix}.shm.json`), JSON.stringify(shm, null, 2));
	}
	let missingAlleles = null;
	if (config.pipeline.missingAlleles.enabled) {
		if (!lineages) throw new Error("Missing-allele screening requires lineage assignment.");
		process.stderr.write("Screening for possible missing V alleles…\n");
		const accumulator = new MissingAlleleAccumulator(config.pipeline.missingAlleles);
		rows.forEach((row, ordinal) => {
			if (activeMask[ordinal]) accumulator.add(row, ordinal, lineages.assignments[ordinal] ?? 0);
		});
		const validator = accumulator.prepareValidation(references.V);
		rows.forEach((row, ordinal) => {
			if (activeMask[ordinal]) validator.add(row, ordinal, lineages.assignments[ordinal] ?? 0);
		});
		missingAlleles = validator.finish();
		await writeFile(join(outputDirectory, `${prefix}.missing-v-alleles.json`), JSON.stringify(missingAlleles, null, 2));
	}
	addHeaders(headers, ["swig_retained"]);
	rows.forEach((row, ordinal) => {
		row.swig_retained = activeMask[ordinal] ? "T" : "F";
	});
	await writeRowsFile(join(outputDirectory, `${prefix}.processed.airr.tsv`), headers, rows, (_, ordinal) => Boolean(activeMask[ordinal]));
	let lineageStudy = null;
	if (config.output.writeLineageStudy && lineages) {
		process.stderr.write("Writing lazy lineage-study bundle…\n");
		lineageStudy = await writeLineageStudy(join(outputDirectory, `${prefix}.lineages.airr.tsv`), join(outputDirectory, `${prefix}.swig-lineage-study.json.gz`), headers, rows, activeMask, lineages, config, references, shm);
	}
	const retained = activeMask.reduce((sum, value) => sum + (value ? 1 : 0), 0);
	const summary = {
		application: "swig-cli",
		version: VERSION,
		completedAt: (/* @__PURE__ */ new Date()).toISOString(),
		inputRecords,
		eligibleRecords,
		annotatedRecords,
		retainedRecords: retained,
		lineages: lineages?.lineageCount ?? 0,
		references: { metadataPreparation: {
			enabled: config.references.prepareMetadata,
			segments: loadedReferences.preparation
		} },
		preprocessing: {
			subsample: config.preprocessing.subsample,
			fastqFilter: {
				options: config.preprocessing.fastqFilter,
				stats: fastqFilterStats
			}
		},
		collapse: collapseResult ? {
			mode: collapseResult.mode,
			inputRecords: collapseResult.inputRecords,
			inputAbundance: collapseResult.inputAbundance,
			uniqueRecords: collapseResult.uniqueRecords,
			collapsedRecords: collapseResult.collapsedRecords,
			warnings: collapseResult.warnings
		} : null,
		chimera: chimeraSummary,
		alleleRefinement: alleleResult ? compactAlleleResult(alleleResult) : null,
		shm: shm ? {
			analyzedRecords: shm.analyzedRecords,
			analyzedAbundance: shm.analyzedAbundance,
			skippedRecords: shm.skippedRecords,
			metric: shm.metric
		} : null,
		missingAlleles: missingAlleles ? {
			candidates: missingAlleles.candidates.length,
			warnings: missingAlleles.warnings
		} : null,
		lineageStudy: lineageStudy ? {
			manifest: `${prefix}.swig-lineage-study.json.gz`,
			airr: `${prefix}.lineages.airr.tsv`,
			records: lineageStudy.linkedAirr.records
		} : null
	};
	await writeFile(join(outputDirectory, `${prefix}.summary.json`), JSON.stringify(summary, null, 2));
	await writeFile(join(outputDirectory, `${prefix}.resolved-config.json`), JSON.stringify(config, null, 2));
	process.stdout.write(`${JSON.stringify(summary, null, 2)}\n`);
}
function vdjTuning(options, callingProfile) {
	if (![
		"-min_D_match",
		"-min_J_length",
		"-num_alignments_D",
		"-num_alignments_J",
		"-D_penalty",
		"-J_penalty"
	].some((key) => options[key] !== void 0)) return void 0;
	const agreement = callingProfile !== "truth_optimized";
	const minD = options["-min_D_match"] === void 0 ? agreement ? 5 : 6 : parseIntegerOption(options["-min_D_match"], "-min_D_match", { minimum: 5 });
	const minJ = options["-min_J_length"] === void 0 ? 10 : parseIntegerOption(options["-min_J_length"], "-min_J_length", { minimum: 0 });
	const topD = options["-num_alignments_D"] === void 0 ? agreement ? 3 : 2 : parseIntegerOption(options["-num_alignments_D"], "-num_alignments_D", {
		minimum: 1,
		allowZero: false
	});
	const topJ = options["-num_alignments_J"] === void 0 ? 2 : parseIntegerOption(options["-num_alignments_J"], "-num_alignments_J", {
		minimum: 1,
		allowZero: false
	});
	const dMismatch = options["-D_penalty"] === void 0 ? agreement ? -4 : -3 : parseFiniteOption(options["-D_penalty"], "-D_penalty");
	const jMismatch = options["-J_penalty"] === void 0 ? agreement ? -4 : -3 : parseFiniteOption(options["-J_penalty"], "-J_penalty");
	if (options["-D_penalty"] !== void 0 && (!Number.isInteger(dMismatch) || dMismatch <= -5 || dMismatch >= 0)) throw new Error("-D_penalty must be an integer greater than -5 and less than 0.");
	if (options["-J_penalty"] !== void 0 && (!Number.isInteger(jMismatch) || jMismatch <= -4 || jMismatch >= 0)) throw new Error("-J_penalty must be an integer greater than -4 and less than 0.");
	return {
		dMatch: 2,
		dMismatch,
		dGapOpen: agreement ? -11 : -13,
		dGapExtend: -1,
		topD,
		minDMatch: minD,
		jMatch: 2,
		jMismatch,
		jGapOpen: agreement ? -13 : -17,
		jGapExtend: agreement ? -1 : -2,
		topJ,
		minJLength: Math.max(1, minJ)
	};
}
async function runVdj(rawArgs, assets) {
	const options = parseVdjArguments(rawArgs);
	if (options.help) {
		process.stdout.write(`${vdjUsage()}\n`);
		return;
	}
	if (options.version) {
		process.stdout.write(`${VERSION}\n`);
		return;
	}
	const outputValue = options["-out"];
	if (!outputValue) throw new Error("--vdj requires -out (or --out) so AIRR rows can be streamed to an explicit destination.");
	if (String(options["-outfmt"] ?? "19").trim() !== "19") throw new Error("swig-cli --vdj currently emits only AIRR rearrangement format; use -outfmt 19.");
	if (String(options["-domain_system"] ?? "imgt").toLowerCase() !== "imgt") throw new Error("swig-cli --vdj supports only -domain_system imgt; Kabat coordinates must not be labeled as IMGT/AIRR annotations.");
	const strand = {
		both: 0,
		plus: 1,
		minus: 2
	}[String(options["-strand"] ?? "both").toLowerCase()];
	if (strand === void 0) throw new Error("-strand must be both, plus, or minus.");
	const minimumIdentity = options["--minimum-identity"] === void 0 ? .6 : parseFiniteOption(options["--minimum-identity"], "--minimum-identity");
	if (minimumIdentity < 0 || minimumIdentity > 1) throw new Error("--minimum-identity must be between 0 and 1.");
	const threadValue = options["--workers"] ?? options["-num_threads"];
	let workers = threadValue === void 0 ? Math.max(1, Math.min(4, availableParallelism())) : parseIntegerOption(threadValue, options["--workers"] !== void 0 ? "--workers" : "-num_threads", { minimum: 0 });
	if (options["--workers"] === void 0 && threadValue !== void 0 && workers === 0) throw new Error("-num_threads must be at least 1; use --workers 0 for automatic selection.");
	if (workers === 0) workers = Math.max(1, Math.min(8, availableParallelism()));
	const batchRecords = (options["--batch-records"] === void 0 ? 0 : parseIntegerOption(options["--batch-records"], "--batch-records", { minimum: 0 })) || automaticBatchRecords(workers);
	const assigner = String(options["--assigner"] ?? "riat_mp");
	if (![
		"standard",
		"riat_mp",
		"aer",
		"aer_robust"
	].includes(assigner)) throw new Error("--assigner must be standard, riat_mp, aer, or aer_robust.");
	const callingProfile = String(options["--calling-profile"] ?? "truth_optimized");
	if (![
		"truth_optimized",
		"igblast_compatible",
		"igblast_balanced"
	].includes(callingProfile)) throw new Error("--calling-profile must be truth_optimized, igblast_compatible, or igblast_balanced.");
	const queryValue = String(options["-query"] ?? "-");
	const queryPath = queryValue === "-" ? "-" : resolve(queryValue);
	const outputPath = String(outputValue) === "-" ? "-" : resolve(String(outputValue));
	if (outputPath !== "-") await mkdir(dirname(outputPath), { recursive: true });
	const prepared = await prepareVdjReferences(options, assets);
	const tuning = vdjTuning(options, callingProfile);
	const pool = new WasmPool(workers, workerInitialization(assets.wasmPath, prepared.references, callingProfile, assigner, tuning));
	await pool.start();
	const output = outputPath === "-" ? process.stdout : createCliWriteStream(outputPath);
	let outputHeader = null;
	let records = 0;
	let completed = false;
	try {
		if (outputPath !== "-") await once(output, "open");
		const assignerLabel = assigner === "riat_mp" ? "RIAT-MP" : assigner === "aer" ? "AER" : assigner === "aer_robust" ? "AER-R (experimental)" : "standard SwiftIG";
		process.stderr.write(`Streaming SwiftIG V(D)J assignments (${prepared.mode}; ${workers} worker${workers === 1 ? "" : "s"}; ${batchRecords.toLocaleString()} records/batch; ${assignerLabel}) to ${outputPath}.\n`);
		const state = {
			inputRecords: 0,
			eligibleRecords: 0,
			selectedRecords: 0,
			fastqFilter: emptyFastqQualityFilterStats(false, false)
		};
		const pending = [];
		const consume = async (item) => {
			const result = await item.promise;
			if (outputHeader === null) {
				outputHeader = result.header;
				await writeChunk(output, `${outputHeader}\n`);
			} else if (result.header !== outputHeader) throw new Error("SwiftIG changed the AIRR schema between V(D)J batches.");
			await writeChunk(output, trimAuxiliaryFwr4(result.header, result.body, prepared.fwr4EndOffsets, prepared.jLengths));
			records += result.count;
		};
		const preprocessing = {
			fastqFilter: {
				...DEFAULT_CLI_CONFIG.preprocessing.fastqFilter,
				trim3Prime: { ...DEFAULT_CLI_CONFIG.preprocessing.fastqFilter.trim3Prime }
			},
			subsample: {
				...DEFAULT_CLI_CONFIG.preprocessing.subsample,
				enabled: false
			}
		};
		const input = {
			path: queryPath,
			format: "fasta"
		};
		for await (const batch of sequenceBatches(input, batchRecords, preprocessing, "reannotate", 0, state)) {
			pending.push({ promise: pool.run(workerAnnotation(batch.text, batch.count, 1, minimumIdentity, strand, { mode: "off" })) });
			if (pending.length >= Math.max(2, workers * 2)) await consume(pending.shift());
		}
		while (pending.length) await consume(pending.shift());
		if (outputPath !== "-") await finishWritable(output);
		completed = true;
	} finally {
		if (!completed && outputPath !== "-") output.destroy();
		await pool.close();
	}
	process.stderr.write(`Completed ${records.toLocaleString()} streaming V(D)J assignment${records === 1 ? "" : "s"}.\n`);
}
async function runCli(assets = defaultCliAssets()) {
	const args = process.argv.slice(2);
	if (args.includes("--precompute_aux") || args.includes("--precompute-aux")) {
		await runPrepareReference(args, assets);
		return;
	}
	if (args.includes("--vdj")) {
		await runVdj(args, assets);
		return;
	}
	const command = args[0] && !args[0].startsWith("-") ? args[0] : "run";
	const rest = command === args[0] ? args.slice(1) : args;
	if (command === "prepare-reference") {
		await runPrepareReference(rest, assets);
		return;
	}
	if (hasFlag(args, "--help") || command === "help") {
		process.stdout.write(`${usage()}\n`);
		return;
	}
	if (hasFlag(args, "--version") || command === "version") {
		process.stdout.write(`${VERSION}\n`);
		return;
	}
	if (hasFlag(args, "--notices") || command === "notices") {
		process.stdout.write(`${thirdPartyNotices()}\n`);
		return;
	}
	if (command === "init") {
		const target = positional(rest)[0] ?? "swig.config.json";
		const example = normalizeCliConfig({
			...DEFAULT_CLI_CONFIG,
			inputs: [{
				path: "reads.fastq.gz",
				sampleId: "sample-1",
				subjectId: "donor-1"
			}]
		});
		writeFileSync(target, `${JSON.stringify(example, null, 2)}\n`, { flag: "wx" });
		process.stdout.write(`Created ${target}\n`);
		return;
	}
	if (command !== "run") throw new Error(`Unknown command ${command}.\n\n${usage()}`);
	const configPath = argumentValue(rest, "--config");
	let base = process.cwd();
	let raw = {};
	if (configPath) {
		const absolute = resolve(configPath);
		base = dirname(absolute);
		raw = JSON.parse(await readFile(absolute, "utf8"));
	}
	const inputs = positional(rest);
	if (inputs.length) raw.inputs = inputs.map((path, index) => ({
		path,
		datasetId: index ? void 0 : argumentValue(rest, "--dataset"),
		sampleId: index ? void 0 : argumentValue(rest, "--sample"),
		subjectId: index ? void 0 : argumentValue(rest, "--donor")
	}));
	const commandOutput = argumentValue(rest, "--out");
	const configuredOutput = typeof raw.output?.directory === "string" && raw.output.directory.trim();
	if (!commandOutput && !configuredOutput) throw new Error("The pipeline CLI requires an explicit output directory in output.directory or --out so results can be streamed to disk.");
	if (commandOutput) raw.output = {
		...raw.output ?? {},
		directory: commandOutput
	};
	const commandWorkers = argumentValue(rest, "--workers");
	if (commandWorkers !== void 0) raw.annotation = {
		...raw.annotation ?? {},
		workers: parseIntegerOption(commandWorkers, "--workers", { minimum: 0 })
	};
	const config = normalizeCliConfig(raw);
	if (![
		"standard",
		"riat_mp",
		"aer",
		"aer_robust"
	].includes(config.annotation.assignerStrategy)) throw new Error("annotation.assignerStrategy must be standard, riat_mp, aer, or aer_robust.");
	if (config.annotation.workers === 0) config.annotation.workers = Math.max(1, Math.min(8, availableParallelism()));
	config.inputs = config.inputs.map((input) => ({
		...input,
		path: typeof input.inline === "string" ? input.path : resolveFrom(base, input.path)
	}));
	await runPipeline(config, base, assets);
}
//#endregion
//#region cli-src/swig-cli-node.mjs
runCli().catch((error) => {
	process.stderr.write(`swig-cli: ${error instanceof Error ? error.message : String(error)}\n`);
	process.exitCode = 1;
});
//#endregion
