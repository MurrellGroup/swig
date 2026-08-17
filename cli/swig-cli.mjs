#!/usr/bin/env node
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
//#region src/study-design.ts
const DEFAULT_PIPELINE_PLAN = {
	enabled: false,
	collapse: {
		enabled: true,
		mode: "exact",
		key: "trimmed",
		scope: "sample",
		unresolvedPolicy: "discard",
		respectConstantCall: true
	},
	chimera: {
		enabled: false,
		segment: "V",
		model: "auto",
		posteriorThreshold: .95,
		retainUnevaluated: true,
		msaSource: "selected",
		uploadedMsa: "",
		uploadedMsaName: ""
	},
	selection: {
		enabled: false,
		datasetId: "",
		sampleId: "",
		subjectId: "",
		cohort: "",
		timepoint: "",
		compartment: "",
		locus: "",
		vCall: "",
		vCallIncludeAmbiguous: false,
		jCall: "",
		jCallIncludeAmbiguous: false,
		cdr3Nt: "",
		cdr3Aa: "",
		productive: "any",
		hasCdr3: "any",
		doubleD: "any"
	},
	alleleRefinement: {
		enabled: false,
		model: "dirichlet",
		scope: "subject",
		segments: ["V", "J"],
		weighting: "unique",
		baselineNeighbourOdds: .01,
		shmLeakageSensitivity: 1,
		reassignmentPolicy: "confidence",
		applyMinimumPosterior: .8
	},
	lineage: {
		enabled: true,
		scope: "subject",
		identity: .85,
		resolution: "gene",
		ambiguity: "overlap",
		productiveOnly: true
	},
	shm: {
		enabled: true,
		metric: "vNtRate"
	},
	missingAlleles: { enabled: false }
};
function datasetScopeKey(record, scope = "global") {
	return `${scope}:${datasetScopeValue(record, scope)}`;
}
function datasetScopeValue(record, scope = "global") {
	if (scope === "global") return "complete study";
	if (scope === "dataset") return record.datasetId || "legacy";
	if (scope === "sample") return record.sampleId || record.datasetId || "legacy";
	if (scope === "subject") return record.subjectId || record.sampleId || record.datasetId || "legacy";
	return record.cohort || record.subjectId || record.sampleId || record.datasetId || "legacy";
}
const METADATA_FIELDS = [
	"swig_dataset_id",
	"sample_id",
	"subject_id",
	"swig_cohort",
	"swig_timepoint",
	"swig_compartment",
	"swig_source_sequence_id"
];
function normalizedLines(body) {
	return (typeof body === "string" ? body : new TextDecoder().decode(body)).split("\n").map((line) => line.replace(/\r$/, "")).filter(Boolean);
}
function replaceOrAppend(values, positions, field, value) {
	const position = positions.get(field);
	if (position === void 0) values.push(value);
	else values[position] = value;
}
/** Adds study metadata to a SwiftIG AIRR batch without materializing the input file. */
function annotateAirrBatch(headerLine, body, dataset) {
	const headers = headerLine.replace(/\r$/, "").split("	");
	const positions = new Map(headers.map((field, index) => [field, index]));
	for (const field of METADATA_FIELDS) if (!positions.has(field)) headers.push(field);
	const sourceSequencePosition = positions.get("sequence_id");
	const bodyText = normalizedLines(body).map((line) => {
		const values = line.split("	");
		const sourceSequenceId = sourceSequencePosition === void 0 ? "" : values[sourceSequencePosition] ?? "";
		if (sourceSequencePosition !== void 0) values[sourceSequencePosition] = `${dataset.datasetId}::${sourceSequenceId || "record"}`;
		replaceOrAppend(values, positions, "swig_dataset_id", dataset.datasetId);
		replaceOrAppend(values, positions, "sample_id", dataset.sampleId);
		replaceOrAppend(values, positions, "subject_id", dataset.subjectId);
		replaceOrAppend(values, positions, "swig_cohort", dataset.cohort);
		replaceOrAppend(values, positions, "swig_timepoint", dataset.timepoint ?? "");
		replaceOrAppend(values, positions, "swig_compartment", dataset.compartment ?? "");
		replaceOrAppend(values, positions, "swig_source_sequence_id", sourceSequenceId);
		return values.join("	");
	}).join("\n");
	return {
		header: headers.join("	"),
		body: bodyText ? `${bodyText}\n` : ""
	};
}
/** Mirrors metadata and collision-safe IDs into the sparse double-D table. */
function annotateDoubleDBatch(headerLine, body, dataset) {
	const headers = headerLine.replace(/\r$/, "").split("	");
	const positions = new Map(headers.map((field, index) => [field, index]));
	for (const field of METADATA_FIELDS.slice(0, 6)) if (!positions.has(field)) headers.push(field);
	const sequencePosition = positions.get("sequence_id");
	const bodyText = normalizedLines(body).map((line) => {
		const values = line.split("	");
		if (sequencePosition !== void 0) values[sequencePosition] = `${dataset.datasetId}::${values[sequencePosition] || "record"}`;
		replaceOrAppend(values, positions, "swig_dataset_id", dataset.datasetId);
		replaceOrAppend(values, positions, "sample_id", dataset.sampleId);
		replaceOrAppend(values, positions, "subject_id", dataset.subjectId);
		replaceOrAppend(values, positions, "swig_cohort", dataset.cohort);
		replaceOrAppend(values, positions, "swig_timepoint", dataset.timepoint ?? "");
		replaceOrAppend(values, positions, "swig_compartment", dataset.compartment ?? "");
		return values.join("	");
	}).join("\n");
	return {
		header: headers.join("	"),
		body: bodyText ? `${bodyText}\n` : ""
	};
}
function stableDatasetSeed(baseSeed, datasetIndex) {
	return Math.trunc(baseSeed) + Math.imul(datasetIndex + 1, 1831565813) | 0;
}
//#endregion
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
//#region src/post-analysis-core.ts
const HASH_SEEDS = [
	2166136261,
	2654435769,
	2246822507,
	3266489909,
	668265263,
	374761393,
	3550635116,
	4251993797
];
function mix32(value) {
	value ^= value >>> 16;
	value = Math.imul(value, 2146121005);
	value ^= value >>> 15;
	value = Math.imul(value, 2221713035);
	value ^= value >>> 16;
	return value >>> 0;
}
function hashSequence(value, seed = HASH_SEEDS[0]) {
	let hash = seed >>> 0;
	for (let index = 0; index < value.length; index += 1) {
		hash ^= value.charCodeAt(index);
		hash = Math.imul(hash, 16777619);
	}
	return mix32(hash);
}
function sequenceFingerprint(value) {
	const normalized = normalizeNt(value);
	return `${normalized.length}:${[
		0,
		2,
		4,
		6
	].map((index) => hashSequence(normalized, HASH_SEEDS[index]).toString(36)).join(":")}`;
}
function normalizeNt(value) {
	return value.toUpperCase().replaceAll("U", "T").replace(/[^ACGTN]/g, "");
}
function normalizeCall(call, resolution) {
	const clean = call.trim().toUpperCase();
	return resolution === "gene" ? clean.replace(/\*.*$/, "") : clean;
}
function callSet(value, resolution, ambiguity) {
	const values = value.split(",").map((call) => normalizeCall(call, resolution)).filter(Boolean);
	return ambiguity === "top" ? values.slice(0, 1) : [...new Set(values)].sort();
}
function callsCompatible(left, right, resolution, ambiguity) {
	const a = callSet(left, resolution, ambiguity);
	const b = callSet(right, resolution, ambiguity);
	if (!a.length || !b.length) return false;
	if (ambiguity === "strict") return a.length === b.length && a.every((value, index) => value === b[index]);
	const values = new Set(a);
	return b.some((value) => values.has(value));
}
function hammingDistanceWithin(left, right, maximum, ambiguousN = true) {
	if (left.length !== right.length) return maximum + 1;
	let distance = 0;
	for (let index = 0; index < left.length; index += 1) if (left[index] !== right[index] && (!ambiguousN || left[index] !== "N" && right[index] !== "N")) {
		distance += 1;
		if (distance > maximum) return distance;
	}
	return distance;
}
function constantCallPartition(record, enabled) {
	if (!enabled) return "";
	return `\u0000C:${callSet(record.cCall ?? "", "gene", "top")[0] ?? "__C_UNASSIGNED__"}`;
}
function dedupKey(record, key, scope, respectConstantCall) {
	const prefix = `${datasetScopeKey(record, scope)}${constantCallPartition(record, respectConstantCall)}\u0000`;
	if (key === "sequence") return `${prefix}${record.sequenceFingerprint}`;
	if (key === "trimmed") return `${prefix}${record.trimmedFingerprint}`;
	if (key === "cdr3") return `${prefix}${record.locus}\u0000${record.cdr3Nt}`;
	return `${prefix}${record.locus}\u0000${record.vCall}\u0000${record.jCall}\u0000${record.cdr3Nt}`;
}
function largestCountGroups(counts, limit = 100) {
	const heap = [];
	const worse = (left, right) => counts[left] < counts[right] || counts[left] === counts[right] && left > right;
	const rise = (start) => {
		let index = start;
		while (index > 0) {
			const parent = Math.floor((index - 1) / 2);
			if (!worse(heap[index], heap[parent])) break;
			[heap[index], heap[parent]] = [heap[parent], heap[index]];
			index = parent;
		}
	};
	const sink = () => {
		let index = 0;
		while (true) {
			const left = index * 2 + 1;
			const right = left + 1;
			let next = index;
			if (left < heap.length && worse(heap[left], heap[next])) next = left;
			if (right < heap.length && worse(heap[right], heap[next])) next = right;
			if (next === index) break;
			[heap[index], heap[next]] = [heap[next], heap[index]];
			index = next;
		}
	};
	for (let ordinal = 0; ordinal < counts.length; ordinal += 1) {
		if (counts[ordinal] <= 1) continue;
		if (heap.length < limit) {
			heap.push(ordinal);
			rise(heap.length - 1);
		} else if (worse(heap[0], ordinal)) {
			heap[0] = ordinal;
			sink();
		}
	}
	return heap.sort((left, right) => counts[right] - counts[left] || left - right).map((ordinal) => ({
		ordinal,
		count: counts[ordinal]
	}));
}
function hasUsableDedupKey(record, key) {
	if (key === "sequence") return Number(record.sequenceFingerprint.split(":", 1)[0]) > 0;
	if (key === "trimmed") return Number(record.trimmedFingerprint.split(":", 1)[0]) > 0;
	if (key === "cdr3") return Boolean(record.locus && record.cdr3Nt);
	return Boolean(record.locus && record.vCall && record.jCall && record.cdr3Nt);
}
function deduplicate(records, key, unresolvedPolicy = "discard", scope = "global", respectConstantCall = true) {
	const representatives = new Int32Array(records.length);
	representatives.fill(-1);
	const counts = new Uint32Array(records.length);
	const seen = /* @__PURE__ */ new Map();
	let inputAbundance = 0;
	let unresolvedRecords = 0;
	let retainedUnresolved = 0;
	for (let index = 0; index < records.length; index += 1) {
		const weight = Math.max(1, Math.floor(records[index].inputCount ?? 1));
		inputAbundance += weight;
		if (!hasUsableDedupKey(records[index], key)) {
			unresolvedRecords += 1;
			if (unresolvedPolicy === "retain") {
				representatives[index] = index;
				counts[index] = weight;
				retainedUnresolved += 1;
			}
			continue;
		}
		const value = dedupKey(records[index], key, scope, respectConstantCall);
		const previous = seen.get(value);
		if (previous === void 0) {
			seen.set(value, index);
			representatives[index] = index;
			counts[index] = weight;
		} else {
			representatives[index] = previous;
			counts[previous] += weight;
		}
	}
	const largestGroups = largestCountGroups(counts);
	const uniqueRecords = seen.size + retainedUnresolved;
	return {
		mode: "exact",
		key,
		algorithm: "Exact key collapse",
		inputRecords: records.length,
		inputAbundance,
		uniqueRecords,
		collapsedRecords: records.length - uniqueRecords,
		representatives,
		counts,
		largestGroups,
		partitions: 1,
		candidateComparisons: 0,
		indelMergedVariants: 0,
		substitutionMergedVariants: 0,
		excludedAmbiguous: 0,
		unresolvedRecords,
		warnings: unresolvedRecords ? [`${unresolvedRecords.toLocaleString()} records without a usable ${key} key were ${unresolvedPolicy === "retain" ? "retained unchanged" : "discarded from the downstream representative set"}.`] : []
	};
}
/**
* Append-only two-bit DNA storage. Denoising may need the complete trimmed VDJ
* sequence, but retaining one JavaScript string per unique read is prohibitive
* for large repertoires. Fixed chunks also avoid a transient 2x allocation
* when an expanding typed array is copied.
*/
var PackedDnaArena = class {
	chunkWords = 262144;
	chunks = [];
	used = [];
	append(sequence) {
		const words = Math.ceil(sequence.length / 16);
		let chunk = this.chunks.length - 1;
		if (chunk < 0 || this.used[chunk] + words > this.chunkWords) {
			chunk = this.chunks.length;
			this.chunks.push(new Uint32Array(Math.max(this.chunkWords, words)));
			this.used.push(0);
		}
		const offset = this.used[chunk];
		this.used[chunk] += words;
		const target = this.chunks[chunk];
		for (let index = 0; index < sequence.length; index += 1) {
			const code = sequence[index] === "C" ? 1 : sequence[index] === "G" ? 2 : sequence[index] === "T" ? 3 : 0;
			target[offset + (index >>> 4)] |= code << (index & 15) * 2;
		}
		return {
			chunk,
			offset,
			length: sequence.length
		};
	}
	base(location, index) {
		const code = this.chunks[location.chunk][location.offset + (index >>> 4)] >>> (index & 15) * 2 & 3;
		return code === 0 ? "A" : code === 1 ? "C" : code === 2 ? "G" : "T";
	}
	equals(location, sequence) {
		if (location.length !== sequence.length) return false;
		for (let index = 0; index < sequence.length; index += 1) if (this.base(location, index) !== sequence[index]) return false;
		return true;
	}
	decode(location) {
		let sequence = "";
		for (let start = 0; start < location.length; start += 1024) {
			const values = [];
			const end = Math.min(location.length, start + 1024);
			for (let index = start; index < end; index += 1) values.push(this.base(location, index));
			sequence += values.join("");
		}
		return sequence;
	}
};
function denoisePartition(record, options) {
	const v = callSet(record.vCall, options.callResolution, options.ambiguity);
	const j = callSet(record.jCall, options.callResolution, options.ambiguity);
	if (!v.length || !j.length) return null;
	return `${datasetScopeKey(record, options.scope ?? "global")}${constantCallPartition(record, options.respectConstantCall ?? true)}\u0000${record.locus || "?"}\u0000${v.join("+")}\u0000${j.join("+")}`;
}
function kmerProfile(sequence, blockCount, k = 6) {
	const space = 1 << 2 * k;
	const counts = new Uint16Array(space);
	const touched = [];
	const mask = space - 1;
	let code = 0;
	for (let index = 0; index < sequence.length; index += 1) {
		const base = sequence[index] === "C" ? 1 : sequence[index] === "G" ? 2 : sequence[index] === "T" ? 3 : 0;
		code = (code << 2 | base) & mask;
		if (index < k - 1) continue;
		if (!counts[code]) touched.push(code);
		if (counts[code] < 65535) counts[code] += 1;
	}
	touched.sort((left, right) => left - right);
	const codes = Uint16Array.from(touched);
	const sparseCounts = Uint16Array.from(touched, (value) => counts[value]);
	const hashes = new Uint32Array(blockCount);
	for (let block = 0; block < blockCount; block += 1) hashes[block] = mix32(2654435769 ^ block);
	for (let index = 0; index < codes.length; index += 1) {
		const block = codes[index] % blockCount;
		hashes[block] = mix32(hashes[block] ^ mix32(Math.imul(codes[index] + 1, 2246822507) ^ sparseCounts[index]));
	}
	return {
		codes,
		counts: sparseCounts,
		hashes
	};
}
function kmerSquaredDistance(left, right, ceiling = Number.POSITIVE_INFINITY) {
	let a = 0;
	let b = 0;
	let distance = 0;
	while (a < left.codes.length || b < right.codes.length) {
		if (b >= right.codes.length || a < left.codes.length && left.codes[a] < right.codes[b]) {
			distance += left.counts[a] ** 2;
			a += 1;
		} else if (a >= left.codes.length || right.codes[b] < left.codes[a]) {
			distance += right.counts[b] ** 2;
			b += 1;
		} else {
			distance += (left.counts[a] - right.counts[b]) ** 2;
			a += 1;
			b += 1;
		}
		if (distance > ceiling) return distance;
	}
	return distance;
}
/** Exact metric index for FAD's final nearest-centroid assignment. */
function buildKmerVpTree(points, profiles, compared) {
	if (!points.length) return null;
	const point = points[points.length - 1];
	if (points.length === 1) return {
		point,
		radius: 0,
		inside: null,
		outside: null
	};
	const distances = points.slice(0, -1).map((candidate) => {
		compared();
		return {
			candidate,
			distance: Math.sqrt(kmerSquaredDistance(profiles.get(point), profiles.get(candidate)))
		};
	}).sort((left, right) => left.distance - right.distance || left.candidate - right.candidate);
	const middle = Math.floor(distances.length / 2);
	return {
		point,
		radius: distances[middle]?.distance ?? 0,
		inside: buildKmerVpTree(distances.slice(0, middle).map((value) => value.candidate), profiles, compared),
		outside: buildKmerVpTree(distances.slice(middle).map((value) => value.candidate), profiles, compared)
	};
}
function nearestKmerPoint(tree, query, profiles, abundance, compared) {
	let best = tree.point;
	let bestDistance = Number.POSITIVE_INFINITY;
	const visit = (node) => {
		if (!node) return;
		compared();
		const distance = Math.sqrt(kmerSquaredDistance(profiles.get(query), profiles.get(node.point)));
		if (distance < bestDistance || distance === bestDistance && abundance(node.point) > abundance(best)) {
			best = node.point;
			bestDistance = distance;
		}
		if (!node.inside && !node.outside) return;
		if (distance < node.radius) {
			visit(node.inside);
			if (distance + bestDistance >= node.radius) visit(node.outside);
		} else {
			visit(node.outside);
			if (distance - bestDistance <= node.radius) visit(node.inside);
		}
	};
	visit(tree);
	return best;
}
function logGamma$1(value) {
	const coefficients = [
		.9999999999998099,
		676.5203681218851,
		-1259.1392167224028,
		771.3234287776531,
		-176.6150291621406,
		12.50734327868691,
		-.1385710952657201,
		9984369578019572e-21,
		1.505632735149312e-7
	];
	if (value < .5) return Math.log(Math.PI) - Math.log(Math.sin(Math.PI * value)) - logGamma$1(1 - value);
	let x = coefficients[0];
	const z = value - 1;
	for (let index = 1; index < coefficients.length; index += 1) x += coefficients[index] / (z + index);
	const t = z + 7.5;
	return .5 * Math.log(2 * Math.PI) + (z + .5) * Math.log(t) - t + Math.log(x);
}
function regularizedGammaP(shape, value) {
	if (!(shape > 0) || value < 0 || Number.isNaN(value)) return NaN;
	if (value === 0) return 0;
	const logScale = -value + shape * Math.log(value) - logGamma$1(shape);
	if (value < shape + 1) {
		let term = 1 / shape;
		let sum = term;
		let denominator = shape;
		for (let iteration = 0; iteration < 1e4; iteration += 1) {
			denominator += 1;
			term *= value / denominator;
			sum += term;
			if (Math.abs(term) < Math.abs(sum) * 1e-14) break;
		}
		return Math.min(1, Math.max(0, sum * Math.exp(logScale)));
	}
	let b = value + 1 - shape;
	const floor = 1e-300;
	let c = 1 / floor;
	let d = 1 / Math.max(floor, b);
	let fraction = d;
	for (let iteration = 1; iteration < 1e4; iteration += 1) {
		const coefficient = -iteration * (iteration - shape);
		b += 2;
		d = coefficient * d + b;
		if (Math.abs(d) < floor) d = floor;
		c = b + coefficient / c;
		if (Math.abs(c) < floor) c = floor;
		d = 1 / d;
		const delta = d * c;
		fraction *= delta;
		if (Math.abs(delta - 1) < 1e-14) break;
	}
	const q = Math.exp(logScale) * fraction;
	return Math.min(1, Math.max(0, 1 - q));
}
/** P(Poisson(lambda) > observed), matching Distributions.jl ccdf. */
function poissonStrictUpperTail(observed, lambda) {
	if (!(lambda > 0)) return 0;
	if (observed < 0) return 1;
	return regularizedGammaP(Math.floor(observed) + 1, lambda);
}
function alternativeCount(length, distance) {
	if (distance <= 0) return 1;
	let combinations = 1;
	for (let index = 1; index <= distance; index += 1) combinations *= (length - distance + index) / index;
	return Math.max(1, combinations * 3 ** distance);
}
function sequenceBlocks(sequence, count) {
	const result = [];
	for (let index = 0; index < count; index += 1) {
		const start = Math.floor(index * sequence.length / count);
		const end = Math.floor((index + 1) * sequence.length / count);
		result.push(sequence.slice(start, end));
	}
	return result;
}
function indexedEditSegments(sequence, count) {
	const result = [];
	for (let index = 0; index < count; index += 1) {
		const start = Math.floor(index * sequence.length / count);
		const end = Math.floor((index + 1) * sequence.length / count);
		result.push({
			index,
			start,
			length: end - start,
			value: sequence.slice(start, end)
		});
	}
	return result;
}
/**
* Allocation-bounded Ukkonen-style edit profiler. The narrow diagonal band is
* O(length × maximum), and tie-breaking chooses fewer indels when two optimal
* paths exist so the aggressive indel rule is not triggered by an arbitrary
* alignment of equal-length sequences.
*/
function createBoundedEditProfiler(maximum) {
	const width = maximum * 2 + 3;
	const infinity = maximum + 1;
	let previousCost = new Int16Array(width);
	let previousSubstitutions = new Int16Array(width);
	let previousInsertions = new Int16Array(width);
	let previousDeletions = new Int16Array(width);
	let currentCost = new Int16Array(width);
	let currentSubstitutions = new Int16Array(width);
	let currentInsertions = new Int16Array(width);
	let currentDeletions = new Int16Array(width);
	return (parent, child) => {
		if (Math.abs(parent.length - child.length) > maximum) return null;
		previousCost.fill(infinity);
		previousSubstitutions.fill(0);
		previousInsertions.fill(0);
		previousDeletions.fill(0);
		for (let column = 0; column <= Math.min(child.length, maximum); column += 1) {
			const offset = column + maximum + 1;
			previousCost[offset] = column;
			previousInsertions[offset] = column;
		}
		for (let row = 1; row <= parent.length; row += 1) {
			currentCost.fill(infinity);
			currentSubstitutions.fill(0);
			currentInsertions.fill(0);
			currentDeletions.fill(0);
			const firstColumn = Math.max(0, row - maximum);
			const lastColumn = Math.min(child.length, row + maximum);
			for (let column = firstColumn; column <= lastColumn; column += 1) {
				const offset = column - row + maximum + 1;
				let bestCost = infinity;
				let bestSubstitutions = 0;
				let bestInsertions = 0;
				let bestDeletions = 0;
				const consider = (cost, substitutions, insertions, deletions) => {
					const indels = insertions + deletions;
					const bestIndels = bestInsertions + bestDeletions;
					if (cost < bestCost || cost === bestCost && indels < bestIndels || cost === bestCost && indels === bestIndels && substitutions < bestSubstitutions || cost === bestCost && indels === bestIndels && substitutions === bestSubstitutions && deletions < bestDeletions) {
						bestCost = cost;
						bestSubstitutions = substitutions;
						bestInsertions = insertions;
						bestDeletions = deletions;
					}
				};
				if (column > 0 && previousCost[offset] <= maximum) {
					const mismatch = parent[row - 1] === child[column - 1] ? 0 : 1;
					consider(previousCost[offset] + mismatch, previousSubstitutions[offset] + mismatch, previousInsertions[offset], previousDeletions[offset]);
				}
				if (previousCost[offset + 1] <= maximum) consider(previousCost[offset + 1] + 1, previousSubstitutions[offset + 1], previousInsertions[offset + 1], previousDeletions[offset + 1] + 1);
				if (column > 0 && currentCost[offset - 1] <= maximum) consider(currentCost[offset - 1] + 1, currentSubstitutions[offset - 1], currentInsertions[offset - 1] + 1, currentDeletions[offset - 1]);
				currentCost[offset] = bestCost;
				currentSubstitutions[offset] = bestSubstitutions;
				currentInsertions[offset] = bestInsertions;
				currentDeletions[offset] = bestDeletions;
			}
			[previousCost, currentCost] = [currentCost, previousCost];
			[previousSubstitutions, currentSubstitutions] = [currentSubstitutions, previousSubstitutions];
			[previousInsertions, currentInsertions] = [currentInsertions, previousInsertions];
			[previousDeletions, currentDeletions] = [currentDeletions, previousDeletions];
		}
		const offset = child.length - parent.length + maximum + 1;
		const distance = previousCost[offset];
		if (distance > maximum) return null;
		return {
			distance,
			substitutions: previousSubstitutions[offset],
			insertions: previousInsertions[offset],
			deletions: previousDeletions[offset]
		};
	};
}
/**
* Streaming builder used by the post-analysis worker. Exact variants are
* dereplicated while batches are scanned, sequences live in a compact 2-bit
* arena, and only one V/J partition's temporary neighbor index exists at a
* time during finalization.
*/
var DenoiseAccumulator = class {
	arena = new PackedDnaArena();
	variants = [];
	variantsByFingerprint = /* @__PURE__ */ new Map();
	variantByOrdinal;
	processed;
	standalone = [];
	retainedAmbiguous = /* @__PURE__ */ new Map();
	excludedAmbiguous = 0;
	unresolvedRecords = 0;
	records;
	options;
	constructor(records, options) {
		this.records = records;
		this.options = options;
		this.variantByOrdinal = new Int32Array(records.length);
		this.variantByOrdinal.fill(-1);
		this.processed = new Uint8Array(records.length);
		if (!(options.errorRate > 0 && options.errorRate < 1)) throw new Error("The denoising error rate must be between 0 and 1.");
		if (!(options.alpha > 0 && options.alpha < 1)) throw new Error("The denoising alpha must be between 0 and 1.");
		if (options.fadNeighborThreshold < 0) throw new Error("The FAD neighbor threshold cannot be negative.");
		if (options.maximumHammingDistance < 1 || options.maximumHammingDistance > 4) throw new Error("The conservative Hamming radius must be from 1 to 4.");
		if (!Number.isInteger(options.maximumEditDistance) || options.maximumEditDistance < 1 || options.maximumEditDistance > 2) throw new Error("The indel-aware edit radius must be 1 or 2.");
		if (!(options.minimumIndelParentRatio > 1)) throw new Error("The indel parent:child abundance ratio must be greater than 1.");
	}
	add(ordinal, rawSequence) {
		if (ordinal < 0 || ordinal >= this.records.length || this.processed[ordinal]) return;
		this.processed[ordinal] = 1;
		const record = this.records[ordinal];
		const sequence = normalizeNt(rawSequence);
		const partition = denoisePartition(record, this.options);
		if (!sequence || !partition) {
			if (this.options.unresolvedPolicy === "retain") this.standalone.push(ordinal);
			this.unresolvedRecords += 1;
			return;
		}
		if (sequence.includes("N")) {
			if (this.options.ambiguousPolicy === "exclude") {
				this.excludedAmbiguous += 1;
				return;
			}
			const key = `${partition}\u0000${sequence}`;
			const existing = this.retainedAmbiguous.get(key);
			const weight = Math.max(1, Math.floor(record.inputCount ?? 1));
			if (existing) {
				existing.members.push(ordinal);
				existing.count += weight;
			} else this.retainedAmbiguous.set(key, {
				representative: ordinal,
				members: [ordinal],
				count: weight
			});
			return;
		}
		const fingerprint = `${partition}\u0000${sequenceFingerprint(sequence)}`;
		const candidates = this.variantsByFingerprint.get(fingerprint) ?? [];
		let variantIndex = candidates.find((candidate) => this.arena.equals(this.variants[candidate].location, sequence));
		const weight = Math.max(1, Math.floor(record.inputCount ?? 1));
		if (variantIndex === void 0) {
			variantIndex = this.variants.length;
			this.variants.push({
				location: this.arena.append(sequence),
				partition,
				representative: ordinal,
				count: weight,
				target: variantIndex
			});
			candidates.push(variantIndex);
			this.variantsByFingerprint.set(fingerprint, candidates);
		} else this.variants[variantIndex].count += weight;
		this.variantByOrdinal[ordinal] = variantIndex;
	}
	finish(onProgress) {
		for (let ordinal = 0; ordinal < this.records.length; ordinal += 1) if (!this.processed[ordinal]) {
			if (this.options.unresolvedPolicy === "retain") this.standalone.push(ordinal);
			this.unresolvedRecords += 1;
		}
		const partitions = /* @__PURE__ */ new Map();
		this.variants.forEach((variant, index) => {
			const values = partitions.get(variant.partition);
			if (values) values.push(index);
			else partitions.set(variant.partition, [index]);
		});
		let candidateComparisons = 0;
		let truncated = 0;
		let indelMergedVariants = 0;
		let substitutionMergedVariants = 0;
		const variantWork = this.options.mode === "fad" ? this.variants.length + this.variants.reduce((count, variant) => count + (variant.count >= this.options.minimumParentCount ? 1 : 0), 0) : this.variants.length;
		let processedVariantWork = 0;
		const progressStride = Math.max(1, Math.floor(Math.max(1, variantWork) / 500));
		const variantProgress = () => {
			processedVariantWork += 1;
			if (processedVariantWork === variantWork || processedVariantWork % progressStride === 0) onProgress?.(processedVariantWork, Math.max(1, variantWork), "variants");
		};
		onProgress?.(0, Math.max(1, variantWork), "variants");
		for (const group of partitions.values()) {
			const result = this.options.mode === "fad" ? this.processFadPartition(group, variantProgress) : this.options.mode === "indel" ? this.processIndelPartition(group, variantProgress) : this.processConservativePartition(group, variantProgress);
			candidateComparisons += result.comparisons;
			truncated += result.truncated;
			indelMergedVariants += result.indelMergedVariants;
			substitutionMergedVariants += result.substitutionMergedVariants;
		}
		const representatives = new Int32Array(this.records.length);
		representatives.fill(-1);
		const counts = new Uint32Array(this.records.length);
		const finalizeTotal = Math.max(1, this.records.length * 3 + this.variants.length + this.standalone.length);
		let finalized = 0;
		const finalizeStride = Math.max(1, Math.floor(finalizeTotal / 500));
		const reportFinalize = (increment = 1) => {
			finalized += increment;
			if (finalized >= finalizeTotal || finalized % finalizeStride === 0) onProgress?.(Math.min(finalized, finalizeTotal), finalizeTotal, "finalize");
		};
		onProgress?.(0, finalizeTotal, "finalize");
		for (let ordinal = 0; ordinal < this.records.length; ordinal += 1) {
			const variantIndex = this.variantByOrdinal[ordinal];
			if (variantIndex >= 0) representatives[ordinal] = this.variants[this.variants[variantIndex].target].representative;
			reportFinalize();
		}
		for (const variant of this.variants) {
			const target = this.variants[variant.target];
			counts[target.representative] += variant.count;
			reportFinalize();
		}
		for (const group of this.retainedAmbiguous.values()) {
			counts[group.representative] = group.count;
			group.members.forEach((ordinal) => {
				representatives[ordinal] = group.representative;
			});
		}
		for (const ordinal of this.standalone) {
			counts[ordinal] = Math.max(1, Math.floor(this.records[ordinal].inputCount ?? 1));
			representatives[ordinal] = ordinal;
			reportFinalize();
		}
		let uniqueRecords = 0;
		for (const count of counts) {
			if (count > 0) uniqueRecords += 1;
			reportFinalize();
		}
		const largestGroups = largestCountGroups(counts);
		const warnings = [];
		if (this.excludedAmbiguous) warnings.push(`${this.excludedAmbiguous.toLocaleString()} records containing ambiguous nucleotide symbols were excluded to match the selected policy.`);
		if (this.unresolvedRecords) warnings.push(`${this.unresolvedRecords.toLocaleString()} records without a usable trimmed sequence or both V/J calls were ${this.options.unresolvedPolicy === "retain" ? "retained unchanged" : "discarded from the downstream representative set"}.`);
		if (truncated) warnings.push(`${truncated.toLocaleString()} variants reached the candidate cap; increase it before treating this denoising result as complete.`);
		let inputAbundance = 0;
		for (const record of this.records) {
			inputAbundance += Math.max(1, Math.floor(record.inputCount ?? 1));
			reportFinalize();
		}
		return {
			mode: this.options.mode,
			key: "trimmed",
			algorithm: this.options.mode === "fad" ? `FAD-compatible corrected 6-mer / method ${this.options.fadMethod}` : this.options.mode === "indel" ? "Indel-aware bounded edit model" : "Conservative exact-neighbor error model",
			inputRecords: this.records.length,
			inputAbundance,
			uniqueRecords,
			collapsedRecords: this.records.length - uniqueRecords,
			representatives,
			counts,
			largestGroups,
			partitions: partitions.size,
			candidateComparisons,
			indelMergedVariants,
			substitutionMergedVariants,
			excludedAmbiguous: this.excludedAmbiguous,
			unresolvedRecords: this.unresolvedRecords,
			warnings
		};
	}
	processFadPartition(group, onVariant) {
		const maximumSquared = Math.max(0, Math.floor(12 * this.options.fadNeighborThreshold + 1e-9));
		const blockCount = Math.max(1, maximumSquared + 1);
		const profiles = /* @__PURE__ */ new Map();
		const sequences = /* @__PURE__ */ new Map();
		for (const index of group) {
			const sequence = this.arena.decode(this.variants[index].location);
			sequences.set(index, sequence);
			profiles.set(index, kmerProfile(sequence, blockCount));
		}
		const ordered = [...group].sort((left, right) => this.variants[right].count - this.variants[left].count || this.variants[left].representative - this.variants[right].representative);
		const accepted = [];
		const index = /* @__PURE__ */ new Map();
		const candidates = /* @__PURE__ */ new Set();
		let comparisons = 0;
		let truncated = 0;
		const addAccepted = (variantIndex) => {
			accepted.push(variantIndex);
			profiles.get(variantIndex).hashes.forEach((hash, block) => {
				const key = `${block}:${hash}`;
				const values = index.get(key);
				if (values) values.push(variantIndex);
				else index.set(key, [variantIndex]);
			});
		};
		for (const variantIndex of ordered.filter((value) => this.variants[value].count >= this.options.minimumParentCount)) {
			const profile = profiles.get(variantIndex);
			candidates.clear();
			profile.hashes.forEach((hash, block) => {
				if (candidates.size >= this.options.maxCandidatesPerVariant) return;
				for (const candidate of index.get(`${block}:${hash}`) ?? []) {
					if (candidates.size >= this.options.maxCandidatesPerVariant) break;
					candidates.add(candidate);
				}
			});
			if (candidates.size >= this.options.maxCandidatesPerVariant) truncated += 1;
			const neighbors = [];
			for (const candidate of candidates) {
				comparisons += 1;
				const distance = kmerSquaredDistance(profile, profiles.get(candidate), maximumSquared);
				if (distance <= maximumSquared) neighbors.push({
					index: candidate,
					distance
				});
			}
			if (!neighbors.length) {
				addAccepted(variantIndex);
				onVariant?.();
				continue;
			}
			neighbors.sort((left, right) => this.variants[right.index].count - this.variants[left.index].count || left.distance - right.distance || this.variants[left.index].representative - this.variants[right.index].representative);
			const parent = neighbors[0].index;
			const child = this.variants[variantIndex];
			const lambda = this.variants[parent].count / Math.max(Number.MIN_VALUE, this.options.expectedZeroErrorFraction) * this.options.errorRate;
			const adjusted = Math.min(1, poissonStrictUpperTail(child.count, lambda) * (sequences.get(variantIndex)?.length ?? 1));
			if (this.options.fadMethod === 2 && adjusted < this.options.alpha) addAccepted(variantIndex);
			else child.target = parent;
			onVariant?.();
		}
		if (!accepted.length && ordered.length) addAccepted(ordered[0]);
		const acceptedSet = new Set(accepted);
		const vpTree = buildKmerVpTree(accepted, profiles, () => {
			comparisons += 1;
		});
		for (const variantIndex of ordered) {
			if (acceptedSet.has(variantIndex)) {
				this.variants[variantIndex].target = variantIndex;
				onVariant?.();
				continue;
			}
			this.variants[variantIndex].target = vpTree ? nearestKmerPoint(vpTree, variantIndex, profiles, (point) => this.variants[point].count, () => {
				comparisons += 1;
			}) : variantIndex;
			onVariant?.();
		}
		return {
			comparisons,
			truncated,
			indelMergedVariants: 0,
			substitutionMergedVariants: 0
		};
	}
	processIndelPartition(group, onVariant) {
		const distanceLimit = this.options.maximumEditDistance;
		const blockCount = distanceLimit + 1;
		const sequences = new Map(group.map((index) => [index, this.arena.decode(this.variants[index].location)]));
		const ordered = [...group].sort((left, right) => this.variants[right].count - this.variants[left].count || this.variants[left].representative - this.variants[right].representative);
		const parentIndex = /* @__PURE__ */ new Map();
		const shortParentsByLength = /* @__PURE__ */ new Map();
		const profileEdit = createBoundedEditProfiler(distanceLimit);
		let comparisons = 0;
		let truncated = 0;
		let indelMergedVariants = 0;
		let substitutionMergedVariants = 0;
		const candidates = /* @__PURE__ */ new Set();
		const addParent = (variantIndex) => {
			if (this.variants[variantIndex].count < this.options.minimumParentCount) return;
			const sequence = sequences.get(variantIndex);
			const segments = indexedEditSegments(sequence, blockCount);
			if (segments.some((segment) => segment.length === 0)) {
				const values = shortParentsByLength.get(sequence.length);
				if (values) values.push(variantIndex);
				else shortParentsByLength.set(sequence.length, [variantIndex]);
				return;
			}
			for (const segment of segments) {
				const key = `${sequence.length}:${segment.index}:${segment.value}`;
				const values = parentIndex.get(key);
				if (values) values.push(variantIndex);
				else parentIndex.set(key, [variantIndex]);
			}
		};
		for (const variantIndex of ordered) {
			const child = this.variants[variantIndex];
			const sequence = sequences.get(variantIndex);
			candidates.clear();
			let capped = false;
			const addCandidates = (values) => {
				for (const candidate of values) {
					if (candidates.has(candidate)) continue;
					if (candidates.size >= this.options.maxCandidatesPerVariant) {
						capped = true;
						return;
					}
					candidates.add(candidate);
				}
			};
			const minimumParentLength = Math.max(1, sequence.length - distanceLimit);
			const maximumParentLength = sequence.length + distanceLimit;
			for (let parentLength = minimumParentLength; parentLength <= maximumParentLength && !capped; parentLength += 1) {
				addCandidates(shortParentsByLength.get(parentLength) ?? []);
				for (let segmentIndex = 0; segmentIndex < blockCount && !capped; segmentIndex += 1) {
					const parentStart = Math.floor(segmentIndex * parentLength / blockCount);
					const segmentLength = Math.floor((segmentIndex + 1) * parentLength / blockCount) - parentStart;
					if (!segmentLength) continue;
					const firstStart = Math.max(0, parentStart - distanceLimit);
					const lastStart = Math.min(sequence.length - segmentLength, parentStart + distanceLimit);
					for (let queryStart = firstStart; queryStart <= lastStart && !capped; queryStart += 1) addCandidates(parentIndex.get(`${parentLength}:${segmentIndex}:${sequence.slice(queryStart, queryStart + segmentLength)}`) ?? []);
				}
			}
			if (capped) truncated += 1;
			let best = -1;
			let bestProfile = null;
			for (const candidate of candidates) {
				const parent = this.variants[candidate];
				if (parent.count <= child.count) continue;
				comparisons += 1;
				const profile = profileEdit(sequences.get(candidate), sequence);
				if (!profile || profile.distance < 1) continue;
				const indels = profile.insertions + profile.deletions;
				let plausible = false;
				if (indels > 0) plausible = parent.count / child.count >= this.options.minimumIndelParentRatio;
				else {
					const exactErrorProbability = (this.options.errorRate / 3) ** profile.substitutions * (1 - this.options.errorRate) ** Math.max(0, sequence.length - profile.substitutions);
					const lambda = parent.count * exactErrorProbability;
					plausible = Math.min(1, poissonStrictUpperTail(child.count, lambda) * alternativeCount(sequence.length, profile.substitutions)) >= this.options.alpha;
				}
				if (!plausible) continue;
				const bestParent = best >= 0 ? this.variants[best] : null;
				if (!bestProfile || profile.distance < bestProfile.distance || profile.distance === bestProfile.distance && profile.substitutions < bestProfile.substitutions || profile.distance === bestProfile.distance && profile.substitutions === bestProfile.substitutions && parent.count > (bestParent?.count ?? -1) || profile.distance === bestProfile.distance && profile.substitutions === bestProfile.substitutions && parent.count === bestParent?.count && parent.representative < (bestParent?.representative ?? Number.POSITIVE_INFINITY)) {
					best = candidate;
					bestProfile = profile;
				}
			}
			if (best >= 0 && bestProfile) {
				child.target = best;
				if (bestProfile.insertions + bestProfile.deletions > 0) indelMergedVariants += 1;
				else substitutionMergedVariants += 1;
			} else {
				child.target = variantIndex;
				addParent(variantIndex);
			}
			onVariant?.();
		}
		return {
			comparisons,
			truncated,
			indelMergedVariants,
			substitutionMergedVariants
		};
	}
	processConservativePartition(group, onVariant) {
		const distanceLimit = this.options.maximumHammingDistance;
		const blockCount = distanceLimit + 1;
		const sequences = new Map(group.map((index) => [index, this.arena.decode(this.variants[index].location)]));
		const ordered = [...group].sort((left, right) => this.variants[right].count - this.variants[left].count || this.variants[left].representative - this.variants[right].representative);
		const parentIndex = /* @__PURE__ */ new Map();
		let comparisons = 0;
		let truncated = 0;
		let substitutionMergedVariants = 0;
		const candidates = /* @__PURE__ */ new Set();
		const addParent = (variantIndex) => {
			if (this.variants[variantIndex].count < this.options.minimumParentCount) return;
			const sequence = sequences.get(variantIndex);
			sequenceBlocks(sequence, blockCount).forEach((block, blockIndex) => {
				const key = `${sequence.length}:${blockIndex}:${block}`;
				const values = parentIndex.get(key);
				if (values) values.push(variantIndex);
				else parentIndex.set(key, [variantIndex]);
			});
		};
		for (const variantIndex of ordered) {
			const child = this.variants[variantIndex];
			const sequence = sequences.get(variantIndex);
			candidates.clear();
			sequenceBlocks(sequence, blockCount).forEach((block, blockIndex) => {
				if (candidates.size >= this.options.maxCandidatesPerVariant) return;
				for (const candidate of parentIndex.get(`${sequence.length}:${blockIndex}:${block}`) ?? []) {
					if (candidates.size >= this.options.maxCandidatesPerVariant) break;
					candidates.add(candidate);
				}
			});
			if (candidates.size >= this.options.maxCandidatesPerVariant) truncated += 1;
			let best = -1;
			let bestLambda = -1;
			let bestDistance = Number.POSITIVE_INFINITY;
			for (const candidate of candidates) {
				comparisons += 1;
				const distance = hammingDistanceWithin(sequence, sequences.get(candidate), distanceLimit, false);
				if (distance < 1 || distance > distanceLimit) continue;
				const lambda = this.variants[candidate].count * ((this.options.errorRate / 3) ** distance * (1 - this.options.errorRate) ** Math.max(0, sequence.length - distance));
				if (Math.min(1, poissonStrictUpperTail(child.count, lambda) * alternativeCount(sequence.length, distance)) >= this.options.alpha && (lambda > bestLambda || lambda === bestLambda && distance < bestDistance)) {
					best = candidate;
					bestLambda = lambda;
					bestDistance = distance;
				}
			}
			if (best >= 0) {
				child.target = best;
				substitutionMergedVariants += 1;
			} else {
				child.target = variantIndex;
				addParent(variantIndex);
			}
			onVariant?.();
		}
		return {
			comparisons,
			truncated,
			indelMergedVariants: 0,
			substitutionMergedVariants
		};
	}
};
var UnionFind = class {
	parent;
	rank;
	constructor(size) {
		this.parent = new Int32Array(size);
		this.rank = new Uint8Array(size);
		for (let index = 0; index < size; index += 1) this.parent[index] = index;
	}
	find(value) {
		let root = value;
		while (this.parent[root] !== root) root = this.parent[root];
		while (this.parent[value] !== value) {
			const next = this.parent[value];
			this.parent[value] = root;
			value = next;
		}
		return root;
	}
	union(left, right) {
		let a = this.find(left);
		let b = this.find(right);
		if (a === b) return;
		if (this.rank[a] < this.rank[b]) [a, b] = [b, a];
		this.parent[b] = a;
		if (this.rank[a] === this.rank[b]) this.rank[a] += 1;
	}
};
function blocks(value, count) {
	const result = [];
	for (let index = 0; index < count; index += 1) {
		const start = Math.floor(index * value.length / count);
		const end = Math.floor((index + 1) * value.length / count);
		result.push({
			index,
			value: value.slice(start, end)
		});
	}
	return result;
}
function compatibleRecords(left, right, options) {
	if (datasetScopeKey(left, options.scope ?? "global") !== datasetScopeKey(right, options.scope ?? "global")) return false;
	if (options.requireSameLocus && left.locus !== right.locus) return false;
	return callsCompatible(left.vCall, right.vCall, options.callResolution, options.ambiguity) && callsCompatible(left.jCall, right.jCall, options.callResolution, options.ambiguity);
}
function recordIndexTokens(record, options) {
	const v = callSet(record.vCall, options.callResolution, options.ambiguity);
	const j = callSet(record.jCall, options.callResolution, options.ambiguity);
	if (!v.length || !j.length) return [];
	const scope = `${datasetScopeKey(record, options.scope ?? "global")}\u0001`;
	if (options.ambiguity === "strict") return [`${scope}${v.join("+")}\u0001${j.join("+")}`];
	return v.flatMap((vCall) => j.map((jCall) => `${scope}${vCall}\u0001${jCall}`));
}
function assignLineages(records, options, dedup, activeMask, doubleDMask) {
	const union = new UnionFind(records.length);
	const bucket = /* @__PURE__ */ new Map();
	const exact = /* @__PURE__ */ new Map();
	let candidateComparisons = 0;
	let truncatedCandidates = 0;
	let assignedRecords = 0;
	let activeAbundance = 0;
	for (let index = 0; index < records.length; index += 1) {
		const record = records[index];
		const weight = activeMask && !activeMask[index] ? 0 : dedup ? dedup.counts[index] : Math.max(1, Math.floor(records[index].inputCount ?? 1));
		if (!weight) continue;
		activeAbundance += weight;
		const cdr3 = normalizeNt(record.cdr3Nt);
		const tokens = recordIndexTokens(record, options);
		if (!cdr3 || !tokens.length || options.productiveOnly && !record.productive) continue;
		assignedRecords += weight;
		const distanceLimit = Math.floor((1 - options.identity) * cdr3.length + 1e-9);
		const blockCount = Math.max(1, Math.min(cdr3.length, distanceLimit + 1));
		const prefix = `${options.requireSameLocus ? record.locus : "*"}\u0000${cdr3.length}\u0000`;
		const exactSignature = `${prefix}${tokens.join("|")}\u0000${cdr3}`;
		const exactRepresentative = exact.get(exactSignature);
		if (exactRepresentative !== void 0) {
			union.union(index, exactRepresentative);
			continue;
		}
		exact.set(exactSignature, index);
		const candidates = /* @__PURE__ */ new Set();
		for (const token of tokens) {
			for (const block of blocks(cdr3, blockCount)) {
				const key = `${prefix}${token}\u0000${block.index}\u0000${block.value}`;
				for (const candidate of bucket.get(key) ?? []) {
					candidates.add(candidate);
					if (candidates.size >= options.maxCandidateComparisons) break;
				}
				if (candidates.size >= options.maxCandidateComparisons) break;
			}
			if (candidates.size >= options.maxCandidateComparisons) break;
		}
		if (candidates.size >= options.maxCandidateComparisons) truncatedCandidates += 1;
		for (const candidate of candidates) {
			candidateComparisons += 1;
			const other = records[candidate];
			if (!compatibleRecords(record, other, options)) continue;
			if (hammingDistanceWithin(cdr3, normalizeNt(other.cdr3Nt), distanceLimit, false) <= distanceLimit) union.union(index, candidate);
		}
		for (const token of tokens) for (const block of blocks(cdr3, blockCount)) {
			const key = `${prefix}${token}\u0000${block.index}\u0000${block.value}`;
			const values = bucket.get(key);
			if (values) values.push(index);
			else bucket.set(key, [index]);
		}
	}
	const eligible = new Uint8Array(records.length);
	const uniqueByRoot = new Uint32Array(records.length);
	const abundanceByRoot = new Float64Array(records.length);
	const representativeByRoot = new Int32Array(records.length);
	representativeByRoot.fill(-1);
	for (let index = 0; index < records.length; index += 1) {
		const weight = activeMask && !activeMask[index] ? 0 : dedup ? dedup.counts[index] : Math.max(1, Math.floor(records[index].inputCount ?? 1));
		if (!weight) continue;
		const record = records[index];
		if (!record.cdr3Nt || !recordIndexTokens(record, options).length || options.productiveOnly && !record.productive) continue;
		const root = union.find(index);
		eligible[index] = 1;
		uniqueByRoot[root] += 1;
		abundanceByRoot[root] += weight;
		if (representativeByRoot[root] < 0) representativeByRoot[root] = index;
	}
	const assignments = new Int32Array(records.length);
	const rootIds = new Int32Array(records.length);
	const top = [];
	const summaryLimit = 1e4;
	const sizeBins = new Uint32Array(7);
	const vUsageMap = /* @__PURE__ */ new Map();
	const jUsageMap = /* @__PURE__ */ new Map();
	let lineageCount = 0;
	const smaller = (left, right) => left.abundance < right.abundance || left.abundance === right.abundance && left.uniqueMembers < right.uniqueMembers || left.abundance === right.abundance && left.uniqueMembers === right.uniqueMembers && left.id > right.id;
	const retainTop = (summary) => {
		if (top.length < summaryLimit) {
			top.push(summary);
			let index = top.length - 1;
			while (index > 0) {
				const parent = Math.floor((index - 1) / 2);
				if (!smaller(top[index], top[parent])) break;
				[top[index], top[parent]] = [top[parent], top[index]];
				index = parent;
			}
			return;
		}
		if (!smaller(top[0], summary)) return;
		top[0] = summary;
		let index = 0;
		while (true) {
			const left = index * 2 + 1;
			const right = left + 1;
			let smallest = index;
			if (left < top.length && smaller(top[left], top[smallest])) smallest = left;
			if (right < top.length && smaller(top[right], top[smallest])) smallest = right;
			if (smallest === index) break;
			[top[index], top[smallest]] = [top[smallest], top[index]];
			index = smallest;
		}
	};
	const addUsage = (map, call, abundance) => {
		const value = map.get(call);
		if (value) {
			value.lineages += 1;
			value.abundance += abundance;
		} else map.set(call, {
			lineages: 1,
			abundance
		});
	};
	for (let root = 0; root < records.length; root += 1) {
		if (!uniqueByRoot[root]) continue;
		const id = ++lineageCount;
		rootIds[root] = id;
		const representativeOrdinal = representativeByRoot[root];
		const representative = records[representativeOrdinal];
		const abundance = abundanceByRoot[root];
		const bin = abundance === 1 ? 0 : abundance <= 3 ? 1 : abundance <= 9 ? 2 : abundance <= 24 ? 3 : abundance <= 99 ? 4 : abundance <= 499 ? 5 : 6;
		sizeBins[bin] += 1;
		const vCalls = callSet(representative.vCall, options.callResolution, options.ambiguity);
		const jCalls = callSet(representative.jCall, options.callResolution, options.ambiguity);
		vCalls.forEach((call) => addUsage(vUsageMap, call, abundance));
		jCalls.forEach((call) => addUsage(jUsageMap, call, abundance));
		retainTop({
			id,
			representativeOrdinal,
			uniqueMembers: uniqueByRoot[root],
			abundance,
			locus: representative.locus,
			vCalls,
			jCalls,
			cdr3Length: representative.cdr3Nt.length,
			studyScope: options.scope ?? "global",
			studyGroup: datasetScopeValue(representative, options.scope ?? "global"),
			sampleIds: [],
			sampleCounts: [],
			subjectIds: [],
			timepoints: [],
			compartments: [],
			doubleDPositiveMembers: 0,
			doubleDPositiveAbundance: 0
		});
	}
	for (let index = 0; index < records.length; index += 1) if (eligible[index]) assignments[index] = rootIds[union.find(index)];
	if (dedup) for (let index = 0; index < dedup.representatives.length; index += 1) {
		const representative = dedup.representatives[index];
		assignments[index] = representative >= 0 ? assignments[representative] : 0;
	}
	const summaryById = new Map(top.map((summary) => [summary.id, summary]));
	const samplesById = /* @__PURE__ */ new Map();
	const sampleCountsById = /* @__PURE__ */ new Map();
	const subjectsById = /* @__PURE__ */ new Map();
	const timepointsById = /* @__PURE__ */ new Map();
	const compartmentsById = /* @__PURE__ */ new Map();
	const addValue = (map, id, value) => {
		if (!value) return;
		const values = map.get(id);
		if (values) values.add(value);
		else map.set(id, new Set([value]));
	};
	for (let index = 0; index < records.length; index += 1) {
		const representative = dedup ? dedup.representatives[index] : index;
		if (representative < 0 || activeMask && !activeMask[representative]) continue;
		const sourceWeight = Math.max(1, Math.floor(records[index].inputCount ?? 1));
		const id = assignments[representative];
		if (!id || !summaryById.has(id)) continue;
		const record = records[index];
		addValue(samplesById, id, record.sampleId);
		if (record.sampleId) {
			let bySample = sampleCountsById.get(id);
			if (!bySample) {
				bySample = /* @__PURE__ */ new Map();
				sampleCountsById.set(id, bySample);
			}
			const counts = bySample.get(record.sampleId) ?? {
				representatives: /* @__PURE__ */ new Set(),
				abundance: 0
			};
			counts.representatives.add(representative);
			counts.abundance += sourceWeight;
			bySample.set(record.sampleId, counts);
		}
		addValue(subjectsById, id, record.subjectId);
		addValue(timepointsById, id, record.timepoint);
		addValue(compartmentsById, id, record.compartment);
	}
	for (let index = 0; index < records.length; index += 1) {
		const activeWeight = dedup ? dedup.counts[index] : Math.max(1, Math.floor(records[index].inputCount ?? 1));
		if (!activeWeight || activeMask && !activeMask[index] || !doubleDMask?.[index]) continue;
		const id = assignments[index];
		if (summaryById.has(id)) {
			const summary = summaryById.get(id);
			summary.doubleDPositiveMembers = (summary.doubleDPositiveMembers ?? 0) + 1;
			summary.doubleDPositiveAbundance = (summary.doubleDPositiveAbundance ?? 0) + activeWeight;
		}
	}
	for (const summary of top) {
		summary.sampleIds = [...samplesById.get(summary.id) ?? []].sort((a, b) => a.localeCompare(b, void 0, { numeric: true }));
		summary.sampleCounts = [...sampleCountsById.get(summary.id) ?? []].map(([sampleId, counts]) => ({
			sampleId,
			uniqueMembers: counts.representatives.size,
			abundance: counts.abundance
		})).sort((left, right) => right.abundance - left.abundance || left.sampleId.localeCompare(right.sampleId, void 0, { numeric: true }));
		summary.subjectIds = [...subjectsById.get(summary.id) ?? []].sort((a, b) => a.localeCompare(b, void 0, { numeric: true }));
		summary.timepoints = [...timepointsById.get(summary.id) ?? []].sort((a, b) => a.localeCompare(b, void 0, { numeric: true }));
		summary.compartments = [...compartmentsById.get(summary.id) ?? []].sort((a, b) => a.localeCompare(b, void 0, { numeric: true }));
	}
	const usage = (map) => [...map.entries()].map(([call, value]) => ({
		call,
		...value
	})).sort((a, b) => b.abundance - a.abundance || b.lineages - a.lineages || a.call.localeCompare(b.call));
	return {
		assignments,
		summaries: top.sort((a, b) => b.abundance - a.abundance || b.uniqueMembers - a.uniqueMembers || a.id - b.id),
		lineageCount,
		sizeHistogram: [
			"1",
			"2–3",
			"4–9",
			"10–24",
			"25–99",
			"100–499",
			"500+"
		].map((label, index) => ({
			label,
			count: sizeBins[index]
		})),
		vUsage: usage(vUsageMap),
		jUsage: usage(jUsageMap),
		assignedRecords,
		unassignedRecords: activeAbundance - assignedRecords,
		candidateComparisons,
		truncatedCandidates
	};
}
function parseFasta$2(text, aligned = false) {
	const records = [];
	let name = "";
	let sequence = [];
	const commit = () => {
		if (!name) return;
		const joined = sequence.join("").toUpperCase().replaceAll(".", "-").replace(/\s/g, "");
		const normalized = aligned ? joined.replace(/[^ACGTNRYKMSWBDHV-]/g, "N") : normalizeNt(joined);
		records.push({
			name,
			sequence: normalized
		});
	};
	for (const line of text.split(/\r?\n/)) if (line.startsWith(">")) {
		commit();
		name = line.slice(1).trim().split(/\s+/, 1)[0] ?? "";
		sequence = [];
	} else if (name) sequence.push(line);
	commit();
	return records;
}
function nameAliases(name) {
	const gene = name.split("|").filter(Boolean).find((field) => /(?:IGH|IGK|IGL|TRA|TRB|TRD|TRG)[VJ]/i.test(field));
	return [...new Set([
		name,
		gene ?? "",
		name.split(/\s+/, 1)[0] ?? ""
	].filter(Boolean))];
}
function prepareReferenceMsa(text) {
	const records = parseFasta$2(text, true);
	if (records.length < 2) throw new Error("CHMMAIRRa requires at least two aligned reference sequences.");
	const length = records[0].sequence.length;
	if (!length || records.some((record) => record.sequence.length !== length)) throw new Error("Reference MSA records must all have the same aligned length.");
	const lookup = /* @__PURE__ */ new Map();
	records.forEach((record, index) => nameAliases(record.name).forEach((alias) => lookup.set(alias, index)));
	return {
		names: records.map((record) => record.name),
		sequences: records.map((record) => record.sequence),
		lookup,
		length
	};
}
function topCallName(value) {
	return value.split(",", 1)[0]?.trim() ?? "";
}
function threadSequenceToMsa(sequenceAlignment, germlineAlignment, call, msa) {
	const callName = topCallName(call);
	let referenceIndex = msa.lookup.get(callName);
	if (referenceIndex === void 0) {
		const gene = callName.replace(/\*.*$/, "");
		referenceIndex = [...msa.lookup.entries()].find(([name]) => name.replace(/\*.*$/, "") === gene)?.[1];
	}
	if (referenceIndex === void 0) throw new Error(`Reference MSA is missing ${callName || "the assigned allele"}.`);
	const gappedFullReference = msa.sequences[referenceIndex];
	const degappedReference = gappedFullReference.replaceAll("-", "");
	const localGermline = germlineAlignment.toUpperCase().replaceAll(".", "-");
	const localQuery = sequenceAlignment.toUpperCase().replaceAll(".", "-").replace(/[^ACGT-]/g, "N");
	const degappedLocal = localGermline.replaceAll("-", "");
	const matchStart = degappedReference.indexOf(degappedLocal);
	if (matchStart < 0) throw new Error(`The ${callName} local germline alignment cannot be located in its MSA reference.`);
	const withoutInsertions = [];
	for (let index = 0; index < localGermline.length; index += 1) if (localGermline[index] !== "-") withoutInsertions.push(localQuery[index] ?? "N");
	const query = [..."N".repeat(matchStart), ...withoutInsertions];
	const threaded = new Array(msa.length).fill("-");
	let position = 0;
	for (let index = 0; index < msa.length; index += 1) {
		if (position >= query.length) break;
		if (gappedFullReference[index] !== "-") threaded[index] = query[position++] ?? "N";
	}
	return threaded.join("");
}
function chmmairraDistanceFromReference(query, germline) {
	const normalizedQuery = query.toUpperCase().replaceAll(".", "-");
	const normalizedGermline = germline.toUpperCase().replaceAll(".", "-");
	let mismatches = 0;
	for (let index = 0; index < Math.min(normalizedQuery.length, normalizedGermline.length); index += 1) if (normalizedQuery[index] !== normalizedGermline[index]) mismatches += 1;
	return mismatches;
}
function emission(reference, observation, mutation) {
	if (observation === "-" || observation === "N") return 1;
	return reference === observation ? 1 - mutation : mutation / 3;
}
function approximateEmissions(msa, observation, mutations) {
	const values = new Float64Array(msa.sequences.length * msa.length);
	for (let state = 0; state < msa.sequences.length; state += 1) for (let site = 0; site < msa.length; site += 1) values[state * msa.length + site] = emission(msa.sequences[state][site], observation[site], mutations[state]);
	return values;
}
function estimateMutationProbabilities(msa, observation, prior, baseMutation) {
	const states = msa.sequences.length;
	const length = msa.length;
	const mutations = new Float64Array(states);
	mutations.fill(baseMutation);
	if (states < 2 || length < 1) return mutations;
	const emissions = approximateEmissions(msa, observation, mutations);
	const alpha = new Float64Array(states * length);
	const beta = new Float64Array(states * length);
	const scale = new Float64Array(length);
	scale[0] = 1;
	const switchProbability = prior / length;
	const same = 1 - switchProbability;
	const different = switchProbability / (states - 1);
	for (let state = 0; state < states; state += 1) alpha[state * length] = emissions[state * length] / states;
	for (let site = 0; site < length - 1; site += 1) {
		let total = 0;
		for (let state = 0; state < states; state += 1) total += alpha[state * length + site];
		let nextTotal = 0;
		for (let state = 0; state < states; state += 1) {
			const value = ((total - alpha[state * length + site]) * different + alpha[state * length + site] * same) * emissions[state * length + site + 1];
			alpha[state * length + site + 1] = value;
			nextTotal += value;
		}
		scale[site + 1] = nextTotal ? 1 / nextTotal : 1;
		for (let state = 0; state < states; state += 1) alpha[state * length + site + 1] *= scale[site + 1];
	}
	for (let state = 0; state < states; state += 1) beta[state * length + length - 1] = scale[length - 1];
	for (let site = length - 2; site >= 0; site -= 1) {
		let total = 0;
		for (let state = 0; state < states; state += 1) total += beta[state * length + site + 1] * emissions[state * length + site + 1];
		for (let state = 0; state < states; state += 1) {
			const own = beta[state * length + site + 1] * emissions[state * length + site + 1];
			beta[state * length + site] = (different * (total - own) + same * own) * scale[site];
		}
	}
	const mutated = new Float64Array(states);
	const sameCount = new Float64Array(states);
	mutated.fill(2);
	sameCount.fill(10);
	for (let site = 0; site < length; site += 1) {
		if (observation[site] === "-" || observation[site] === "N") continue;
		let normalization = 0;
		for (let state = 0; state < states; state += 1) normalization += alpha[state * length + site] * beta[state * length + site];
		if (!normalization) continue;
		for (let state = 0; state < states; state += 1) {
			const posterior = alpha[state * length + site] * beta[state * length + site] / normalization;
			if (observation[site] === msa.sequences[state][site]) sameCount[state] += posterior;
			else mutated[state] += posterior;
		}
	}
	for (let state = 0; state < states; state += 1) mutations[state] = mutated[state] / (mutated[state] + sameCount[state]);
	return mutations;
}
function approximateProbability(msa, observation, prior, mutations) {
	const states = msa.sequences.length;
	if (states < 2) return 0;
	const emissions = approximateEmissions(msa, observation, mutations);
	const nonChimeric = new Float64Array(states);
	const chimeric = new Float64Array(states);
	const switchProbability = prior / msa.length;
	const same = 1 - switchProbability;
	const different = switchProbability / (states - 1);
	for (let state = 0; state < states; state += 1) nonChimeric[state] = emissions[state * msa.length] / states;
	for (let site = 0; site < msa.length - 1; site += 1) {
		let total = 0;
		for (let state = 0; state < states; state += 1) total += nonChimeric[state] + chimeric[state];
		let nextTotal = 0;
		for (let state = 0; state < states; state += 1) {
			const emissionValue = emissions[state * msa.length + site + 1];
			const own = nonChimeric[state] + chimeric[state];
			chimeric[state] = ((total - own) * different + chimeric[state] * same) * emissionValue;
			nonChimeric[state] = nonChimeric[state] * same * emissionValue;
			nextTotal += chimeric[state] + nonChimeric[state];
		}
		const scaling = nextTotal ? 1 / nextTotal : 1;
		for (let state = 0; state < states; state += 1) {
			chimeric[state] *= scaling;
			nonChimeric[state] *= scaling;
		}
	}
	let chimericTotal = 0;
	let total = 0;
	for (let state = 0; state < states; state += 1) {
		chimericTotal += chimeric[state];
		total += chimeric[state] + nonChimeric[state];
	}
	return total ? chimericTotal / total : 0;
}
function fullProbability(msa, observation, options) {
	const references = msa.sequences.length;
	const rates = options.mutationRates.length ? options.mutationRates : [.005];
	const rateCount = rates.length;
	const states = references * rateCount;
	if (references < 2) return 0;
	const nonChimeric = new Float64Array(states);
	const chimeric = new Float64Array(states);
	const nextNonChimeric = new Float64Array(states);
	const nextChimeric = new Float64Array(states);
	const switchProbability = options.priorProbability / msa.length;
	const self = 1 - switchProbability - options.mutationSwitchProbability;
	const differentReference = switchProbability / ((references - 1) * rateCount);
	const differentMutation = rateCount === 1 ? 0 : options.mutationSwitchProbability / (rateCount - 1);
	const emissionAt = (state, site) => {
		const reference = Math.floor(state / rateCount);
		const rate = state % rateCount;
		return emission(msa.sequences[reference][site], observation[site], rates[rate]);
	};
	for (let state = 0; state < states; state += 1) nonChimeric[state] = emissionAt(state, 0) / states;
	for (let site = 0; site < msa.length - 1; site += 1) {
		let total = 0;
		for (let state = 0; state < states; state += 1) total += nonChimeric[state] + chimeric[state];
		let nextTotal = 0;
		for (let reference = 0; reference < references; reference += 1) {
			let nonReference = 0;
			let chimReference = 0;
			for (let rate = 0; rate < rateCount; rate += 1) {
				const state = reference * rateCount + rate;
				nonReference += nonChimeric[state];
				chimReference += chimeric[state];
			}
			for (let rate = 0; rate < rateCount; rate += 1) {
				const state = reference * rateCount + rate;
				const emit = emissionAt(state, site + 1);
				nextNonChimeric[state] = (nonChimeric[state] * self + (nonReference - nonChimeric[state]) * differentMutation) * emit;
				nextChimeric[state] = (chimeric[state] * self + (chimReference - chimeric[state]) * differentMutation + (total - chimReference - nonReference) * differentReference) * emit;
				nextTotal += nextNonChimeric[state] + nextChimeric[state];
			}
		}
		const scaling = nextTotal ? 1 / nextTotal : 1;
		for (let state = 0; state < states; state += 1) {
			nonChimeric[state] = nextNonChimeric[state] * scaling;
			chimeric[state] = nextChimeric[state] * scaling;
		}
	}
	let chimericTotal = 0;
	let total = 0;
	for (let state = 0; state < states; state += 1) {
		chimericTotal += chimeric[state];
		total += chimeric[state] + nonChimeric[state];
	}
	return total ? chimericTotal / total : 0;
}
function fullViterbi(msa, observation, options) {
	const references = msa.sequences.length;
	const rates = options.mutationRates.length ? options.mutationRates : [.005];
	const rateCount = rates.length;
	const states = references * rateCount;
	const length = msa.length;
	if (references < 2) return {
		startingReference: msa.names[0] ?? "",
		recombinations: [],
		referencePath: options.tracePath ? new Int32Array(length) : void 0
	};
	const previous = new Float64Array(states);
	const next = new Float64Array(states);
	const from = new Int32Array(states * length);
	const switchProbability = options.priorProbability / length;
	const logSelf = Math.log(1 - switchProbability - options.mutationSwitchProbability);
	const logDifferentReference = Math.log(switchProbability / ((references - 1) * rateCount));
	const logDifferentMutation = rateCount === 1 ? Number.NEGATIVE_INFINITY : Math.log(options.mutationSwitchProbability / (rateCount - 1));
	const emissionAt = (state, site) => emission(msa.sequences[Math.floor(state / rateCount)][site], observation[site], rates[state % rateCount]);
	for (let state = 0; state < states; state += 1) {
		previous[state] = Math.log(emissionAt(state, 0) / states);
		from[state * length] = state;
	}
	for (let site = 0; site < length - 1; site += 1) {
		const bestInReference = new Int32Array(references);
		for (let reference = 0; reference < references; reference += 1) {
			let best = reference * rateCount;
			for (let rate = 1; rate < rateCount; rate += 1) {
				const state = reference * rateCount + rate;
				if (previous[state] > previous[best]) best = state;
			}
			bestInReference[reference] = best;
		}
		let bestReference = 0;
		let secondReference = references > 1 ? 1 : 0;
		if (previous[bestInReference[secondReference]] > previous[bestInReference[bestReference]]) [bestReference, secondReference] = [secondReference, bestReference];
		for (let reference = 2; reference < references; reference += 1) if (previous[bestInReference[reference]] > previous[bestInReference[bestReference]]) {
			secondReference = bestReference;
			bestReference = reference;
		} else if (previous[bestInReference[reference]] > previous[bestInReference[secondReference]]) secondReference = reference;
		for (let reference = 0; reference < references; reference += 1) {
			const outside = bestInReference[reference === bestReference ? secondReference : bestReference];
			for (let rate = 0; rate < rateCount; rate += 1) {
				const state = reference * rateCount + rate;
				let selected = state;
				let score = previous[state] + logSelf;
				const referenceScore = previous[outside] + logDifferentReference;
				if (referenceScore > score) {
					selected = outside;
					score = referenceScore;
				}
				if (rateCount > 1) {
					let otherRate = reference * rateCount + (rate === 0 ? 1 : 0);
					for (let candidateRate = 0; candidateRate < rateCount; candidateRate += 1) {
						const candidate = reference * rateCount + candidateRate;
						if (candidate !== state && previous[candidate] > previous[otherRate]) otherRate = candidate;
					}
					const mutationScore = previous[otherRate] + logDifferentMutation;
					if (mutationScore > score) {
						selected = otherRate;
						score = mutationScore;
					}
				}
				from[state * length + site + 1] = selected;
				next[state] = score + Math.log(emissionAt(state, site + 1));
			}
		}
		previous.set(next);
	}
	let current = 0;
	for (let state = 1; state < states; state += 1) if (previous[state] > previous[current]) current = state;
	const recombinations = [];
	const referencePath = options.tracePath ? new Int32Array(length) : void 0;
	if (referencePath) referencePath[length - 1] = Math.floor(current / rateCount);
	for (let site = length - 1; site >= 1; site -= 1) {
		const parent = from[current * length + site];
		const leftReference = Math.floor(parent / rateCount);
		const rightReference = Math.floor(current / rateCount);
		if (leftReference !== rightReference) recombinations.push({
			position: site + 1,
			left: msa.names[leftReference],
			right: msa.names[rightReference]
		});
		current = parent;
		if (referencePath) referencePath[site - 1] = Math.floor(current / rateCount);
	}
	return {
		startingReference: msa.names[Math.floor(current / rateCount)],
		recombinations: recombinations.reverse(),
		referencePath
	};
}
function approximateViterbi(msa, observation, prior, mutations, tracePath = false) {
	const states = msa.sequences.length;
	const length = msa.length;
	const emissions = approximateEmissions(msa, observation, mutations);
	const previous = new Float64Array(states);
	const next = new Float64Array(states);
	const from = new Int32Array(states * length);
	const same = Math.log(1 - prior / length);
	const different = Math.log(prior / length / (states - 1));
	for (let state = 0; state < states; state += 1) {
		previous[state] = Math.log(emissions[state * length] / states);
		from[state * length] = state;
	}
	for (let site = 0; site < length - 1; site += 1) {
		let best = 0;
		for (let state = 1; state < states; state += 1) if (previous[state] > previous[best]) best = state;
		for (let state = 0; state < states; state += 1) {
			const stay = previous[state] + same;
			const change = previous[best] + different;
			const remain = state === best || stay > change;
			from[state * length + site + 1] = remain ? state : best;
			next[state] = (remain ? stay : change) + Math.log(emissions[state * length + site + 1]);
		}
		previous.set(next);
	}
	let current = 0;
	for (let state = 1; state < states; state += 1) if (previous[state] > previous[current]) current = state;
	const reverse = [];
	const referencePath = tracePath ? new Int32Array(length) : void 0;
	if (referencePath) referencePath[length - 1] = current;
	for (let site = length - 1; site >= 1; site -= 1) {
		const parent = from[current * length + site];
		if (parent !== current) reverse.push({
			position: site + 1,
			left: msa.names[parent],
			right: msa.names[current]
		});
		current = parent;
		if (referencePath) referencePath[site - 1] = current;
	}
	return {
		startingReference: msa.names[current],
		recombinations: reverse.reverse(),
		referencePath
	};
}
function runChmm(msa, threadedObservation, localSequenceAlignment, localGermlineAlignment, options) {
	if (threadedObservation.length !== msa.length) throw new Error("Threaded query and reference MSA lengths differ.");
	const dfr = chmmairraDistanceFromReference(localSequenceAlignment, localGermlineAlignment);
	if (options.method === "DB") {
		const path = options.detailed ? fullViterbi(msa, threadedObservation, options) : {
			startingReference: "",
			recombinations: []
		};
		return {
			probability: fullProbability(msa, threadedObservation, options),
			dfr,
			...path
		};
	}
	const mutations = estimateMutationProbabilities(msa, threadedObservation, options.priorProbability, options.baseMutationProbability);
	return {
		probability: approximateProbability(msa, threadedObservation, options.priorProbability, mutations),
		dfr,
		...options.detailed ? approximateViterbi(msa, threadedObservation, options.priorProbability, mutations, options.tracePath) : {
			startingReference: "",
			recombinations: []
		}
	};
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
		batchRecords: 2e3,
		callingProfile: "truth_optimized",
		assignerStrategy: "aer",
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
	annotation.batchRecords = Math.max(1, Math.floor(finite(annotation.batchRecords, 2e3)));
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
const SOURCE_TRANSFERRED_IMGT = 3;
const SOURCE_VALIDATED_J_MOTIF = 4;
const SOURCE_PROVIDED = 5;
const SOURCE_TRANSFERRED_J = 6;
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
function nearestTemplates(query, templates, limit = 12) {
	const queryKmers = kmerSet(query);
	return templates.map((template) => {
		const templateKmers = kmerSet(template[1]);
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
function templateCandidateTiers(queryName, query, templates, eligible = hasRegionMetadata) {
	const delineated = templates.filter((template) => eligible(template[2]));
	const alleleName = canonicalAllele(queryName);
	const geneName = canonicalGene(queryName);
	const exactAllele = delineated.filter(([name]) => canonicalAllele(name) === alleleName);
	const exactGene = delineated.filter(([name]) => canonicalGene(name) === geneName);
	const preferred = exactAllele.length ? exactAllele : exactGene;
	const preferredKeys = new Set(preferred.map(([name, sequence]) => `${name}\u0000${sequence}`));
	return [preferred, nearestTemplates(query, delineated.filter(([name, sequence]) => !preferredKeys.has(`${name}\u0000${sequence}`)))].filter((tier) => tier.length);
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
function transferMetadata(name, sequence, templates) {
	for (const candidates of templateCandidateTiers(name, sequence, templates)) {
		let selected;
		for (const template of candidates) {
			const alignment = globalAlignment(sequence, template[1]);
			const named = canonicalGene(template[0]) === canonicalGene(name);
			if (alignment.identity < (named ? .8 : .72)) continue;
			const templateMetadata = template[2];
			const templateBounds = templateMetadata.slice(2, 12);
			const mapped = mapReferenceCoordinates(alignment, template[1].length);
			const bounds = templateBounds.map((boundary) => mapped[boundary]);
			if (bounds.some((value, index) => value < 0 || value > sequence.length || index && value < bounds[index - 1])) continue;
			let valid = true;
			for (let index = 0; index < bounds.length; index += 2) if (bounds[index + 1] <= bounds[index]) valid = false;
			if (!valid) continue;
			const templateFrame = templateMetadata[0] >= 0 ? templateMetadata[0] : templateBounds[0] % 3;
			const frame = positiveModulo(bounds[0] + templateFrame - templateBounds[0], 3);
			const anchorEnd = nearestFrameCysEnd(sequence, bounds[9], frame, 24);
			if (!anchorEnd || anchorEnd <= bounds[8]) continue;
			bounds[9] = anchorEnd;
			const metadata = compactMetadata(frame, -1, bounds, SOURCE_TRANSFERRED_IMGT);
			if (!validateMetadata(metadata, sequence.length, "V")) continue;
			if (!selected || alignment.identity > selected.identity) selected = {
				metadata,
				identity: alignment.identity
			};
		}
		if (selected) return selected.metadata;
	}
}
function transferJMetadata(name, sequence, templates) {
	for (const candidates of templateCandidateTiers(name, sequence, templates, hasJMetadata)) {
		let selected;
		for (const template of candidates) {
			const alignment = globalAlignment(sequence, template[1]);
			const named = canonicalGene(template[0]) === canonicalGene(name);
			if (alignment.identity < (named ? .75 : .68)) continue;
			const referenceAnchor = template[2][1] + 1;
			if (referenceAnchor < 0 || referenceAnchor + 6 > template[1].length) continue;
			const mapped = mapReferenceCoordinates(alignment, template[1].length);
			const anchor = mapped[referenceAnchor];
			const anchorEnd = mapped[referenceAnchor + 6];
			if (anchor < 0 || anchorEnd - anchor !== 6 || !isJAnchor(sequence, anchor)) continue;
			const transferred = compactMetadata(anchor % 3, anchor - 1, EMPTY_BOUNDS, SOURCE_TRANSFERRED_J);
			if (!validateMetadata(transferred, sequence.length, "J")) continue;
			if (!selected || alignment.identity > selected.identity) selected = {
				metadata: transferred,
				identity: alignment.identity
			};
		}
		if (selected) return selected.metadata;
	}
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
			const normalized = normalizeIndexSequence(input.rawSequence);
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
function preprocessGermlineFasta(text, segment, templates = [], allowedLoci = LOCI) {
	const inputRecords = parseFasta(text);
	const records = [];
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
		const normalized = normalizeSequence(input.rawSequence);
		if (!normalized.sequence) throw new Error(`${name} has an empty nucleotide sequence.`);
		ambiguousBases += normalized.ambiguous;
		const fields = parseImgtFields(input.header);
		let metadata = metadataFromHeader(input.header);
		if (metadata && !validateMetadata(metadata, normalized.sequence.length, segment)) throw new Error(`${name} contains invalid SWIGMETA coordinates.`);
		if (metadata) metadata = [...metadata.slice(0, 12), metadata[12] || SOURCE_PROVIDED];
		if (!metadata && segment === "V") metadata = imgtVMetadata(input.rawSequence, fields);
		if (!metadata && segment === "V" && templates.length) metadata = transferMetadata(name, normalized.sequence, templates);
		if (!metadata && segment === "J" && templates.length) metadata = transferJMetadata(name, normalized.sequence, templates);
		if (!metadata && segment === "J") metadata = jMetadata(normalized.sequence, fields);
		if (!metadata && segment === "D") {
			const frame = imgtFrame(fields);
			if (frame >= 0) metadata = compactMetadata(frame, -1, EMPTY_BOUNDS, SOURCE_IMGT_GAPPED);
		}
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
		warnings
	};
}
/**
* Apply the same progressively broadened template search used by the browser
* worker. Existing valid SWIGMETA is retained, while records still lacking V
* or J metadata are offered to each successive template tier.
*/
function preprocessGermlineFastaAcrossTiers(text, segment, templateTiers, allowedLoci = LOCI) {
	let report;
	for (const templates of templateTiers.length ? templateTiers : [[]]) {
		report = preprocessGermlineFasta(report?.fasta ?? text, segment, [...templates], allowedLoci);
		if (segment !== "V" && segment !== "J") break;
		if (report.annotated === report.count) break;
	}
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
const VERSION = "0.33.0";
const CLI_DIRECTORY = dirname(fileURLToPath(import.meta.url));
function defaultCliAssets() {
	const directory = join(CLI_DIRECTORY, "assets");
	return {
		wasmPath: join(directory, "swiftig.wasm"),
		referencePackPath: join(directory, "imgt-reference-pack.json.gz")
	};
}
function usage() {
	return `swig-cli ${VERSION}\n\nRun a complete non-phylogenetic Swig pipeline:\n  swig-cli run reads.fastq.gz --out swig-output\n  swig-cli run --config swig.config.json [--out DIRECTORY] [--workers N]\n\nRun only streaming V(D)J assignment (AIRR outfmt 19):\n  swig-cli --vdj -query reads.fasta -germline_db_V V.fasta -germline_db_D D.fasta \\\n    -germline_db_J J.fasta -out calls.airr.tsv\n\nCreate an editable config:\n  swig-cli init swig.config.json\n\nSingle-input metadata options:\n  --sample SAMPLE_ID  --donor SUBJECT_ID  --dataset DATASET_ID\n\nSamples with the same subjectId/--donor are treated as the same donor.\nLineage phylogenetics is intentionally not run by swig-cli.`;
}
function vdjUsage() {
	return `swig-cli ${VERSION} --vdj\n\nLow-overhead, streaming SwiftIG V(D)J assignment with IgBLAST-style option names.\n\nRequired:\n  -germline_db_V FASTA  -germline_db_J FASTA  -out AIRR_TSV\n\nInput and optional references:\n  -query FASTA            Query FASTA or '-' for stdin (default '-')\n  -germline_db_D FASTA    D germline FASTA\n  -c_region_db FASTA      Constant-region FASTA\n\nAnnotation modes (default: assignments only; CDR/FWR fields remain empty):\n  -custom_internal_data FILE  IgBLAST V .ndm.imgt data (1-based inclusive intervals)\n  -auxiliary_data FILE        IgBLAST J .aux data (0-based frame/CDR3 stop)\n  -d_frame_data FILE          IgBLAST D frame-one starts\n  --swigannots                Infer/validate metadata as in Swig Web\n\nExecution:\n  -num_threads N          Worker count; --workers N overrides it\n  --workers N             Exact workers, or 0 for automatic\n  --batch-records N       Records per bounded WASM batch (default 2000)\n  -strand both|plus|minus -outfmt 19 -organism NAME -ig_seqtype Ig|TCR\n\nThe germline options take FASTA files, not makeblastdb binary prefixes. Output is SwiftIG AIRR,\nnot IgBLAST pairwise/tabular formatting. The output path is mandatory and is written incrementally.`;
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
function parseVdjArguments(rawArgs) {
	const aliases = new Map([
		["--query", "-query"],
		["--output", "-out"],
		["--out", "-out"],
		["--germline-db-v", "-germline_db_V"],
		["--germline-db-d", "-germline_db_D"],
		["--germline-db-j", "-germline_db_J"],
		["--c-region-db", "-c_region_db"]
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
		"-J_penalty"
	]);
	const flags = new Set(["--swigannots", "-show_translation"]);
	const options = {};
	for (let index = 0; index < rawArgs.length; index += 1) {
		let token = rawArgs[index];
		if (token === "--vdj") continue;
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
		if (!valued.has(token)) throw new Error(`Unsupported --vdj option ${token}.\n\n${vdjUsage()}`);
		const value = inline !== void 0 ? inline : rawArgs[++index];
		if (value === void 0) throw new Error(`${token} requires a value.`);
		options[token] = value;
	}
	return options;
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
	const raw = input.path === "-" ? process.stdin : createReadStream(input.path, range ? {
		start: range.start,
		end: range.end - 1
	} : void 0);
	return createInterface({
		input: /\.gz$/i.test(input.path) || range ? raw.pipe(createGunzip()) : raw,
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
		this.pending = /* @__PURE__ */ new Map();
		this.nextId = 1;
		this.nextWorker = 0;
	}
	async start() {
		for (let index = 0; index < this.size; index += 1) {
			const webWorker = Boolean(process.versions.bun && globalThis.Worker);
			const worker = webWorker ? new globalThis.Worker(new URL("./swig-worker.js", import.meta.url)) : new Worker(new URL("./swig-worker.mjs", import.meta.url));
			const receive = (message) => {
				const pending = this.pending.get(message.id);
				if (!pending) return;
				this.pending.delete(message.id);
				if (message.error) pending.reject(new Error(message.error));
				else pending.resolve(message.result);
			};
			const fail = (error) => {
				for (const [id, pending] of this.pending) if (pending.worker === worker) {
					this.pending.delete(id);
					pending.reject(error instanceof Error ? error : new Error(error?.message ?? String(error)));
				}
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
	}
	request(worker, message) {
		const id = this.nextId++;
		return new Promise((resolvePromise, reject) => {
			this.pending.set(id, {
				resolve: resolvePromise,
				reject,
				worker
			});
			worker.postMessage({
				id,
				...message
			});
		});
	}
	run(message) {
		const worker = this.workers[this.nextWorker++ % this.workers.length];
		return this.request(worker, {
			type: "annotate",
			...message
		});
	}
	async close() {
		for (const worker of this.workers) await worker.terminate();
		this.workers = [];
	}
};
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
async function prepareVdjReferences(options, assets) {
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
			const report = preprocessGermlineFastaAcrossTiers(raw[segment], segment, germlineTemplateTiers(pack, species, scope, segment), allowedLoci);
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
		const report = preprocessGermlineFastaAcrossTiers(raw.V, "V", germlineTemplateTiers(pack, species, scope, "V"), allowedLoci);
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
function runChimera(rows, activeMask, config, scope, headers, references) {
	const selectedSegment = config.segment.toUpperCase();
	const msaText = config.uploadedMsa?.trim() || (config.msaSource === "selected" ? references[selectedSegment]?.trim() : "");
	if (!msaText) throw new Error("CLI chimera filtering requires either pipeline.chimera.uploadedMsa or an aligned selected-segment reference.");
	let msa;
	try {
		msa = prepareReferenceMsa(msaText);
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
	const cache = /* @__PURE__ */ new Map();
	let evaluated = 0, flagged = 0, unevaluated = 0;
	addHeaders(headers, ["swig_chimera_probability", "swig_chimera_status"]);
	if (config.detailed) addHeaders(headers, ["swig_chimera_starting_reference", "swig_chimera_recombinations"]);
	for (let ordinal = 0; ordinal < rows.length; ordinal += 1) {
		if (!activeMask[ordinal]) continue;
		const row = rows[ordinal], call = row[`${segment}_call`] ?? "", sequence = row[`${segment}_sequence_alignment`] ?? "", germline = row[`${segment}_germline_alignment`] ?? "";
		if (!call || !sequence || !germline) {
			row.swig_chimera_status = "missing_alignment";
			unevaluated += 1;
			if (!config.retainUnevaluated) activeMask[ordinal] = 0;
			continue;
		}
		if (chmmairraDistanceFromReference(sequence, germline) < config.minimumDfr) {
			row.swig_chimera_status = "low_dfr";
			unevaluated += 1;
			if (!config.retainUnevaluated) activeMask[ordinal] = 0;
			continue;
		}
		const key = `${call}\0${sequence}\0${germline}`;
		try {
			let result = cache.get(key);
			if (result === void 0) {
				result = runChmm(msa, threadSequenceToMsa(sequence, germline, call, msa), sequence, germline, options);
				cache.set(key, result);
			}
			row.swig_chimera_probability = String(result.probability);
			row.swig_chimera_status = "evaluated";
			evaluated += 1;
			if (config.detailed) {
				row.swig_chimera_starting_reference = result.startingReference;
				row.swig_chimera_recombinations = result.recombinations.map((event) => `${event.left}->${event.right}@${event.position}`).join(";");
			}
			if (result.probability >= config.posteriorThreshold) {
				activeMask[ordinal] = 0;
				flagged += 1;
			}
		} catch (error) {
			row.swig_chimera_status = "error";
			unevaluated += 1;
			if (!config.retainUnevaluated) activeMask[ordinal] = 0;
		}
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
async function finishWritable(stream) {
	stream.end();
	await once(stream, "finish");
}
async function writeRowsFile(path, headers, rows, include = () => true) {
	const stream = createWriteStream(path);
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
	const stream = createWriteStream(path);
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
	const shmSummaries = shm?.lineages.flatMap((group) => {
		const match = /^Lineage\s+(\d+)$/.exec(group.label);
		return match ? [{
			lineageId: Number(match[1]),
			mean: group.mean,
			p95: group.p95 ?? 0
		}] : [];
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
	const pool = config.inputs.some((input) => detectFormat(input.path, input.format) !== "airr" || config.annotation.airrMode === "reannotate") ? new WasmPool(config.annotation.workers, {
		wasmPath: assets.wasmPath,
		references,
		callingProfile: config.annotation.callingProfile,
		assignerStrategy: config.annotation.assignerStrategy
	}) : null;
	if (pool) await pool.start();
	const annotatedStream = config.output.writeAnnotatedAirr ? createWriteStream(join(outputDirectory, `${prefix}.annotated.airr.tsv`)) : null;
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
			for await (const batch of sequenceBatches(input, config.annotation.batchRecords, config.preprocessing, config.annotation.airrMode, datasetIndex, preprocessingState)) {
				const promise = batch.format === "airr" && config.annotation.airrMode === "preserve" ? Promise.resolve({
					direct: true,
					header: batch.header,
					body: batch.body,
					doubleDHeader: "",
					doubleDBody: ""
				}) : pool.run({
					text: batch.text,
					count: batch.count,
					format: batch.format === "fasta" ? 1 : batch.format === "fastq" ? 2 : 3,
					minimumIdentity: config.annotation.minimumIdentity,
					strand: config.annotation.strand,
					doubleD: config.annotation.doubleD
				});
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
		if (config.pipeline.collapse.mode === "exact") collapseResult = deduplicate(records, config.pipeline.collapse.key, config.pipeline.collapse.unresolvedPolicy, config.pipeline.collapse.scope, config.pipeline.collapse.respectConstantCall);
		else {
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
			collapseResult = accumulator.finish();
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
		chimeraSummary = runChimera(rows, activeMask, config.pipeline.chimera, config.references.scope, headers, references);
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
	const batchRecords = options["--batch-records"] === void 0 ? 2e3 : parseIntegerOption(options["--batch-records"], "--batch-records", {
		minimum: 1,
		allowZero: false
	});
	const threadValue = options["--workers"] ?? options["-num_threads"];
	let workers = threadValue === void 0 ? Math.max(1, Math.min(4, availableParallelism())) : parseIntegerOption(threadValue, options["--workers"] !== void 0 ? "--workers" : "-num_threads", { minimum: 0 });
	if (options["--workers"] === void 0 && threadValue !== void 0 && workers === 0) throw new Error("-num_threads must be at least 1; use --workers 0 for automatic selection.");
	if (workers === 0) workers = Math.max(1, Math.min(8, availableParallelism()));
	const assigner = String(options["--assigner"] ?? "aer");
	if (![
		"standard",
		"riat_mp",
		"aer"
	].includes(assigner)) throw new Error("--assigner must be standard, riat_mp, or aer.");
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
	const pool = new WasmPool(workers, {
		wasmPath: assets.wasmPath,
		references: prepared.references,
		callingProfile,
		assignerStrategy: assigner,
		tuning
	});
	await pool.start();
	const output = outputPath === "-" ? process.stdout : createWriteStream(outputPath);
	let outputHeader = null;
	let records = 0;
	let completed = false;
	try {
		if (outputPath !== "-") await once(output, "open");
		process.stderr.write(`Streaming SwiftIG V(D)J assignments (${prepared.mode}; ${workers} worker${workers === 1 ? "" : "s"}) to ${outputPath}.\n`);
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
			pending.push({ promise: pool.run({
				text: batch.text,
				count: batch.count,
				format: 1,
				minimumIdentity,
				strand,
				doubleD: { mode: "off" }
			}) });
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
	if (args.includes("--vdj")) {
		await runVdj(args, assets);
		return;
	}
	const command = args[0] && !args[0].startsWith("-") ? args[0] : "run";
	const rest = command === args[0] ? args.slice(1) : args;
	if (hasFlag(args, "--help") || command === "help") {
		process.stdout.write(`${usage()}\n`);
		return;
	}
	if (hasFlag(args, "--version") || command === "version") {
		process.stdout.write(`${VERSION}\n`);
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
