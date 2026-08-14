import assert from "node:assert/strict";
import fs from "node:fs";
import test from "node:test";
import zlib from "node:zlib";

import { preprocessGermlineFasta, type MetadataAllele } from "../src/germline-preprocess.ts";
import { composeReferenceOverrides, referenceCellKey, segmentAppliesToLocus } from "../src/reference-composition.ts";
import {
  collectionsFor,
  collectionsForDatabase,
  databaseOptionsFor,
  databasesForCell,
  DEFAULT_DATABASE_ID,
  preferredDatabaseIdFor,
  REFERENCE_DATABASES,
} from "../src/reference-catalog.ts";

const pack = JSON.parse(zlib.gunzipSync(fs.readFileSync(new URL("../public/references/imgt-202632-7-swig-0.7.json.gz", import.meta.url))).toString());
const human = pack.species.find((entry: { name: string }) => entry.name === "Homo sapiens");
const humanIghV = human.loci.IGH.V as MetadataAllele[];
const humanIghJ = human.loci.IGH.J as MetadataAllele[];

function imgtTemplateTiers(speciesName: string, segment: "V" | "J"): MetadataAllele[][] {
  const selected = pack.species.find((entry: { name: string }) => entry.name === speciesName);
  assert.ok(selected);
  const baseTaxon = speciesName.split("_", 1)[0];
  const genus = baseTaxon.split(" ", 1)[0];
  return [
    [selected],
    pack.species.filter((entry: { name: string }) => entry.name !== speciesName && entry.name.split("_", 1)[0] === baseTaxon),
    pack.species.filter((entry: { name: string }) => entry.name.split("_", 1)[0] !== baseTaxon && entry.name.startsWith(`${genus} `)),
    pack.species.filter((entry: { name: string }) => !entry.name.startsWith(`${genus} `)),
  ].map((group) => group.flatMap((entry: { loci: { IGH?: Record<string, MetadataAllele[]> } }) => entry.loci.IGH?.[segment] ?? [])).filter((group) => group.length);
}

function firstAnnotated(alleles: MetadataAllele[], segment: "V" | "J"): MetadataAllele {
  const allele = alleles.find((entry) => segment === "V"
    ? entry[2]?.slice(2, 12).every((value) => value >= 0)
    : Boolean(entry[2] && entry[2]![0] >= 0 && entry[2]![1] >= 0));
  assert.ok(allele);
  return allele;
}

test("V delineation transfers from the nearest annotated relative without using functionality labels", () => {
  const template = firstAnnotated(humanIghV, "V");
  const changed = `${template[1].slice(0, 145)}${template[1][145] === "A" ? "C" : "A"}${template[1].slice(146)}`;
  const report = preprocessGermlineFasta(`>${template[0]}_S9999 pseudogene\n${changed}\n`, "V", humanIghV, ["IGH"]);
  assert.equal(report.annotated, 1);
  assert.equal(report.transferred, 1);
  assert.match(report.fasta, /SWIGMETA=[^\n]*,3\n/);
});

test("V coordinate transfer retains boundaries across a sequence-level frameshift", () => {
  const template = firstAnnotated(humanIghV, "V");
  const changed = `${template[1].slice(0, 130)}${template[1].slice(131)}`;
  const report = preprocessGermlineFasta(`>${template[0]}_S9998\n${changed}\n`, "V", humanIghV, ["IGH"]);
  assert.equal(report.annotated, 1);
  assert.equal(report.transferred, 1);
});

test("a V sequence that does not span the mapped IMGT intervals is not assigned fabricated boundaries", () => {
  const template = firstAnnotated(humanIghV, "V");
  const report = preprocessGermlineFasta(`>${template[0]}_S9997\n${template[1].slice(0, 90)}\n`, "V", humanIghV, ["IGH"]);
  assert.equal(report.annotated, 0);
  assert.equal(report.unannotated, 1);
  assert.doesNotMatch(report.fasta, /SWIGMETA=/);
});

test("J homology resolves multiple motif candidates but still verifies the mapped F/W-G anchor", () => {
  const template = firstAnnotated(humanIghJ, "J");
  const anchor = template[2]![1] + 1;
  assert.ok(anchor > 6);
  const changed = `TTCGGT${template[1].slice(6)}`;
  assert.notEqual(changed.slice(anchor, anchor + 6), "TTCGGT");
  const report = preprocessGermlineFasta(`>${template[0]}_S9999\n${changed}\n`, "J", humanIghJ, ["IGH"]);
  assert.equal(report.annotated, 1);
  assert.equal(report.transferred, 1);
  assert.match(report.fasta, /SWIGMETA=[^\n]*,6\n/);
});

test("germline preprocessing rejects a locus mismatch", () => {
  assert.throws(
    () => preprocessGermlineFasta(">TRBV1*01\nACGTACGTACGT\n", "V", humanIghV, ["IGH"]),
    /outside the selected IGH search space/,
  );
});

test("the KI catalog exposes human IGH/TCR and macaque IGH collections by matching taxon and locus", () => {
  assert.deepEqual(collectionsFor("Homo sapiens", "IGH").map((entry) => entry.id), ["kiarva-human-igh"]);
  for (const locus of ["TRA", "TRB", "TRD", "TRG"] as const) {
    const collections = collectionsFor("Homo sapiens", locus);
    assert.equal(collections.length, 1);
    assert.equal(collections[0].locus, locus);
    assert.ok(collections[0].segments.V);
    assert.ok(collections[0].segments.J);
  }
  assert.equal(collectionsFor("Macaca mulatta_AG07107", "IGH")[0]?.id, "kimdb-macaca_mulatta");
  assert.equal(collectionsFor("Macaca fascicularis", "IGH")[0]?.id, "kimdb-macaca_fascicularis");
  assert.deepEqual(collectionsFor("Mus musculus_C57BL/6", "IGH"), []);
});

test("database choices prefer the compatible KI collection and retain IMGT as fallback", () => {
  const cat = databaseOptionsFor("Felis catus_Abyssinian", "202632-7");
  assert.deepEqual(cat.map((option) => option.id), [DEFAULT_DATABASE_ID]);
  assert.match(cat[0].label, /^IMGT\/GENE-DB 202632-7/);
  assert.deepEqual(databaseOptionsFor("Homo sapiens", "202632-7", "BCR").map((option) => option.id), ["kiarva-human-igh", DEFAULT_DATABASE_ID, "ki-human-tcr"]);
  assert.deepEqual(databaseOptionsFor("Homo sapiens", "202632-7", "TCR").map((option) => option.id), ["ki-human-tcr", DEFAULT_DATABASE_ID, "kiarva-human-igh"]);
  assert.deepEqual(databaseOptionsFor("Macaca mulatta_AG07107", "202632-7", "BCR").map((option) => option.id), ["kimdb-macaca_mulatta", DEFAULT_DATABASE_ID]);
  assert.equal(preferredDatabaseIdFor("Homo sapiens","BCR"),"kiarva-human-igh");
  assert.equal(preferredDatabaseIdFor("Homo sapiens","TCR"),"ki-human-tcr");
  assert.equal(preferredDatabaseIdFor("Macaca fascicularis_Cynomolgus","IGH"),"kimdb-macaca_fascicularis");
  assert.equal(preferredDatabaseIdFor("Felis catus_Abyssinian","BCR"),DEFAULT_DATABASE_ID);
});

test("published databases compose by receptor, locus, and segment", () => {
  const kiarva = REFERENCE_DATABASES.find((database) => database.id === "kiarva-human-igh");
  const tcr = REFERENCE_DATABASES.find((database) => database.id === "ki-human-tcr");
  assert.ok(kiarva);
  assert.ok(tcr);
  assert.deepEqual(collectionsForDatabase(kiarva, "BCR").map((collection) => collection.locus), ["IGH"]);
  assert.deepEqual(collectionsForDatabase(tcr, "TCR").map((collection) => collection.locus), ["TRA", "TRB", "TRD", "TRG"]);
  assert.deepEqual(collectionsForDatabase(tcr, "TRB").map((collection) => collection.locus), ["TRB"]);
  assert.deepEqual(databasesForCell("Homo sapiens", "IGH", "V").map((database) => database.id), ["kiarva-human-igh"]);
  assert.deepEqual(databasesForCell("Homo sapiens", "IGK", "V"), []);
  assert.deepEqual(databasesForCell("Homo sapiens", "TRB", "D").map((database) => database.id), ["ki-human-tcr"]);
});

test("a locus-specific replacement retains IMGT records from the other combined-chain loci", () => {
  const customIghV = ">IGHV_CUSTOM*01\nACGTACGTACGT\n";
  const composed = composeReferenceOverrides(
    ["IGH", "IGK", "IGL"],
    { [referenceCellKey("IGH", "V")]: { text: customIghV } },
    (locus, segment) => (human.loci[locus]?.[segment] ?? []).map((allele: MetadataAllele) => `>${allele[0]}\n${allele[1]}\n`).join(""),
  );
  assert.match(composed.V ?? "", /^>IGHV_CUSTOM\*01/m);
  assert.match(composed.V ?? "", new RegExp(`^>${human.loci.IGK.V[0][0].replace(/[.*+?^${}()|[\]\\]/g, "\\$&")}(?:\\s|$)`, "m"));
  assert.match(composed.V ?? "", new RegExp(`^>${human.loci.IGL.V[0][0].replace(/[.*+?^${}()|[\]\\]/g, "\\$&")}(?:\\s|$)`, "m"));
  assert.equal(composed.D, undefined);
  assert.equal(segmentAppliesToLocus("IGK", "D"), false);
  assert.equal(segmentAppliesToLocus("IGH", "D"), true);
});

test("per-cell allele exclusions alter the initial composed assignment FASTA only in that locus and segment", () => {
  const excludedIghV = human.loci.IGH.V[0][0];
  const retainedIghV = human.loci.IGH.V[1][0];
  const representativeIgkV = human.loci.IGK.V[0][0];
  const representativeIghJ = human.loci.IGH.J[0][0];
  const composed = composeReferenceOverrides(
    ["IGH", "IGK", "IGL"],
    {},
    (locus, segment) => (human.loci[locus]?.[segment] ?? []).map((allele: MetadataAllele) => `>${allele[0]}\n${allele[1]}\n`).join(""),
    { [referenceCellKey("IGH", "V")]: [excludedIghV] },
  );
  assert.ok(composed.V);
  assert.doesNotMatch(composed.V!, new RegExp(`^>${excludedIghV.replace(/[.*+?^${}()|[\]\\]/g, "\\$&")}(?:\\s|$)`, "m"));
  assert.match(composed.V!, new RegExp(`^>${retainedIghV.replace(/[.*+?^${}()|[\]\\]/g, "\\$&")}(?:\\s|$)`, "m"));
  assert.match(composed.V!, new RegExp(`^>${representativeIgkV.replace(/[.*+?^${}()|[\]\\]/g, "\\$&")}(?:\\s|$)`, "m"));
  assert.equal(composed.J, undefined);
  const baselineVCount = human.loci.IGH.V.length + human.loci.IGK.V.length + human.loci.IGL.V.length;
  assert.equal((composed.V!.match(/^>/gm) ?? []).length, baselineVCount - 1);
  const baselineJ = (human.loci.IGH.J ?? []).map((allele: MetadataAllele) => `>${allele[0]}\n${allele[1]}\n`).join("");
  assert.match(baselineJ, new RegExp(`^>${representativeIghJ.replace(/[.*+?^${}()|[\]\\]/g, "\\$&")}(?:\\s|$)`, "m"));
});

test("bundled KIMDB files are intact FASTA exports", () => {
  for (const species of ["Macaca_mulatta", "Macaca_fascicularis"]) {
    for (const segment of ["V", "D", "J"]) {
      const text = fs.readFileSync(new URL(`../public/references/kimdb-1.1/${species}/IGH/${segment}.fasta`, import.meta.url), "utf8");
      assert.ok(text.startsWith(">"));
      assert.ok((text.match(/^>/gm) ?? []).length >= (segment === "V" ? 600 : 14));
    }
  }
});

test("the built IMGT pack includes assembled heavy-chain and TCR constant references", () => {
  const rhesus = pack.species.find((entry: { name: string }) => entry.name === "Macaca mulatta_AG07107");
  const cynomolgus = pack.species.find((entry: { name: string }) => entry.name === "Macaca fascicularis");
  assert.ok(human.loci.IGH.C.length > 50);
  for (const prefix of ["IGHM", "IGHD", "IGHG", "IGHA", "IGHE"]) {
    assert.ok(human.loci.IGH.C.some(([name]: MetadataAllele) => name.startsWith(prefix)), prefix);
  }
  assert.ok(human.loci.TRA.C.length > 0);
  assert.ok(human.loci.TRB.C.length > 0);
  assert.ok(rhesus.loci.IGH.C.length > 0);
  assert.ok(cynomolgus.loci.IGH.C.length > 0);
  assert.ok((rhesus.loci.IGH.C as MetadataAllele[]).every((allele) => allele[1].length >= 300));
});

test("every bundled KIMDB macaque V and J record receives validated in-browser annotation", { timeout: 30_000 }, () => {
  for (const [speciesName, slug] of [["Macaca mulatta_AG07107", "Macaca_mulatta"], ["Macaca fascicularis", "Macaca_fascicularis"]]) {
    for (const segment of ["V", "J"] as const) {
      let text = fs.readFileSync(new URL(`../public/references/kimdb-1.1/${slug}/IGH/${segment}.fasta`, import.meta.url), "utf8");
      let report;
      for (const templates of imgtTemplateTiers(speciesName, segment)) {
        report = preprocessGermlineFasta(report?.fasta ?? text, segment, templates, ["IGH"]);
        if (report.annotated === report.count) break;
      }
      assert.ok(report);
      assert.equal(report.annotated, report.count, `${speciesName} ${segment}: ${report.warnings.join("; ")}`);
      assert.equal((report.fasta.match(/SWIGMETA=/g) ?? []).length, report.count);
    }
  }
});
