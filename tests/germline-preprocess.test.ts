import assert from "node:assert/strict";
import fs from "node:fs";
import test from "node:test";
import zlib from "node:zlib";

import { preprocessGermlineFasta, type MetadataAllele } from "../src/germline-preprocess.ts";
import { collectionsFor, databaseOptionsFor, DEFAULT_DATABASE_ID } from "../src/reference-catalog.ts";

const pack = JSON.parse(zlib.gunzipSync(fs.readFileSync(new URL("../public/references/imgt-202632-7.json.gz", import.meta.url))).toString());
const human = pack.species.find((entry: { name: string }) => entry.name === "Homo sapiens");
const humanIghV = human.loci.IGH.V as MetadataAllele[];
const humanIghJ = human.loci.IGH.J as MetadataAllele[];

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

test("database choices always default to IMGT and never expose a mismatched KI collection", () => {
  const cat = databaseOptionsFor("Felis catus_Abyssinian", "IGH", "202632-7");
  assert.deepEqual(cat.map((option) => option.id), [DEFAULT_DATABASE_ID]);
  assert.match(cat[0].label, /^IMGT\/GENE-DB 202632-7/);
  assert.deepEqual(databaseOptionsFor("Homo sapiens", null, "202632-7").map((option) => option.id), [DEFAULT_DATABASE_ID]);
  assert.deepEqual(databaseOptionsFor("Homo sapiens", "IGH", "202632-7").map((option) => option.id), [DEFAULT_DATABASE_ID, "kiarva-human-igh"]);
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
