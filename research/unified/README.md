# Joint-model stress tests

From the source root, with a Node version supporting TypeScript stripping:

```
node --experimental-strip-types research/unified/stress.ts 0 201 2 mixed 3000 null.json
node --experimental-strip-types research/unified/stress.ts .01 202 1 mixed 3000 rare.json
node --experimental-strip-types research/unified/stress.ts .01 203 2 memory 3000 memory.json
node --experimental-strip-types research/unified/stress.ts .05 204 2 mixed 3000 strong.json
```

Arguments: true allele frequency, deterministic seed, number of changed sites (1 or 2), repertoire regime, target-gene lineages, output JSON. The generator is deliberately not the fitted HS5F model. It includes uniform background mutation, exceptional hotspots, adjacent coordinated changes and sequencing error. See the 0.38.9 benchmark for failures as well as successes. These runs do not include assignment/denoising/lineage construction; they test discovery on independent simulated lineages with known truth.

The implemented likelihood lives in `src/unified-germline.ts` and `src/shm-model/unified-{mixture,kernel}.ts`. Mechanistic and UI invariants are tested in `tests/unified-germline.test.ts`. The other personalized-germline method is retained unchanged.
