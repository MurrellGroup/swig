# Integration and editorial notes

## Contents and entry point

`index.md` is the reader-facing entry point. The numbered chapters form a continuous introduction; `13-analysis-guide.md` provides navigation by analysis. `glossary.md` and `references.md` are supporting pages. This file contains integration notes and need not appear in the reader-facing navigation.

Keep the Markdown files together and retain the `assets/` subdirectory. Links between pages and to figures are relative. External concept links point to Wikipedia; supporting literature and implementation references are collected in `references.md`.

The pages use ordinary Markdown headings, links, images, fenced text blocks, and pipe tables. A renderer with table support is needed. They use no MDX, JavaScript, theme-specific components, footnote extension, or mathematical typesetting dependency. Internal section links use conventional GitHub-style heading fragments; a site with different heading-ID rules should adapt those links or configure matching IDs.

## Repository integration

The collection lives in `public/immunology/`. Vite serves that directory in development and copies it into the production site. The browser reader is `index.html`; its chapter URLs use `?page=02-vdj-recombination`, with ordinary fragments for sections. All resource URLs are relative, so the reader works at both `/immunology/` and a project path such as `/swig/immunology/`.

The application links to the guide and analysis index from `src/main.tsx`, using Vite's `BASE_URL`. Both links open a separate tab so the current analysis remains open. The repository README also links directly to the Markdown entry point.

Edit the Markdown and SVGs in place. The reader fetches the selected chapter from this directory; there is no second generated copy of its text. `guide-pages.js` defines the chapter order and permitted filenames. Add new reader-facing pages there and update incoming links.

`guide-markdown.js` handles the syntax used by this collection: headings, paragraphs, inline links and images, emphasis, code spans and fences, simple lists, horizontal rules, and pipe tables. It escapes raw HTML and restricts image paths to the local SVG assets. Extend its tests before introducing additional Markdown syntax. The reader requires JavaScript; its source link and GitHub fallback provide access to the Markdown separately.

Run `npm run test:docs` after editing. It checks every page's headings, internal links, section anchors, SVG accessibility, navigation coverage, base-path handling, and renderer fixtures. The command also runs as part of `npm test`. No additional runtime packages, external scripts, analytics, or font downloads are used by the reader.

## Figures

The four SVGs are original, self-contained vector diagrams, with editable text, `viewBox`, a white background, and embedded accessibility titles and descriptions. The Markdown image tags include independent alt text and captions. No font files, external images, scripts, or network resources are required. Text uses a system sans-serif fallback stack, with monospace for nucleotide sequences.

For responsive display, use the site's usual equivalent of `max-width: 100%; height: auto`. The diagrams contain detailed labels and should be openable at full size on small screens. A renderer that sanitizes inline SVG can still display these files through ordinary Markdown image links, subject to the site's asset policy.

All drawn sequences, tree branches, and numerical examples are synthetic. The figure files are:

- `assets/01-vdj-recombination.svg`
- `assets/02-b-cell-response.svg`
- `assets/03-lineage-and-shm.svg`
- `assets/04-uca-inference.svg`

## Scope of the software review

The guide was checked against the public Swig repository and method documentation on **12 September 2026**. The inspected `package.json` and methods index identified **version 0.38.10**. This was a documentation review, rather than an end-to-end test of the deployed application. The method links in the bibliography follow `main`; they may change after this review.

The analysis guide includes the separate repertoire allele-pooling, missing-V, personalized expressed V-set, and joint inherited/SHM analyses. It also distinguishes the standard germline-rooted display and parsimony mapping from phylogenetic UCA inference. Future edits should preserve those distinctions as the implementation changes.

Only a small number of defaults are quoted. The 85% CDR3 identity threshold in Chapter 9 and the three-observation minimum in Chapter 11 describe the inspected implementation, rather than universal biological recommendations. Readers should consult the active method specification and settings for their installed version.

## Biological scope

The foundational explanation mainly concerns human and mouse immunology. The guide calls out important species and locus dependencies, particularly in recombination priors, SHM models, and TCR interpretation. The vaccine chapter uses dated research examples and makes no claim to survey the latest clinical-trial status.

Wikipedia links serve as conceptual navigation. The bibliography supplies foundational teaching sources, primary studies, standards, and implementation specifications. Research precedents cited alongside Swig are identified as precedents, rather than interchangeable implementations.

## Updating the collection

After changing a filename or heading, update incoming relative links. After changing a scientific example, check both its prose and figure. For the lineage figure, verify every edge against the displayed sequences and preserve the distinction between the UCA and the sampled MRCA. For the probability figure, preserve the dependencies between codon positions.

When software behaviour changes, update the relevant chapter, analysis-guide row, source reference, and review version together. New analyses should receive an explanation of the biological quantity being estimated and the observations that can support it.
