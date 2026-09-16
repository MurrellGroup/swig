/** Reader order. Filenames are also the allow-list for same-origin Markdown loads. */
export const pages = [
  { slug: 'index', title: 'Guide overview', group: 'Start here' },
  { slug: '01-recognition-and-receptors', title: '1. Recognition and receptors', group: 'Foundations' },
  { slug: '02-vdj-recombination', title: '2. Building a receptor', group: 'Foundations' },
  { slug: '03-b-cell-responses', title: '3. B-cell responses', group: 'Foundations' },
  { slug: '04-t-cell-responses', title: '4. T-cell responses', group: 'Foundations' },
  { slug: '05-samples-sequences-and-counts', title: '5. Samples, sequences, and counts', group: 'Sequence analysis' },
  { slug: '06-vdj-assignment', title: '6. Assigning V, D, J, and C', group: 'Sequence analysis' },
  { slug: '07-errors-and-chimeras', title: '7. Errors, duplicates, and chimeras', group: 'Sequence analysis' },
  { slug: '08-germline-variation', title: '8. Inherited variation', group: 'Sequence analysis' },
  { slug: '09-repertoires-and-clones', title: '9. Repertoires and clonal families', group: 'Sequence analysis' },
  { slug: '10-shm-alignments-and-trees', title: '10. Mutations, alignments, and trees', group: 'Sequence analysis' },
  { slug: '11-unmutated-common-ancestors', title: '11. Unmutated common ancestors', group: 'Sequence analysis' },
  { slug: '12-vaccine-design', title: '12. Antibody histories and vaccine design', group: 'Sequence analysis' },
  { slug: '13-analysis-guide', title: 'Find background by analysis', group: 'Reference' },
  { slug: 'glossary', title: 'Glossary', group: 'Reference' },
  { slug: 'references', title: 'References', group: 'Reference' },
];

export function findPage(slug = 'index') {
  return pages.find((page) => page.slug === slug) || null;
}

export function pageUrl(slug = 'index') {
  return slug === 'index' ? './' : `?page=${encodeURIComponent(slug)}`;
}
