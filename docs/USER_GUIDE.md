# SCN1A2A-Atlas — user guide

SCN1A2A-Atlas is a manually curated collection of functional-genomics datasets addressing **SCN1A** and **SCN2A**, the paralogous genes encoding the voltage-gated sodium channels NaV1.1 and NaV1.2. Every dataset was selected through a systematic review and reprocessed from raw reads through a single quantification pipeline, so that differences between studies reflect biology and design rather than the processing convention of each laboratory. The allele mechanism of every record was curated from the full text of the source publication, so that studies can be grouped by what was actually perturbed rather than by how the accession was labelled.

The application requires no programming and no prior download of data. Gene identifiers are resolved through one-to-one orthologue groups, so a single query returns human and mouse records together.

**Application:** <https://iza-mamede.shinyapps.io/SCN1A2A-Atlas/>
**Scope:** 18 datasets (10 SCN1A, 8 SCN2A; 8 human, 10 mouse) · version [N] · updated [date]

---

## Contents

- [Starting questions and corresponding modules](#starting-questions-and-corresponding-modules)
- [Recommended workflow](#recommended-workflow)
- [Glossary](#glossary)
- [Module descriptions](#module-descriptions)
- [Worked examples](#worked-examples)
- [Interpretation and limitations](#interpretation-and-limitations)
- [Troubleshooting](#troubleshooting)
- [Methods summary and citation](#methods-summary-and-citation)

---

## Starting questions and corresponding modules

The three most frequent queries and the modules in which they are answered.

| Question | Module |
|---|---|
| Change in gene-level expression between perturbation and control | Gene DE |
| Identification of the affected isoform | DET |
| Biological processes altered in more than one study | Pathways |

## Recommended workflow

1. **Define the context.** Select organism, channel and gene in the sidebar. Orthologue mapping returns human and mouse records for the same query.
2. **Read the main panel.** Each module carries a note on the interpretation of the plot displayed.
3. **Verify the source.** Consult the Datasets module for the sample size, the curated allele mechanism and the tissue of each study. The allele mechanism determines what a measurement at the perturbed locus can report, and no result should be interpreted without it.

---

## Glossary

| Term | Definition |
|---|---|
| Allele mechanism | The curated description of what was actually perturbed in a study: gene trap, engineered protein-null, patient nonsense or frameshift, knock-in missense, poison-exon knock-in, regulatory-interval deletion, or conditional knockout. Curated from the full text of each source publication rather than from accession metadata. |
| Contrast | The comparison tested in a study: perturbation versus the matched control arm designated by the original authors. Contrasts are defined within a study and are never constructed across studies. |
| Control / WT | The unperturbed arm of a study, used as the reference for its contrast. |
| DEG | Differentially expressed gene. The summarised expression of the gene differs between the two arms of a contrast. |
| DET | Differentially expressed transcript. A specific isoform differs between arms, which may occur while the gene-level total remains unchanged. |
| het / hom / KO | Perturbation of one allele (heterozygous), of both alleles (homozygous), or abolition of gene function (knockout). |
| Isoform (transcript) | One of several RNA species produced from the same gene by alternative exon usage. Identifiers such as `SCN2A-204` denote individual transcripts, which may differ in function. |
| log2FC | Log2 of the ratio between perturbation and control. Zero indicates no change; +1 a doubling; −1 a halving. Negative values indicate lower abundance in the perturbed arm. |
| NaV1.1 / NaV1.2 | The voltage-gated sodium channels encoded by SCN1A and SCN2A. NaV1.1 predominates in fast-spiking parvalbumin-expressing GABAergic interneurons; NaV1.2 in excitatory glutamatergic neurons. |
| Orthologue mapping | One-to-one correspondence between human and mouse genes, applied so that a single query returns records from both species together. Mapping status is recorded in the data dictionary. |
| padj / qvalue | The p-value adjusted for multiple testing. Lower values indicate greater statistical confidence; the conventional threshold is 0.05. |
| Poison exon | An exon whose inclusion introduces a premature termination codon, directing the transcript to nonsense-mediated decay rather than translation. A regulatory mechanism of established relevance in SCN1A and SCN2A. |
| Haploinsufficiency | A single functional copy of the gene is insufficient for normal function, producing a disease phenotype. |
| Bulk RNA-seq | Sequencing that measures the average expression across all cells in a sample. |
| Single-cell / single-nucleus | Sequencing that measures expression in individual cells or nuclei, allowing changes to be attributed to specific cell types. |
| iPSC neuron | A neuron differentiated in vitro from an induced pluripotent stem cell line, either patient-derived or engineered. |
| t-SNE / UMAP | Dimensionality-reduction methods that project high-dimensional data into two dimensions. Proximity indicates similarity of profile; axis values have no biological interpretation. |
| Pathway | A set of genes associated with a common biological function or process. |

---

## Module descriptions

### Datasets — curated master table

Filters the collection by gene, organism, assay, mechanism and processing status. Each row is a study. The `open` link resolves to the source publication and the `SRA/GEO` link to the deposited raw data. The gene-model map places each allele on the exon structure of the canonical transcripts, with gene-wide perturbations drawn as a track beneath the point variants.

> **Interpretation.** The mechanism column is the entry point to the collection. Studies described in the literature as loss of function differ in whether a transcript from the perturbed allele survives to be counted, in which cells the perturbation acts, and in whether transcript abundance is expected to change at all.

> **Caveat.** Several records carry a curated mechanism and no differential expression results, because the deposited material did not permit reprocessing from raw reads. The processing status column identifies these records.

### Gene DE — gene-level differential expression

Returns the gene-level result for the selected gene across all contrasts at once. Each bar corresponds to one contrast, defined as perturbation versus the matched control arm designated by the original authors.

> **Interpretation.** Bars to the left indicate lower abundance in the perturbed arm; bars to the right, higher abundance. Opaque bars meet the significance threshold; translucent bars do not.

> **Caveat.** The expected result at the perturbed locus depends on the allele class. A knock-in missense allele substitutes a single residue and is not expected to alter abundance, so a measured effect of zero is the correct result rather than a failed experiment. A knockout-first gene trap leaves residual message, so the magnitude at the perturbed locus understates the functional loss. A heterozygous deletion assayed in whole brain is diluted by every cell type that does not express the gene.

### DET — transcript-level differential expression

Resolves changes that gene-level analysis does not detect, since the total abundance of a gene may remain constant while the proportion between its isoforms shifts. Testing was performed with swish, which propagates inferential replicate uncertainty and is therefore appropriate where two transcripts of the same gene compete for the same reads. Transcripts were grouped with Terminus before testing, so that transcripts the expectation-maximisation step cannot separate are collapsed into a single testable unit rather than reported as independent features.

> **Interpretation.** In the transcript-level embedding, examine whether samples sharing an allele mechanism group together. Grouping indicates an isoform profile characteristic of that mechanism; dispersion indicates no defined signature.

### t-SNE — sample-level structure

Projects all samples in two dimensions and allows the colouring to be switched between gene, organism, dataset and mutation, in order to establish which variable accounts for the observed structure.

> **Interpretation.** Begin by colouring the embedding by dataset. This module is descriptive: samples were quantified through a common pipeline, but no cross-study batch correction was applied, and contrasts are never constructed across studies. Structure that follows the dataset annotation reflects study of origin rather than biology.

> **Caveat.** t-SNE does not preserve global distances. Proximity between points indicates similarity, but the distance between well-separated groups is not interpretable.

### Single-cell — cell-type resolution

Displays the processed droplet single-cell and single-nucleus embeddings, quantified with alevin, annotated to cell types and harmonised to a common schema with a standardised per-cell channel-expression column. Colouring by cell type displays the annotation; colouring by channel expression displays SCN1A or SCN2A abundance.

> **Interpretation.** Consult the per-type cell counts. Cell types represented by few cells do not support conclusions, and the compartment in which a channel is expressed constrains what a whole-tissue measurement can show.

### Studies — per-study view

Presents the complete result of a single study rather than the genes of interest alone, combining gene-level and transcript-level volcano plots with enrichment summaries for the selected contrast.

> **Interpretation.** In the volcano plot, the horizontal axis represents the magnitude of change and the vertical axis the statistical confidence. Highlighted points meet both thresholds.

### Pathways — processes recurring across studies

Identifies pathways altered in more than one study. The minimum-datasets control sets the number of studies in which a pathway must appear before it is displayed.

> **Caveat.** Enrichment modules are presented as hypothesis-generating summaries rather than as evidence of direct regulation. Recurrence across studies does not imply that the change occurred in the same direction; the direction of effect should be verified per study in the Studies module.

---

## Worked examples

**Gene-level response to SCN1A haploinsufficiency in mouse models.** Open the Gene DE module, select organism *Mus musculus* and channel SCN1A, and set the padj threshold to 0.05. Count the opaque bars with negative values. Return to the Datasets module to identify the allele mechanism of each study, since the expected magnitude at the perturbed locus differs between a gene trap, a heterozygous deletion and a knock-in missense allele.

**Isoform-level changes at poison exons.** Open the DET module, select the gene of interest and examine the bar plot. Identify the transcripts with the largest change and verify, in the Datasets module, whether the study design addressed that mechanism.

**Cross-species comparison for a candidate gene.** Query the gene once in the Gene DE module. Orthologue mapping returns human and mouse records together, allowing the response to SCN1A haploinsufficiency in mouse and to SCN2A loss of function in human to be read side by side.

**Biological processes shared between studies.** Open the Pathways module and set the minimum number of datasets to 3. Record the pathways and the corresponding studies, then verify the direction of effect for each in the Studies module.

---

## Interpretation and limitations

### Expected concordance between studies

Liao and colleagues perturbed SCN2A by three independent strategies in the same system, with the same differentiation protocol and in the same laboratory, and reported expression profiles correlating at 0.202 between the two CRISPR edits and at 0.343 and 0.335 between CRISPR interference and each edit, with the number of differentially expressed genes ranging from 223 to 4,587. Two of those three perturbations correspond to datasets in this collection. Concordance between studies that differ additionally in species, tissue, developmental stage and allele class should therefore be interpreted against 0.2 to 0.35 rather than against zero, and low convergence between such studies is the expected outcome.

### Scope of the collection

- **No gain-of-function perturbation is represented.** The records model loss of function, haploinsufficiency, a hypomorphic allele and regulatory changes. Studies addressing the gain-versus-loss axis have been published, but none identified deposited raw sequencing data that could be reprocessed.
- **Coverage is stated as of the search date.** Candidate datasets were identified between 2 and 6 May 2025. The master table is maintained as a living record, re-queried with the same string at each release, with the query date stamped in the data dictionary. Exhaustiveness is not claimed.
- **Author-processed outputs were not aggregated.** Records deposited only as count matrices are curated for mechanism but carry no differential expression results.
- **The resource does not establish causality** and is not intended for diagnostic or therapeutic decisions.
- **The resource is not a meta-analysis.** Presenting contrasts alongside one another is not equivalent to combining them statistically.

---

## Troubleshooting

| Observation | Probable cause and action |
|---|---|
| Slow initial load | The server instance was idle. Wait approximately 20 seconds and reload the page. |
| A dataset appears in Datasets but returns nothing in Gene DE or DET | Several records carry a curated allele mechanism without differential expression results, because the deposited material did not permit reprocessing from raw reads. This is expected behaviour, not an error. |
| Empty plot | The padj threshold is too restrictive, or the gene was not quantified in the selected dataset. Raise the threshold or change the organism. |
| Gene not found in the search field | The gene was not detected in that dataset. Not every study quantifies every gene. |
| Plot truncated at the margin | Insufficient window width. Collapse the sidebar or use a larger display. |
| Overlapping labels | Too many categories displayed simultaneously. Apply additional filters. |

---

## Methods summary and citation

The harmonised objects consumed by the application are produced from pipeline outputs by a single reproducible build script, so that no value displayed in the interface is entered by hand.

- **Quantification.** Bulk RNA sequencing reads quantified with Salmon v1.10.9 in selective-alignment mode with Gibbs posterior sampling enabled (100 replicates), against GENCODE v49 (human) and GENCODE vM39 (mouse).
- **Import and testing.** Transcript-level estimates imported with tximport and tximeta for provenance-aware import, summarised to gene level, and tested with DESeq2. Transcript-level differential expression computed with swish (fishpond), with transcripts grouped by Terminus before testing.
- **Single-cell data.** Droplet single-cell and single-nucleus data quantified with alevin, annotated to cell types and harmonised to a common schema.
- **Contrasts.** Defined per study as perturbation versus matched control, using the control arm designated by the original authors, and never constructed across studies. No cross-study batch correction was applied.
- **Implementation.** R with Shiny; tables rendered with DT and figures with ggplot2.

**Citation:** [full reference]
**Code and data dictionary:** <https://github.com/iza-mcac/SCN1A2A-Atlas>
**Reporting errors:** [issue tracker address]

> Bracketed fields are to be completed before release.
