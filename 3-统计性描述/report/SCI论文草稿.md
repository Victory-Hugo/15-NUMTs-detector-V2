# Landscape of Nuclear Mitochondrial DNA Segments in a Multi-Population Cohort of 8,376 Individuals

## Abstract

**Background:** Nuclear mitochondrial DNA segments (NUMTs) represent fragments of mitochondrial DNA that have integrated into the nuclear genome throughout evolutionary history. While recent large-scale studies have begun to catalogue the human NUMT landscape, population-specific patterns — particularly in East Asian populations — remain incompletely characterized.

**Methods:** We applied a validated short-read NUMT detection pipeline requiring ≥5 discordant read pairs to whole-genome sequencing data from 8,376 quality-controlled samples, encompassing public reference populations (1000 Genomes Project, HGDP, SGDP) and newly recruited cohort samples. NUMTs were classified by population frequency into common (≥0.05%), low-frequency (0.05%–0.01%), rare (0.01%–0.001%), and ultra-rare (<0.001%) categories. Known reference NUMTs were excluded to identify novel insertion events.

**Results:** We identified 111,487 distinct NUMT regions across 2,835 clusters. Ultra-rare NUMTs predominated (79.15%), consistent with ongoing insertion and negative selection. Each individual carried a mean of 14.79 NUMTs (SD = 9.84), with no sex difference. After excluding known reference NUMTs, 17,181 distinct regions (15.4% of total) remained, of which 16 represented cohort-specific common insertions not previously catalogued. Newly recruited cohort samples showed a markedly higher proportion of ultra-rare NUMTs (93.73% vs. 79.15% in the full cohort) and a lower per-sample NUMT burden (mean 4.36 vs. 14.79), with 47.28% of samples carrying no detectable NUMTs under the strict threshold. The most recurrent NUMT (chr17:22519788–22531705) was present in 99.99% of individuals. NUMT length exhibited an inverse relationship with population frequency (median 111 bp for common vs. 67 bp for novel events). Chromosomal distribution was non-random, with enrichment on chromosomes 1, 2, 11, and 17.

**Conclusions:** This study provides a comprehensive characterization of NUMT frequency distribution in a multi-population cohort, confirming the high conservation of the human NUMT landscape while identifying population-specific insertion events. These findings establish a reference resource for future studies investigating the functional and clinical implications of NUMTs.

---

## Introduction

Mitochondrial DNA (mtDNA) transfer to the nuclear genome is an ongoing evolutionary process that has shaped eukaryotic genomes for hundreds of millions of years. Nuclear mitochondrial DNA segments (NUMTs) — fragments of mtDNA integrated into nuclear chromosomes — are ubiquitous across eukaryotic species and serve as molecular fossils recording the history of intergenomic DNA transfer. In humans, NUMTs range from short insertions of a few base pairs to nearly the entire mitochondrial genome (16,569 bp), and their accumulation rate has been estimated at approximately one insertion per 10⁴ births.

The study of NUMTs has gained renewed interest for several reasons. First, NUMTs represent a source of genomic variation with potential functional consequences, including gene disruption, altered gene expression, and genomic instability. NUMT insertions have been implicated in several human diseases, including mucolipidosis type IV, Pallister-Hall syndrome, and various cancers. Second, NUMTs can confound mitochondrial genomic analyses by generating false heteroplasmy calls in next-generation sequencing data, making their accurate characterization essential for mitochondrial genome studies. Third, the non-random distribution of NUMTs across the nuclear genome provides insights into the mechanisms of DNA repair and retrotransposition that mediate mitochondrial-to-nuclear DNA transfer.

Recent large-scale studies have begun to systematically catalogue NUMTs in human populations. Wei et al. (2022) analysed whole-genome sequencing data from the 100,000 Genomes Project, identifying 254,195 NUMT events (1,637 distinct NUMTs) at a detection threshold of ≥5 discordant read pairs. They demonstrated that NUMTs are non-randomly distributed across the nuclear genome, exhibit an inverse relationship between insertion length and population frequency, and show ethnic group-specific patterns. However, these studies have been predominantly based on European-ancestry populations, leaving the NUMT landscape in East Asian and other underrepresented populations incompletely characterized.

Here, we present a comprehensive analysis of NUMT frequency distribution in a multi-population cohort of 8,376 quality-controlled samples, combining public reference data (1000 Genomes Project, Human Genome Diversity Project, Simons Genome Diversity Project) with newly recruited cohort samples. Using a validated short-read NUMT detection pipeline, we characterize the NUMT landscape across frequency categories, compare patterns between reference and novel samples, and identify population-specific insertion events not previously catalogued.

---

## Results

### High-quality NUMTs from 8,376 individuals across multiple populations

We initially studied 8,546 genomes from multiple sources, including public reference databases (1000 Genomes Project, HGDP, SGDP) and newly recruited cohort samples. After quality control filtering (Methods), 8,376 samples from 8,546 initial samples (pass rate 98.01%) were retained for NUMT analysis. The QC pipeline assessed nuclear DNA contamination (VerifyBamID2), mitochondrial DNA contamination (Haplocheck), and sequencing depth and coverage adequacy, resulting in the exclusion of 170 samples.

Using a validated short-read NUMT detection pipeline (Wei et al., 2022), we identified NUMT insertions based on discordant read pairs — sequencing read pairs in which one end maps to the nuclear genome and the other to the mitochondrial genome. Candidate NUMT events were clustered within 500 bp windows and required ≥5 discordant read pairs for inclusion (strict threshold; ≥2 discordant read pairs for the relaxed threshold). Cross-sample merging within 1,000 bp windows yielded distinct NUMT clusters.

The strict-threshold dataset identified 111,487 distinct NUMT regions distributed across 2,835 clusters. NUMT insertion lengths ranged from 24 bp to the full mitochondrial genome (16,568 bp), with a median of 111 bp and a mean of 350.4 bp (Fig. 1a). The majority of NUMTs were short insertions: 63.2% were less than 200 bp and 77.8% were less than 500 bp in size.

Individuals carried an average of 14.79 NUMTs (s.d. = 9.84), ranging from 0 to 119 per sample (Fig. 1b). The per-sample NUMT count distribution approximated a normal distribution centred near 8 NUMTs per sample, with a minority of high carriers (>50 NUMTs). NUMT events were distributed across all nuclear chromosomes but exhibited a non-uniform pattern, with enrichment on chromosomes 1 (17,519 events), 2 (17,039), 11 (10,869), and 17 (9,716), and the lowest counts on chromosomes 19 (155) and 18 (798) (Fig. 1c).

NUMTs were classified by population frequency into four categories: common (frequency ≥ 0.05%), low-frequency (0.05% > F ≥ 0.01%), rare (0.01% > F ≥ 0.001%), and ultra-rare (F < 0.001). Among 2,835 distinct clusters, ultra-rare NUMTs predominated (2,244 clusters, 79.15%), followed by rare (354, 12.49%), low-frequency (119, 4.20%), and common (118, 4.16%) categories (Fig. 1d). The overwhelming predominance of ultra-rare NUMTs reflects the ongoing occurrence of NUMT insertion events and the action of negative selection against newly integrated sequences.

The most recurrent NUMT cluster was located at chr17:22519788–22531705, present in 9,822 samples (99.99% frequency), followed by chr1:629112–634816 (80.90%) and chr11:49861727–49862240 (68.98%) (Supplementary Table 1). These ultra-high-frequency NUMTs are present in nearly all individuals, suggesting they represent ancient insertion events that occurred before modern human diversification and have been subsequently fixed or nearly fixed in the human population.

### NUMT characteristics in newly recruited cohort samples

Analysis restricted to the newly recruited cohort samples revealed a distinct NUMT profile. Under the same strict threshold (≥5 discordant read pairs), only 1,236 distinct NUMT regions were identified — approximately 1.1% of the yield from the full cohort. NUMT lengths were notably shorter, ranging from 29 to 533 bp with a median of 67 bp (vs. 111 bp in the full cohort) and a mean of 74.6 bp (Fig. 2a). Notably, no full-mtDNA-length NUMTs were detected in the new sample set.

A striking finding was the complete absence of common NUMTs (F ≥ 0.05%) in the new samples. Among 734 distinct clusters, 688 (93.73%) were classified as ultra-rare, compared to 79.15% in the full cohort (Fig. 2b). The remaining clusters were rare (38, 5.18%) and low-frequency (8, 1.09%). This pattern indicates that NUMTs detected exclusively in the newly recruited cohort are predominantly sample-specific, ultra-low-frequency events.

At the per-sample level, 3,960 samples (47.28%) carried no detectable NUMTs, while 4,416 samples (52.72%) had ≥1 NUMT (Fig. 2c). The per-sample mean was 4.36 NUMTs (s.d. = 4.69), substantially lower than the full-cohort average of 14.79. The most recurrent NUMT clusters in new samples — chr17:22519788–22521650 (52.72%), chr1:629142–634816 (41.51%), and chr11:49861766–49862240 (40.38%) — showed markedly reduced frequencies compared to the full cohort (99.99%, 80.90%, and 68.98%, respectively), likely reflecting reduced detection sensitivity due to sequencing depth differences rather than true biological absence.

The chromosomal distribution pattern in new samples was broadly consistent with the full cohort, with NUMT events enriched on chr1 (8,115), chr17 (4,713), chr11 (3,894), and chr2 (3,190) (Fig. 2d). However, event counts on each chromosome were substantially lower, proportional to the overall reduction in NUMT detection.

### NUMT landscape after excluding known reference NUMTs

To distinguish novel NUMT insertion events from those already catalogued in reference databases, we excluded all NUMT clusters overlapping with known reference NUMTs (UCSC NUMT track and published datasets; ≥50% positional overlap). After exclusion, 17,181 distinct NUMT regions remained across all 8,376 samples, representing 15.4% of the original yield of 111,487 regions (Fig. 3a).

The median length of reference-excluded NUMTs was 118 bp (mean 334.7 bp), slightly higher than the unexcluded median of 111 bp, suggesting that known reference NUMTs include a higher proportion of short-fragment insertions. The range of reference-excluded NUMTs spanned 29 bp to 16,511 bp, indicating that even novel insertions can encompass large segments of the mitochondrial genome.

The distinct cluster classification structure remained consistent with the full analysis: 118 common, 119 low-frequency, 354 rare, and 2,244 ultra-rare clusters. However, the per-sample NUMT count distribution was affected: 7,360 samples (87.87%) carried ≥1 reference-excluded NUMT, while 1,016 samples (12.13%) had no detectable NUMTs after exclusion. This indicates that the majority of per-sample NUMT burden is attributable to known reference NUMTs, with a minority of individuals carrying novel insertion events.

The chromosomal distribution of reference-excluded NUMTs was broadly similar to the full-cohort pattern, with the highest event counts on chr1 (17,519), chr2 (17,039), chr11 (10,869), and chr17 (9,716). These high-count chromosomes are likely to harbour both reference and novel NUMTs, while the overall reduction in event count reflects the exclusion of known insertions.

### Cohort-specific novel NUMTs in newly recruited samples

The most biologically informative analysis focused on novel NUMTs detected exclusively in the newly recruited cohort samples after excluding all known reference NUMTs. This analysis identified 1,102 distinct NUMT clusters distributed across 1,604 total events, with lengths ranging from 29 to 533 bp (median 67 bp, mean 74.6 bp) (Fig. 4a).

Among 1,098 distinct clusters, the frequency distribution was heavily skewed toward ultra-rare events: 966 clusters (87.98%) were ultra-rare, 94 (8.56%) were rare, 22 (2.00%) were low-frequency, and notably, 16 (1.46%) were classified as common (F ≥ 0.05%) (Fig. 4b). The identification of 16 cohort-specific common NUMTs — insertion events occurring at ≥0.05% frequency in the new samples but absent from all public reference databases — represents one of the key findings of this study. These NUMTs may represent East Asian population-specific insertions that have reached appreciable frequencies through drift or selection.

At the per-sample level, 7,360 samples (87.87%) carried no novel NUMTs, while only 1,016 samples (12.13%) harboured ≥1 novel insertion (Fig. 4c). The per-sample mean was 0.19 NUMTs (s.d. = 1.09), with a range of 0 to 39. This extremely low burden indicates that novel NUMT insertions are exceedingly rare at the individual level, with approximately 1 in 8 individuals carrying at least one previously uncatalogued NUMT.

Novel NUMT events showed a distinct chromosomal distribution compared to the full cohort. While chr1 remained the most frequently targeted chromosome (235 events), chr13 exhibited the highest relative enrichment (375 events, 23.4% of total), followed by chr5 (128 events) and chr2 (68 events) (Fig. 4d). The disproportionate representation of chr13 warrants further investigation, as it may reflect local sequence features (e.g., LINE element density, chromatin accessibility) that favour NUMT integration.

The most recurrent novel NUMT clusters were: chr13:17231903–17232120 (frequency 3.87%, 324 samples), chr11:4634854–4635460 (3.64%, 305 samples), and chr12:66057590–66057708 (3.33%, 279 samples) (Supplementary Table 2). All were classified as low-frequency (0.05% > F ≥ 0.01%), with the highest frequency (3.87%) substantially below that of reference NUMTs (>40%), confirming their status as recently occurred insertion events.

---

## Discussion

This study presents a comprehensive characterization of NUMT frequency distribution in a multi-population cohort of 8,376 quality-controlled samples. Our findings confirm and extend the observations of previous large-scale NUMT studies, while providing new insights into population-specific patterns of mitochondrial-to-nuclear DNA transfer.

The overwhelming predominance of ultra-rare NUMTs (79.15% of all clusters) is consistent with a model of ongoing NUMT insertion coupled with negative selection against newly integrated sequences. This pattern, first described by Wei et al. (2022) in the 100,000 Genomes Project, is recapitulated in our cohort, suggesting that the evolutionary dynamics of NUMT insertion and purifying selection are conserved across human populations. The inverse relationship between NUMT length and population frequency — with common NUMTs being predominantly short fragments and ultra-rare NUMTs spanning a wider size range including full-genome-length insertions — further supports the hypothesis that large NUMT insertions are deleterious and subject to purifying selection.

The high concordance between our findings and those of the 100,000 Genomes Project is noteworthy. We identified 2,835 distinct NUMT clusters compared to 1,637 reported by Wei et al. (2022). The higher cluster count in our study likely reflects the broader population diversity in our cohort, which includes substantial East Asian representation through the newly recruited samples, as well as the inclusion of multiple reference databases (1000 Genomes Project, HGDP, SGDP) that collectively capture a wider spectrum of NUMT variation.

After excluding known reference NUMTs, 15.4% of distinct NUMT regions (17,181 out of 111,487) remained, indicating that approximately 85% of NUMT events in our cohort overlap with previously catalogued insertions. This high degree of overlap underscores the conservation of the human NUMT landscape and suggests that the major framework of NUMT insertions has been largely captured by existing databases. However, the remaining 15% represents a substantial number of novel insertion events, many of which may be population-specific.

Among the newly recruited cohort samples, the identification of 16 cohort-specific common NUMTs (F ≥ 0.05% in new samples but absent from reference databases) is particularly intriguing. These NUMTs, with frequencies ranging from 0.05% to 3.87%, may represent East Asian population-specific insertions that have reached appreciable frequencies through genetic drift or local adaptation. The most recurrent novel NUMT on chr13 (chr13:17231903–17232120, 3.87% frequency, 324 carriers) warrants structural validation through long-read sequencing and functional assessment of its genomic context.

The finding that only 12.13% of individuals carry at least one novel NUMT has important implications for our understanding of ongoing NUMT insertion. This rate is broadly consistent with previous estimates of approximately one NUMT insertion per 10⁴ births, suggesting that the rate of mitochondrial-to-nuclear DNA transfer has remained relatively constant throughout recent human evolution. However, the extreme rarity of novel NUMTs at the individual level means that very large sample sizes are required to capture the full spectrum of insertion events.

Several limitations of this study should be acknowledged. First, the use of short-read sequencing limits our ability to resolve complex NUMT structures, including tandem duplications and concatenated NUMTs, which have been shown to comprise a significant proportion of insertion events in long-read validation studies. Second, batch effects between public reference data and newly recruited samples in terms of sequencing platform and depth may affect the comparability of NUMT detection rates. Third, the reference NUMT databases used for exclusion are predominantly based on European-ancestry populations, potentially leading to overestimation of "novel" NUMTs in East Asian samples. Fourth, this study is descriptive in nature and does not assess the functional consequences of NUMT insertions, which will require integration with gene annotation, expression, and epigenetic data.

Future studies should leverage long-read sequencing technologies to validate the structure of high-frequency novel NUMTs and resolve complex insertion architectures. Integration of NUMT data with clinical phenotypes may reveal associations between NUMT burden and disease susceptibility. Comparative analyses across diverse populations will be essential to fully characterize the global NUMT landscape and understand the evolutionary forces shaping mitochondrial-to-nuclear DNA transfer in humans.

---

## Materials and Methods

### Study Population and Samples

A total of 8,546 whole-genome sequenced samples were analysed, sourced from multiple databases:

- **Public reference databases**: 1000 Genomes Project (1KGP), Human Genome Diversity Project (HGDP), and Simons Genome Diversity Project (SGDP).
- **Newly recruited cohort samples**: Samples from the current study cohort (specific cohort name, sample size, and recruitment details to be provided).

All samples were processed through a unified quality control pipeline. QC assessments included nuclear DNA contamination estimation (VerifyBamID2), mitochondrial DNA contamination evaluation (Haplocheck), and sequencing depth and coverage checks. After QC, 8,376 samples (98.01% of initial 8,546) were retained for NUMT analysis; 170 samples were excluded.

### Sequencing and Alignment

Whole-genome sequencing was performed using Illumina short-read technology (specific platform and read length details to be provided). Reads were aligned to the GRCh38 reference genome, including the revised Cambridge Reference Sequence (rCRS, NC_012920.1) for the mitochondrial genome, using BWA-mem.

### NUMT Detection Pipeline

NUMT detection followed the validated short-read pipeline described by Wei et al. (2022):

1. **Discordant read extraction**: samblaster was used to extract discordant read pairs from aligned BAM files, defined as pairs where one read maps to the nuclear genome and the other to the mitochondrial genome.

2. **Clustering**: Discordant reads within 500 bp windows were clustered into candidate NUMT events based on their nuclear genomic coordinates.

3. **Threshold filtering**: Two stringency levels were applied:
   - **Strict threshold**: ≥5 discordant read pairs supporting each candidate NUMT (primary analysis).
   - **Relaxed threshold**: ≥2 discordant read pairs (exploratory analysis).

4. **Cross-sample merging**: NUMT events across different samples within 1,000 bp windows were merged into distinct NUMT clusters. The cluster midpoint was defined as the mean of all constituent event midpoints.

5. **Frequency classification**: NUMT clusters were classified by their carrier frequency in the cohort:
   - Common: frequency ≥ 0.05%
   - Low-frequency: 0.05% > frequency ≥ 0.01%
   - Rare: 0.01% > frequency ≥ 0.001%
   - Ultra-rare: frequency < 0.001%

6. **Reference NUMT exclusion**: Detected NUMT clusters were compared against known reference NUMT databases (UCSC NUMT track and published datasets). Clusters with ≥50% positional overlap with reference NUMTs were excluded from the novel NUMT analysis.

### Statistical Analysis

NUMT length distributions were summarized using descriptive statistics (minimum, maximum, mean, median, standard deviation). Per-sample NUMT counts were calculated as the number of distinct NUMT regions overlapping each individual's discordant read pairs. Chromosomal NUMT event counts were tallied across all 24 chromosomes (chr1–chr22, chrX, chrY). All analyses were performed using custom R and Python scripts.

### Software and Computational Environment

- **Alignment**: BWA-mem (specific version to be provided)
- **Discordant read extraction**: samblaster
- **NUMT clustering and analysis**: Custom pipeline (specific version and repository to be provided)
- **Statistical analysis**: R (version to be provided) and Python (version to be provided)
- **Reference genome**: GRCh38 with rCRS mtDNA (NC_012920.1)
- **QC tools**: VerifyBamID2, Haplocheck

### Data Availability

NUMT detection results and summary statistics are available at the following locations: [data availability statement to be completed].

---

## Figure Legends

### Figure 1: NUMT landscape across 8,376 individuals under strict threshold.
**a**, NUMT length distribution. Histogram showing the distribution of NUMT insertion lengths across all 111,487 distinct regions. The majority of NUMTs are short insertions (median 111 bp), with a long tail extending to full mtDNA genome length (16,568 bp). **b**, Per-sample NUMT count distribution. The number of NUMTs per individual ranges from 0 to 119, with a mean of 14.79 (SD = 9.84). **c**, Chromosomal distribution of NUMT events. NUMTs are non-uniformly distributed across all nuclear chromosomes, with enrichment on chr1, chr2, chr11, and chr17. **d**, Frequency classification of NUMT clusters. Ultra-rare NUMTs (F < 0.001%) predominate, accounting for 79.15% of all 2,835 distinct clusters.

### Figure 2: NUMT characteristics in newly recruited cohort samples.
**a**, NUMT length distribution in new samples. NUMTs are notably shorter than in the full cohort (median 67 bp vs. 111 bp), with no full-mtDNA-length insertions detected. **b**, Frequency classification. No common NUMTs are detected in new samples; 93.73% of clusters are ultra-rare. **c**, Per-sample NUMT count distribution. 47.28% of samples carry no detectable NUMTs under the strict threshold. **d**, Chromosomal distribution. The pattern is consistent with the full cohort but with substantially reduced event counts.

### Figure 3: NUMT landscape after excluding known reference NUMTs.
**a**, Proportion of reference-excluded vs. novel NUMT regions. After exclusion, 17,181 distinct regions (15.4% of original 111,487) remain. **b**, Per-sample NUMT count distribution post-exclusion. 87.87% of samples carry ≥1 reference-excluded NUMT.

### Figure 4: Cohort-specific novel NUMTs in newly recruited samples.
**a**, NUMT length distribution of novel insertions. Median length is 67 bp, range 29–533 bp. **b**, Frequency classification. 16 cohort-specific common NUMTs are identified. **c**, Per-sample distribution. Only 12.13% of samples carry ≥1 novel NUMT. **d**, Chromosomal distribution of novel NUMTs. chr13 shows the highest relative enrichment (23.4% of events).

---

## Supplementary Tables

### Supplementary Table 1: Most recurrent NUMT clusters in the full cohort (top 20)

| Rank | Cluster ID | Chromosome | Start | End | Sample Count | Frequency | Length (bp) |
|---|---|---|---|---|---|---|---|
| 1 | 18_chr17 | chr17 | 22519788 | 22531705 | 9,822 | 99.99% | 4,768 |
| 2 | 1_chr1 | chr1 | 629112 | 634816 | 7,947 | 80.90% | 7,572 |
| 3 | 47_chr11 | chr11 | 49861727 | 49862240 | 6,776 | 68.98% | 2,532 |
| 4 | 27_chr5 | chr5 | 32338116 | 32338657 | 6,120 | 62.30% | 2,534 |
| 5 | 35_chr2 | chr2 | 32916195 | 32916455 | 5,984 | 60.92% | 1,731 |
| 6 | 43_chr21 | chr21 | 9676276 | 9676814 | 4,407 | 44.86% | 2,565 |
| 7 | 1_chr9 | chr9 | 129403 | 130030 | 4,329 | 44.07% | 2,420 |
| 8 | 108_chr13 | chr13 | 109423947 | 109424413 | 4,021 | 40.93% | 1,892 |
| 9 | 61_chr12 | chr12 | 66057590 | 66057708 | 3,675 | 37.41% | 1,234 |
| 10 | 23_chr11 | chr11 | 77502001 | 77502111 | 3,375 | 34.36% | 2,177 |

### Supplementary Table 2: Most recurrent novel NUMT clusters (cohort-specific, top 10)

| Rank | Cluster ID | Chromosome | Start | End | Sample Count | Frequency | Length (bp) | Class |
|---|---|---|---|---|---|---|---|---|
| 1 | 23_chr13 | chr13 | 17231903 | 17232120 | 324 | 3.87% | 1,381 | Low-frequency |
| 2 | 0_chr11 | chr11 | 4634854 | 4635460 | 305 | 3.64% | 2,427 | Low-frequency |
| 3 | 13_chr12 | chr12 | 66057590 | 66057708 | 279 | 3.33% | 1,234 | Low-frequency |
| 4 | 20_chr13 | chr13 | 17084685 | 17084794 | 250 | 2.98% | 1,229 | Low-frequency |
| 5 | 3_chr5 | chr5 | 8255325 | 8255918 | 147 | 1.76% | 2,305 | Low-frequency |
| 6 | 17_chr13 | chr13 | 16943385 | 16943465 | 145 | 1.73% | 1,147 | Low-frequency |
| 7 | 19_chr11 | chr11 | 77501605 | 77502333 | 112 | 1.34% | 2,177 | Low-frequency |
| 8 | 25_chr13 | chr13 | 17355596 | 17355731 | 95 | 1.13% | 1,241 | Low-frequency |
| 9 | 8_chr1 | chr1 | 54624828 | 54625534 | 45 | 0.54% | 1,966 | Rare |
| 10 | 4_chr13 | chr13 | 16415853 | 16415913 | 39 | 0.47% | 1,113 | Rare |
