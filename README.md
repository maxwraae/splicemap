# splicemap

Annotate splicing regulatory elements on genomic DNA sequences.

## What It Does

Given a genomic sequence (GenBank format) and an mRNA transcript accession, splicemap discovers exon boundaries by aligning the transcript against the genomic DNA, annotates the resulting introns, and maps all splicing regulatory elements: 5' and 3' splice sites, branch points, polypyrimidine tracts, exonic splicing enhancers (ESE), and exonic splicing silencers (ESS). Output is written back to the GenBank file — viewable in SnapGene or UGENE — with a companion markdown report.

## Quick Start

```bash
git clone https://github.com/maxwraae/splicemap.git
cd splicemap
pip install -r requirements.txt
python splicemap.py splicemap examples/MECP2_CS.gb -t NM_004992.4
```

## Example Output

![MECP2 gene with all splicing annotations](images/splicemap-overview.png)
*Full MECP2 gene (76 kb) showing exons, introns, and ESE/ESS sites annotated by splicemap.*

![Exon 2 detail showing splice signals and ESE sites](images/splicemap-exon2-detail.png)
*Zoomed into Exon 2: branch point (orange), PPT (gold), 3'SS, and ESE sites (SRSF1/2/5/6) in distinct colors.*

```
Found 3 intron(s)
  intron 1: 116-5411 (5,296 bp)
  intron 2: 5536-65161 (59,626 bp)
  intron 3: 65512-66267 (756 bp)

Splice Map: MECP2_CS (76,145 bp, linear)
============================================================

Intron                 Length     5'SS     3'SS   BPS(z)   PPT%
-------------------- -------- -------- -------- -------- ------
intron 1               5,296    -15.5      0.1      6.6    75%
intron 2              59,626    -15.6     -4.7      6.2    63%
intron 3                 756     10.1     12.4      6.6    80%

Exon                 Length  ESEfinder  hnRNP  ESRseq+  ESRseq-
-------------------- ------  ---------  -----  -------  -------
exon_us_intron_1       114         25      0       63       13
exon_ds_intron_1       124         22      0       24       37
exon_ds_intron_2       351         57      1      118       69
exon_ds_intron_3      9878       1436     55     2232     3125
```

The annotated GenBank file can be opened in SnapGene or UGENE to visualize all annotations with color-coded features.

## What Gets Annotated

| Element | Method | Color |
|---------|--------|-------|
| 5' Splice Site (5'SS) | MaxEntScan | Royal Blue |
| 3' Splice Site (3'SS) | MaxEntScan | Crimson |
| Branch Point (BPS) | BPP | Orange |
| Polypyrimidine Tract (PPT) | Pyrimidine density | Gold |
| SRSF1 binding (ESE) | ESEfinder | Violet |
| SRSF2 binding (ESE) | ESEfinder | Pink |
| SRSF5 binding (ESE) | ESEfinder | Amber |
| SRSF6 binding (ESE) | ESEfinder | Emerald |
| ESRseq enhancer (ESE) | ESRseq hexamer lookup | Green |
| ESRseq silencer (ESS) | ESRseq hexamer lookup | Orange |
| hnRNP A1 binding (ESS) | Motif matching | Red |
| hnRNP H binding (ESS) | G-run detection | Dark Red |

Splicemap reads exon annotations from the input GenBank file (it does not annotate exons or introns itself). Download your gene as a RefSeqGene from [NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/) to get a file with exon annotations already included.

## Methods and references

**Splice sites** scored with [MaxEntScan](https://pubmed.ncbi.nlm.nih.gov/15285897/) (Yeo & Burge 2004). Above 6 is strong, 3-6 moderate, below 3 weak. Non-canonical dinucleotides (not GT...AG) are flagged.

**Branch points** predicted with [BPP](https://github.com/zhqingit/BPP) (position weight matrix on verified human branch points). z-score above 2 is strong.

**ESEfinder** predicts binding sites for 4 SR proteins (SRSF1, SRSF2, SRSF5, SRSF6) using in vitro SELEX-derived matrices ([Cartegni et al. 2003](https://pubmed.ncbi.nlm.nih.gov/12824367/)). Tells you which protein likely binds where. Does not cover other SR proteins (SRSF3, SRSF7, Tra2-beta, RBFOX). ~44% accuracy on known splicing mutations.

**ESRseq** looks up every hexamer against experimentally measured splicing activity scores ([Ke et al. 2011](https://pubmed.ncbi.nlm.nih.gov/21659425/)). Each hexamer was tested in a minigene assay and scored by RNA-seq. Positive = promotes inclusion, negative = promotes skipping. Captures the net effect of all regulatory proteins, not just four. ~83% accuracy on known splicing mutations. Does not tell you which protein is responsible.

**hnRNP motifs** detect two silencer proteins by pattern matching: hnRNP A1 (Burd & Dreyfuss 1994) and hnRNP H G-runs (Caputi & Bhatt 2003). For broader silencer coverage, use ESRseq negative scores.

### Limitations

- Scans flat sequence only. RNA secondary structure is not considered.
- No positional weighting. ESEs near splice sites matter more than those in the exon center.
- No combinatorial effects between adjacent elements.
- No cell-type specificity. Splicing regulation varies between tissues.
- Silencer protein coverage is limited to hnRNP A1 and H by motif. ESRseq provides broader but protein-anonymous coverage.

## Commands

### Reading and inspection

| Command | Description |
|---------|-------------|
| `read <file>` | Parse .gb/.fasta, show clean summary |
| `features <file>` | List all annotations as a table |
| `seq <file> <start> <end>` | Extract a sequence region (1-based) |
| `translate <file> <start> <end>` | Extract a region and translate to protein |
| `search <file> <sequence>` | Find all occurrences of a DNA motif (both strands) |
| `orfs <file>` | Find ORFs |
| `sites <file>` | Find restriction sites |
| `current` | Show files currently open in UGENE |
| `open <file>` | Open file in default viewer (UGENE/SnapGene) |

### Annotation

| Command | Description |
|---------|-------------|
| `splicemap <file> -t <accession>` | Full splice map: discover exons, annotate all splicing elements |
| `exons <file> --transcript <accession>` | Find and optionally annotate exon boundaries from mRNA alignment |
| `annotate <file> <start> <end> <label>` | Add a feature annotation |
| `annotate-seq <file> <sequence> <label>` | Find a sequence and annotate its position |
| `splice-signals <file>` | Annotate 5'SS, 3'SS, BPS, PPT on detected or specified introns |
| `branchpoint <file>` | Predict branch point locations in introns |
| `remove <file> <label>` | Remove an annotation by label |

### Sequence editing

| Command | Description |
|---------|-------------|
| `insert <file> <position> <sequence>` | Insert a DNA sequence, shifting downstream features |
| `delete <file> <start> <end>` | Delete a region, shifting downstream features |
| `replace <file> <start> <end> <sequence>` | Replace a region with a new sequence |
| `revcomp <file>` | Reverse complement the sequence |

### Analysis

| Command | Description |
|---------|-------------|
| `diff <file1> <file2>` | Compare two constructs |
| `blast <file>` | Run a remote NCBI BLAST search |
| `stitch <file> [labels...]` | Extract annotated regions, stitch together, optionally translate |
| `check <file>` | Preflight validation of a construct |
| `gibson <file> --enzymes E1,E2 --insert SEQ` | Design Gibson assembly eBlocks |
| `varmap <file> <variants_csv>` | Visualize variant positions mapped onto a sequence |
| `export <file> <format>` | Convert between formats: fasta, genbank, tab (feature TSV) |

## Dependencies

```
Python 3.8+
pip install -r requirements.txt  # biopython, pydna
```

Branch point prediction tools (BPP, SVM-BPfinder) are automatically downloaded on first use. No manual setup required.

## External Services

- **NCBI Entrez**: Used to fetch mRNA transcripts for exon discovery. Requires internet connection.
- **NCBI BLAST**: Optional remote BLAST searches via the `blast` command.
- **IDT API**: Optional complexity screening for synthesis orders. Set `IDT_API_KEY` environment variable to enable.

## License

GPLv3. See LICENSE file.
