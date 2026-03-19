# splicemap

Annotate splicing regulatory elements on genomic DNA sequences. Given a GenBank file with exon annotations (or an mRNA accession to discover them), splicemap maps splice sites, branch points, polypyrimidine tracts, and exonic splicing enhancers and silencers. Annotations are color-coded and written back to the GenBank file for viewing in SnapGene or UGENE. A markdown report is generated alongside.

## Quick Start

```bash
git clone https://github.com/maxwraae/splicemap.git
cd splicemap
pip install -r requirements.txt
python splicemap.py splicemap examples/MECP2_CS.gb -t NM_004992.4
```

## What You See

Annotations are named by what they are, not which tool found them. Colors group by function. Shades indicate relative confidence. Double-click any annotation in SnapGene to see the tool, score, motif, and other details.

| Annotation | Color | What it is |
|------------|-------|------------|
| Branch Point | Orange (dark to light by rank) | Candidate branch point adenosine |
| Polypyrimidine Tract | Amber | Pyrimidine-rich region between BPS and 3'SS |
| 5' Splice Site | Teal | Donor splice site (GT) |
| 3' Splice Site | Teal (lighter) | Acceptor splice site (AG) |
| Splice Enhancer (blue) | Blue (shades) | SR protein binding site (ESEfinder) |
| Splice Enhancer (green) | Green | Hexamer with positive splicing activity (ESRseq) |
| Splice Silencer (red) | Red (shades) | Hexamer or motif that suppresses exon inclusion |

## Example Output

```
Splice Map: MECP2_CS (76,145 bp, linear)
============================================================

Intron                 Length     5'SS     3'SS   BPS(z)   PPT%  U-run
-------------------- -------- -------- -------- -------- ------ ------
intron 1               5,296      7.9     10.8      6.5    80%      4
intron 2              59,626     10.9      5.3      6.1    67%      3
intron 3                 756     10.1     12.4      6.6    80%      3

Exon                 Length  ESEfinder  hnRNP  ESRseq+  ESRseq-
-------------------- ------  ---------  -----  -------  -------
exon_us_intron_1       114         25      0       63       13
exon_ds_intron_1       124         22      0       24       37
exon_ds_intron_2       351         57      1      118       69
exon_ds_intron_3      9878       1436     55     2232     3125
```

## Methods

### Splice sites

[MaxEntScan](https://pubmed.ncbi.nlm.nih.gov/15285897/) (Yeo & Burge 2004). 5'SS scored on a 9-mer (3 exonic + 6 intronic), 3'SS on a 23-mer (20 intronic + 3 exonic). Log-odds scores. Above 6 is strong, 3-6 moderate, below 3 weak. Non-canonical dinucleotides are flagged.

### Branch points

[BPP](https://github.com/zhqingit/BPP) (PWM trained on verified human branch points) and [SVM-BPfinder](https://github.com/comprna/SVM-BPfinder-3M) (SVM classifier). Both run independently on each intron. Up to 4 candidates shown, ranked by score across both tools. Darker orange = higher confidence.

Branch point prediction is roughly 75-80% accurate. No tool reliably identifies the correct branch point across all intron contexts.

### Polypyrimidine tract

Defined as the region between the top branch point candidate and the 3'SS AG. Reports length, pyrimidine percentage, and longest uninterrupted U-run. The U-run is the most informative single feature for U2AF65 binding (crystal structures show its two RRM domains each grab 4-5 uridines).

No validated computational model for U2AF65 binding affinity exists. PPT is scored by composition, which is standard practice. The PPT window depends on the branch point prediction.

### Exonic splicing enhancers and silencers

Two methods, measuring different things.

**ESEfinder** ([Cartegni et al. 2003](https://pubmed.ncbi.nlm.nih.gov/12824367/)). Position weight matrices from SELEX experiments for four SR proteins: SRSF1, SRSF2, SRSF5, SRSF6. Tells you which protein binds where. Only covers 4 of ~12 SR proteins. In vitro binding preference does not always match in vivo function. ~44% accuracy on known splicing mutations.

**ESRseq** ([Ke et al. 2011](https://pubmed.ncbi.nlm.nih.gov/21659425/)). All 4,096 possible hexamers tested in a minigene assay and scored by RNA-seq. Positive score = promotes exon inclusion (enhancer). Negative = promotes skipping (silencer). Captures the combined effect of all proteins that bind a given sequence. Does not identify which protein is responsible. ~83% accuracy on known splicing mutations.

**hnRNP motifs.** Pattern matching for hnRNP A1 ([Burd & Dreyfuss 1994](https://pubmed.ncbi.nlm.nih.gov/7520568/)) and hnRNP H G-runs ([Caputi & Bhatt 2003](https://pubmed.ncbi.nlm.nih.gov/12554860/)).

### Limitations

- Flat sequence only. No RNA secondary structure.
- No positional weighting (ESEs near splice sites matter more than those mid-exon).
- No combinatorial effects between adjacent elements.
- No cell-type or tissue specificity.
- Branch point prediction accuracy is inherently limited. PPT analysis depends on it.

## Commands

### Reading and inspection

| Command | Description |
|---------|-------------|
| `read <file>` | Parse .gb/.fasta, show summary |
| `features <file>` | List all annotations |
| `seq <file> <start> <end>` | Extract sequence (1-based) |
| `translate <file> <start> <end>` | Translate a region |
| `search <file> <sequence>` | Find motif occurrences (both strands) |
| `orfs <file>` | Find open reading frames |
| `sites <file>` | Find restriction sites |
| `open <file>` | Open in default viewer |

### Annotation

| Command | Description |
|---------|-------------|
| `splicemap <file> -t <accession>` | Full splice map |
| `exons <file> -t <accession>` | Find and annotate exon boundaries |
| `annotate <file> <start> <end> <label>` | Add a feature |
| `annotate-seq <file> <sequence> <label>` | Find and annotate a sequence |
| `splice-signals <file>` | Annotate splice signals on detected introns |
| `branchpoint <file>` | Predict branch points |
| `remove <file> <label>` | Remove an annotation |

### Sequence editing

| Command | Description |
|---------|-------------|
| `insert <file> <pos> <seq>` | Insert sequence, shift features |
| `delete <file> <start> <end>` | Delete region, shift features |
| `replace <file> <start> <end> <seq>` | Replace region |
| `revcomp <file>` | Reverse complement |

### Analysis

| Command | Description |
|---------|-------------|
| `diff <file1> <file2>` | Compare two constructs |
| `blast <file>` | Remote NCBI BLAST |
| `stitch <file> [labels...]` | Stitch regions, optionally translate |
| `check <file>` | Preflight validation |
| `gibson <file> --enzymes E1,E2 --insert SEQ` | Design Gibson assembly |
| `varmap <file> <variants_csv>` | Map variant positions |
| `export <file> <format>` | Convert format (fasta, genbank, tab) |

## Dependencies

```
Python 3.8+
pip install -r requirements.txt  # biopython, pydna
```

Branch point tools (BPP, SVM-BPfinder) are downloaded automatically on first use.

## License

GPLv3
