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

Exon                 Length  ESE sites  ESS sites
-------------------- ------  ---------  ---------
exon_us_intron_1       115         25          0
exon_ds_intron_1       124         24          0
exon_ds_intron_2       350         57          1
exon_ds_intron_3      9878       1436         55
```

The annotated GenBank file can be opened in SnapGene or UGENE to visualize all annotations with color-coded features.

## What Gets Annotated

| Element | Method | Color |
|---------|--------|-------|
| Exons | mRNA transcript alignment | Steel Blue |
| Introns | Derived from exon gaps | Grey |
| 5' Splice Site (5'SS) | MaxEntScan (Yeo & Burge 2004) | Royal Blue |
| 3' Splice Site (3'SS) | MaxEntScan (Yeo & Burge 2004) | Crimson |
| Branch Point (BPS) | BPP (position weight matrix) | Orange |
| Polypyrimidine Tract (PPT) | Pyrimidine density scoring | Gold |
| SRSF1 (ESE) | ESEfinder matrix | Violet |
| SRSF2 (ESE) | ESEfinder matrix | Pink |
| SRSF5 (ESE) | ESEfinder matrix | Amber |
| SRSF6 (ESE) | ESEfinder matrix | Emerald |
| hnRNP A1 (ESS) | Motif matching (TAGG[GT][TA]) | Red |
| hnRNP H (ESS) | G-run detection (GGGG+) | Dark Red |

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
