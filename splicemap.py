#!/usr/bin/env python3
"""
splicemap - Splice site annotation and analysis for GenBank files.

Commands:
  current                           - Show files open in UGENE
  read <file>                       - Parse .gb/.fasta, show clean summary
  features <file>                   - List all annotations as a table
  seq <file> <start> <end>          - Extract a sequence region (1-based)
  translate <file> <start> <end>    - Extract a region and translate to protein
  annotate <file> <start> <end> <label> [--type TYPE]
                                    - Add a feature annotation, write back
  remove <file> <label>             - Remove annotation by label, write back
  orfs <file> [--min-length N]      - Find ORFs via ugenecl
  sites <file> [--enzymes E1,E2]    - Find restriction sites
  revcomp <file> [--out FILE]       - Reverse complement, write to file
  diff <file1> <file2>              - Compare two constructs
  search <file> <sequence> [--fwd] [--rc]
                                    - Find all occurrences of a DNA motif (both strands)
  annotate-seq <file> <sequence> <label> [--antisense] [--type TYPE] [--color COLOR] [--note NOTE]
                                    - Find a sequence and annotate its position
  open <file>                       - Open file in default viewer (UGENE/SnapGene)
  insert <file> <position> <sequence> [--label LABEL] [--type TYPE]
                                    - Insert a DNA sequence at a given position, shift features, write back
  delete <file> <start> <end>       - Delete a DNA region, shift downstream features, write back
  replace <file> <start> <end> <sequence> [--label LABEL] [--type TYPE]
                                    - Replace a DNA region with a new sequence, adjusting features
  export <file> <format> [--out FILE]
                                    - Convert between formats: fasta, genbank (or gb), tab (feature TSV)
  blast <file> [--db DATABASE] [--program PROGRAM] [--evalue EVALUE] [--max-hits N] [--region START END]
                                    - Run a remote NCBI BLAST search and display parsed results
  exons <file> --transcript ACCESSION [--annotate]
                                    - Find exon boundaries by aligning mRNA against genomic sequence
  stitch <file> [labels...] [--exons] [--translate]
                                    - Extract regions by annotation, stitch together, optionally translate
  gibson <file> --enzymes E1,E2 --insert SEQ [--label NAME]
                                    - Design Gibson assembly eBlocks from backbone + insert
  check <file> [--backbone FILE] [--enzymes E1,E2]
                                    - Preflight validation of constructs
  branchpoint <file> [--intron LABEL | --region S E | --seq SEQ]
                                    - Predict branch point locations in introns
  splice-signals <file> [--intron START END]
                                    - Annotate splice signals (5'SS, 3'SS, BPS, PPT) on an intron.
                                      If --intron not provided, introns are auto-detected from exon annotations.
  varmap <file> <variants_csv> [--region S E] [--intron LABEL]
                                    - Visualize variant positions mapped onto a DNA sequence
"""

import argparse
import csv
import os
import shutil
import subprocess
import sys
from io import StringIO

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from pydna.genbankfixer import gbtext_clean

UGENECL = os.environ.get("UGENECL_PATH", "/Applications/Unipro UGENE.app/Contents/MacOS/ugenecl")

# Relative adaptiveness (fraction of max frequency for each amino acid)
# Source: Kazusa Codon Usage Database, Homo sapiens
HUMAN_CODON_FREQ = {
    'TTT': 0.45, 'TTC': 0.55, 'TTA': 0.07, 'TTG': 0.13,
    'CTT': 0.13, 'CTC': 0.20, 'CTA': 0.07, 'CTG': 0.41,
    'ATT': 0.36, 'ATC': 0.48, 'ATA': 0.16, 'ATG': 1.00,
    'GTT': 0.18, 'GTC': 0.24, 'GTA': 0.11, 'GTG': 0.47,
    'TAT': 0.43, 'TAC': 0.57, 'TAA': 0.28, 'TAG': 0.20,
    'CAT': 0.41, 'CAC': 0.59, 'CAA': 0.25, 'CAG': 0.75,
    'AAT': 0.46, 'AAC': 0.54, 'AAA': 0.42, 'AAG': 0.58,
    'GAT': 0.46, 'GAC': 0.54, 'GAA': 0.42, 'GAG': 0.58,
    'TCT': 0.15, 'TCC': 0.22, 'TCA': 0.12, 'TCG': 0.06,
    'CCT': 0.28, 'CCC': 0.33, 'CCA': 0.27, 'CCG': 0.11,
    'ACT': 0.24, 'ACC': 0.36, 'ACA': 0.28, 'ACG': 0.12,
    'GCT': 0.26, 'GCC': 0.40, 'GCA': 0.23, 'GCG': 0.11,
    'TGT': 0.45, 'TGC': 0.55, 'TGA': 0.52, 'TGG': 1.00,
    'CGT': 0.08, 'CGC': 0.19, 'CGA': 0.11, 'CGG': 0.21,
    'AGT': 0.15, 'AGC': 0.24, 'AGA': 0.20, 'AGG': 0.20,
    'GGT': 0.16, 'GGC': 0.34, 'GGA': 0.25, 'GGG': 0.25,
}

# ── ESEfinder matrices (Cartegni et al. 2003, NAR) ─────────────────────────
# Position weight matrices for SR protein binding site detection.
# Scores are log-odds vs background. Sliding window sum >= threshold = hit.
# Source: ESEfinder 3.0 (esefinder.ahc.umn.edu)
ESEFINDER_MATRICES = {
    'SRSF1': {
        'alt_name': 'SF2/ASF',
        'matrix': [
            {'A': -1.14, 'C': 1.37, 'G': -0.21, 'T': -1.58},
            {'A': 0.62,  'C': -1.10, 'G': 0.17,  'T': -0.50},
            {'A': -1.58, 'C': 0.73, 'G': 0.48,  'T': -1.58},
            {'A': 1.32,  'C': 0.33, 'G': -1.58, 'T': -1.13},
            {'A': -1.58, 'C': 0.94, 'G': 0.33,  'T': -1.58},
            {'A': -1.58, 'C': -1.58, 'G': 0.99, 'T': -1.13},
            {'A': 0.62,  'C': -1.58, 'G': -0.11, 'T': 0.27},
        ],
        'threshold': 1.956,
        'color': '#1D4ED8',
    },
    'SRSF2': {
        'alt_name': 'SC35',
        'matrix': [
            {'A': -0.88, 'C': -1.16, 'G': 0.87, 'T': -1.18},
            {'A': 0.09,  'C': -1.58, 'G': 0.45, 'T': -0.20},
            {'A': -0.06, 'C': 0.95,  'G': -1.36, 'T': 0.38},
            {'A': -1.58, 'C': 1.11,  'G': -1.58, 'T': 0.88},
            {'A': 0.09,  'C': 0.56,  'G': -0.33, 'T': -0.20},
            {'A': -0.41, 'C': 0.86,  'G': -0.05, 'T': -0.86},
            {'A': -0.06, 'C': 0.32,  'G': -1.36, 'T': 0.96},
            {'A': 0.23,  'C': -1.58, 'G': 0.68,  'T': -1.58},
        ],
        'threshold': 2.383,
        'color': '#2563EB',
    },
    'SRSF5': {
        'alt_name': 'SRp40',
        'matrix': [
            {'A': -0.13, 'C': 0.56, 'G': -1.58, 'T': 0.92},
            {'A': -1.58, 'C': 0.68, 'G': -0.14, 'T': 0.37},
            {'A': 1.28,  'C': -1.12, 'G': -1.33, 'T': 0.23},
            {'A': -0.33, 'C': 1.24, 'G': -0.48, 'T': -1.14},
            {'A': 0.97,  'C': -0.77, 'G': -1.58, 'T': 0.72},
            {'A': -0.13, 'C': 0.13, 'G': 0.44,  'T': -1.58},
            {'A': -1.58, 'C': -0.05, 'G': 0.80,  'T': -1.58},
        ],
        'threshold': 2.670,
        'color': '#3B82F6',
    },
    'SRSF6': {
        'alt_name': 'SRp55',
        'matrix': [
            {'A': -0.66, 'C': 0.39, 'G': -1.58, 'T': 1.22},
            {'A': 0.11,  'C': -1.58, 'G': 0.72,  'T': -1.58},
            {'A': -0.66, 'C': 1.48, 'G': -1.58, 'T': -0.07},
            {'A': 0.11,  'C': -1.58, 'G': 0.72,  'T': -1.58},
            {'A': -1.58, 'C': -1.58, 'G': 0.21,  'T': 1.02},
            {'A': 0.61,  'C': 0.98, 'G': -0.79,  'T': -1.58},
        ],
        'threshold': 2.676,
        'color': '#60A5FA',
    },
}

# ── Splicemap color palette ────────────────────────────────────────────────
SPLICEMAP_COLORS = {
    # Splice sites: teal (distinct from enhancer blue/green)
    '5SS':      '#0D9488',  # teal (donor)
    '3SS':      '#14B8A6',  # lighter teal (acceptor)
    # Branch machinery: orange family
    'BPS':      '#C2410C',  # dark orange (branch point)
    'PPT':      '#F59E0B',  # amber (polypyrimidine tract)
    # ESEfinder: blue family (enhancers, method 1)
    'SRSF1':    '#1D4ED8',  # dark blue
    'SRSF2':    '#2563EB',  # medium blue
    'SRSF5':    '#3B82F6',  # blue
    'SRSF6':    '#60A5FA',  # light blue
    # ESRseq: green (enhancers, method 2) / red (silencers)
    'ESRseq_ESE': '#16A34A',  # green
    'ESRseq_ESS': '#DC2626',  # red
    # hnRNP motifs: red family (silencers)
    'hnRNP_A1': '#EF4444',  # red
    'hnRNP_H':  '#B91C1C',  # dark red
}

# ── hnRNP binding motifs (ESSs) ────────────────────────────────────────────
# Simple regex patterns for exonic splicing silencer detection.
# hnRNP A1: UAGGGU/A consensus (Burd & Dreyfuss 1994)
# hnRNP H: G-runs of 4+ (Caputi & Bhatt 2003)
HNRNP_MOTIFS = {
    'hnRNP_A1': {
        'pattern': r'TAGG[GT][TA]',
        'color': '#EF4444',
        'description': 'hnRNP A1 binding site (ESS)',
    },
    'hnRNP_H': {
        'pattern': r'G{4,}',
        'color': '#B91C1C',
        'description': 'hnRNP H binding site (G-run ESS)',
    },
}

# ── ESRseq hexamer lookup (Ke et al., Genome Research 2011) ───────────────────
# Measured splicing activity of all hexamers in a minigene assay.
# Positive score = promotes exon inclusion (ESE). Negative = promotes skipping (ESS).
# 83% accuracy on known splicing mutations vs ESEfinder's 44%.
_ESRSEQ_ESE = None  # dict: hexamer -> score (positive)
_ESRSEQ_ESS = None  # dict: hexamer -> score (negative)

def _load_esrseq():
    """Load ESRseq hexamer tables (lazy, once)."""
    global _ESRSEQ_ESE, _ESRSEQ_ESS
    if _ESRSEQ_ESE is not None:
        return
    tools_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "tools", "esrseq")
    _ESRSEQ_ESE = {}
    _ESRSEQ_ESS = {}
    for fname, target in [("ESE.txt", _ESRSEQ_ESE), ("ESS.txt", _ESRSEQ_ESS)]:
        path = os.path.join(tools_dir, fname)
        if os.path.exists(path):
            with open(path) as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) == 2:
                        target[parts[0].upper()] = float(parts[1])


def _find_esrseq_sites(exon_seq):
    """Scan exon sequence for ESRseq hexamers (ESE and ESS).

    Returns:
        list of dicts with 'type' ('ESE'/'ESS'), 'start_0', 'end_0', 'score', 'motif'
    """
    _load_esrseq()
    seq = str(exon_seq).upper()
    hits = []
    for i in range(len(seq) - 5):
        hexamer = seq[i:i + 6]
        if hexamer in _ESRSEQ_ESE:
            hits.append({
                'type': 'ESE',
                'start_0': i,
                'end_0': i + 6,
                'score': _ESRSEQ_ESE[hexamer],
                'motif': hexamer,
            })
        elif hexamer in _ESRSEQ_ESS:
            hits.append({
                'type': 'ESS',
                'start_0': i,
                'end_0': i + 6,
                'score': _ESRSEQ_ESS[hexamer],
                'motif': hexamer,
            })
    return hits


ESRSEQ_COLORS = {
    'ESE': '#22C55E',   # green (enhancer, method 2)
    'ESS': '#DC2626',   # red (silencer)
}

# Lazy-loaded maxentpy for splice site scoring
_MAXENT_MATRIX5 = None
_MAXENT_MATRIX3 = None

def _load_maxent():
    """Load MaxEntScan scoring matrices (lazy, once)."""
    global _MAXENT_MATRIX5, _MAXENT_MATRIX3
    if _MAXENT_MATRIX5 is None:
        tools_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'tools')
        import sys
        if tools_dir not in sys.path:
            sys.path.insert(0, tools_dir)
        from maxentpy.maxent import load_matrix5, load_matrix3
        _MAXENT_MATRIX5 = load_matrix5()
        _MAXENT_MATRIX3 = load_matrix3()
    return _MAXENT_MATRIX5, _MAXENT_MATRIX3


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def read_genbank(path):
    """Read a GenBank file, auto-repairing malformed files via pydna."""
    path = os.path.expanduser(path)
    if not os.path.exists(path):
        print(f"Error: file not found: {path}", file=sys.stderr)
        sys.exit(1)

    ext = os.path.splitext(path)[1].lower()

    if ext in (".gb", ".gbk", ".genbank"):
        with open(path, "r") as fh:
            raw = fh.read()
        try:
            cleaned = gbtext_clean(raw)
            record = SeqIO.read(StringIO(cleaned), "genbank")
        except Exception:
            # Fallback: try reading raw
            try:
                record = SeqIO.read(path, "genbank")
            except Exception as e:
                print(f"Error parsing {path}: {e}", file=sys.stderr)
                sys.exit(1)
    elif ext == ".dna":
        try:
            from snapgene_reader import snapgene_file_to_dict
            d = snapgene_file_to_dict(path)
            seq_str = d.get("seq", "")
            record = SeqRecord(Seq(seq_str), id=os.path.splitext(os.path.basename(path))[0],
                               name=os.path.splitext(os.path.basename(path))[0],
                               description="")
            # Topology
            is_circular = d.get("is_circular", False)
            record.annotations["topology"] = "circular" if is_circular else "linear"
            record.annotations["molecule_type"] = "DNA"
            # Features
            for feat_dict in d.get("features", []):
                strand_raw = feat_dict.get("strand", "+")
                if strand_raw == "+":
                    strand = 1
                elif strand_raw == "-":
                    strand = -1
                else:
                    strand = 0
                start = feat_dict.get("start", 0)
                end = feat_dict.get("end", 0)
                feat_type = feat_dict.get("type", "misc_feature")
                name = feat_dict.get("name", "")
                qualifiers = dict(feat_dict.get("qualifiers", {})) if feat_dict.get("qualifiers") else {}
                if name and "label" not in qualifiers:
                    qualifiers["label"] = [name]
                location = FeatureLocation(start, end, strand=strand)
                sf = SeqFeature(location, type=feat_type, qualifiers=qualifiers)
                record.features.append(sf)
        except Exception as e:
            print(f"Error parsing SnapGene .dna file {path}: {e}", file=sys.stderr)
            sys.exit(1)
    elif ext in (".fa", ".fasta", ".fna"):
        try:
            record = SeqIO.read(path, "fasta")
        except Exception as e:
            print(f"Error parsing {path}: {e}", file=sys.stderr)
            sys.exit(1)
    else:
        # Try genbank first, then fasta
        try:
            with open(path, "r") as fh:
                raw = fh.read()
            cleaned = gbtext_clean(raw)
            record = SeqIO.read(StringIO(cleaned), "genbank")
        except Exception:
            try:
                record = SeqIO.read(path, "fasta")
            except Exception as e:
                print(f"Error parsing {path}: {e}", file=sys.stderr)
                sys.exit(1)

    return record


def write_genbank(record, path):
    """Write a GenBank record, saving a .bak first. Ensures required annotations."""
    path = os.path.expanduser(path)

    # Sanitize name/id for GenBank LOCUS line (no spaces allowed)
    record.name = record.name.replace(" ", "_")[:16]
    record.id = record.id.replace(" ", "_")[:16]

    # Ensure molecule_type is set (required by BioPython GenBank writer)
    if "molecule_type" not in record.annotations:
        record.annotations["molecule_type"] = "DNA"

    # Preserve topology if not set
    if "topology" not in record.annotations:
        record.annotations["topology"] = "linear"

    # Save backup
    bak_path = path + ".bak"
    if os.path.exists(path):
        shutil.copy2(path, bak_path)

    with open(path, "w") as fh:
        SeqIO.write(record, fh, "genbank")


def open_in_viewer(path):
    """Open a file in the default application (SnapGene/UGENE)."""
    path = os.path.expanduser(path)
    if sys.platform == "darwin":
        subprocess.run(["open", path], check=False)
    elif sys.platform == "linux":
        subprocess.run(["xdg-open", path], check=False, stderr=subprocess.DEVNULL)


def get_label(feature):
    """Extract label from a feature's qualifiers."""
    for key in ("label", "gene", "product", "note", "name"):
        val = feature.qualifiers.get(key)
        if val:
            return val[0] if isinstance(val, list) else val
    return ""


def _truncate_seq(seq, max_len=10):
    """Truncate a sequence string for display in labels."""
    if len(seq) <= max_len:
        return seq
    half = max_len // 2
    return seq[:half] + "..." + seq[-half:]


def _add_edit_mark(record, start_0based, end_0based, label, note):
    """Add a bright red CDS annotation marking an edit region.

    Parameters are 0-based half-open coordinates (BioPython convention).
    """
    # Ensure we have at least a 1bp annotation
    if end_0based <= start_0based:
        end_0based = start_0based + 1
    # Clamp to sequence length
    seq_len = len(record.seq)
    if end_0based > seq_len:
        end_0based = seq_len
    if start_0based >= seq_len:
        start_0based = max(0, seq_len - 1)

    edit_feature = SeqFeature(
        FeatureLocation(start_0based, end_0based, strand=1),
        type="CDS",
        qualifiers={
            "label": [label],
            "ApEinfo_fwdcolor": ["#FF0000"],
            "ApEinfo_revcolor": ["#FF0000"],
            "note": [note],
        }
    )
    record.features.append(edit_feature)


def strand_str(strand):
    if strand == 1:
        return "+"
    elif strand == -1:
        return "-"
    return "."


# ---------------------------------------------------------------------------
# Commands
# ---------------------------------------------------------------------------

def cmd_current(args):
    """Show files open in UGENE via lsof."""
    try:
        result = subprocess.run(
            ["lsof", "-c", "UGENE", "-F", "n"],
            capture_output=True, text=True
        )
        lines = result.stdout.splitlines()
        files = []
        for line in lines:
            if line.startswith("n") and len(line) > 1:
                f = line[1:]
                if any(f.endswith(ext) for ext in (".gb", ".gbk", ".fa", ".fasta", ".dna", ".snapgene")):
                    files.append(f)
        if files:
            print("Open in UGENE:")
            for f in files:
                print(f"  {f}")
        else:
            print("No GenBank/FASTA files detected open in UGENE.")
    except FileNotFoundError:
        print("lsof not available.", file=sys.stderr)
        sys.exit(1)


def cmd_read(args):
    """Parse a .gb/.fasta file and print a clean summary."""
    record = read_genbank(args.file)

    name = record.name or record.id or os.path.basename(args.file)
    length = len(record.seq)
    topology = record.annotations.get("topology", "linear")
    mol_type = record.annotations.get("molecule_type", "DNA")
    n_features = len(record.features)

    # Filter out the implicit 'source' feature from display count
    display_features = [f for f in record.features if f.type != "source"]

    print(f"Name:     {name}")
    print(f"Length:   {length} bp")
    print(f"Topology: {topology}")
    print(f"Molecule: {mol_type}")
    print(f"Features ({len(display_features)}):")

    for feat in display_features:
        label = get_label(feat)
        start = int(feat.location.start)
        end = int(feat.location.end)
        strand = strand_str(feat.location.strand)
        feat_type = feat.type.ljust(12)
        loc = f"[{start}:{end}]".ljust(16)
        print(f"  {feat_type}  {loc}  {strand}  {label}")


def cmd_features(args):
    """List all annotations as a table."""
    record = read_genbank(args.file)

    # Header
    print(f"{'#':<4}  {'Type':<14}  {'Start':>7}  {'End':>7}  {'Strand':>6}  {'Size':>7}  Label")
    print("-" * 80)

    for i, feat in enumerate(record.features):
        if feat.type == "source":
            continue
        label = get_label(feat)
        start = int(feat.location.start) + 1  # Convert to 1-based
        end = int(feat.location.end)            # Already inclusive in 1-based
        length = end - start + 1
        strand = strand_str(feat.location.strand)
        print(f"{i:<4}  {feat.type:<14}  {start:>7}  {end:>7}  {strand:>6}  {length:>5} bp  {label}")


def cmd_seq(args):
    """Extract and display a sequence region (1-based, inclusive)."""
    record = read_genbank(args.file)

    start = args.start - 1  # Convert to 0-based
    end = args.end

    seq_len = len(record.seq)
    if start < 0 or end > seq_len or start >= end:
        print(f"Error: coordinates out of range (sequence is {seq_len} bp)", file=sys.stderr)
        sys.exit(1)

    subseq = str(record.seq[start:end])
    name = record.name or record.id

    print(f">{name}:{args.start}-{args.end} ({end - start} bp)")

    # Print in 60-char lines
    for i in range(0, len(subseq), 60):
        print(subseq[i:i+60])


def cmd_translate(args):
    """Extract a DNA region and translate to protein (1-based, inclusive)."""
    record = read_genbank(args.file)

    start = args.start - 1  # Convert to 0-based
    end = args.end

    seq_len = len(record.seq)
    if start < 0 or end > seq_len or start >= end:
        print(f"Error: coordinates out of range (sequence is {seq_len} bp)", file=sys.stderr)
        sys.exit(1)

    subseq = record.seq[start:end]

    # Apply reverse complement if requested
    if args.reverse:
        subseq = subseq.reverse_complement()

    # Apply frame offset (1-based frame: 1, 2, 3 → offset 0, 1, 2)
    frame_offset = args.frame - 1
    coding_seq = subseq[frame_offset:]

    region_len = len(subseq)
    coding_len = len(coding_seq)
    remainder = coding_len % 3
    n_codons = coding_len // 3
    n_aa = n_codons  # before stop codon trimming

    strand_label = "RC" if args.reverse else "+"
    frame_label = f"{strand_label}{args.frame}"

    print(f"Region: {args.start}-{args.end} ({region_len} bp)")
    if args.frame != 1 or args.reverse:
        effective_len = coding_len - remainder
        print(f"Frame: {frame_label} (offset {frame_offset}, translating {effective_len} bp = {n_codons} codons)")
    else:
        effective_len = coding_len - remainder
        print(f"Frame: {frame_label} ({effective_len} bp = {n_codons} codons)")

    if remainder != 0:
        print(f"Warning: {remainder} bp incomplete codon at end (ignored)")

    # Translate
    protein = coding_seq[:n_codons * 3].translate(table=1)
    protein_str = str(protein)

    # Find stop codons
    stop_positions = [i + 1 for i, aa in enumerate(protein_str) if aa == "*"]

    aa_count = len(protein_str)
    print(f"Length: {aa_count} aa")
    print()

    # Print DNA
    dna_str = str(subseq)
    print("DNA:")
    for i in range(0, len(dna_str), 60):
        print(dna_str[i:i+60])
    print()

    # Print protein
    print("Protein:")
    for i in range(0, len(protein_str), 60):
        print(protein_str[i:i+60])

    # Stop codon report
    if stop_positions:
        stops_str = ", ".join(f"position {p}" for p in stop_positions)
        print(f"\nStop codons: {stops_str}")
    else:
        print("\nNo stop codons.")


def cmd_annotate(args):
    """Add a feature annotation and write back to file."""
    record = read_genbank(args.file)

    start = args.start - 1  # Convert to 0-based
    end = args.end
    seq_len = len(record.seq)

    if start < 0 or end > seq_len or start >= end:
        print(f"Error: coordinates out of range (sequence is {seq_len} bp)", file=sys.stderr)
        sys.exit(1)

    feat_type = args.type if args.type else "misc_feature"

    qualifiers = {
        "label": [args.label]
    }
    if args.color:
        qualifiers["ApEinfo_fwdcolor"] = [args.color]
        qualifiers["ApEinfo_revcolor"] = [args.color]
    if args.note:
        qualifiers["note"] = [args.note]

    new_feature = SeqFeature(
        FeatureLocation(start, end, strand=1),
        type=feat_type,
        qualifiers=qualifiers
    )

    record.features.append(new_feature)

    write_genbank(record, args.file)
    print(f"Added {feat_type} '{args.label}' at [{start}:{end}] to {args.file}")
    print(f"Backup saved to {args.file}.bak")
    open_in_viewer(args.file)


def cmd_annotate_seq(args):
    """Find a sequence in a GenBank file and annotate its position."""
    record = read_genbank(args.file)

    query = args.sequence.upper()

    # If --antisense, reverse complement the input sequence before searching
    complement = str.maketrans("ACGT", "TGCA")
    if args.antisense:
        query = query.translate(complement)[::-1]

    seq = str(record.seq).upper()
    seq_len = len(seq)
    rc_query = query.translate(complement)[::-1]
    context_bp = 10

    matches = []

    # Search forward strand
    pos = seq.find(query)
    while pos != -1:
        matches.append((pos, 1))  # (0-based start, strand)
        pos = seq.find(query, pos + 1)

    # Search reverse complement strand (skip if palindrome)
    if rc_query != query:
        pos = seq.find(rc_query)
        while pos != -1:
            matches.append((pos, -1))
            pos = seq.find(rc_query, pos + 1)

    if len(matches) == 0:
        print(f"No match found for {args.sequence}", file=sys.stderr)
        sys.exit(1)

    if len(matches) > 1:
        # Print all match positions with context, like cmd_search
        name = record.name or record.id or os.path.basename(args.file)
        print(f"Multiple matches found for: {query} ({len(query)} bp) in {name}\n")
        print(f"  {'#':>2}  {'Position':>8}  {'Strand':>6}  Context")
        for i, (idx0, strand) in enumerate(sorted(matches, key=lambda x: x[0]), 1):
            pos1 = idx0 + 1
            strand_label = "+" if strand == 1 else "-"
            ctx_start = max(0, idx0 - context_bp)
            ctx_end = min(seq_len, idx0 + len(query) + context_bp)
            left = seq[ctx_start:idx0]
            middle = seq[idx0:idx0 + len(query)]
            right = seq[idx0 + len(query):ctx_end]
            left_display = ("..." + left) if ctx_start > 0 else left
            right_display = (right + "...") if ctx_end < seq_len else right
            context_str = left_display + middle + right_display
            print(f"  {i:>2}  {pos1:>8}  {strand_label:>6}  {context_str}")
        print(f"\nMultiple matches found. Use 'bio annotate <file> <start> <end> <label>' with specific coordinates.")
        sys.exit(1)

    # Exactly one match
    idx0, strand = matches[0]
    start = idx0          # 0-based start (BioPython FeatureLocation convention)
    end = idx0 + len(query)  # 0-based exclusive end

    feat_type = args.type if args.type else "misc_feature"

    qualifiers = {
        "label": [args.label]
    }
    if args.color:
        qualifiers["ApEinfo_fwdcolor"] = [args.color]
        qualifiers["ApEinfo_revcolor"] = [args.color]
    if args.note:
        qualifiers["note"] = [args.note]

    new_feature = SeqFeature(
        FeatureLocation(start, end, strand=strand),
        type=feat_type,
        qualifiers=qualifiers
    )

    record.features.append(new_feature)

    write_genbank(record, args.file)
    strand_label = "+" if strand == 1 else "-"
    print(f"Added {feat_type} '{args.label}' at [{start + 1}:{end}] (strand {strand_label}) to {args.file}")
    print(f"Backup saved to {args.file}.bak")
    open_in_viewer(args.file)


def cmd_remove(args):
    """Remove annotation by label and write back to file."""
    record = read_genbank(args.file)

    before = len(record.features)
    record.features = [
        f for f in record.features
        if get_label(f) != args.label
    ]
    after = len(record.features)
    removed = before - after

    if removed == 0:
        print(f"No features with label '{args.label}' found.")
        sys.exit(1)

    write_genbank(record, args.file)
    print(f"Removed {removed} feature(s) with label '{args.label}' from {args.file}")
    print(f"Backup saved to {args.file}.bak")
    open_in_viewer(args.file)


def cmd_orfs(args):
    """Find ORFs using ugenecl."""
    if not os.path.exists(UGENECL):
        print(f"Error: ugenecl not found at {UGENECL}", file=sys.stderr)
        sys.exit(1)

    min_len = args.min_length if args.min_length else 100

    cmd = [
        UGENECL,
        "--task=find-orfs",
        f"--in={os.path.expanduser(args.file)}",
        f"--min-length={min_len}"
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True)
        print(result.stdout)
        if result.stderr:
            print(result.stderr, file=sys.stderr)
    except Exception as e:
        print(f"Error running ugenecl: {e}", file=sys.stderr)
        sys.exit(1)


def cmd_sites(args):
    """Find restriction sites in a sequence."""
    from Bio.Restriction import RestrictionBatch, Analysis, AllEnzymes

    record = read_genbank(args.file)

    if args.enzymes:
        enzyme_names = [e.strip() for e in args.enzymes.split(",")]
        try:
            rb = RestrictionBatch(enzyme_names)
        except Exception as e:
            print(f"Error loading enzymes: {e}", file=sys.stderr)
            sys.exit(1)
    else:
        # Default to common cloning enzymes
        common = ["EcoRI", "BamHI", "HindIII", "NotI", "XhoI", "SalI", "NcoI",
                  "NheI", "XbaI", "SpeI", "PstI", "SacI", "KpnI", "AvaI",
                  "ClaI", "EcoRV", "SmaI"]
        rb = RestrictionBatch(common)

    analysis = Analysis(rb, record.seq, linear=False)
    results = analysis.full()

    name = record.name or record.id
    print(f"Restriction sites in {name} ({len(record.seq)} bp, circular):\n")

    found = {enz: sites for enz, sites in results.items() if sites}
    not_found = [str(enz) for enz, sites in results.items() if not sites]

    if found:
        print(f"{'Enzyme':<12}  {'Cuts':>5}  Sites")
        print("-" * 50)
        for enz, sites in sorted(found.items(), key=lambda x: str(x[0])):
            sites_str = ", ".join(str(s) for s in sites)
            print(f"{str(enz):<12}  {len(sites):>5}  {sites_str}")
    else:
        print("None of the queried enzymes cut this sequence.")

    if not_found and args.enzymes:
        # Only show "not found" when user specified enzymes explicitly
        print(f"\nNo cuts: {', '.join(not_found)}")


def cmd_revcomp(args):
    """Reverse complement a sequence and write to file."""
    record = read_genbank(args.file)

    rc_record = record.reverse_complement(id=True, name=True, description=True, annotations=True, features=True)

    if args.out:
        out_path = os.path.expanduser(args.out)
    else:
        base, ext = os.path.splitext(os.path.expanduser(args.file))
        out_path = base + "_rc" + ext

    # Ensure molecule_type for write
    if "molecule_type" not in rc_record.annotations:
        rc_record.annotations["molecule_type"] = "DNA"

    write_genbank(rc_record, out_path)
    print(f"Reverse complement written to {out_path}")
    open_in_viewer(out_path)


def cmd_diff(args):
    """Compare two constructs — sequence and feature differences."""
    record1 = read_genbank(args.file1)
    record2 = read_genbank(args.file2)

    name1 = record1.name or os.path.basename(args.file1)
    name2 = record2.name or os.path.basename(args.file2)

    print(f"Diff: {name1}  vs  {name2}\n")

    # Length comparison
    len1, len2 = len(record1.seq), len(record2.seq)
    if len1 != len2:
        print(f"Length:  {name1} = {len1} bp,  {name2} = {len2} bp  (delta: {len2 - len1:+d})")
    else:
        print(f"Length:  {len1} bp (identical)")

    # Topology
    topo1 = record1.annotations.get("topology", "linear")
    topo2 = record2.annotations.get("topology", "linear")
    if topo1 != topo2:
        print(f"Topology: {name1}={topo1}, {name2}={topo2}")
    else:
        print(f"Topology: {topo1} (same)")

    # Sequence identity
    min_len = min(len1, len2)
    mismatches = sum(1 for a, b in zip(str(record1.seq[:min_len]), str(record2.seq[:min_len])) if a != b)
    pct = 100.0 * (min_len - mismatches) / min_len if min_len else 0
    print(f"Identity: {pct:.2f}%  ({mismatches} mismatches in first {min_len} bp)")

    # Feature comparison
    labels1 = {get_label(f): f for f in record1.features if f.type != "source" and get_label(f)}
    labels2 = {get_label(f): f for f in record2.features if f.type != "source" and get_label(f)}

    only_in_1 = sorted(set(labels1) - set(labels2))
    only_in_2 = sorted(set(labels2) - set(labels1))
    in_both = sorted(set(labels1) & set(labels2))

    print(f"\nFeatures:")
    print(f"  {name1}: {len(labels1)} features")
    print(f"  {name2}: {len(labels2)} features")

    if only_in_1:
        print(f"\nOnly in {name1}:")
        for lbl in only_in_1:
            f = labels1[lbl]
            start, end = int(f.location.start), int(f.location.end)
            print(f"  - {f.type:<14} [{start}:{end}]  {lbl}")

    if only_in_2:
        print(f"\nOnly in {name2}:")
        for lbl in only_in_2:
            f = labels2[lbl]
            start, end = int(f.location.start), int(f.location.end)
            print(f"  + {f.type:<14} [{start}:{end}]  {lbl}")

    if in_both:
        moved = []
        for lbl in in_both:
            f1, f2 = labels1[lbl], labels2[lbl]
            s1, e1 = int(f1.location.start), int(f1.location.end)
            s2, e2 = int(f2.location.start), int(f2.location.end)
            if s1 != s2 or e1 != e2:
                moved.append((lbl, f1, f2, s1, e1, s2, e2))

        if moved:
            print(f"\nShifted features (same label, different coordinates):")
            for lbl, f1, f2, s1, e1, s2, e2 in moved:
                print(f"  ~ {lbl:<20}  [{s1}:{e1}] → [{s2}:{e2}]")

    if not only_in_1 and not only_in_2 and not moved:
        print("\n  Features are identical.")


def cmd_search(args):
    """Search for a DNA sequence motif in both strands (exact match, case-insensitive)."""
    record = read_genbank(args.file)

    query = args.sequence.upper()
    seq = str(record.seq).upper()
    seq_len = len(seq)
    name = record.name or record.id or os.path.basename(args.file)
    context_bp = 10

    # Reverse complement of query
    complement = str.maketrans("ACGT", "TGCA")
    rc_query = query.translate(complement)[::-1]

    # Determine which strands to search
    search_fwd = not args.rc   # search forward unless --rc only
    search_rev = not args.fwd  # search reverse unless --fwd only

    matches = []

    if search_fwd:
        pos = seq.find(query)
        while pos != -1:
            matches.append((pos + 1, "+", pos))  # (1-based pos, strand, 0-based idx)
            pos = seq.find(query, pos + 1)

    if search_rev and rc_query != query:
        pos = seq.find(rc_query)
        while pos != -1:
            matches.append((pos + 1, "-", pos))
            pos = seq.find(rc_query, pos + 1)
    elif search_rev and rc_query == query:
        # Palindrome: already captured all occurrences on fwd; mark them as both strands
        pass

    # Sort by position
    matches.sort(key=lambda x: x[0])

    strand_label = ""
    if args.fwd and not args.rc:
        strand_label = " (forward strand only)"
    elif args.rc and not args.fwd:
        strand_label = " (reverse complement strand only)"

    print(f"Searching for: {query} ({len(query)} bp) in {name}{strand_label}")
    print()

    if not matches:
        print(f"No matches found.")
        return

    print(f"Found {len(matches)} match(es):\n")
    print(f"  {'#':>2}  {'Position':>8}  {'Strand':>6}  Context")

    for i, (pos1, strand, idx0) in enumerate(matches, 1):
        ctx_start = max(0, idx0 - context_bp)
        ctx_end = min(seq_len, idx0 + len(query) + context_bp)
        left = seq[ctx_start:idx0]
        middle = seq[idx0:idx0 + len(query)]
        right = seq[idx0 + len(query):ctx_end]

        if ctx_start > 0:
            left_display = "..." + left
        else:
            left_display = left

        if ctx_end < seq_len:
            right_display = right + "..."
        else:
            right_display = right

        context_str = left_display + middle + right_display
        print(f"  {i:>2}  {pos1:>8}  {strand:>6}  {context_str}")


def cmd_insert(args):
    """Insert a DNA sequence at a given position and write back to file."""
    import re
    from Bio.Seq import Seq

    record = read_genbank(args.file)

    # Validate input sequence: only ACGT allowed (case-insensitive)
    seq_input = args.sequence.upper()
    if not re.fullmatch(r"[ACGT]+", seq_input):
        print(f"Error: sequence contains invalid characters (only ACGT allowed): {args.sequence}", file=sys.stderr)
        sys.exit(1)

    seq_len = len(record.seq)
    position_1based = args.position

    # Validate position: 1-based, must be in range [1, seq_len + 1]
    if position_1based < 1 or position_1based > seq_len + 1:
        print(f"Error: position {position_1based} out of range (valid: 1 to {seq_len + 1})", file=sys.stderr)
        sys.exit(1)

    # Convert to 0-based index for slicing
    insert_idx = position_1based - 1
    insert_len = len(seq_input)

    # Build the new sequence by direct string manipulation
    old_seq_str = str(record.seq)
    new_seq_str = old_seq_str[:insert_idx] + seq_input + old_seq_str[insert_idx:]
    record.seq = Seq(new_seq_str)

    # Shift all features whose start >= insert_idx (downstream features)
    for feat in record.features:
        loc = feat.location
        # Handle simple FeatureLocation
        try:
            from Bio.SeqFeature import CompoundLocation
            if isinstance(loc, CompoundLocation):
                # Shift each part of a compound location
                new_parts = []
                for part in loc.parts:
                    new_start = int(part.start) + insert_len if int(part.start) >= insert_idx else int(part.start)
                    new_end = int(part.end) + insert_len if int(part.end) > insert_idx else int(part.end)
                    new_parts.append(FeatureLocation(new_start, new_end, strand=part.strand))
                feat.location = CompoundLocation(new_parts, operator=loc.operator)
            else:
                new_start = int(loc.start) + insert_len if int(loc.start) >= insert_idx else int(loc.start)
                new_end = int(loc.end) + insert_len if int(loc.end) > insert_idx else int(loc.end)
                feat.location = FeatureLocation(new_start, new_end, strand=loc.strand)
        except Exception:
            # Fallback for simple locations if import or isinstance fails
            new_start = int(loc.start) + insert_len if int(loc.start) >= insert_idx else int(loc.start)
            new_end = int(loc.end) + insert_len if int(loc.end) > insert_idx else int(loc.end)
            feat.location = FeatureLocation(new_start, new_end, strand=loc.strand)

    # Optionally annotate the inserted region
    label = args.label if args.label else None
    feat_type = args.type if args.type else "misc_feature"

    if label:
        new_feature = SeqFeature(
            FeatureLocation(insert_idx, insert_idx + insert_len, strand=1),
            type=feat_type,
            qualifiers={"label": [label]}
        )
        record.features.append(new_feature)

    # Auto-mark the edit with a red CDS annotation
    if not args.no_mark:
        mark_label = f"INSERT: {insert_len}bp"
        seq_display = _truncate_seq(seq_input)
        mark_note = (f"Inserted {insert_len}bp at position {position_1based}. "
                     f"Sequence: {seq_input}")
        _add_edit_mark(record, insert_idx, insert_idx + insert_len, mark_label, mark_note)

    name = record.name or record.id or os.path.basename(args.file)
    new_length = len(record.seq)

    write_genbank(record, args.file)

    if label:
        print(f'Inserted {insert_len} bp at position {position_1based} in {name} (labeled: "{label}", type: {feat_type})')
    else:
        print(f"Inserted {insert_len} bp at position {position_1based} in {name}")
    print(f"New length: {new_length} bp")
    print(f"Backup saved to {args.file}.bak")

    open_in_viewer(args.file)


def cmd_delete(args):
    """Delete a DNA region and shift downstream feature coordinates."""
    from Bio.Seq import Seq

    record = read_genbank(args.file)

    seq_len = len(record.seq)
    start_1based = args.start
    end_1based = args.end

    # Validate coordinates (1-based, inclusive, matching seq/annotate conventions)
    if start_1based < 1 or end_1based > seq_len:
        print(f"Error: coordinates out of range (sequence is {seq_len} bp)", file=sys.stderr)
        sys.exit(1)
    if start_1based >= end_1based:
        print(f"Error: start must be less than end", file=sys.stderr)
        sys.exit(1)

    # Convert to 0-based half-open interval for slicing
    del_start = start_1based - 1   # 0-based inclusive start
    del_end = end_1based            # 0-based exclusive end (== 1-based inclusive end)
    del_len = del_end - del_start

    # Categorize features relative to deleted region
    kept_features = []
    removed_features = []
    warned_features = []

    for feat in record.features:
        loc = feat.location
        try:
            from Bio.SeqFeature import CompoundLocation
            is_compound = isinstance(loc, CompoundLocation)
        except Exception:
            is_compound = False

        if is_compound:
            # For compound locations, check against overall span
            feat_start = int(loc.start)
            feat_end = int(loc.end)
        else:
            feat_start = int(loc.start)
            feat_end = int(loc.end)

        # Entirely within the deleted region — remove
        if feat_start >= del_start and feat_end <= del_end:
            removed_features.append(feat)
            continue

        # Partially overlaps — warn and keep (coordinates may be invalid)
        if feat_start < del_end and feat_end > del_start:
            # Overlaps but not entirely contained
            warned_features.append(feat)
            kept_features.append(feat)
            continue

        # Entirely upstream — unchanged
        if feat_end <= del_start:
            kept_features.append(feat)
            continue

        # Entirely downstream — shift left by del_len
        if is_compound:
            new_parts = []
            for part in loc.parts:
                new_start = int(part.start) - del_len
                new_end = int(part.end) - del_len
                new_parts.append(FeatureLocation(new_start, new_end, strand=part.strand))
            feat.location = CompoundLocation(new_parts, operator=loc.operator)
        else:
            new_start = int(loc.start) - del_len
            new_end = int(loc.end) - del_len
            feat.location = FeatureLocation(new_start, new_end, strand=loc.strand)

        kept_features.append(feat)

    record.features = kept_features

    # Build new sequence
    old_seq_str = str(record.seq)
    new_seq_str = old_seq_str[:del_start] + old_seq_str[del_end:]
    record.seq = Seq(new_seq_str)

    # Auto-mark the deletion junction with a red CDS annotation
    if not args.no_mark:
        # 2bp window around the junction point (del_start in the new sequence)
        mark_start = max(0, del_start - 1)
        mark_end = min(len(record.seq), del_start + 1)
        mark_label = f"DELETED: {del_len}bp"
        mark_note = (f"Deleted {del_len}bp from positions {start_1based}-{end_1based}")
        _add_edit_mark(record, mark_start, mark_end, mark_label, mark_note)

    name = record.name or record.id or os.path.basename(args.file)
    new_length = len(record.seq)

    write_genbank(record, args.file)

    print(f"Deleted {del_len} bp from position {start_1based}-{end_1based} in {name}")

    if removed_features:
        print("Removed features within deleted region:")
        for feat in removed_features:
            label = get_label(feat)
            fs = int(feat.location.start)
            fe = int(feat.location.end)
            print(f"  - {feat.type} \"{label}\" [{fs}:{fe}]")

    if warned_features:
        print("Warning: features partially overlapping the deleted region (kept, may be invalid):")
        for feat in warned_features:
            label = get_label(feat)
            fs = int(feat.location.start)
            fe = int(feat.location.end)
            print(f"  ! {feat.type} \"{label}\" [{fs}:{fe}]")

    print(f"New length: {new_length} bp")
    print(f"Backup saved to {args.file}.bak")

    open_in_viewer(args.file)


def cmd_replace(args):
    """Replace a DNA region with a new sequence, adjusting features accordingly."""
    import re
    from Bio.Seq import Seq

    record = read_genbank(args.file)

    # Validate input sequence: only ACGT allowed (case-insensitive)
    seq_input = args.sequence.upper()
    if not re.fullmatch(r"[ACGT]+", seq_input):
        print(f"Error: sequence contains invalid characters (only ACGT allowed): {args.sequence}", file=sys.stderr)
        sys.exit(1)

    seq_len = len(record.seq)
    start_1based = args.start
    end_1based = args.end

    # Validate coordinates (1-based, inclusive)
    if start_1based < 1 or end_1based > seq_len:
        print(f"Error: coordinates out of range (sequence is {seq_len} bp)", file=sys.stderr)
        sys.exit(1)
    if start_1based >= end_1based:
        print(f"Error: start must be less than end", file=sys.stderr)
        sys.exit(1)

    # Convert to 0-based half-open interval for slicing
    del_start = start_1based - 1   # 0-based inclusive start
    del_end = end_1based            # 0-based exclusive end (== 1-based inclusive end)
    del_len = del_end - del_start
    new_len = len(seq_input)
    net_shift = new_len - del_len

    # Categorize features relative to replaced region
    kept_features = []
    removed_features = []
    warned_features = []

    for feat in record.features:
        loc = feat.location
        try:
            from Bio.SeqFeature import CompoundLocation
            is_compound = isinstance(loc, CompoundLocation)
        except Exception:
            is_compound = False

        feat_start = int(loc.start)
        feat_end = int(loc.end)

        # Entirely within the replaced region — remove
        if feat_start >= del_start and feat_end <= del_end:
            removed_features.append(feat)
            continue

        # Partially overlaps — warn and keep (coordinates may be invalid)
        if feat_start < del_end and feat_end > del_start:
            warned_features.append(feat)
            kept_features.append(feat)
            continue

        # Entirely upstream — unchanged
        if feat_end <= del_start:
            kept_features.append(feat)
            continue

        # Entirely downstream — shift by net_shift
        if is_compound:
            new_parts = []
            for part in loc.parts:
                new_start = int(part.start) + net_shift
                new_end = int(part.end) + net_shift
                new_parts.append(FeatureLocation(new_start, new_end, strand=part.strand))
            feat.location = CompoundLocation(new_parts, operator=loc.operator)
        else:
            new_start = int(loc.start) + net_shift
            new_end = int(loc.end) + net_shift
            feat.location = FeatureLocation(new_start, new_end, strand=loc.strand)

        kept_features.append(feat)

    record.features = kept_features

    # Build new sequence: replace old region with new sequence
    old_seq_str = str(record.seq)
    new_seq_str = old_seq_str[:del_start] + seq_input + old_seq_str[del_end:]
    record.seq = Seq(new_seq_str)

    # Optionally annotate the new region
    label = args.label if args.label else None
    feat_type = args.type if args.type else "misc_feature"

    if label:
        new_feature = SeqFeature(
            FeatureLocation(del_start, del_start + new_len, strand=1),
            type=feat_type,
            qualifiers={"label": [label]}
        )
        record.features.append(new_feature)

    # Auto-mark the edit with a red CDS annotation
    if not args.no_mark:
        old_seq_display = _truncate_seq(old_seq_str[del_start:del_end])
        new_seq_display = _truncate_seq(seq_input)
        mark_label = f"EDIT: {old_seq_display}\u2192{new_seq_display}"
        mark_note = (f"Replaced {del_len}bp at positions {start_1based}-{end_1based}. "
                     f"Old: {old_seq_str[del_start:del_end]}. New: {seq_input}.")
        _add_edit_mark(record, del_start, del_start + new_len, mark_label, mark_note)

    name = record.name or record.id or os.path.basename(args.file)
    new_length = len(record.seq)

    write_genbank(record, args.file)

    print(f"Replaced {del_len} bp ({start_1based}-{end_1based}) with {new_len} bp in {name}")

    if removed_features:
        print("Removed features within replaced region:")
        for feat in removed_features:
            lbl = get_label(feat)
            fs = int(feat.location.start)
            fe = int(feat.location.end)
            print(f"  - {feat.type} \"{lbl}\" [{fs}:{fe}]")

    if warned_features:
        print("Warning: features partially overlapping the replaced region (kept, may be invalid):")
        for feat in warned_features:
            lbl = get_label(feat)
            fs = int(feat.location.start)
            fe = int(feat.location.end)
            print(f"  ! {feat.type} \"{lbl}\" [{fs}:{fe}]")

    net_str = f"{net_shift:+d}"
    print(f"Net change: {net_str} bp")
    print(f"New length: {new_length} bp")
    print(f"Backup saved to {args.file}.bak")

    open_in_viewer(args.file)


def cmd_export(args):
    """Convert a sequence file between formats (fasta, genbank, tab)."""
    import csv

    record = read_genbank(args.file)
    name = record.name or record.id or os.path.splitext(os.path.basename(args.file))[0]
    in_path = os.path.expanduser(args.file)
    fmt = args.format.lower()

    # Normalise alias
    if fmt == "gb":
        fmt = "genbank"

    if fmt not in ("fasta", "genbank", "tab"):
        print(f"Error: unsupported format '{args.format}'. Choose: fasta, genbank, tab", file=sys.stderr)
        sys.exit(1)

    # Determine output path
    if args.out:
        out_path = os.path.expanduser(args.out)
    else:
        base = os.path.splitext(in_path)[0]
        ext_map = {"fasta": ".fasta", "genbank": ".gb", "tab": ".tsv"}
        out_path = base + ext_map[fmt]

    # Write output
    if fmt == "fasta":
        with open(out_path, "w") as fh:
            SeqIO.write(record, fh, "fasta")
        print(f"Exported {name} to FASTA: {out_path}")
        open_in_viewer(out_path)

    elif fmt == "genbank":
        write_genbank(record, out_path)
        print(f"Exported {name} to GenBank: {out_path}")
        open_in_viewer(out_path)

    elif fmt == "tab":
        features = [f for f in record.features if f.type != "source"]
        with open(out_path, "w", newline="") as fh:
            writer = csv.writer(fh, delimiter="\t")
            writer.writerow(["type", "start", "end", "strand", "label"])
            for feat in features:
                writer.writerow([
                    feat.type,
                    int(feat.location.start),
                    int(feat.location.end),
                    strand_str(feat.location.strand),
                    get_label(feat),
                ])
        print(f"Exported {name} features to TSV: {out_path} ({len(features)} features)")


def cmd_open(args):
    """Open a file in the default viewer."""
    path = os.path.expanduser(args.file)
    if not os.path.exists(path):
        print(f"Error: file not found: {path}", file=sys.stderr)
        sys.exit(1)
    subprocess.run(["open", path], check=False)
    print(f"Opened {path}")


def cmd_blast(args):
    """Run a remote NCBI BLAST search and display parsed results."""
    from Bio.Blast import NCBIWWW, NCBIXML

    record = read_genbank(args.file)

    program = args.program
    database = args.db
    evalue = args.evalue
    max_hits = args.max_hits

    name = record.name or record.id or os.path.basename(args.file)
    seq_str = str(record.seq)

    # Handle --region (1-based inclusive)
    if args.region:
        start_1based, end_1based = args.region
        seq_len = len(seq_str)
        if start_1based < 1 or end_1based > seq_len or start_1based >= end_1based:
            print(f"Error: region {start_1based}-{end_1based} out of range (sequence is {seq_len} bp)", file=sys.stderr)
            sys.exit(1)
        seq_str = seq_str[start_1based - 1:end_1based]
        query_label = f"{name}:{start_1based}-{end_1based} ({len(seq_str)} bp)"
    else:
        query_label = f"{name} ({len(seq_str)} bp)"

    print(f"Submitting to NCBI BLAST...")
    print(f"  Program:  {program}")
    print(f"  Database: {database}")
    print(f"  Query:    {query_label}")
    print(f"  E-value:  {evalue}")
    print()
    print("Waiting for results... (this may take 30-60 seconds)")

    try:
        result_handle = NCBIWWW.qblast(
            program,
            database,
            seq_str,
            expect=evalue,
            hitlist_size=max_hits,
        )
    except Exception as e:
        print(f"\nError: BLAST request failed: {e}", file=sys.stderr)
        sys.exit(1)

    try:
        blast_record = NCBIXML.read(result_handle)
    except Exception as e:
        print(f"\nError: Failed to parse BLAST results: {e}", file=sys.stderr)
        sys.exit(1)

    alignments = blast_record.alignments
    n_hits = len(alignments)

    print()
    print(f"Results ({n_hits} hits):")
    print()

    if n_hits == 0:
        print("  No significant hits found.")
        return

    # Header
    print(f"  {'#':>3}  {'Score':>8}  {'E-value':>10}  {'Identity':>9}  {'Accession':<18}  Description")

    for i, alignment in enumerate(alignments, 1):
        if not alignment.hsps:
            continue
        hsp = alignment.hsps[0]  # Best HSP for this alignment

        score = int(hsp.score)
        expect = hsp.expect
        identity_pct = 100.0 * hsp.identities / hsp.align_length if hsp.align_length else 0.0

        accession = alignment.accession or ""
        description = alignment.hit_def or ""
        # Truncate description to 50 chars
        if len(description) > 50:
            description = description[:47] + "..."

        # Format e-value: use scientific notation for very small values
        if expect == 0.0:
            expect_str = "0.0"
        elif expect < 1e-4:
            expect_str = f"{expect:.2e}"
        else:
            expect_str = f"{expect:.4f}"

        print(f"  {i:>3}  {score:>8}  {expect_str:>10}  {identity_pct:>8.1f}%  {accession:<18}  {description}")


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def cmd_exons(args):
    """Find exon boundaries by aligning an mRNA transcript against a genomic sequence."""
    from Bio import Entrez

    Entrez.email = args.email

    # Fetch mRNA from NCBI
    accession = args.transcript
    print(f"Fetching {accession} from NCBI...")
    try:
        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
        mrna_record = SeqIO.read(handle, "fasta")
        handle.close()
    except Exception as e:
        print(f"Error fetching {accession}: {e}", file=sys.stderr)
        sys.exit(1)

    mrna = str(mrna_record.seq)
    print(f"mRNA: {len(mrna)} bp")

    # Read genomic file
    record = read_genbank(args.file)
    genomic = str(record.seq)
    name = record.name or record.id or os.path.basename(args.file)
    print(f"Genomic: {len(genomic)} bp")
    print()

    # Align mRNA against genomic sequence to find exon boundaries
    exons = []
    mrna_pos = 0
    min_match = 20

    while mrna_pos < len(mrna):
        query = mrna[mrna_pos:mrna_pos + 30]
        if len(query) < 20:
            break
        gpos = genomic.upper().find(query.upper())
        if gpos == -1:
            mrna_pos += 1
            continue
        # Extend match
        match_len = 0
        while (mrna_pos + match_len < len(mrna) and
               gpos + match_len < len(genomic) and
               mrna[mrna_pos + match_len].upper() == genomic[gpos + match_len].upper()):
            match_len += 1
        if match_len >= min_match:
            exons.append({
                "exon": len(exons) + 1,
                "mrna_start": mrna_pos + 1,
                "mrna_end": mrna_pos + match_len,
                "genomic_start": gpos + 1,
                "genomic_end": gpos + match_len,
                "length": match_len,
            })
            mrna_pos += match_len
        else:
            mrna_pos += 1

    if not exons:
        print("No exons found. The mRNA may not match this genomic sequence.")
        return

    print(f"Found {len(exons)} exon(s):\n")
    print(f"  {'#':>3}  {'mRNA start':>10}  {'mRNA end':>8}  {'Genomic start':>13}  {'Genomic end':>11}  {'Length':>8}")
    for e in exons:
        print(f"  {e['exon']:>3}  {e['mrna_start']:>10}  {e['mrna_end']:>8}  {e['genomic_start']:>13}  {e['genomic_end']:>11}  {e['length']:>5} bp")

    if not args.annotate:
        print(f"\nUse --annotate to add exon features to the file.")
        return

    # Add exon annotations with distinct colors for SnapGene
    exon_palette = [
        "#4682B4",  # steel blue
        "#2E8B57",  # sea green
        "#DAA520",  # goldenrod
        "#CD5C5C",  # indian red
        "#8A2BE2",  # blue violet
        "#FF8C00",  # dark orange
        "#20B2AA",  # light sea green
        "#DC143C",  # crimson
    ]
    for e in exons:
        color = exon_palette[(e["exon"] - 1) % len(exon_palette)]
        feature = SeqFeature(
            FeatureLocation(e["genomic_start"] - 1, e["genomic_end"], strand=1),
            type="exon",
            qualifiers={
                "label": [f"Exon {e['exon']}"],
                "ApEinfo_fwdcolor": [color],
                "ApEinfo_revcolor": [color],
            }
        )
        record.features.append(feature)

    write_genbank(record, args.file)
    print(f"\nAdded {len(exons)} exon annotations to {os.path.basename(args.file)}")


# ---------------------------------------------------------------------------
# stitch – extract regions by annotation, stitch, optionally translate
# ---------------------------------------------------------------------------

def cmd_stitch(args):
    record = read_genbank(args.file)
    basename = os.path.basename(args.file)

    # Collect matching features
    matched = []
    if args.exons:
        for f in record.features:
            if f.type == "exon":
                label = f.qualifiers.get("label", [""])[0] or f.qualifiers.get("gene", [""])[0] or f"exon@{int(f.location.start)}"
                matched.append((label, f))
    elif args.labels:
        label_set = set(args.labels)
        for f in record.features:
            fl = f.qualifiers.get("label", [""])[0]
            if fl in label_set:
                matched.append((fl, f))
                label_set.discard(fl)
        if label_set:
            print(f"Warning: labels not found: {', '.join(sorted(label_set))}", file=sys.stderr)
    else:
        print("Error: provide annotation labels or use --exons", file=sys.stderr)
        sys.exit(1)

    if not matched:
        print("No matching features found.", file=sys.stderr)
        sys.exit(1)

    # Sort by genomic start position
    matched.sort(key=lambda x: int(x[1].location.start))

    print(f"Stitching {len(matched)} regions from {basename}:\n")
    print(f"  {'#':>3}  {'Label':<14}  {'Position':<16}  {'Length':>8}")

    segments = []
    for i, (label, feat) in enumerate(matched, 1):
        start = int(feat.location.start)
        end = int(feat.location.end)
        seg_seq = str(record.seq[start:end])
        segments.append(seg_seq)
        length = end - start
        print(f"  {i:>3}  {label:<14}  {start + 1}-{end:<12}  {length:>6} bp")

    stitched = "".join(segments)
    total_len = len(stitched)
    print(f"\nStitched: {total_len} bp")

    if args.translate:
        # Find ATG that starts the longest ORF
        upper_stitched = stitched.upper()
        best_atg = -1
        best_prot = ""
        pos = 0
        while True:
            atg_pos = upper_stitched.find("ATG", pos)
            if atg_pos == -1:
                break
            coding = stitched[atg_pos:]
            prot = str(Seq(coding).translate(to_stop=True))
            if len(prot) > len(best_prot):
                best_atg = atg_pos
                best_prot = prot
            pos = atg_pos + 1

        if best_atg == -1:
            print("\nNo ATG found in stitched sequence.", file=sys.stderr)
            sys.exit(1)

        atg_pos = best_atg
        protein = best_prot

        print(f"\nCDS: ATG found at position {atg_pos + 1} of stitched sequence")
        print("Translating...")
        prot_len = len(protein)

        # Calculate MW
        mw_kda = None
        try:
            from Bio.SeqUtils import molecular_weight
            mw = molecular_weight(Seq(protein), seq_type="protein")
            mw_kda = mw / 1000.0
        except Exception:
            mw_kda = prot_len * 110.0 / 1000.0

        print(f"\nProtein: {prot_len} aa")
        print(f"Predicted MW: ~{mw_kda:.1f} kDa")

        # Build a short name from the file
        name_base = os.path.splitext(basename)[0]
        header = f"{name_base}_protein ({prot_len} aa, {mw_kda:.1f} kDa)"
        print(f"\n>{header}")
        # Print protein in 60-char lines
        for j in range(0, prot_len, 60):
            print(protein[j:j+60])
    else:
        # Print stitched DNA in FASTA format
        name_base = os.path.splitext(basename)[0]
        print(f"\n>{name_base}_stitched ({total_len} bp)")
        for j in range(0, total_len, 60):
            print(stitched[j:j+60])


def find_enzyme_sites(seq, enzyme_names):
    """Find restriction sites for named enzymes in a sequence. Returns dict: name → [0-based positions]."""
    from Bio.Restriction import RestrictionBatch, Analysis

    try:
        rb = RestrictionBatch(enzyme_names)
    except Exception as e:
        print(f"Error loading enzymes: {e}", file=sys.stderr)
        sys.exit(1)

    analysis = Analysis(rb, seq, linear=False)
    results = analysis.full()

    out = {}
    for enz, sites in results.items():
        # BioPython Restriction returns 1-based positions; convert to 0-based
        out[str(enz)] = [s - 1 for s in sites]
    return out


def calculate_overlap(backbone_seq, cut_pos, direction, tm_target, preserve_site_seq=None):
    """Calculate a Gibson overlap from backbone sequence around a cut position.

    Args:
        backbone_seq: full backbone sequence (Bio.Seq.Seq or str)
        cut_pos: 0-based cut position
        direction: "upstream" (5' of cut) or "downstream" (3' of cut)
        tm_target: target melting temperature in °C
        preserve_site_seq: if given, ensure this sequence is included in the overlap

    Returns:
        (overlap_seq_str, tm_float)
    """
    from Bio.SeqUtils.MeltingTemp import Tm_NN, DNA_NN3

    seq_str = str(backbone_seq)
    seq_len = len(seq_str)
    min_overlap = 20
    max_overlap = 40

    if direction == "upstream":
        # Extend backward from cut_pos
        for length in range(min_overlap, max_overlap + 1):
            start = cut_pos - length
            if start < 0:
                # Wrap around for circular
                overlap = seq_str[start % seq_len:] + seq_str[:cut_pos]
            else:
                overlap = seq_str[start:cut_pos]
            tm = Tm_NN(Seq(overlap), nn_table=DNA_NN3)
            if tm >= tm_target:
                # If preserving site, check we include it
                if preserve_site_seq and preserve_site_seq.upper() not in overlap.upper():
                    continue
                return overlap, round(tm, 1)
        # Return max length if target not reached
        start = cut_pos - max_overlap
        if start < 0:
            overlap = seq_str[start % seq_len:] + seq_str[:cut_pos]
        else:
            overlap = seq_str[start:cut_pos]
        tm = Tm_NN(Seq(overlap), nn_table=DNA_NN3)
        return overlap, round(tm, 1)
    else:
        # Extend forward from cut_pos
        for length in range(min_overlap, max_overlap + 1):
            end = cut_pos + length
            if end > seq_len:
                overlap = seq_str[cut_pos:] + seq_str[:end % seq_len]
            else:
                overlap = seq_str[cut_pos:end]
            tm = Tm_NN(Seq(overlap), nn_table=DNA_NN3)
            if tm >= tm_target:
                if preserve_site_seq and preserve_site_seq.upper() not in overlap.upper():
                    continue
                return overlap, round(tm, 1)
        # Return max length if target not reached
        end = cut_pos + max_overlap
        if end > seq_len:
            overlap = seq_str[cut_pos:] + seq_str[:end % seq_len]
        else:
            overlap = seq_str[cut_pos:end]
        tm = Tm_NN(Seq(overlap), nn_table=DNA_NN3)
        return overlap, round(tm, 1)


def stitch_features(record, feature_labels):
    """Extract features by label from a record and concatenate their sequences in order.

    Returns the combined sequence as a string.
    """
    segments = []
    found_labels = set()
    for label in feature_labels:
        matched = False
        for feat in record.features:
            if get_label(feat) == label:
                seg = str(feat.extract(record.seq))
                segments.append(seg)
                found_labels.add(label)
                matched = True
                break
        if not matched:
            print(f"Error: feature label '{label}' not found in record", file=sys.stderr)
            sys.exit(1)
    return "".join(segments)


# ---------------------------------------------------------------------------
# Preflight check helpers
# ---------------------------------------------------------------------------

CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

KNOWN_TAGS = {
    'FLAG': 'DYKDDDDK',
    'HA': 'YPYDVPDYA',
    'MYC': 'EQKLISEEDL',
    'V5': 'GKPIPNPLLGLD',
    'His6': 'HHHHHH',
}


def check_reading_frame(seq_str):
    """Analyze reading frame of a DNA sequence. Returns list of (status, category, message)."""
    results = []
    seq_upper = seq_str.upper()

    # Find all ATGs
    atg_positions = []
    pos = 0
    while True:
        pos = seq_upper.find('ATG', pos)
        if pos == -1:
            break
        atg_positions.append(pos)
        pos += 1

    if not atg_positions:
        results.append(("error", "READING FRAME", "No ATG found in sequence"))
        return results

    # Primary CDS from first ATG
    atg_pos = atg_positions[0]
    results.append(("ok", "READING FRAME", f"ATG at position {atg_pos + 1}"))

    cds = seq_upper[atg_pos:]

    # Check divisible by 3
    if len(cds) % 3 != 0:
        results.append(("warn", "READING FRAME", f"CDS length ({len(cds)} bp) not divisible by 3"))
    else:
        results.append(("ok", "READING FRAME", f"CDS length: {len(cds)} bp ({len(cds) // 3} codons)"))

    # Translate and check for stops
    protein = []
    premature_stop = False
    stop_pos = None
    for i in range(0, len(cds) - 2, 3):
        codon = cds[i:i+3]
        if len(codon) < 3:
            break
        aa = CODON_TABLE.get(codon, 'X')
        if aa == '*':
            if i + 3 < len(cds):
                premature_stop = True
                stop_pos = i // 3 + 1
                protein.append(aa)
                break
            else:
                # Terminal stop
                protein.append(aa)
                break
        protein.append(aa)

    protein_str = ''.join(protein)

    if premature_stop:
        results.append(("error", "READING FRAME", f"Premature stop codon at codon {stop_pos}"))
    else:
        results.append(("ok", "READING FRAME", "No premature stop codons"))

    # Check terminal stop
    if protein_str.endswith('*'):
        results.append(("ok", "READING FRAME", "Terminal stop codon present"))
    else:
        results.append(("warn", "READING FRAME", "No terminal stop codon"))

    # Protein length (excluding stop)
    prot_no_stop = protein_str.rstrip('*')
    results.append(("ok", "READING FRAME", f"Protein: {len(prot_no_stop)} aa"))

    # Check for known tags
    for tag_name, tag_seq in KNOWN_TAGS.items():
        if tag_seq in prot_no_stop:
            tag_pos = prot_no_stop.index(tag_seq) + 1
            results.append(("ok", "READING FRAME", f"{tag_name} tag detected at aa {tag_pos}"))

    return results


def check_kozak(seq_str, atg_pos):
    """Check Kozak consensus around an ATG. Returns (status, message)."""
    seq_upper = seq_str.upper()

    # Need at least 3 bases before ATG and 1 after ATG+G
    # Kozak: (gcc)RccATGG where R = A or G at -3, G at +4
    has_minus3 = False
    has_plus4 = False
    context = ""

    # Build context string
    start = max(0, atg_pos - 6)
    end = min(len(seq_upper), atg_pos + 7)
    context = seq_upper[start:end]

    # Check -3 position (purine: A or G)
    if atg_pos >= 3:
        minus3 = seq_upper[atg_pos - 3]
        has_minus3 = minus3 in ('A', 'G')

    # Check +4 position (G after ATG)
    plus4_pos = atg_pos + 3
    if plus4_pos < len(seq_upper):
        has_plus4 = seq_upper[plus4_pos] == 'G'

    if has_minus3 and has_plus4:
        return ("ok", f"Kozak: {context} — strong (-3 purine, +4 G)")
    elif has_minus3 or has_plus4:
        detail = "-3 purine" if has_minus3 else "+4 G"
        missing = "+4 not G" if has_minus3 else "-3 is pyrimidine"
        return ("warn", f"Kozak: {context} — moderate ({detail}, {missing})")
    else:
        return ("warn", f"Kozak: {context} — weak (-3 is pyrimidine, +4 not G)")


def check_codon_usage(protein_seq, dna_seq):
    """Check for rare human codons. Returns (rare_count, total_codons, rare_list)."""
    dna_upper = dna_seq.upper()
    rare_codons = []
    total = 0

    for i in range(0, len(dna_upper) - 2, 3):
        codon = dna_upper[i:i+3]
        if len(codon) < 3:
            break
        aa = CODON_TABLE.get(codon, 'X')
        if aa == '*':
            break
        total += 1
        freq = HUMAN_CODON_FREQ.get(codon, 0)
        if freq < 0.10:
            rare_codons.append((i // 3 + 1, codon, aa, freq))

    return len(rare_codons), total, rare_codons


def check_intron(seq_str):
    """Search for GT...AG intron-like patterns. Returns list of (status, category, message)."""
    results = []
    seq_upper = seq_str.upper()
    min_intron = 60

    # Find all GT positions
    intron_candidates = []
    gt_pos = 0
    while True:
        gt_pos = seq_upper.find('GT', gt_pos)
        if gt_pos == -1:
            break
        # Look for AG at least min_intron bp downstream
        ag_pos = gt_pos + min_intron
        while ag_pos < len(seq_upper) - 1:
            ag_pos = seq_upper.find('AG', ag_pos)
            if ag_pos == -1:
                break
            intron_len = ag_pos + 2 - gt_pos
            # Check extended donor: GT[A/G]AG pattern (GTAAG or GTGAG)
            donor_extended = False
            if gt_pos + 4 < len(seq_upper):
                donor_context = seq_upper[gt_pos:gt_pos + 5]
                if donor_context in ('GTAAG', 'GTGAG'):
                    donor_extended = True

            # Check polypyrimidine tract (8+ C/T in a row within 40bp upstream of AG)
            ppt_found = False
            ppt_region_start = max(gt_pos + 2, ag_pos - 40)
            ppt_region = seq_upper[ppt_region_start:ag_pos]
            run = 0
            max_run = 0
            for base in ppt_region:
                if base in ('C', 'T'):
                    run += 1
                    max_run = max(max_run, run)
                else:
                    run = 0
            if max_run >= 8:
                ppt_found = True

            if donor_extended and ppt_found:
                intron_candidates.append((gt_pos, ag_pos, intron_len, max_run))

            ag_pos += 1
        gt_pos += 1

    if not intron_candidates:
        return results

    # Report the best candidate (longest polypyrimidine tract)
    intron_candidates.sort(key=lambda x: x[3], reverse=True)
    best = intron_candidates[0]
    gt_p, ag_p, i_len, ppt_run = best

    results.append(("ok", "INTRON", f"Candidate intron: {gt_p + 1}..{ag_p + 2} ({i_len} bp)"))

    donor = seq_upper[gt_p:gt_p + 5] if gt_p + 5 <= len(seq_upper) else seq_upper[gt_p:]
    results.append(("ok", "INTRON", f"5' splice donor: {donor}"))

    acceptor = seq_upper[ag_p - 2:ag_p + 2] if ag_p >= 2 else seq_upper[:ag_p + 2]
    results.append(("ok", "INTRON", f"3' splice acceptor: ...{acceptor}"))

    results.append(("ok", "INTRON", f"Polypyrimidine tract: {ppt_run} bp run"))

    if i_len >= 60:
        results.append(("ok", "INTRON", f"Intron length: {i_len} bp (≥60 bp)"))
    else:
        results.append(("warn", "INTRON", f"Intron length: {i_len} bp (<60 bp, unusually short)"))

    # Frame maintenance: check that exon1 + exon2 maintains reading frame
    exon1_len = gt_p  # everything before GT
    exon2_len = len(seq_upper) - (ag_p + 2)  # everything after AG
    intron_bases = i_len
    if intron_bases % 3 == 0:
        results.append(("ok", "INTRON", "Intron removal maintains reading frame"))
    else:
        results.append(("warn", "INTRON", f"Intron removal shifts reading frame by {intron_bases % 3} bases"))

    return results


def check_idt_ordering(seq_str):
    """Check IDT ordering constraints. Returns list of (status, category, message)."""
    results = []
    seq_upper = seq_str.upper()
    seq_len = len(seq_upper)

    # Length check
    if 300 <= seq_len <= 1500:
        results.append(("ok", "IDT ORDERING", f"Length: {seq_len} bp (eBlock range 300-1500)"))
    elif 125 <= seq_len <= 3000:
        results.append(("ok", "IDT ORDERING", f"Length: {seq_len} bp (gBlock range 125-3000)"))
    elif seq_len < 125:
        results.append(("error", "IDT ORDERING", f"Length: {seq_len} bp (too short, min 125 bp for gBlock)"))
    else:
        results.append(("error", "IDT ORDERING", f"Length: {seq_len} bp (exceeds gBlock max 3000 bp)"))

    # Overall GC content
    gc_count = seq_upper.count('G') + seq_upper.count('C')
    gc_pct = gc_count / seq_len * 100 if seq_len > 0 else 0
    if gc_pct < 25:
        results.append(("warn", "IDT ORDERING", f"GC content: {gc_pct:.1f}% (low, <25%)"))
    elif gc_pct > 65:
        results.append(("warn", "IDT ORDERING", f"GC content: {gc_pct:.1f}% (high, >65%)"))
    else:
        results.append(("ok", "IDT ORDERING", f"GC content: {gc_pct:.1f}%"))

    # 50bp sliding window GC
    extreme_windows = 0
    worst_gc = None
    worst_pos = 0
    for i in range(0, seq_len - 49):
        window = seq_upper[i:i + 50]
        w_gc = (window.count('G') + window.count('C')) / 50 * 100
        if w_gc < 20 or w_gc > 80:
            extreme_windows += 1
            if worst_gc is None or abs(w_gc - 50) > abs(worst_gc - 50):
                worst_gc = w_gc
                worst_pos = i + 1
    if extreme_windows > 0:
        results.append(("warn", "IDT ORDERING", f"GC windows: {extreme_windows} extreme 50bp windows (worst: {worst_gc:.0f}% at pos {worst_pos})"))
    else:
        results.append(("ok", "IDT ORDERING", "GC windows: no extreme 50bp windows"))

    # Homopolymer runs
    max_run = 0
    max_run_base = ''
    max_run_pos = 0
    run = 1
    for i in range(1, seq_len):
        if seq_upper[i] == seq_upper[i - 1]:
            run += 1
        else:
            if run > max_run:
                max_run = run
                max_run_base = seq_upper[i - 1]
                max_run_pos = i - run + 1
            run = 1
    if run > max_run:
        max_run = run
        max_run_base = seq_upper[-1]
        max_run_pos = seq_len - run + 1

    if max_run > 8:
        results.append(("warn", "IDT ORDERING", f"Homopolymer: {max_run}x {max_run_base} at pos {max_run_pos}"))
    else:
        results.append(("ok", "IDT ORDERING", f"Homopolymers: longest run {max_run} bp"))

    # Direct repeats >20bp
    direct_repeats = 0
    seen_kmers = {}
    for i in range(seq_len - 19):
        kmer = seq_upper[i:i + 20]
        if kmer in seen_kmers:
            direct_repeats += 1
        else:
            seen_kmers[kmer] = i
    if direct_repeats > 0:
        results.append(("warn", "IDT ORDERING", f"Direct repeats: {direct_repeats} repeated 20bp sequences"))
    else:
        results.append(("ok", "IDT ORDERING", "Direct repeats: none >20bp"))

    # Inverted repeats >15bp
    inverted_repeats = 0
    complement = str.maketrans('ACGT', 'TGCA')
    seen_15mers = set()
    for i in range(seq_len - 14):
        kmer = seq_upper[i:i + 15]
        seen_15mers.add(kmer)

    for i in range(seq_len - 14):
        kmer = seq_upper[i:i + 15]
        rc = kmer[::-1].translate(complement)
        if rc in seen_15mers and rc != kmer:
            inverted_repeats += 1

    # Each pair is counted twice, so divide
    inverted_repeats = inverted_repeats // 2
    if inverted_repeats > 0:
        results.append(("warn", "IDT ORDERING", f"Inverted repeats: {inverted_repeats} pairs of 15bp inverted repeats"))
    else:
        results.append(("ok", "IDT ORDERING", "Inverted repeats: none >15bp"))

    return results


def try_idt_api(seq_str):
    """Try IDT complexity screening API. Returns list of issues or None."""
    import json

    api_key = os.environ.get("IDT_API_KEY")
    if not api_key:
        secrets_path = os.path.expanduser("~/.config/splicemap/secrets.env")
        if os.path.exists(secrets_path):
            with open(secrets_path) as f:
                for line in f:
                    line = line.strip()
                    if line.startswith('IDT_API_KEY='):
                        api_key = line.split('=', 1)[1].strip().strip('"').strip("'")
                        break

    if not api_key:
        return None

    try:
        import urllib.request
        import urllib.error

        url = "https://www.idtdna.com/api/complexities/screengBlockSequences"
        data = json.dumps([{"Name": "check", "Sequence": seq_str}]).encode('utf-8')
        req = urllib.request.Request(url, data=data, method='POST')
        req.add_header('Content-Type', 'application/json')
        req.add_header('Authorization', f'Bearer {api_key}')

        with urllib.request.urlopen(req, timeout=15) as resp:
            result = json.loads(resp.read().decode('utf-8'))

        issues = []
        if isinstance(result, list):
            for item in result:
                if isinstance(item, list):
                    for issue in item:
                        if isinstance(issue, dict):
                            name = issue.get('Name', 'Unknown')
                            desc = issue.get('Description', '')
                            issues.append(f"{name}: {desc}")
                elif isinstance(item, dict):
                    name = item.get('Name', 'Unknown')
                    desc = item.get('Description', '')
                    issues.append(f"{name}: {desc}")

        return issues if issues else []
    except Exception:
        return None


def cmd_check(args):
    """Preflight validation of construct files."""
    record = read_genbank(args.file)
    seq_str = str(record.seq).upper()
    seq_len = len(seq_str)
    basename = os.path.basename(args.file)

    # Collect all results: list of (category, status, message)
    all_results = []

    # ── Gibson Assembly checks ──────────────────────────────────────────
    if args.backbone and args.enzymes:
        backbone_record = read_genbank(args.backbone)
        backbone_seq = backbone_record.seq
        enzyme_names = [e.strip() for e in args.enzymes.split(",")]

        # Check enzyme sites in backbone
        sites = find_enzyme_sites(backbone_seq, enzyme_names)
        for enz_name in enzyme_names:
            enz_sites = sites.get(enz_name, [])
            if len(enz_sites) == 1:
                all_results.append(("GIBSON ASSEMBLY", "ok", f"{enz_name}: 1 site in backbone (single cutter)"))
            elif len(enz_sites) == 0:
                all_results.append(("GIBSON ASSEMBLY", "error", f"{enz_name}: no sites in backbone"))
            else:
                all_results.append(("GIBSON ASSEMBLY", "error", f"{enz_name}: {len(enz_sites)} sites in backbone (not single cutter)"))

        # Check for enzyme sites in the construct/insert
        insert_sites = find_enzyme_sites(Seq(seq_str), enzyme_names)
        for enz_name in enzyme_names:
            i_sites = insert_sites.get(enz_name, [])
            if i_sites:
                all_results.append(("GIBSON ASSEMBLY", "warn", f"{enz_name}: {len(i_sites)} site(s) in construct"))
            else:
                all_results.append(("GIBSON ASSEMBLY", "ok", f"{enz_name}: no sites in construct"))

        # Look for overlap annotations
        from Bio.SeqUtils.MeltingTemp import Tm_NN, DNA_NN3
        for feat in record.features:
            label = get_label(feat)
            if label and ("overlap" in label.lower()):
                overlap_seq = str(feat.extract(record.seq))
                try:
                    tm = Tm_NN(Seq(overlap_seq), nn_table=DNA_NN3)
                    tm = round(tm, 1)
                    if 55 <= tm <= 65:
                        all_results.append(("GIBSON ASSEMBLY", "ok", f"{label}: {len(overlap_seq)} bp, Tm={tm}°C"))
                    elif tm < 55:
                        all_results.append(("GIBSON ASSEMBLY", "warn", f"{label}: Tm={tm}°C (low, target 55-65°C)"))
                    else:
                        all_results.append(("GIBSON ASSEMBLY", "warn", f"{label}: Tm={tm}°C (high, target 55-65°C)"))
                except Exception:
                    all_results.append(("GIBSON ASSEMBLY", "warn", f"{label}: could not calculate Tm"))

    # ── Reading frame checks ────────────────────────────────────────────
    rf_results = check_reading_frame(seq_str)
    all_results.extend([(cat, st, msg) for st, cat, msg in rf_results])

    # ── Expression checks ───────────────────────────────────────────────
    # Find primary ATG
    atg_pos = seq_str.find('ATG')
    if atg_pos >= 0:
        # Kozak
        kozak_status, kozak_msg = check_kozak(seq_str, atg_pos)
        all_results.append(("EXPRESSION", kozak_status, kozak_msg))

        # Upstream ATGs (in 5' overlap/UTR region before primary CDS)
        if atg_pos == 0:
            all_results.append(("EXPRESSION", "ok", "No upstream ATGs (CDS starts at position 1)"))
        else:
            all_results.append(("EXPRESSION", "ok", f"No upstream ATGs (primary ATG is first in sequence)"))

        # Codon usage
        cds = seq_str[atg_pos:]
        # Translate to get protein
        protein = []
        for i in range(0, len(cds) - 2, 3):
            codon = cds[i:i+3]
            aa = CODON_TABLE.get(codon, 'X')
            if aa == '*':
                break
            protein.append(aa)
        protein_str = ''.join(protein)
        rare_count, total_codons, rare_list = check_codon_usage(protein_str, cds)

        if rare_count == 0:
            all_results.append(("EXPRESSION", "ok", f"Codon usage: no rare codons ({total_codons} codons checked)"))
        elif rare_count <= 3:
            all_results.append(("EXPRESSION", "warn", f"Codon usage: {rare_count} rare codon(s) of {total_codons}"))
            for pos, codon, aa, freq in rare_list:
                all_results.append(("EXPRESSION", "warn", f"  codon {pos}: {codon} ({aa}) freq={freq:.2f}"))
        else:
            all_results.append(("EXPRESSION", "warn", f"Codon usage: {rare_count} rare codons of {total_codons}"))
            for pos, codon, aa, freq in rare_list[:5]:
                all_results.append(("EXPRESSION", "warn", f"  codon {pos}: {codon} ({aa}) freq={freq:.2f}"))
            if len(rare_list) > 5:
                all_results.append(("EXPRESSION", "warn", f"  ... and {len(rare_list) - 5} more"))

    # ── Intron checks ───────────────────────────────────────────────────
    intron_results = check_intron(seq_str)
    if intron_results:
        all_results.extend([(cat, st, msg) for st, cat, msg in intron_results])
    # We'll note "no intron" in the report formatting below

    # ── IDT ordering checks ─────────────────────────────────────────────
    idt_results = check_idt_ordering(seq_str)
    all_results.extend([(cat, st, msg) for st, cat, msg in idt_results])

    # IDT API check
    api_issues = try_idt_api(seq_str)
    if api_issues is not None:
        if not api_issues:
            all_results.append(("IDT ORDERING", "ok", "IDT API: no complexity issues"))
        else:
            for issue in api_issues:
                all_results.append(("IDT ORDERING", "warn", f"IDT API: {issue}"))
    else:
        all_results.append(("IDT ORDERING", "warn", "IDT API: no API key configured (heuristics only)"))

    # ── Print formatted report ──────────────────────────────────────────
    status_sym = {"ok": "\u2713", "warn": "\u26A0", "error": "\u2717"}

    print(f"PREFLIGHT: {basename} ({seq_len} bp)")
    print("\u2550" * 40)

    # Group by category, preserving order
    categories_seen = []
    cat_results = {}
    for cat, st, msg in all_results:
        if cat not in cat_results:
            cat_results[cat] = []
            categories_seen.append(cat)
        cat_results[cat].append((st, msg))

    # Ensure we always show these categories in order
    standard_order = ["GIBSON ASSEMBLY", "READING FRAME", "EXPRESSION", "INTRON", "IDT ORDERING"]
    ordered_cats = []
    for cat in standard_order:
        if cat in cat_results:
            ordered_cats.append(cat)
    # Add any extras
    for cat in categories_seen:
        if cat not in ordered_cats:
            ordered_cats.append(cat)

    warn_count = 0
    error_count = 0

    for cat in standard_order:
        print(f"\n{cat}")
        if cat not in cat_results:
            if cat == "GIBSON ASSEMBLY":
                if not (args.backbone and args.enzymes):
                    print("  (no backbone/enzymes provided, skipped)")
            elif cat == "INTRON":
                print("  (no intron detected)")
            continue
        for st, msg in cat_results[cat]:
            sym = status_sym.get(st, "?")
            print(f"  {sym} {msg}")
            if st == "warn":
                warn_count += 1
            elif st == "error":
                error_count += 1

    # Also print any extra categories not in standard_order
    for cat in ordered_cats:
        if cat not in standard_order:
            print(f"\n{cat}")
            for st, msg in cat_results[cat]:
                sym = status_sym.get(st, "?")
                print(f"  {sym} {msg}")
                if st == "warn":
                    warn_count += 1
                elif st == "error":
                    error_count += 1

    # Summary
    if error_count > 0:
        verdict = "NEEDS ATTENTION"
    else:
        verdict = "READY TO ORDER"

    print(f"\nSUMMARY: {warn_count} warnings, {error_count} errors \u2014 {verdict}")


def cmd_gibson(args):
    """Design Gibson Assembly eBlocks from a backbone, restriction enzymes, and insert."""
    import Bio.Restriction

    record = read_genbank(args.file)
    backbone_seq = record.seq
    backbone_name = record.name or record.id
    warnings = []

    # Parse enzyme names
    enzyme_names = [e.strip() for e in args.enzymes.split(",")]
    if len(enzyme_names) != 2:
        print("Error: exactly 2 enzymes required (comma-separated)", file=sys.stderr)
        sys.exit(1)

    # Get insert sequence
    insert_seq = None
    insert_features = []  # list of (label, type, seq_len) for annotation

    if args.insert:
        insert_seq = args.insert.upper().strip()
        if not all(c in "ACGT" for c in insert_seq):
            print("Error: insert sequence must contain only ACGT characters", file=sys.stderr)
            sys.exit(1)
    elif args.insert_from:
        insert_record = read_genbank(args.insert_from)
        if args.features:
            feature_labels = [f.strip() for f in args.features.split(",")]
            # Collect feature info for annotation before stitching
            offset = 0
            for label in feature_labels:
                for feat in insert_record.features:
                    if get_label(feat) == label:
                        seg_len = len(feat.extract(insert_record.seq))
                        insert_features.append((label, feat.type, offset, seg_len))
                        offset += seg_len
                        break
            insert_seq = stitch_features(insert_record, feature_labels)
        else:
            # Use the full sequence from the insert file
            insert_seq = str(insert_record.seq).upper()
    else:
        print("Error: provide --insert or --insert-from", file=sys.stderr)
        sys.exit(1)

    # Find enzyme sites in backbone
    sites = find_enzyme_sites(backbone_seq, enzyme_names)

    for enz_name in enzyme_names:
        enz_sites = sites.get(enz_name, [])
        if len(enz_sites) == 0:
            print(f"Error: {enz_name} does not cut the backbone", file=sys.stderr)
            sys.exit(1)
        elif len(enz_sites) > 1:
            print(f"Error: {enz_name} cuts {len(enz_sites)} times (must be single cutter)", file=sys.stderr)
            sys.exit(1)

    cut1 = sites[enzyme_names[0]][0]
    cut2 = sites[enzyme_names[1]][0]

    # Sort so cut_upstream < cut_downstream
    cut_upstream = min(cut1, cut2)
    cut_downstream = max(cut1, cut2)

    # Check insert for enzyme sites
    insert_sites = find_enzyme_sites(Seq(insert_seq), enzyme_names)
    for enz_name in enzyme_names:
        if insert_sites.get(enz_name, []):
            msg = f"Warning: insert contains {enz_name} site(s)"
            warnings.append(msg)
            print(msg, file=sys.stderr)

    # Get enzyme recognition sequences for --preserve-sites
    preserve_5p = None
    preserve_3p = None
    if args.preserve_sites:
        enz1_name = enzyme_names[0] if sites[enzyme_names[0]][0] == cut_upstream else enzyme_names[1]
        enz2_name = enzyme_names[1] if enz1_name == enzyme_names[0] else enzyme_names[0]
        enz1_obj = getattr(Bio.Restriction, enz1_name)
        enz2_obj = getattr(Bio.Restriction, enz2_name)
        preserve_5p = str(enz1_obj.site)
        preserve_3p = str(enz2_obj.site)

    # Calculate overlaps
    overlap_5p, tm_5p = calculate_overlap(backbone_seq, cut_upstream, "upstream", args.tm_target, preserve_5p)
    overlap_3p, tm_3p = calculate_overlap(backbone_seq, cut_downstream, "downstream", args.tm_target, preserve_3p)

    # Build eBlock sequence
    eblock_seq_str = overlap_5p + insert_seq + overlap_3p
    eblock_record = SeqRecord(
        Seq(eblock_seq_str),
        id=args.label.replace(" ", "_")[:16],
        name=args.label.replace(" ", "_")[:16],
        description=f"Gibson eBlock: {args.label}",
    )
    eblock_record.annotations["molecule_type"] = "DNA"
    eblock_record.annotations["topology"] = "linear"

    # Annotate eBlock
    pos = 0
    # 5' overlap
    eblock_record.features.append(SeqFeature(
        FeatureLocation(pos, pos + len(overlap_5p)),
        type="misc_feature",
        qualifiers={"label": [f"5' overlap (Tm={tm_5p}\u00b0C)"]}
    ))
    pos += len(overlap_5p)

    # Insert features
    if insert_features:
        for label, feat_type, feat_offset, feat_len in insert_features:
            eblock_record.features.append(SeqFeature(
                FeatureLocation(pos + feat_offset, pos + feat_offset + feat_len),
                type=feat_type,
                qualifiers={"label": [label]}
            ))
    else:
        eblock_record.features.append(SeqFeature(
            FeatureLocation(pos, pos + len(insert_seq)),
            type="CDS",
            qualifiers={"label": [args.label]}
        ))
    pos += len(insert_seq)

    # 3' overlap
    eblock_record.features.append(SeqFeature(
        FeatureLocation(pos, pos + len(overlap_3p)),
        type="misc_feature",
        qualifiers={"label": [f"3' overlap (Tm={tm_3p}\u00b0C)"]}
    ))

    # Build assembled plasmid
    bb_str = str(backbone_seq)
    assembled_seq_str = bb_str[:cut_upstream] + insert_seq + bb_str[cut_downstream:]
    assembled_record = SeqRecord(
        Seq(assembled_seq_str),
        id=record.id,
        name=record.name,
        description=f"Assembled: {args.label} in {backbone_name}",
    )
    assembled_record.annotations["molecule_type"] = "DNA"
    assembled_record.annotations["topology"] = record.annotations.get("topology", "circular")

    # Adjust backbone features for assembled plasmid
    delta = len(insert_seq) - (cut_downstream - cut_upstream)
    for feat in record.features:
        start = int(feat.location.start)
        end = int(feat.location.end)
        # Skip features entirely within the replaced region
        if start >= cut_upstream and end <= cut_downstream:
            continue
        # Shift features downstream of the cut
        if start >= cut_downstream:
            new_start = start + delta
            new_end = end + delta
        elif end <= cut_upstream:
            new_start = start
            new_end = end
        else:
            # Partial overlap — skip with warning
            label = get_label(feat)
            msg = f"Warning: feature '{label}' overlaps cut region, skipped in assembly"
            warnings.append(msg)
            continue
        assembled_record.features.append(SeqFeature(
            FeatureLocation(new_start, new_end, strand=feat.location.strand),
            type=feat.type,
            qualifiers=dict(feat.qualifiers)
        ))

    # Add insert annotation to assembled plasmid
    assembled_record.features.append(SeqFeature(
        FeatureLocation(cut_upstream, cut_upstream + len(insert_seq)),
        type="CDS",
        qualifiers={"label": [args.label]}
    ))

    # Output directory
    if args.out_dir:
        out_dir = os.path.expanduser(args.out_dir)
    else:
        out_dir = os.path.dirname(os.path.expanduser(args.file)) or "."
    os.makedirs(out_dir, exist_ok=True)

    label_safe = args.label.replace(" ", "_")
    eblock_path = os.path.join(out_dir, f"{label_safe}_eBlock.gb")
    assembled_path = os.path.join(out_dir, f"{label_safe}_assembled.gb")

    write_genbank(eblock_record, eblock_path)
    write_genbank(assembled_record, assembled_path)

    open_in_viewer(eblock_path)
    open_in_viewer(assembled_path)

    # Print summary
    print(f"Gibson Assembly Design: {args.label}")
    print(f"  Backbone:       {backbone_name} ({len(backbone_seq)} bp)")
    print(f"  Enzymes:        {enzyme_names[0]} (pos {cut_upstream + 1}), {enzyme_names[1]} (pos {cut_downstream + 1})")
    print(f"  Insert:         {len(insert_seq)} bp")
    print(f"  5' overlap:     {len(overlap_5p)} bp (Tm={tm_5p}\u00b0C)")
    print(f"  3' overlap:     {len(overlap_3p)} bp (Tm={tm_3p}\u00b0C)")
    print(f"  eBlock:         {len(eblock_seq_str)} bp → {eblock_path}")
    print(f"  Assembled:      {len(assembled_seq_str)} bp → {assembled_path}")

    if warnings:
        print(f"\nWarnings ({len(warnings)}):")
        for w in warnings:
            print(f"  {w}")

    print("\nNote: run 'bio.py sites' on the assembled plasmid to verify restriction sites.")


# ---------------------------------------------------------------------------
# Branch Point Prediction
# ---------------------------------------------------------------------------

TOOLS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "tools")
BPP_DIR = os.path.join(TOOLS_DIR, "BPP")
SVM_DIR = os.path.join(TOOLS_DIR, "SVM-BPfinder")


def _ensure_branchpoint_tools():
    """Ensure BPP and SVM-BPfinder are installed. Clone from GitHub if missing."""
    bpp_script = os.path.join(BPP_DIR, "BP_PPT.py")
    svm_script = os.path.join(SVM_DIR, "svm_bpfinder.py")

    if not os.path.exists(bpp_script):
        print("Installing BPP...", file=sys.stderr)
        os.makedirs(BPP_DIR, exist_ok=True)
        tmp_dir = os.path.join(TOOLS_DIR, "_bpp_clone")
        subprocess.run(["git", "clone", "https://github.com/zhqingit/BPP.git", tmp_dir],
                        capture_output=True, check=True)
        # BPP is Python 2 — we need to convert BP_PPT.py to Python 3
        src = os.path.join(tmp_dir, "BP_PPT.py")
        with open(src) as f:
            code = f.read()
        # Apply Python 2 → 3 fixes
        code = code.replace('#!/usr/bin/python', '#!/usr/bin/python3')
        import re
        # Fix print statements to print() calls
        code = re.sub(r'^(\s*)print\s+"([^"]*)"', r'\1print("\2")', code, flags=re.MULTILINE)
        code = re.sub(r'^(\s*)print\s+(\w)', r'\1print(\2', code, flags=re.MULTILINE)
        # Fix unclosed print( from the regex above - add ) at end of print lines
        lines = code.split('\n')
        fixed = []
        for line in lines:
            stripped = line.lstrip()
            if stripped.startswith('print(') and not stripped.endswith(')') and not stripped.endswith('\\'):
                # Count parens
                opens = stripped.count('(')
                closes = stripped.count(')')
                if opens > closes:
                    line = line + ')' * (opens - closes)
            fixed.append(line)
        code = '\n'.join(fixed)
        with open(bpp_script, 'w') as f:
            f.write(code)
        # Copy supporting files
        import shutil as _shutil
        demo_src = os.path.join(tmp_dir, "demo")
        demo_dst = os.path.join(BPP_DIR, "demo")
        if os.path.isdir(demo_src):
            _shutil.copytree(demo_src, demo_dst, dirs_exist_ok=True)
        lic = os.path.join(tmp_dir, "LICENSE")
        if os.path.exists(lic):
            _shutil.copy2(lic, BPP_DIR)
        _shutil.rmtree(tmp_dir, ignore_errors=True)
        print("BPP installed.", file=sys.stderr)

    if not os.path.exists(svm_script):
        print("Installing SVM-BPfinder...", file=sys.stderr)
        os.makedirs(SVM_DIR, exist_ok=True)
        tmp_dir = os.path.join(TOOLS_DIR, "_svm_clone")
        subprocess.run(["git", "clone", "https://github.com/comprna/SVM-BPfinder-3M.git", tmp_dir],
                        capture_output=True, check=True)
        import shutil as _shutil
        for item in ["svm_bpfinder.py", "SCRIPTS", "MODELS", "LICENSE"]:
            src = os.path.join(tmp_dir, item)
            dst = os.path.join(SVM_DIR, item)
            if os.path.isdir(src):
                _shutil.copytree(src, dst, dirs_exist_ok=True)
            elif os.path.isfile(src):
                _shutil.copy2(src, dst)
        # Ensure executables
        for p in [os.path.join(SVM_DIR, "svm_bpfinder.py"),
                  os.path.join(SVM_DIR, "SCRIPTS", "svm_classify"),
                  os.path.join(SVM_DIR, "SCRIPTS", "svm_getfeat.py")]:
            if os.path.exists(p):
                os.chmod(p, 0o755)
        _shutil.rmtree(tmp_dir, ignore_errors=True)
        print("SVM-BPfinder installed.", file=sys.stderr)

    return bpp_script, svm_script


def _run_bpp(intron_seq, n_results=10):
    """Run BPP on an intron sequence. Returns list of dicts with predictions."""
    import tempfile
    bpp_script, _ = _ensure_branchpoint_tools()
    pwm = os.path.join(BPP_DIR, "demo", "pwmBP_human.txt")
    ppt = os.path.join(BPP_DIR, "demo", "scPPT_human.txt")

    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
        f.write(">query\n")
        f.write(intron_seq.upper() + "\n")
        fa_path = f.name

    try:
        result = subprocess.run(
            ["python3", bpp_script, "-b", pwm, "-p", ppt, "-i", fa_path, "-r", str(n_results)],
            capture_output=True, text=True, timeout=60
        )
        lines = result.stdout.strip().split('\n')
        predictions = []
        for line in lines:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split('\t')
            if len(parts) >= 9:
                predictions.append({
                    'motif_7mer': parts[1],
                    'dist_to_3ss': int(parts[2]),
                    'sc_bps': float(parts[3]),
                    'sc_ppt': float(parts[4]),
                    'sc': float(parts[5]),
                    'zsc_bps': float(parts[6]),
                    'zsc_ppt': float(parts[7]),
                    'zsc': float(parts[8]),
                })
        return predictions
    except Exception as e:
        print(f"Warning: BPP failed: {e}", file=sys.stderr)
        return []
    finally:
        os.unlink(fa_path)


def _run_svm_bpfinder(intron_seq, scan_length=None):
    """Run SVM-BPfinder on an intron sequence. Returns list of dicts with predictions."""
    import tempfile
    _, svm_script = _ensure_branchpoint_tools()

    if scan_length is None:
        scan_length = min(len(intron_seq), 500)

    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False, dir=SVM_DIR) as f:
        f.write(">query\n")
        f.write(intron_seq.upper() + "\n")
        fa_path = f.name

    old_cwd = os.getcwd()
    try:
        os.chdir(SVM_DIR)
        result = subprocess.run(
            ["python3", svm_script, "-i", fa_path, "-s", "Hsap",
             "-l", str(scan_length), "-d", "10"],
            capture_output=True, text=True, timeout=60
        )
        lines = result.stdout.strip().split('\n')
        predictions = []
        for line in lines:
            if line.startswith('seq_id') or not line.strip():
                continue
            parts = line.split('\t')
            if len(parts) >= 10:
                predictions.append({
                    'agez': int(parts[1]),
                    'dist_to_3ss': int(parts[2]),
                    'motif_9mer': parts[3],
                    'bp_scr': float(parts[4]),
                    'y_cont': float(parts[5]),
                    'ppt_off': int(parts[6]),
                    'ppt_len': int(parts[7]),
                    'ppt_scr': float(parts[8]) if parts[8] != '' else 0.0,
                    'svm_scr': float(parts[9]),
                })
        return predictions
    except Exception as e:
        print(f"Warning: SVM-BPfinder failed: {e}", file=sys.stderr)
        return []
    finally:
        os.chdir(old_cwd)
        os.unlink(fa_path)


def _merge_predictions(bpp_preds, svm_preds, intron_len):
    """Merge BPP and SVM-BPfinder predictions. Match within 2nt tolerance."""
    # Build unified list. Key = distance to 3'SS (position from end).
    merged = {}

    for bp in bpp_preds:
        dist = bp['dist_to_3ss']
        # intron position (0-based from start) = intron_len - dist - 4
        # BPP's bp_pos is distance from the 3' end, specifically from the last AG
        # The branch point A is at position 5 of the 7-mer (index 4, 0-based)
        # Intron pos (1-based) of BP A = intron_len - dist + 1
        # Actually from the output: bp_pos=20 means 20 nt from 3'SS
        intron_pos_1based = intron_len - dist
        merged[dist] = {
            'dist_to_3ss': dist,
            'intron_pos': intron_pos_1based,
            'bpp_zscore': bp['zsc'],
            'bpp_motif': bp['motif_7mer'],
            'svm_score': None,
            'svm_motif': None,
        }

    for sv in svm_preds:
        dist = sv['dist_to_3ss']
        # Find if there's a matching BPP prediction within 2nt
        matched_key = None
        for existing_dist in list(merged.keys()):
            if abs(existing_dist - dist) <= 2:
                matched_key = existing_dist
                break

        if matched_key is not None:
            merged[matched_key]['svm_score'] = sv['svm_scr']
            merged[matched_key]['svm_motif'] = sv['motif_9mer']
        else:
            intron_pos_1based = intron_len - dist
            merged[dist] = {
                'dist_to_3ss': dist,
                'intron_pos': intron_pos_1based,
                'bpp_zscore': None,
                'bpp_motif': None,
                'svm_score': sv['svm_scr'],
                'svm_motif': sv['motif_9mer'],
            }

    # Score consensus
    results = list(merged.values())
    for r in results:
        has_bpp = r['bpp_zscore'] is not None
        has_svm = r['svm_score'] is not None
        bpp_high = has_bpp and r['bpp_zscore'] > 1.0
        svm_high = has_svm and r['svm_score'] > 0.0

        if has_bpp and has_svm and bpp_high and svm_high:
            r['consensus'] = 3  # ★★★
        elif has_bpp and has_svm:
            r['consensus'] = 2  # ★★☆
        else:
            r['consensus'] = 1  # ★☆☆

    # Sort by consensus desc, then by BPP z-score desc (if available), then SVM score desc
    def sort_key(r):
        bpp = r['bpp_zscore'] if r['bpp_zscore'] is not None else -999
        svm = r['svm_score'] if r['svm_score'] is not None else -999
        return (-r['consensus'], -bpp, -svm)

    results.sort(key=sort_key)
    return results


def _consensus_stars(n):
    """Return star string for consensus level."""
    if n == 3:
        return "\u2605\u2605\u2605"
    elif n == 2:
        return "\u2605\u2605\u2606"
    else:
        return "\u2605\u2606\u2606"


def _extract_intron_sequence(args):
    """Extract intron sequence from args. Returns (seq, name, intron_label, genomic_start, genomic_end, file_dir)."""
    if hasattr(args, 'seq') and args.seq:
        seq = args.seq.upper().replace('\n', '').replace(' ', '')
        return seq, "query", "raw sequence", None, None, os.getcwd()

    record = read_genbank(args.file)
    file_dir = os.path.dirname(os.path.abspath(os.path.expanduser(args.file)))
    name = record.name or record.id or os.path.splitext(os.path.basename(args.file))[0]

    if hasattr(args, 'intron') and args.intron:
        # Find feature by label
        for feat in record.features:
            label = get_label(feat)
            if label and args.intron.lower() in label.lower():
                start = int(feat.location.start)  # 0-based
                end = int(feat.location.end)       # 0-based exclusive
                seq = str(record.seq[start:end]).upper()
                return seq, name, label, start + 1, end, file_dir
        print(f"Error: no feature matching '{args.intron}' found", file=sys.stderr)
        sys.exit(1)

    elif hasattr(args, 'region') and args.region:
        start_1 = args.region[0]
        end_1 = args.region[1]
        start_0 = start_1 - 1
        seq = str(record.seq[start_0:end_1]).upper()
        return seq, name, f"region {start_1}-{end_1}", start_1, end_1, file_dir

    else:
        print("Error: specify --intron, --region, or --seq", file=sys.stderr)
        sys.exit(1)


def cmd_branchpoint(args):
    """Predict branch point locations in intron sequences using BPP and SVM-BPfinder."""
    intron_seq, name, intron_label, genomic_start, genomic_end, file_dir = _extract_intron_sequence(args)
    intron_len = len(intron_seq)

    # Validate: intron should end with AG (3' splice site)
    three_ss = intron_seq[-2:].upper()
    if three_ss != "AG":
        print(f"Warning: intron does not end with AG (found: {three_ss}). 3'SS detection may be inaccurate.", file=sys.stderr)

    # Run both tools
    bpp_preds = _run_bpp(intron_seq, n_results=10)
    svm_preds = _run_svm_bpfinder(intron_seq, scan_length=min(intron_len, 500))

    bpp_ok = len(bpp_preds) > 0
    svm_ok = len(svm_preds) > 0

    if not bpp_ok and not svm_ok:
        print("Error: both BPP and SVM-BPfinder failed to produce results.", file=sys.stderr)
        sys.exit(1)

    # Merge and rank
    results = _merge_predictions(bpp_preds, svm_preds, intron_len)

    # Limit to top 10
    results = results[:10]

    # ── Terminal output ──────────────────────────────────────────────────
    coord_str = f" (positions {genomic_start}-{genomic_end})" if genomic_start else ""
    print(f"\nBranch Point Prediction: {name} / {intron_label} ({intron_len} bp){coord_str}")
    if not bpp_ok:
        print("  Note: BPP failed — showing SVM-BPfinder results only")
    if not svm_ok:
        print("  Note: SVM-BPfinder failed — showing BPP results only")
    print()

    # Table header
    print(f"  {'#':>2}  {'Position':>8}  {'Score(BPP)':>10}  {'Score(SVM)':>10}  {'Consensus':>9}  {'Dist to 3\'SS':>12}  Motif")

    for i, r in enumerate(results):
        pos_str = str(r['intron_pos'])
        bpp_str = f"{r['bpp_zscore']:.2f}" if r['bpp_zscore'] is not None else "---"
        svm_str = f"{r['svm_score']:.3f}" if r['svm_score'] is not None else "---"
        stars = _consensus_stars(r['consensus'])
        dist_str = f"{r['dist_to_3ss']} nt"
        motif = r['svm_motif'] or r['bpp_motif'] or "---"
        # Capitalize the branch point A in the motif (position 5 of 9-mer or 4 of 7-mer)
        if r['svm_motif'] and len(r['svm_motif']) == 9:
            m = list(r['svm_motif'].lower())
            m[4] = m[4].upper()  # Branch point A is at center of 9-mer
            motif = ''.join(m)
        elif r['bpp_motif'] and len(r['bpp_motif']) == 7:
            m = list(r['bpp_motif'].lower())
            # In BPP 7-mers like TACTAAC, the BP A is typically position 5 (index 4) or 6 (index 5)
            # The canonical is YNYURAY where A is at position 6 (0-based: 5)
            m[5] = m[5].upper()
            motif = ''.join(m)

        print(f"  {i+1:>2}  {pos_str:>8}  {bpp_str:>10}  {svm_str:>10}  {stars:>9}  {dist_str:>12}  {motif}")

    # ── Generate detailed .md output ─────────────────────────────────────
    from datetime import date
    today = date.today().isoformat()

    if hasattr(args, 'file') and args.file:
        base = os.path.splitext(os.path.basename(args.file))[0]
        md_path = os.path.join(file_dir, f"{base}_branchpoint.md")
    else:
        md_path = os.path.join(file_dir, "branchpoint_results.md")

    # Build sequence map showing top predictions near 3' end
    seq_end = intron_seq[-50:] if len(intron_seq) >= 50 else intron_seq
    end_offset = intron_len - len(seq_end)  # 0-based offset of seq_end[0] in intron

    md_lines = []
    md_lines.append(f"# Branch Point Prediction: {name} / {intron_label}")
    md_lines.append("")
    md_lines.append(f"**Date:** {today}")
    if genomic_start and genomic_end:
        md_lines.append(f"**Intron:** positions {genomic_start}-{genomic_end} ({intron_len} bp)")
        md_lines.append(f"**3' splice site:** AG at position {genomic_end - 1}-{genomic_end}")
    else:
        md_lines.append(f"**Intron:** {intron_len} bp")
        md_lines.append(f"**3' splice site:** AG at intron positions {intron_len - 1}-{intron_len}")
    md_lines.append("")

    # Sequence map
    md_lines.append("## Sequence Map")
    md_lines.append("")
    md_lines.append("```")
    # Show the last ~50 nt with annotations
    seq_display = seq_end.lower()
    md_lines.append(f"...{seq_display}")

    # Build annotation line showing top predictions
    anno_line = list(" " * (len(seq_display) + 3))  # +3 for "..."
    label_line = list(" " * (len(seq_display) + 3))
    for i, r in enumerate(results[:3]):
        # Position in the displayed sequence
        intron_pos_0 = r['intron_pos'] - 1  # 0-based in full intron
        display_pos = intron_pos_0 - end_offset + 3  # +3 for "..."
        if 3 <= display_pos < len(anno_line):
            anno_line[display_pos] = "\u25B2"
            tag = f"BP{i+1}"
            for ci, ch in enumerate(tag):
                pos = display_pos + ci
                if pos < len(label_line):
                    label_line[pos] = ch

    # Mark 3'SS
    ss_pos = len(seq_display) + 1  # position of the last two chars
    if ss_pos < len(anno_line):
        anno_line[ss_pos] = "\u25B2"

    md_lines.append(''.join(anno_line).rstrip())
    md_lines.append(''.join(label_line).rstrip())
    md_lines.append("```")
    md_lines.append("")

    # Predictions table
    md_lines.append("## Predictions")
    md_lines.append("")

    header = "| # | Intron Pos | "
    if genomic_start:
        header += "Genomic Pos | "
    header += "BPP z-score | SVM score | Consensus | Dist to 3'SS | Motif |"
    md_lines.append(header)

    sep = "|---|---|"
    if genomic_start:
        sep += "---|"
    sep += "---|---|---|---|---|"
    md_lines.append(sep)

    for i, r in enumerate(results):
        bpp_str = f"{r['bpp_zscore']:.2f}" if r['bpp_zscore'] is not None else "---"
        svm_str = f"{r['svm_score']:.3f}" if r['svm_score'] is not None else "---"
        stars = _consensus_stars(r['consensus'])
        motif = (r['svm_motif'] or r['bpp_motif'] or "---").upper()

        row = f"| {i+1} | {r['intron_pos']} | "
        if genomic_start:
            gpos = genomic_start + r['intron_pos'] - 1
            row += f"{gpos} | "
        row += f"{bpp_str} | {svm_str} | {stars} | {r['dist_to_3ss']} nt | {motif} |"
        md_lines.append(row)

    md_lines.append("")

    # Tools section
    md_lines.append("## Tools Used")
    md_lines.append("")
    md_lines.append("- **BPP** (zhqingit/BPP) — position weight matrix + PPT scoring, z-score ranking")
    md_lines.append("- **SVM-BPfinder** (comprna/SVM-BPfinder-3M) — SVM classifier, human model")
    md_lines.append("")

    if not bpp_ok:
        md_lines.append("*Note: BPP failed to produce results for this sequence.*")
        md_lines.append("")
    if not svm_ok:
        md_lines.append("*Note: SVM-BPfinder failed to produce results for this sequence.*")
        md_lines.append("")

    md_content = '\n'.join(md_lines)
    with open(md_path, 'w') as f:
        f.write(md_content)

    print(f"\nResults saved to: {md_path}")


def cmd_varmap(args):
    """Visualize variant positions mapped onto a DNA sequence."""
    import csv as csv_mod

    record = read_genbank(args.file)
    seq = str(record.seq).upper()
    name = record.name or record.id or os.path.splitext(os.path.basename(args.file))[0]

    # ── Parse CSV ────────────────────────────────────────────────────────
    csv_path = os.path.expanduser(args.variants_csv)
    if not os.path.exists(csv_path):
        print(f"Error: file not found: {csv_path}", file=sys.stderr)
        sys.exit(1)

    variants = []  # list of (variant_name, [positions])
    with open(csv_path, "r") as fh:
        reader = csv_mod.DictReader(fh)
        # Find the "Reverted Positions" column (flexible naming)
        fieldnames = reader.fieldnames or []
        pos_col = None
        for fn in fieldnames:
            if "revert" in fn.lower() and "pos" in fn.lower():
                pos_col = fn
                break
        if pos_col is None:
            # Try exact match
            for fn in fieldnames:
                if fn.strip().lower() == "reverted positions":
                    pos_col = fn
                    break
        if pos_col is None:
            print(f"Error: CSV must have a 'Reverted Positions' column. Found: {fieldnames}", file=sys.stderr)
            sys.exit(1)

        name_col = None
        for fn in fieldnames:
            if fn.strip().lower() == "name":
                name_col = fn
                break
        if name_col is None:
            name_col = fieldnames[1] if len(fieldnames) > 1 else fieldnames[0]

        for row in reader:
            pos_str = row.get(pos_col, "").strip()
            vname = row.get(name_col, "").strip()
            if not pos_str or pos_str.lower() in ("none", "all", "wt", "-", ""):
                continue
            try:
                positions = [int(p.strip()) for p in pos_str.split(",") if p.strip()]
            except ValueError:
                continue
            if positions:
                variants.append((vname, positions))

    if not variants:
        print("No variants with valid positions found in CSV.", file=sys.stderr)
        sys.exit(1)

    # ── Coordinate mapping (intron-relative → genomic) ───────────────────
    intron_start_0 = None  # 0-based genomic start of intron
    intron_label = None
    if args.intron:
        for feat in record.features:
            label = get_label(feat)
            if label and args.intron.lower() in label.lower():
                intron_start_0 = int(feat.location.start)
                intron_end_0 = int(feat.location.end)
                intron_label = label
                break
        if intron_start_0 is None:
            print(f"Error: no feature matching '{args.intron}' found", file=sys.stderr)
            sys.exit(1)

    # Convert variant positions: if --intron, positions are 1-based intron-relative
    # Convert to 0-based genomic positions for sequence indexing
    all_positions = set()
    converted_variants = []
    for vname, positions in variants:
        if args.intron:
            # Intron positions are 1-based relative to intron start
            genomic_positions = [intron_start_0 + (p - 1) for p in positions]
        else:
            # Positions are 1-based genomic
            genomic_positions = [p - 1 for p in positions]
        converted_variants.append((vname, positions, genomic_positions))
        all_positions.update(genomic_positions)

    # ── Determine display window ─────────────────────────────────────────
    if args.region:
        disp_start_0 = args.region[0] - 1  # convert 1-based to 0-based
        disp_end_0 = args.region[1] - 1     # inclusive 0-based
    else:
        # Auto from variant positions + context
        context = 5
        disp_start_0 = max(0, min(all_positions) - context)
        disp_end_0 = min(len(seq) - 1, max(all_positions) + context)

    display_positions = list(range(disp_start_0, disp_end_0 + 1))
    n_pos = len(display_positions)
    spacing = 4  # characters per position

    # ── Determine variable region ────────────────────────────────────────
    var_min_0 = min(all_positions)
    var_max_0 = max(all_positions)

    # ── Render ───────────────────────────────────────────────────────────
    def pos_label(genomic_0):
        """Return the display label for a genomic 0-based position."""
        if args.intron:
            return str(genomic_0 - intron_start_0 + 1)  # 1-based intron position
        else:
            return str(genomic_0 + 1)  # 1-based genomic

    # 1. Position header
    pos_labels = [pos_label(p) for p in display_positions]
    max_label_len = max(len(lbl) for lbl in pos_labels)
    header_label = "Intron pos:" if args.intron else "Position:"

    # Print position numbers row by row (vertically stacked for alignment)
    # Use vertical layout if numbers are wide
    if max_label_len <= 3:
        # Single-row: pad each label to `spacing` width
        line = header_label.ljust(14)
        for lbl in pos_labels:
            line += lbl.rjust(spacing)
        print(line)
    else:
        # Multi-row: stack digits vertically
        padded = [lbl.rjust(max_label_len) for lbl in pos_labels]
        for row in range(max_label_len):
            if row == 0:
                prefix = header_label.ljust(14)
            else:
                prefix = " " * 14
            line = prefix
            for lbl in padded:
                line += lbl[row].rjust(spacing)
            print(line)

    # 2. Sequence line
    seq_line = "Sequence:".ljust(14)
    for p in display_positions:
        base = seq[p] if p < len(seq) else "?"
        seq_line += base.rjust(spacing)
    print(seq_line)

    # 3. Separator showing variable region
    print()
    sep_line = " " * 14
    for p in display_positions:
        if p < var_min_0:
            sep_line += "─" * spacing
        elif p == var_min_0:
            sep_line += "┤" + "─" * (spacing - 1)
        elif p == var_max_0:
            sep_line += "─" * (spacing - 1) + "├"
        elif var_min_0 < p < var_max_0:
            sep_line += "─" * spacing
        else:
            sep_line += "─" * spacing

    # Annotate the separator
    context_left = var_min_0 - disp_start_0
    context_right = disp_end_0 - var_max_0
    var_label = f" VARIABLE REGION ({pos_label(var_min_0)}-{pos_label(var_max_0)}) "
    var_width = (var_max_0 - var_min_0 + 1) * spacing
    if len(var_label) < var_width:
        # Center the label in the variable region
        before = 14 + context_left * spacing
        sep_line = " " * 14
        # Context left
        if context_left > 0:
            ctx = "─" * (context_left * spacing - 4)
            sep_line += ctx + " ctx ┤"
        # Variable region
        remaining = var_width
        if context_left == 0:
            sep_line += "├"
            remaining -= 1
        label_pad = remaining - len(var_label)
        left_pad = label_pad // 2
        right_pad = label_pad - left_pad
        sep_line += "─" * left_pad + var_label + "─" * right_pad
        if context_right > 0:
            sep_line += "├ ctx ─"
            sep_line += "─" * max(0, context_right * spacing - 7) + "→"
        else:
            sep_line += "┤"
    print(sep_line)

    # 4. Variant groups — group by number of reverted positions
    print()
    groups = {}  # count → list of (name, orig_positions, genomic_positions)
    for vname, orig_pos, geno_pos in converted_variants:
        count = len(geno_pos)
        groups.setdefault(count, []).append((vname, orig_pos, geno_pos))

    disp_set = set(display_positions)

    def _render_marker_line(geno_positions, prefix, suffix):
        """Render a line with * at each position, connected by ─── for adjacent positions."""
        gp_set = set(geno_positions)
        gp_sorted = sorted(geno_positions)
        line = prefix.ljust(14)
        for p in display_positions:
            if p in gp_set:
                line += "*".rjust(spacing)
            else:
                line += " " * spacing
        line += suffix
        # Now replace spaces between consecutive * positions with ───
        # Re-render with connectors
        line = prefix.ljust(14)
        for i, p in enumerate(display_positions):
            if p in gp_set:
                # Check if previous displayed position was also in set (connect)
                if i > 0 and display_positions[i - 1] in gp_set:
                    line += "───*"
                else:
                    line += "*".rjust(spacing)
            else:
                line += " " * spacing
        line += suffix
        return line

    for count in sorted(groups.keys()):
        group = groups[count]
        if count == 1:
            group_label = "SINGLE REVERTANTS:"
        elif count == 2:
            group_label = "2-POSITION:"
        elif count == 3:
            group_label = "3-POSITION:"
        else:
            group_label = f"{count}-POSITION:"

        if count == 1:
            # Singles: one line with * at every unique single-variant position
            # No connectors — each * is independent
            all_single_pos = set()
            for vname, orig_pos, geno_pos in group:
                all_single_pos.update(geno_pos)
            suffix = f"       ({len(group)} variants)"
            line = group_label.ljust(14)
            for p in display_positions:
                if p in all_single_pos:
                    line += "*".rjust(spacing)
                else:
                    line += " " * spacing
            line += suffix
            print(line)
            print()
        else:
            # Multi-position: group by unique position set
            unique_patterns = []  # (pos_set_tuple, count_matching)
            seen = set()
            for vname, orig_pos, geno_pos in group:
                key = tuple(sorted(geno_pos))
                if key not in seen:
                    seen.add(key)
                    matching_count = sum(1 for _, _, gp in group if tuple(sorted(gp)) == key)
                    unique_patterns.append((key, matching_count))

            for pos_set, match_count in unique_patterns:
                geno_pos = sorted(pos_set)
                suffix = f"       ({match_count} variant{'s' if match_count != 1 else ''})"
                print(_render_marker_line(geno_pos, group_label, suffix))
                # Use group_label only for first pattern, indent rest
                group_label = ""

                # Position label line below
                if args.intron:
                    orig_positions = [g - intron_start_0 + 1 for g in geno_pos]
                else:
                    orig_positions = [g + 1 for g in geno_pos]
                pos_range = "-".join(str(p) for p in sorted(orig_positions))
                first_idx = display_positions.index(min(geno_pos))
                last_idx = display_positions.index(max(geno_pos))
                center_col = 14 + ((first_idx + last_idx) * spacing) // 2
                label_line = " " * max(0, center_col - len(pos_range) // 2) + pos_range
                print(label_line)
                print()


def _find_flanking_exons(record, intron_start_0, intron_end_0):
    """Find exon features flanking an intron.

    Args:
        record: BioPython SeqRecord
        intron_start_0: 0-based start of the intron
        intron_end_0: 0-based exclusive end of the intron

    Returns:
        (upstream_exon, downstream_exon) where each is a (start_0, end_0) tuple or None
    """
    exon_features = [f for f in record.features if f.type == "exon"]
    # Also include misc_features labeled as exon
    for f in record.features:
        if f.type in ("misc_feature", "regulatory"):
            for key in ("label", "gene", "note", "standard_name"):
                val = f.qualifiers.get(key, [""])
                if isinstance(val, list):
                    val = " ".join(val)
                if "exon" in val.lower():
                    exon_features.append(f)
                    break

    upstream = None
    downstream = None
    for f in exon_features:
        f_start = int(f.location.start)
        f_end = int(f.location.end)
        # Upstream exon ends at intron start
        if f_end == intron_start_0:
            upstream = (f_start, f_end)
        # Downstream exon starts at intron end
        if f_start == intron_end_0:
            downstream = (f_start, f_end)
    return upstream, downstream


def _discover_introns(record, transcript_accession=None, email="user@example.com"):
    """Discover intron boundaries from a GenBank record.

    Tries multiple strategies in order:
    1. If exon-type features exist in the file, use gaps between them (best: NCBI-annotated)
    2. If transcript_accession given, fetch RefSeqGene from NCBI to get exon coordinates
    3. If intron-type features exist, use those directly
    4. Scan misc_feature/regulatory labels for 'intron' keyword, use those
    5. Error with helpful message

    Returns: list of dicts, each with:
        - 'start_0': int, 0-based intron start
        - 'end_0': int, 0-based exclusive intron end
        - 'upstream_exon': (start_0, end_0) tuple or None
        - 'downstream_exon': (start_0, end_0) tuple or None
        - 'label': str, e.g. "intron 1"

    Introns are sorted by start position.
    """
    introns = []

    # ------------------------------------------------------------------
    # Strategy 1: File already has exon features — use them directly.
    # This is the preferred path: download RefSeqGene from NCBI which
    # comes with exon annotations, then splicemap reads them.
    # ------------------------------------------------------------------
    exon_features = sorted(
        [f for f in record.features if f.type == "exon"],
        key=lambda f: int(f.location.start)
    )
    # Filter out any _SM suffixed exons (those are our own annotations)
    exon_features = [f for f in exon_features
                     if not any(str(l).endswith("_SM") for l in f.qualifiers.get("label", []))]
    if len(exon_features) >= 2:
        for i in range(len(exon_features) - 1):
            up = exon_features[i]
            dn = exon_features[i + 1]
            intron_start = int(up.location.end)
            intron_end = int(dn.location.start)
            if intron_end > intron_start:
                introns.append({
                    "start_0": intron_start,
                    "end_0": intron_end,
                    "upstream_exon": (int(up.location.start), int(up.location.end)),
                    "downstream_exon": (int(dn.location.start), int(dn.location.end)),
                    "label": f"intron {i + 1}",
                })
        if introns:
            return sorted(introns, key=lambda x: x["start_0"])

    # ------------------------------------------------------------------
    # Strategy 2: Fetch RefSeqGene from NCBI to get exon coordinates.
    # The transcript accession (e.g. NM_004992.4) is used to look up
    # the gene, then the RefSeqGene record which has exon annotations.
    # Exon sequences are mapped onto the user's genomic sequence.
    # ------------------------------------------------------------------
    if transcript_accession:
        from Bio import Entrez
        from Bio import SeqIO as _SeqIO
        Entrez.email = email

        try:
            # Fetch the mRNA as GenBank to get exon sequences
            handle = Entrez.efetch(db="nucleotide", id=transcript_accession,
                                   rettype="gb", retmode="text")
            mrna_record = _SeqIO.read(handle, "genbank")
            handle.close()

            mrna_seq = str(mrna_record.seq).upper()
            genomic = str(record.seq).upper()

            # Extract exon sequences from mRNA record
            mrna_exons = sorted(
                [f for f in mrna_record.features if f.type == "exon"],
                key=lambda f: int(f.location.start)
            )

            exon_seqs = []
            if mrna_exons:
                for ef in mrna_exons:
                    s, e = int(ef.location.start), int(ef.location.end)
                    exon_seqs.append(mrna_seq[s:e])
            else:
                # No exon features on mRNA: use the whole mRNA as one "exon"
                # and fall through to other strategies
                pass

            # Map each exon sequence onto the user's genomic sequence
            if len(exon_seqs) >= 2:
                mapped_exons = []
                search_start = 0
                for exon_seq in exon_seqs:
                    anchor = exon_seq[:min(30, len(exon_seq))]
                    gpos = genomic.find(anchor, search_start)
                    if gpos == -1:
                        gpos = genomic.find(anchor)
                    if gpos == -1:
                        continue
                    # Extend match
                    match_len = 0
                    while (match_len < len(exon_seq) and
                           gpos + match_len < len(genomic) and
                           exon_seq[match_len] == genomic[gpos + match_len]):
                        match_len += 1
                    mapped_exons.append({
                        "genomic_start_0": gpos,
                        "genomic_end_0": gpos + match_len,
                    })
                    search_start = gpos + match_len

                if len(mapped_exons) >= 2:
                    mapped_exons.sort(key=lambda e: e["genomic_start_0"])
                    for i in range(len(mapped_exons) - 1):
                        up = mapped_exons[i]
                        dn = mapped_exons[i + 1]
                        intron_start = up["genomic_end_0"]
                        intron_end = dn["genomic_start_0"]
                        if intron_end > intron_start:
                            # Snap to nearest GT-AG
                            best_5, best_5_dist = intron_start, 999
                            for off in range(-5, 6):
                                p = intron_start + off
                                if 0 <= p < len(genomic) - 1 and genomic[p:p+2] == 'GT' and abs(off) < best_5_dist:
                                    best_5, best_5_dist = p, abs(off)
                            best_3, best_3_dist = intron_end, 999
                            for off in range(-5, 6):
                                p = intron_end + off
                                if 1 <= p <= len(genomic) and genomic[p-2:p] == 'AG' and abs(off) < best_3_dist:
                                    best_3, best_3_dist = p, abs(off)
                            introns.append({
                                "start_0": best_5,
                                "end_0": best_3,
                                "upstream_exon": (up["genomic_start_0"], best_5),
                                "downstream_exon": (best_3, dn["genomic_end_0"]),
                                "label": f"intron {i + 1}",
                            })
                    if introns:
                        return sorted(introns, key=lambda x: x["start_0"])
        except Exception as e:
            print(f"Warning: could not fetch {transcript_accession}: {e}", file=sys.stderr)

    # ------------------------------------------------------------------
    # Strategy 3: Intron-type features — use directly
    # ------------------------------------------------------------------
    intron_features = sorted(
        [f for f in record.features if f.type == "intron"],
        key=lambda f: int(f.location.start)
    )
    if intron_features:
        for i, f in enumerate(intron_features):
            intron_start = int(f.location.start)
            intron_end = int(f.location.end)
            upstream, downstream = _find_flanking_exons(record, intron_start, intron_end)
            introns.append({
                "start_0": intron_start,
                "end_0": intron_end,
                "upstream_exon": upstream,
                "downstream_exon": downstream,
                "label": f"intron {i + 1}",
            })
        return sorted(introns, key=lambda x: x["start_0"])

    # ------------------------------------------------------------------
    # Strategy 4: misc_feature / regulatory with 'intron' in label
    # ------------------------------------------------------------------
    labeled_introns = []
    for f in record.features:
        if f.type not in ("misc_feature", "regulatory"):
            continue
        for key in ("label", "gene", "note", "standard_name"):
            val = f.qualifiers.get(key, [""])
            if isinstance(val, list):
                val = " ".join(val)
            if "intron" in val.lower():
                labeled_introns.append(f)
                break

    # Filter: require GT...AG boundaries (canonical splice signals) to reduce false positives
    seq_str = str(record.seq).upper()
    verified_introns = []
    for f in labeled_introns:
        s, e = int(f.location.start), int(f.location.end)
        if e - s >= 4 and s < len(seq_str) and e <= len(seq_str):
            dinuc_5 = seq_str[s:s+2]
            dinuc_3 = seq_str[e-2:e]
            if dinuc_5 in ("GT", "GC") and dinuc_3 == "AG":
                verified_introns.append(f)
    labeled_introns = verified_introns

    labeled_introns.sort(key=lambda f: int(f.location.start))
    if labeled_introns:
        for i, f in enumerate(labeled_introns):
            intron_start = int(f.location.start)
            intron_end = int(f.location.end)
            upstream, downstream = _find_flanking_exons(record, intron_start, intron_end)
            # Try to get a meaningful label from qualifiers
            label_text = f"intron {i + 1}"
            for key in ("label", "note", "standard_name"):
                val = f.qualifiers.get(key, [""])
                if isinstance(val, list):
                    val = " ".join(val)
                if val.strip():
                    label_text = val.strip()
                    break
            introns.append({
                "start_0": intron_start,
                "end_0": intron_end,
                "upstream_exon": upstream,
                "downstream_exon": downstream,
                "label": label_text,
            })
        return sorted(introns, key=lambda x: x["start_0"])

    # ------------------------------------------------------------------
    # Strategy 5: Nothing found — error with helpful message
    # ------------------------------------------------------------------
    filename = getattr(record, "_filename", None) or getattr(record, "id", "<unknown>")
    print(
        f"Error: No introns found in {filename}.\n"
        f"  - Use --transcript to provide an mRNA accession (e.g. -t NM_004992.4)\n"
        f"  - Or annotate exons first with: bio.py exons <file> -t <accession> --annotate",
        file=sys.stderr,
    )
    sys.exit(1)


def _score_splice_site_5(record_seq, intron_start_0):
    """Score the 5' splice site (donor) using MaxEntScan.

    Extracts 3 exonic bases before the intron + 6 intronic bases.
    Total window: 9 bases.

    Args:
        record_seq: full sequence as string (or Seq object)
        intron_start_0: 0-based position where the intron starts
                        (first base of the intron, typically G of GT)

    Returns:
        (score, seq_9mer) tuple, or (None, None) if insufficient flanking sequence
    """
    seq = str(record_seq).upper()
    start = intron_start_0 - 3
    end = intron_start_0 + 6
    if intron_start_0 < 3 or end > len(seq):
        return (None, None)
    seq_9mer = seq[start:end]
    matrix5, _matrix3 = _load_maxent()  # also sets up sys.path for maxentpy
    from maxentpy.maxent import score5
    try:
        return (score5(seq_9mer, matrix=matrix5), seq_9mer)
    except Exception:
        return (None, None)


def _score_splice_site_3(record_seq, intron_end_0):
    """Score the 3' splice site (acceptor) using MaxEntScan.

    Extracts 20 intronic bases before the exon + 3 exonic bases.
    Total window: 23 bases.

    Args:
        record_seq: full sequence as string (or Seq object)
        intron_end_0: 0-based exclusive end of intron
                      (= 0-based start of the downstream exon)

    Returns:
        (score, seq_23mer) tuple, or (None, None) if insufficient flanking sequence
    """
    seq = str(record_seq).upper()
    start = intron_end_0 - 20
    end = intron_end_0 + 3
    if intron_end_0 < 20 or end > len(seq):
        return (None, None)
    seq_23mer = seq[start:end]
    _matrix5, matrix3 = _load_maxent()  # also sets up sys.path for maxentpy
    from maxentpy.maxent import score3
    try:
        return (score3(seq_23mer, matrix=matrix3), seq_23mer)
    except Exception:
        return (None, None)


def _find_ese_sites(exon_seq):
    """Scan exon sequence for SR protein binding sites (ESEs) using ESEfinder matrices.

    Slides each ESEfinder matrix across the sequence. At each position, sums the
    log-odds scores for each nucleotide. If the sum >= threshold, it's a hit.

    Args:
        exon_seq: DNA sequence string (uppercase recommended)

    Returns:
        list of dicts, each with:
            'protein': str (e.g. 'SRSF1')
            'start_0': int (0-based position in exon_seq)
            'end_0': int (0-based exclusive end)
            'score': float (PWM score)
            'motif': str (the matching sequence)
    """
    seq = str(exon_seq).upper()
    hits = []
    for protein_name, data in ESEFINDER_MATRICES.items():
        matrix = data['matrix']
        threshold = data['threshold']
        window_len = len(matrix)
        for i in range(len(seq) - window_len + 1):
            kmer = seq[i:i + window_len]
            score = sum(matrix[j].get(kmer[j], -1.58) for j in range(window_len))
            if score >= threshold:
                hits.append({
                    'protein': protein_name,
                    'start_0': i,
                    'end_0': i + window_len,
                    'score': score,
                    'motif': kmer,
                })
    return hits


def _find_ess_sites(exon_seq):
    """Scan exon sequence for hnRNP binding sites (ESSs) using regex patterns.

    Args:
        exon_seq: DNA sequence string

    Returns:
        list of dicts, each with:
            'protein': str (e.g. 'hnRNP_A1')
            'start_0': int
            'end_0': int
            'motif': str
    """
    import re
    seq = str(exon_seq).upper()
    hits = []
    for protein_name, data in HNRNP_MOTIFS.items():
        pattern = data['pattern']
        for match in re.finditer(pattern, seq):
            hits.append({
                'protein': protein_name,
                'start_0': match.start(),
                'end_0': match.end(),
                'motif': match.group(),
            })
    return hits


def _merge_splicing_regions(hits, gap=3):
    """Merge adjacent or overlapping hits for the same protein into regions.

    Two hits are merged if they're for the same protein and their positions
    overlap or are within `gap` bases of each other.

    Args:
        hits: list of dicts with 'protein', 'start_0', 'end_0', 'score' (optional)
        gap: maximum gap between hits to merge (default 3)

    Returns:
        list of dicts with:
            'protein': str
            'start_0': int (start of merged region)
            'end_0': int (end of merged region)
            'top_score': float (highest score in the region, or 0.0 for ESS)
            'hit_count': int (number of individual hits merged)
    """
    from collections import defaultdict
    by_protein = defaultdict(list)
    for hit in hits:
        by_protein[hit['protein']].append(hit)

    merged = []
    for protein_name, protein_hits in by_protein.items():
        sorted_hits = sorted(protein_hits, key=lambda h: h['start_0'])
        current_start = sorted_hits[0]['start_0']
        current_end = sorted_hits[0]['end_0']
        current_top_score = sorted_hits[0].get('score', 0.0)
        current_count = 1

        for hit in sorted_hits[1:]:
            if hit['start_0'] <= current_end + gap:
                # Merge
                current_end = max(current_end, hit['end_0'])
                current_top_score = max(current_top_score, hit.get('score', 0.0))
                current_count += 1
            else:
                merged.append({
                    'protein': protein_name,
                    'start_0': current_start,
                    'end_0': current_end,
                    'top_score': current_top_score,
                    'hit_count': current_count,
                })
                current_start = hit['start_0']
                current_end = hit['end_0']
                current_top_score = hit.get('score', 0.0)
                current_count = 1

        merged.append({
            'protein': protein_name,
            'start_0': current_start,
            'end_0': current_end,
            'top_score': current_top_score,
            'hit_count': current_count,
        })

    return merged


def _clean_splicemap_features(record):
    """Remove all features created by splicemap (those with _SM suffix in label).
    Returns the cleaned record.
    """
    def _has_sm_label(feat):
        labels = feat.qualifiers.get('label', [])
        return any(str(lbl).endswith('_SM') for lbl in labels)

    record.features = [f for f in record.features if not _has_sm_label(f)]
    return record


def _splicemap_annotate(record, introns, skip_ese=False):
    """Annotate all splicing regulatory elements on a GenBank record.

    Args:
        record: BioPython SeqRecord
        introns: list from _discover_introns()
        skip_ese: if True, skip ESE/ESS scanning

    Returns:
        (record, report_data) where report_data is a dict with all scores
        and annotations for the markdown report generator.
    """

    def _confidence_level(value, strong_threshold, moderate_threshold):
        if value is None:
            return 'unknown'
        if value > strong_threshold:
            return 'strong'
        if value > moderate_threshold:
            return 'moderate'
        return 'weak'

    # Step 1: clean existing SM: features
    record = _clean_splicemap_features(record)

    record_seq = str(record.seq).upper()

    report_data = {
        'name': record.name,
        'length': len(record.seq),
        'topology': record.annotations.get('topology', 'linear'),
        'introns': [],
        'exons': [],
        'warnings': [],
    }

    new_features = []

    # Track unique exons: key = (start_0, end_0) -> exon label
    exon_registry = {}

    # Step 2: annotate each intron
    for intron in introns:
        intron_label = intron['label']  # e.g. "intron 1"
        label_key = intron_label.replace(' ', '_')  # e.g. "intron_1"
        start_0 = intron['start_0']
        end_0 = intron['end_0']
        intron_len = end_0 - start_0

        intron_report = {
            'label': intron_label,
            'start_0': start_0,
            'end_0': end_0,
            'length': intron_len,
            'score_5ss': None,
            'seq_5ss': None,
            'score_3ss': None,
            'seq_3ss': None,
            'confidence_5ss': 'unknown',
            'confidence_3ss': 'unknown',
            'bps_score': None,
            'bps_motif': None,
            'bps_dist': None,
            'confidence_bps': 'unknown',
            'ppt_length': None,
            'ppt_pyr_pct': None,
            'confidence_ppt': 'unknown',
        }

        # ── 5'SS ──────────────────────────────────────────────────────────
        score5, seq5 = _score_splice_site_5(record_seq, start_0)
        intron_report['score_5ss'] = score5
        intron_report['seq_5ss'] = seq5
        if score5 is not None:
            conf5 = _confidence_level(score5, 6, 3)
            intron_report['confidence_5ss'] = conf5
            if score5 < 3:
                report_data['warnings'].append(
                    f"Weak 5'SS for {intron_label} (score: {score5:.2f})"
                )
        else:
            conf5 = 'unknown'

        # Check canonical GT dinucleotide
        dinuc5 = record_seq[start_0:start_0 + 2]
        if dinuc5 != 'GT':
            report_data['warnings'].append(
                f"5'SS of {intron_label} is non-canonical ({dinuc5})"
            )

        note5 = f"MaxEntScan: {score5:.2f} ({conf5}). Seq: {seq5}" if score5 is not None else f"Seq: {dinuc5}"
        feat5 = SeqFeature(
            FeatureLocation(start_0, start_0 + 2, strand=1),
            type="regulatory",
            qualifiers={
                "label": [f"5' Splice Site_{label_key}_SM"],
                "ApEinfo_fwdcolor": [SPLICEMAP_COLORS['5SS']],
                "ApEinfo_revcolor": [SPLICEMAP_COLORS['5SS']],
                "note": [note5],
            }
        )
        new_features.append(feat5)

        # ── 3'SS ──────────────────────────────────────────────────────────
        score3, seq3 = _score_splice_site_3(record_seq, end_0)
        intron_report['score_3ss'] = score3
        intron_report['seq_3ss'] = seq3
        if score3 is not None:
            conf3 = _confidence_level(score3, 6, 3)
            intron_report['confidence_3ss'] = conf3
            if score3 < 3:
                report_data['warnings'].append(
                    f"Weak 3'SS for {intron_label} (score: {score3:.2f})"
                )
        else:
            conf3 = 'unknown'

        # Check canonical AG dinucleotide
        dinuc3 = record_seq[end_0 - 2:end_0]
        if dinuc3 != 'AG':
            report_data['warnings'].append(
                f"3'SS of {intron_label} is non-canonical ({dinuc3})"
            )

        note3 = f"MaxEntScan: {score3:.2f} ({conf3}). Seq: {seq3}" if score3 is not None else f"Seq: {dinuc3}"
        feat3 = SeqFeature(
            FeatureLocation(end_0 - 2, end_0, strand=1),
            type="regulatory",
            qualifiers={
                "label": [f"3' Splice Site_{label_key}_SM"],
                "ApEinfo_fwdcolor": [SPLICEMAP_COLORS['3SS']],
                "ApEinfo_revcolor": [SPLICEMAP_COLORS['3SS']],
                "note": [note3],
            }
        )
        new_features.append(feat3)

        # ── BPS ───────────────────────────────────────────────────────────
        ppt_scan_start_intron = None  # 0-based index within intron where PPT scan begins
        intron_report['bps_candidates'] = []
        try:
            _ensure_branchpoint_tools()
            intron_seq = record_seq[start_0:end_0]

            # Run BPP — take top 2
            bpp_preds = _run_bpp(intron_seq, n_results=5)
            bpp_top2 = bpp_preds[:2]

            # Run SVM-BPfinder — take top 2
            try:
                svm_preds = _run_svm_bpfinder(intron_seq)
                svm_top2 = svm_preds[:2]
            except Exception:
                svm_top2 = []

            # Collect all candidates with scores for unified ranking
            all_bps_candidates = []
            for pred in bpp_top2:
                all_bps_candidates.append(('BPP', pred, pred['zsc']))
            for pred in svm_top2:
                all_bps_candidates.append(('SVM', pred, pred['svm_scr']))
            # Sort by score descending (BPP z-score and SVM score are both higher=better)
            all_bps_candidates.sort(key=lambda x: x[2], reverse=True)

            # Rank-based colors: best candidate brightest, fading out
            bps_rank_colors = ['#C2410C', '#E8601A', '#F59E0B', '#FCD34D']

            for overall_rank, (tool, pred, _score) in enumerate(all_bps_candidates[:4]):
                rank_color = bps_rank_colors[overall_rank]

                if tool == 'BPP':
                    dist = pred['dist_to_3ss']
                    bp_a_pos_1 = intron_len - dist
                    motif_start_intron_1 = bp_a_pos_1 - 4
                    motif_start_intron_0 = max(motif_start_intron_1 - 1, 0)
                    motif_end_intron_0 = motif_start_intron_0 + 7
                    motif_len = 7

                    motif_str = pred['motif_7mer']
                    m = list(motif_str.lower())
                    if len(m) >= 6:
                        m[5] = m[5].upper()
                    motif_display = ''.join(m)
                    bps_score = pred['zsc']
                    conf_bps = _confidence_level(bps_score, 2, 0)
                    note_text = f"BPP. z-score: {bps_score:.2f}. Motif: {motif_display}. -{dist}nt from 3'SS"
                else:
                    dist = pred['dist_to_3ss']
                    bp_a_pos_1 = intron_len - dist
                    motif_start_intron_1 = bp_a_pos_1 - 4
                    motif_start_intron_0 = max(motif_start_intron_1 - 1, 0)
                    motif_end_intron_0 = motif_start_intron_0 + 9
                    motif_len = 9

                    motif_display = pred['motif_9mer']
                    bps_score = pred['svm_scr']
                    conf_bps = None
                    note_text = f"SVM-BPfinder. score: {bps_score:.3f}. Motif: {motif_display}. -{dist}nt from 3'SS"

                genomic_bps_start = start_0 + motif_start_intron_0
                genomic_bps_end = start_0 + motif_end_intron_0

                # Store rank-1 as primary BPS for backward compat + PPT anchor
                if overall_rank == 0:
                    intron_report['bps_score'] = bps_score
                    intron_report['bps_motif'] = motif_display
                    intron_report['bps_dist'] = dist
                    intron_report['confidence_bps'] = conf_bps if conf_bps else _confidence_level(bps_score, 0.5, -0.5)
                    ppt_scan_start_intron = motif_end_intron_0

                intron_report['bps_candidates'].append({
                    'tool': tool,
                    'rank': overall_rank + 1,
                    'score': bps_score,
                    'motif': motif_display,
                    'dist': dist,
                    'confidence': conf_bps,
                    'genomic_start': genomic_bps_start,
                    'genomic_end': genomic_bps_end,
                })

                feat_bps = SeqFeature(
                    FeatureLocation(genomic_bps_start, genomic_bps_end, strand=1),
                    type="regulatory",
                    qualifiers={
                        "label": [f"Branch Point {overall_rank+1}_{label_key}_SM"],
                        "ApEinfo_fwdcolor": [rank_color],
                        "ApEinfo_revcolor": [rank_color],
                        "note": [note_text],
                    }
                )
                new_features.append(feat_bps)

            # Fallback: if neither tool returned anything
            if not bpp_top2 and not svm_top2:
                report_data['warnings'].append(f"BPS prediction failed for {intron_label}")
            elif not bpp_top2 and svm_top2:
                # BPP failed — use SVM rank-1 for PPT window anchor
                report_data['warnings'].append(f"BPP failed for {intron_label}, using SVM-BPfinder for PPT anchor")
                pred = svm_top2[0]
                dist = pred['dist_to_3ss']
                bp_a_pos_1 = intron_len - dist
                motif_start_intron_0 = max(bp_a_pos_1 - 5, 0)
                ppt_scan_start_intron = motif_start_intron_0 + 9

        except Exception:
            report_data['warnings'].append(f"BPS prediction failed for {intron_label}")

        # ── PPT ───────────────────────────────────────────────────────────
        if ppt_scan_start_intron is not None:
            ppt_region = record_seq[start_0 + ppt_scan_start_intron:end_0 - 2]
            if ppt_region:
                import re as _re
                ppt_len = len(ppt_region)
                py_count = sum(1 for b in ppt_region if b in ('C', 'T'))
                pyr_pct = int(round(100 * py_count / ppt_len)) if ppt_len > 0 else 0

                # Longest uninterrupted U-run (T in DNA)
                u_runs = _re.findall(r'T+', ppt_region)
                longest_u_run = max(len(r) for r in u_runs) if u_runs else 0

                intron_report['ppt_length'] = ppt_len
                intron_report['ppt_pyr_pct'] = pyr_pct
                intron_report['ppt_longest_u_run'] = longest_u_run
                conf_ppt = _confidence_level(pyr_pct, 80, 60)
                intron_report['confidence_ppt'] = conf_ppt

                if pyr_pct < 60:
                    report_data['warnings'].append(
                        f"Weak PPT for {intron_label} ({pyr_pct}%)"
                    )

                genomic_ppt_start = start_0 + ppt_scan_start_intron
                genomic_ppt_end = end_0 - 2

                feat_ppt = SeqFeature(
                    FeatureLocation(genomic_ppt_start, genomic_ppt_end, strand=1),
                    type="regulatory",
                    qualifiers={
                        "label": [f"Polypyrimidine Tract_{label_key}_SM"],
                        "ApEinfo_fwdcolor": [SPLICEMAP_COLORS['PPT']],
                        "ApEinfo_revcolor": [SPLICEMAP_COLORS['PPT']],
                        "note": [f"{ppt_len} bp, {pyr_pct}% pyrimidine, longest U-run: {longest_u_run} ({conf_ppt})"],
                    }
                )
                new_features.append(feat_ppt)

        report_data['introns'].append(intron_report)

        # Track flanking exons for ESE/ESS scanning
        if not skip_ese:
            up_exon = intron.get('upstream_exon')  # tuple (start_0, end_0) or None
            dn_exon = intron.get('downstream_exon')
            if up_exon:
                key = (up_exon[0], up_exon[1]) if isinstance(up_exon, tuple) else (up_exon.get('start_0'), up_exon.get('end_0'))
                if key not in exon_registry and None not in key:
                    exon_registry[key] = f"exon_us_{label_key}"
            if dn_exon:
                key = (dn_exon[0], dn_exon[1]) if isinstance(dn_exon, tuple) else (dn_exon.get('start_0'), dn_exon.get('end_0'))
                if key not in exon_registry and None not in key:
                    exon_registry[key] = f"exon_ds_{label_key}"

    # Step 3: ESE/ESS annotations
    if not skip_ese:
        for (ex_start_0, ex_end_0), exon_label in exon_registry.items():
            exon_seq = record_seq[ex_start_0:ex_end_0]
            exon_len = ex_end_0 - ex_start_0
            exon_key = exon_label.replace(' ', '_')

            exon_report = {
                'label': exon_label,
                'start_0': ex_start_0,
                'end_0': ex_end_0,
                'length': exon_len,
                'ese_summary': {},
                'ess_summary': {},
                'esrseq_ese_count': 0,
                'esrseq_ess_count': 0,
            }

            # ESE sites
            ese_hits = _find_ese_sites(exon_seq)
            # Group by protein
            ese_by_protein = {}
            for hit in ese_hits:
                p = hit['protein']
                ese_by_protein.setdefault(p, []).append(hit)

            for protein, hits in ese_by_protein.items():
                merged = _merge_splicing_regions(hits, gap=3)
                density = round(len(hits) / exon_len * 100, 1) if exon_len > 0 else 0
                top_score = max(h['top_score'] for h in merged) if merged else 0
                exon_report['ese_summary'][protein] = {
                    'count': len(hits),
                    'top_score': round(top_score, 2),
                    'density': density,
                }

                color = ESEFINDER_MATRICES.get(protein, {}).get('color', '#888888')
                # Annotate top 5 merged regions
                for ri, region in enumerate(merged[:5], 1):
                    genomic_start = ex_start_0 + region['start_0']
                    genomic_end = ex_start_0 + region['end_0']
                    feat_ese = SeqFeature(
                        FeatureLocation(genomic_start, genomic_end, strand=1),
                        type="misc_feature",
                        qualifiers={
                            "label": [f"Splice Enhancer {ri}_{exon_key}_SM"],
                            "ApEinfo_fwdcolor": [color],
                            "ApEinfo_revcolor": [color],
                            "note": [f"ESEfinder: {protein}. score: {region['top_score']:.2f}, {region['hit_count']} hit(s)"],
                        }
                    )
                    new_features.append(feat_ese)

            # ESS sites
            ess_hits = _find_ess_sites(exon_seq)
            ess_by_protein = {}
            for hit in ess_hits:
                p = hit['protein']
                ess_by_protein.setdefault(p, []).append(hit)

            for protein in HNRNP_MOTIFS:
                hits = ess_by_protein.get(protein, [])
                exon_report['ess_summary'][protein] = {'count': len(hits)}
                if not hits:
                    continue
                merged = _merge_splicing_regions(hits, gap=3)
                color = HNRNP_MOTIFS.get(protein, {}).get('color', SPLICEMAP_COLORS.get(protein, '#888888'))
                for ri, region in enumerate(merged[:5], 1):
                    genomic_start = ex_start_0 + region['start_0']
                    genomic_end = ex_start_0 + region['end_0']
                    feat_ess = SeqFeature(
                        FeatureLocation(genomic_start, genomic_end, strand=1),
                        type="misc_feature",
                        qualifiers={
                            "label": [f"Splice Silencer {ri}_{exon_key}_SM"],
                            "ApEinfo_fwdcolor": [color],
                            "ApEinfo_revcolor": [color],
                            "note": [f"ESS: {protein}. {region['hit_count']} hit(s)"],
                        }
                    )
                    new_features.append(feat_ess)

            # ESRseq sites (Ke et al. 2011)
            esrseq_hits = _find_esrseq_sites(exon_seq)
            esrseq_ese = [h for h in esrseq_hits if h['type'] == 'ESE']
            esrseq_ess = [h for h in esrseq_hits if h['type'] == 'ESS']
            exon_report['esrseq_ese_count'] = len(esrseq_ese)
            exon_report['esrseq_ess_count'] = len(esrseq_ess)

            # Annotate top ESRseq ESE regions (merge overlapping, take top 5)
            if esrseq_ese:
                ese_for_merge = [{'protein': 'ESRseq_ESE', 'start_0': h['start_0'],
                                  'end_0': h['end_0'], 'score': h['score']} for h in esrseq_ese]
                merged_ese = _merge_splicing_regions(ese_for_merge, gap=3)
                for ri, region in enumerate(merged_ese[:5], 1):
                    genomic_start = ex_start_0 + region['start_0']
                    genomic_end = ex_start_0 + region['end_0']
                    new_features.append(SeqFeature(
                        FeatureLocation(genomic_start, genomic_end, strand=1),
                        type="misc_feature",
                        qualifiers={
                            "label": [f"Splice Enhancer E{ri}_{exon_key}_SM"],
                            "ApEinfo_fwdcolor": [ESRSEQ_COLORS['ESE']],
                            "ApEinfo_revcolor": [ESRSEQ_COLORS['ESE']],
                            "note": [f"ESRseq ESE. score: {region['top_score']:.3f}, {region['hit_count']} hexamer(s). Ke et al. 2011"],
                        }
                    ))

            # Annotate top ESRseq ESS regions
            if esrseq_ess:
                ess_for_merge = [{'protein': 'ESRseq_ESS', 'start_0': h['start_0'],
                                  'end_0': h['end_0'], 'score': h['score']} for h in esrseq_ess]
                merged_ess = _merge_splicing_regions(ess_for_merge, gap=3)
                for ri, region in enumerate(merged_ess[:5], 1):
                    genomic_start = ex_start_0 + region['start_0']
                    genomic_end = ex_start_0 + region['end_0']
                    new_features.append(SeqFeature(
                        FeatureLocation(genomic_start, genomic_end, strand=1),
                        type="misc_feature",
                        qualifiers={
                            "label": [f"Splice Silencer S{ri}_{exon_key}_SM"],
                            "ApEinfo_fwdcolor": [ESRSEQ_COLORS['ESS']],
                            "ApEinfo_revcolor": [ESRSEQ_COLORS['ESS']],
                            "note": [f"ESRseq ESS. score: {region['top_score']:.3f}, {region['hit_count']} hexamer(s). Ke et al. 2011"],
                        }
                    ))

            report_data['exons'].append(exon_report)

    # Append all new features to the record
    record.features.extend(new_features)

    return record, report_data


def _annotate_single_intron(intron_seq, start_0, end_1, label_suffix=""):
    """
    Annotate splice signals for a single intron.

    Args:
        intron_seq: uppercase DNA sequence of the intron
        start_0: 0-based genomic start of the intron
        end_1: 0-based exclusive genomic end of the intron (= 1-based inclusive end)
        label_suffix: string appended to feature labels (e.g. "_intron1"), or "" for no suffix

    Returns:
        (new_features, checks)
        new_features: list of SeqFeature objects
        checks: list of (signal, position_str, size_str, details_str, status_str) tuples
    """
    intron_len = len(intron_seq)
    new_features = []
    checks = []

    ss5_label = f"5SS{label_suffix}"
    ss3_label = f"3SS{label_suffix}"
    bps_label = f"BPS{label_suffix}"
    ppt_label = f"PPT{label_suffix}"

    # ── 5'SS annotation ───────────────────────────────────────────────────
    ss5_start = start_0        # 0-based
    ss5_end = start_0 + 2      # 0-based exclusive
    ss5_seq = intron_seq[:2]
    new_features.append(SeqFeature(
        FeatureLocation(ss5_start, ss5_end, strand=1),
        type="regulatory",
        qualifiers={
            "label": [ss5_label],
            "ApEinfo_fwdcolor": ["#00CC00"],
            "ApEinfo_revcolor": ["#00CC00"],
        }
    ))
    if ss5_seq == "GT":
        ss5_status = "GT \u2713"
    elif ss5_seq == "GC":
        ss5_status = "GC \u2713 (non-canonical)"
    else:
        ss5_status = f"\u26a0 {ss5_seq}, expected GT"
    checks.append((ss5_label, f"{ss5_start + 1}-{ss5_end}", f"{ss5_end - ss5_start} bp", ss5_seq, ss5_status))

    # ── 3'SS annotation ───────────────────────────────────────────────────
    ss3_start = end_1 - 2      # 0-based
    ss3_end = end_1            # 0-based exclusive
    ss3_seq = intron_seq[-2:]
    new_features.append(SeqFeature(
        FeatureLocation(ss3_start, ss3_end, strand=1),
        type="regulatory",
        qualifiers={
            "label": [ss3_label],
            "ApEinfo_fwdcolor": ["#CC0000"],
            "ApEinfo_revcolor": ["#CC0000"],
        }
    ))
    if ss3_seq == "AG":
        ss3_status = "AG \u2713"
    else:
        ss3_status = f"\u26a0 {ss3_seq}, expected AG"
    checks.append((ss3_label, f"{ss3_start + 1}-{ss3_end}", f"{ss3_end - ss3_start} bp", ss3_seq, ss3_status))

    # ── BPS annotation ────────────────────────────────────────────────────
    bps_genomic_start = None
    bps_genomic_end = None
    ppt_scan_start_intron = None  # 0-based intron index to start PPT scan from
    bps_dist_to_3ss = None

    try:
        _ensure_branchpoint_tools()
        bpp_preds = _run_bpp(intron_seq, n_results=5)

        if bpp_preds:
            top = bpp_preds[0]
            dist = top['dist_to_3ss']
            bp_a_pos_1 = intron_len - dist  # 1-based position of BP A in intron
            motif_start_intron_1 = bp_a_pos_1 - 4  # 1-based start of 7-mer
            motif_start_intron_0 = motif_start_intron_1 - 1  # 0-based

            if motif_start_intron_0 < 0:
                motif_start_intron_0 = 0
            motif_end_intron_0 = motif_start_intron_0 + 7  # 0-based exclusive

            bps_genomic_start = start_0 + motif_start_intron_0
            bps_genomic_end = start_0 + motif_end_intron_0
            bps_dist_to_3ss = dist

            motif_str = top['motif_7mer']
            m = list(motif_str.lower())
            if len(m) >= 6:
                m[5] = m[5].upper()
            motif_display = ''.join(m)

            score = top['zsc']

            new_features.append(SeqFeature(
                FeatureLocation(bps_genomic_start, bps_genomic_end, strand=1),
                type="regulatory",
                qualifiers={
                    "label": [bps_label],
                    "ApEinfo_fwdcolor": ["#FF6600"],
                    "ApEinfo_revcolor": ["#FF6600"],
                    "note": [f"score: {score:.2f}, motif: {motif_display}"],
                }
            ))

            # Validate BPS distance to 3'SS (19-37 nt is canonical)
            if 19 <= dist <= 37:
                bps_dist_status = f"-{dist}nt \u2713"
            else:
                bps_dist_status = f"\u26a0 -{dist}nt from 3'SS (expected 19-37)"

            checks.append((
                bps_label,
                f"{bps_genomic_start + 1}-{bps_genomic_end}",
                f"{bps_genomic_end - bps_genomic_start} bp",
                f"{motif_display} (score: {score:.2f}) {bps_dist_status}",
                bps_dist_status,
            ))

            ppt_scan_start_intron = motif_end_intron_0

        else:
            print("Warning: BPP produced no predictions. Skipping BPS and PPT annotations.", file=sys.stderr)

    except Exception as e:
        print(f"Warning: BPP failed ({e}). Skipping BPS and PPT annotations.", file=sys.stderr)

    # ── PPT annotation ────────────────────────────────────────────────────
    ppt_genomic_start = None
    ppt_genomic_end = None

    if ppt_scan_start_intron is not None:
        ppt_region_intron = intron_seq[ppt_scan_start_intron:intron_len - 2]

        if ppt_region_intron:
            ppt_len = len(ppt_region_intron)
            py_count = sum(1 for b in ppt_region_intron if b in ('C', 'T'))
            pyr_pct = int(round(100 * py_count / ppt_len)) if ppt_len > 0 else 0

            ppt_genomic_start = start_0 + ppt_scan_start_intron
            ppt_genomic_end = start_0 + intron_len - 2

            new_features.append(SeqFeature(
                FeatureLocation(ppt_genomic_start, ppt_genomic_end, strand=1),
                type="regulatory",
                qualifiers={
                    "label": [ppt_label],
                    "ApEinfo_fwdcolor": ["#FFCC00"],
                    "ApEinfo_revcolor": ["#FFCC00"],
                    "note": [f"{ppt_len} bp, {pyr_pct}% pyrimidine (U2AF65 binding region)"],
                }
            ))

            if pyr_pct > 60:
                ppt_status = f"{pyr_pct}% \u2713"
            else:
                ppt_status = f"\u26a0 {pyr_pct}% pyrimidine (weak)"

            checks.append((
                ppt_label,
                f"{ppt_genomic_start + 1}-{ppt_genomic_end}",
                f"{ppt_len} bp",
                f"{pyr_pct}% pyrimidine",
                ppt_status,
            ))

    # ── Order check ───────────────────────────────────────────────────────
    positions = {
        "5SS": ss5_start,
        "BPS": bps_genomic_start,
        "PPT": ppt_genomic_start,
        "3SS": ss3_start,
    }
    order_ok = True
    if positions["BPS"] is not None and positions["5SS"] >= positions["BPS"]:
        order_ok = False
    if positions["PPT"] is not None and positions["BPS"] is not None and positions["BPS"] >= positions["PPT"]:
        order_ok = False
    if positions["PPT"] is not None and positions["PPT"] >= positions["3SS"]:
        order_ok = False
    if positions["BPS"] is not None and positions["BPS"] >= positions["3SS"]:
        order_ok = False
    if not order_ok:
        print("  \u26a0 Feature order check failed: expected 5'SS < BPS < PPT < 3'SS", file=sys.stderr)

    return new_features, checks


def cmd_splice_signals(args):
    """Annotate all four splice signals (5'SS, 3'SS, BPS, PPT) on an intron.

    If --intron is provided, annotates a single intron at the given coordinates.
    If --intron is not provided, auto-detects introns from exon annotations.
    """
    record = read_genbank(args.file)
    seq_len = len(record.seq)

    # ── Build list of (start_1, end_1, label) introns to process ──────────
    if args.intron is not None:
        # Manual mode: single intron, no suffix on labels
        start_1, end_1 = args.intron
        if start_1 < 1 or end_1 > seq_len or start_1 >= end_1:
            print(f"Error: intron coordinates out of range (sequence is {seq_len} bp)", file=sys.stderr)
            sys.exit(1)
        introns_to_process = [(start_1, end_1, "")]
        auto_mode = False
    else:
        # Auto-detect mode: find introns from exon gaps
        exon_features = [f for f in record.features if f.type == "exon"]
        if not exon_features:
            print("Error: no exon features found in file. Use --intron START END to specify coordinates manually.", file=sys.stderr)
            sys.exit(1)

        # Sort exons by start position (0-based)
        exon_features.sort(key=lambda f: int(f.location.start))

        introns_to_process = []
        for i in range(len(exon_features) - 1):
            exon_end_0 = int(exon_features[i].location.end)      # 0-based exclusive end of exon i
            next_exon_start_0 = int(exon_features[i + 1].location.start)  # 0-based start of exon i+1
            if next_exon_start_0 > exon_end_0:
                intron_start_1 = exon_end_0 + 1   # 1-based
                intron_end_1 = next_exon_start_0   # 1-based inclusive (= 0-based exclusive of next exon)
                introns_to_process.append((intron_start_1, intron_end_1, f"_intron{i + 1}"))

        if not introns_to_process:
            print("Error: exon features found but no gaps between them — no introns detected.", file=sys.stderr)
            sys.exit(1)

        print(f"Auto-detected {len(introns_to_process)} intron(s) from {len(exon_features)} exon annotation(s).")
        auto_mode = True

    # ── Annotate each intron ───────────────────────────────────────────────
    all_new_features = []

    for intron_start_1, intron_end_1, label_suffix in introns_to_process:
        start_0 = intron_start_1 - 1
        end_1_exclusive = intron_end_1  # 0-based exclusive = 1-based inclusive end

        intron_seq = str(record.seq[start_0:end_1_exclusive]).upper()
        intron_len = len(intron_seq)
        intron_num = label_suffix.replace("_intron", "") if label_suffix else None

        if auto_mode:
            print(f"\nIntron {intron_num} ({intron_start_1}-{intron_end_1}, {intron_len:,} bp)")
        else:
            print(f"\nIntron ({intron_start_1}-{intron_end_1}, {intron_len:,} bp)")

        new_features, checks = _annotate_single_intron(intron_seq, start_0, end_1_exclusive, label_suffix)
        all_new_features.extend(new_features)

        # Print summary table for this intron
        for signal, position, size, details, status in checks:
            print(f"  {signal:<12}  {position:<18}  {size:<8}  {status}")

    print()

    # ── Write all features and open ───────────────────────────────────────
    for feat in all_new_features:
        record.features.append(feat)

    write_genbank(record, args.file)
    print(f"Added {len(all_new_features)} splice signal annotation(s) to {args.file}")
    print(f"Backup saved to {args.file}.bak")
    open_in_viewer(args.file)


def _generate_splicemap_report(report_data, output_path):
    """Generate a markdown report from splicemap report_data dict."""
    from datetime import date

    name = report_data.get('name', '?')
    length = report_data.get('length', 0)
    topology = report_data.get('topology', 'unknown')
    introns = report_data.get('introns', [])
    exons = report_data.get('exons', [])
    warnings = report_data.get('warnings', [])

    today = date.today().strftime('%Y-%m-%d')

    def fmt_score(score):
        if score is None:
            return 'n/a'
        return f'{score:.1f}'

    def fmt_conf(conf):
        if not conf:
            return '?'
        return conf.capitalize()

    def short_conf(conf):
        if not conf:
            return '?'
        c = conf.lower()
        if c == 'strong':
            return 'strong'
        if c == 'moderate':
            return 'mod'
        if c == 'weak':
            return 'weak'
        return c

    def status_for_intron(intron):
        confs = [
            intron.get('confidence_5ss', ''),
            intron.get('confidence_3ss', ''),
            intron.get('confidence_bps', ''),
            intron.get('confidence_ppt', ''),
        ]
        weak = [c for c in confs if c and c.lower() == 'weak']
        moderate = [c for c in confs if c and c.lower() == 'moderate']
        if weak:
            return 'Weak 5SS/3SS/BPS/PPT'
        if not moderate:
            return 'All strong'
        # List which elements are moderate
        labels = []
        if intron.get('confidence_5ss', '').lower() == 'moderate':
            labels.append("5'SS")
        if intron.get('confidence_3ss', '').lower() == 'moderate':
            labels.append("3'SS")
        if intron.get('confidence_bps', '').lower() == 'moderate':
            labels.append('BPS')
        if intron.get('confidence_ppt', '').lower() == 'moderate':
            labels.append('PPT')
        return 'Moderate: ' + ', '.join(labels)

    lines = []
    lines.append(f'# Splice Map: {name}')
    lines.append('')
    lines.append(f'**Date:** {today}')
    lines.append(f'**Sequence:** {length:,} bp, {topology}')
    lines.append(f'**Introns:** {len(introns)}')
    lines.append('')
    lines.append('---')
    lines.append('')

    # Summary table
    lines.append('## Summary')
    lines.append('')
    lines.append("| Intron | Length | 5'SS | 3'SS | BPS (z) | PPT | U-run | Status |")
    lines.append('|--------|--------|------|------|---------|-----|-------|--------|')
    for intron in introns:
        label = intron.get('label', '?')
        ilen = intron.get('length', 0)
        s5 = f"{fmt_score(intron.get('score_5ss'))} ({short_conf(intron.get('confidence_5ss'))})"
        s3 = f"{fmt_score(intron.get('score_3ss'))} ({short_conf(intron.get('confidence_3ss'))})"
        bps = f"{fmt_score(intron.get('bps_score'))} ({short_conf(intron.get('confidence_bps'))})"
        ppt_pct = intron.get('ppt_pyr_pct')
        ppt_str = f"{ppt_pct}% ({short_conf(intron.get('confidence_ppt'))})" if ppt_pct is not None else 'n/a'
        urun = str(intron.get('ppt_longest_u_run', 'n/a'))
        status = status_for_intron(intron)
        lines.append(f'| {label} | {ilen} | {s5} | {s3} | {bps} | {ppt_str} | {urun} | {status} |')
    lines.append('')
    lines.append('---')
    lines.append('')

    # Per-intron detail sections
    for intron in introns:
        label = intron.get('label', '?')
        start_1 = intron.get('start_0', 0) + 1
        end_1 = intron.get('end_0', 0)
        ilen = intron.get('length', 0)

        lines.append(f'## {label} ({start_1:,}\u2013{end_1:,}, {ilen} bp)')
        lines.append('')

        # Find matching upstream exon (exon whose end_0 == intron start_0)
        up_exon = None
        dn_exon = None
        for exon in exons:
            if exon.get('end_0') == intron.get('start_0'):
                up_exon = exon
            if exon.get('start_0') == intron.get('end_0'):
                dn_exon = exon

        # ASCII diagram
        lines.append('### Splice Signals')
        lines.append('')
        lines.append('```')

        if up_exon:
            exon_len = up_exon.get('length', 0)
            exon_name = up_exon.get('label', 'exon')
            lines.append('     SF1   U2AF65  U2AF35              SR proteins              U1')
            lines.append('      \u25bc      \u25bc       \u25bc                \u25bc\u25bc\u25bc\u25bc\u25bc\u25bc\u25bc\u25bc\u25bc                \u25bc')
            lines.append('\u2500\u2500 BPS \u2500\u2500 PPT \u2500\u2500 AG \u2502 . . E S E . . E S E . . \u2502 GT \u2500\u2500 intron \u2500\u2500')
            lines.append(f'                     \u25c4\u2500\u2500\u2500\u2500\u2500 {exon_name} ({exon_len} bp) \u2500\u2500\u2500\u2500\u2500\u25ba')
            lines.append("                  3'SS                          5'SS")
        elif dn_exon:
            exon_len = dn_exon.get('length', 0)
            exon_name = dn_exon.get('label', 'exon')
            lines.append('     U1                              SR proteins    SF1   U2AF65  U2AF35')
            lines.append('      \u25bc                             \u25bc\u25bc\u25bc\u25bc\u25bc\u25bc\u25bc\u25bc\u25bc              \u25bc      \u25bc       \u25bc')
            lines.append('\u2500\u2500 intron \u2500\u2500 GT \u2502 . . E S E . . E S E . . \u2502 AG \u2500\u2500 PPT \u2500\u2500 BPS \u2500\u2500')
            lines.append(f'                    \u25c4\u2500\u2500\u2500\u2500\u2500 {exon_name} ({exon_len} bp) \u2500\u2500\u2500\u2500\u2500\u25ba')
            lines.append("                  5'SS                          3'SS")
        else:
            lines.append('\u2500\u2500 BPS \u2500\u2500 PPT \u2500\u2500 3\'SS \u2500\u2500 EXON \u2500\u2500 5\'SS \u2500\u2500 intron \u2500\u2500')

        lines.append('```')
        lines.append('')

        # Element details table
        lines.append('### Element Details')
        lines.append('')
        lines.append('| Element | Position | Score | Confidence | Details |')
        lines.append('|---------|----------|-------|------------|---------|')

        # 5'SS
        s5_start = intron.get('start_0', 0) + 1
        s5_end = intron.get('start_0', 0) + 2
        s5_score = intron.get('score_5ss')
        s5_conf = fmt_conf(intron.get('confidence_5ss'))
        s5_seq = intron.get('seq_5ss', '')
        lines.append(f"| 5'SS | {s5_start:,}\u2013{s5_end:,} | {fmt_score(s5_score)} | {s5_conf} | GT, seq: {s5_seq} |")

        # 3'SS
        s3_end = intron.get('end_0', 0)
        s3_start = s3_end - 1
        s3_score = intron.get('score_3ss')
        s3_conf = fmt_conf(intron.get('confidence_3ss'))
        s3_seq = intron.get('seq_3ss', '')
        # Truncate long 3'SS sequences
        if len(s3_seq) > 11:
            s3_seq_disp = '...' + s3_seq[-11:]
        else:
            s3_seq_disp = s3_seq
        lines.append(f"| 3'SS | {s3_start:,}\u2013{s3_end:,} | {fmt_score(s3_score)} | {s3_conf} | AG, seq: {s3_seq_disp} |")

        # BPS — show all candidates (BPP and SVM separately)
        bps_candidates = intron.get('bps_candidates', [])
        # Fallback to legacy fields if bps_candidates is empty (backward compat)
        if bps_candidates:
            for cand in bps_candidates:
                tool = cand.get('tool', '?')
                rank = cand.get('rank', '?')
                score = cand.get('score')
                motif = cand.get('motif', '')
                dist = cand.get('dist')
                conf = fmt_conf(cand.get('confidence')) if cand.get('confidence') else '-'
                g_start = cand.get('genomic_start')
                g_end = cand.get('genomic_end')
                pos_str = f'{g_start + 1:,}\u2013{g_end:,}' if g_start is not None else 'n/a'
                dist_str = f'-{dist}nt from 3\'SS' if dist is not None else ''
                if tool == 'BPP':
                    score_str = f'z={fmt_score(score)}'
                    detail = f'{motif}, {dist_str}'
                else:
                    score_str = f'score={score:.3f}' if score is not None else 'n/a'
                    detail = f'{motif}, {dist_str}'
                lines.append(f'| BPS ({tool} #{rank}) | {pos_str} | {score_str} | {conf} | {detail} |')
        else:
            # Legacy single BPS fallback
            bps_dist = intron.get('bps_dist')
            bps_score = intron.get('bps_score')
            bps_motif = intron.get('bps_motif', '')
            bps_conf = fmt_conf(intron.get('confidence_bps'))
            if bps_dist is not None:
                bps_pos_end = intron.get('end_0', 0) - bps_dist
                bps_pos_start = bps_pos_end - 6
                bps_detail = f'{bps_motif}, -{bps_dist}nt from 3\'SS'
                lines.append(f'| BPS | {bps_pos_start:,}\u2013{bps_pos_end:,} | z={fmt_score(bps_score)} | {bps_conf} | {bps_detail} |')
            else:
                lines.append(f'| BPS | n/a | z={fmt_score(bps_score)} | {bps_conf} | {bps_motif} |')

        # PPT
        bps_dist = intron.get('bps_dist')
        ppt_len = intron.get('ppt_length')
        ppt_pct = intron.get('ppt_pyr_pct')
        ppt_conf = fmt_conf(intron.get('confidence_ppt'))
        longest_u_run = intron.get('ppt_longest_u_run')
        urun_str = f', longest U-run: {longest_u_run}' if longest_u_run is not None else ''
        if bps_dist is not None and ppt_len is not None:
            ppt_end = intron.get('end_0', 0) - bps_dist + 1  # just after BPS region
            ppt_start = ppt_end - ppt_len
            lines.append(f'| PPT | {ppt_start:,}\u2013{ppt_end:,} | {ppt_pct}% pyr | {ppt_conf} | {ppt_len} bp{urun_str} |')
        else:
            lines.append(f'| PPT | n/a | {ppt_pct}% pyr | {ppt_conf} | {ppt_len} bp{urun_str} |')

        lines.append('')

        # ESE/ESS table for flanking exon
        target_exon = up_exon if up_exon else dn_exon
        if target_exon:
            exon_label = target_exon.get('label', 'exon')
            exon_len = target_exon.get('length', 0)
            ese_summary = target_exon.get('ese_summary', {})
            ess_summary = target_exon.get('ess_summary', {})

            if ese_summary or ess_summary:
                lines.append(f'### ESE/ESS in flanking exon ({exon_len} bp)')
                lines.append('')
                lines.append('| Protein | Type | Sites | Top Score | Density/100nt |')
                lines.append('|---------|------|-------|-----------|---------------|')

                for protein, data in sorted(ese_summary.items()):
                    count = data.get('count', 0)
                    if count == 0:
                        continue
                    top = data.get('top_score')
                    density = data.get('density', 0)
                    top_str = f'{top:.2f}' if top is not None else '-'
                    lines.append(f'| {protein} | ESE | {count} | {top_str} | {density:.1f} |')

                for protein, data in sorted(ess_summary.items()):
                    count = data.get('count', 0)
                    if count == 0:
                        continue
                    top = data.get('top_score')
                    density = data.get('density', 0)
                    top_str = f'{top:.2f}' if top is not None else '-'
                    lines.append(f'| {protein} | ESS | {count} | {top_str} | {density:.1f} |')

                lines.append('')

        lines.append('---')
        lines.append('')

    # Warnings
    lines.append('## Warnings')
    lines.append('')
    if warnings:
        for w in warnings:
            lines.append(f'- {w}')
    else:
        lines.append('(none)')
    lines.append('')
    lines.append('---')
    lines.append('')
    lines.append('*Generated by bio.py splicemap*')
    lines.append('')

    with open(output_path, 'w') as f:
        f.write('\n'.join(lines))


def _print_splicemap_summary(report):
    """Print a terminal summary of the splicemap results."""
    name = report.get('name', '?')
    length = report.get('length', 0)
    topology = report.get('topology', 'linear')
    intron_count = len(report.get('introns', []))

    print(f"Splice Map: {name} ({length:,} bp, {topology})")
    print(f"{'='*60}")

    # Intron summary table
    if report['introns']:
        col_5ss = "5'SS"
        col_3ss = "3'SS"
        print(f"\n{'Intron':<20} {'Length':>8} {col_5ss:>8} {col_3ss:>8} {'BPS(z)':>8} {'PPT%':>6} {'U-run':>6}")
        print(f"{'-'*20} {'-'*8} {'-'*8} {'-'*8} {'-'*8} {'-'*6} {'-'*6}")
        for i in report['introns']:
            s5 = f"{i['score_5ss']:.1f}" if i['score_5ss'] is not None else "N/A"
            s3 = f"{i['score_3ss']:.1f}" if i['score_3ss'] is not None else "N/A"
            bps = f"{i.get('bps_score', 0):.1f}" if i.get('bps_score') is not None else "N/A"
            ppt = f"{i.get('ppt_pyr_pct', 0)}%" if i.get('ppt_pyr_pct') is not None else "N/A"
            urun = str(i.get('ppt_longest_u_run', 'N/A'))
            print(f"{i['label']:<20} {i['length']:>7,} {s5:>8} {s3:>8} {bps:>8} {ppt:>6} {urun:>6}")

    # ESE/ESS summary per exon
    if report.get('exons'):
        print(f"\n{'Exon':<20} {'Length':>6}  {'ESEfinder':>9}  {'hnRNP':>5}  {'ESRseq+':>7}  {'ESRseq-':>7}")
        print(f"{'-'*20} {'-'*6}  {'-'*9}  {'-'*5}  {'-'*7}  {'-'*7}")
        for e in report['exons']:
            ese_total = sum(info['count'] for info in e.get('ese_summary', {}).values())
            ess_total = sum(info['count'] for info in e.get('ess_summary', {}).values())
            esrseq_ese = e.get('esrseq_ese_count', 0)
            esrseq_ess = e.get('esrseq_ess_count', 0)
            print(f"{e['label']:<20} {e['length']:>5}  {ese_total:>9}  {ess_total:>5}  {esrseq_ese:>7}  {esrseq_ess:>7}")

    # Warnings
    if report.get('warnings'):
        print(f"\nWarnings:")
        for w in report['warnings']:
            print(f"  \u26a0 {w}")

    # Feature count
    sm_count = len(report.get('introns', [])) * 4  # rough estimate
    print(f"\nAnnotations added (_SM suffix). Re-run to update, --clean to remove.")


def cmd_splicemap(args):
    """Map all splicing regulatory elements on a sequence.

    Discovers introns (from exon annotations, intron features, or transcript alignment),
    then annotates: 5'SS, 3'SS, branch points, PPT, ESE sites, and ESS sites.
    Writes color-coded annotations to the GenBank file and generates a markdown report.
    """
    record = read_genbank(args.file)

    if args.clean:
        record = _clean_splicemap_features(record)
        write_genbank(record, args.file)
        print(f"Removed splicemap annotations from {args.file}")
        open_in_viewer(args.file)
        return

    # Discover introns
    transcript = getattr(args, 'transcript', None)
    introns = _discover_introns(record, transcript_accession=transcript)
    print(f"Found {len(introns)} intron(s)")
    for i in introns:
        print(f"  {i['label']}: {i['start_0']+1}-{i['end_0']} ({i['end_0']-i['start_0']:,} bp)")
    print()

    # Annotate
    skip_ese = getattr(args, 'no_ese', False)
    record, report = _splicemap_annotate(record, introns, skip_ese=skip_ese)

    # Write annotated GenBank
    write_genbank(record, args.file)

    # Print summary table
    _print_splicemap_summary(report)

    # Generate markdown report (for now just print path, report generator comes next)
    report_path = args.file.rsplit('.', 1)[0] + '_splicemap.md'
    _generate_splicemap_report(report, report_path)
    print(f"\nReport: {report_path}")

    # Open file
    open_in_viewer(args.file)


def main():
    parser = argparse.ArgumentParser(
        prog="splicemap",
        description="Molecular biology CLI for GenBank files.",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    subs = parser.add_subparsers(dest="command", metavar="command")
    subs.required = True

    # current
    p = subs.add_parser("current", help="Show files open in UGENE")
    p.set_defaults(func=cmd_current)

    # read
    p = subs.add_parser("read", help="Parse .gb/.fasta, show clean summary")
    p.add_argument("file", help="GenBank or FASTA file")
    p.set_defaults(func=cmd_read)

    # features
    p = subs.add_parser("features", help="List all annotations as a table")
    p.add_argument("file", help="GenBank file")
    p.set_defaults(func=cmd_features)

    # seq
    p = subs.add_parser("seq", help="Extract a sequence region (1-based, inclusive)")
    p.add_argument("file", help="GenBank or FASTA file")
    p.add_argument("start", type=int, help="Start position (1-based)")
    p.add_argument("end", type=int, help="End position (1-based, inclusive)")
    p.set_defaults(func=cmd_seq)

    # translate
    p = subs.add_parser("translate", help="Extract a region and translate to protein (1-based, inclusive)")
    p.add_argument("file", help="GenBank or FASTA file")
    p.add_argument("start", type=int, help="Start position (1-based)")
    p.add_argument("end", type=int, help="End position (1-based, inclusive)")
    p.add_argument("--frame", type=int, choices=[1, 2, 3], default=1, help="Reading frame offset (default: 1)")
    p.add_argument("--reverse", action="store_true", help="Translate the reverse complement of the region")
    p.set_defaults(func=cmd_translate)

    # annotate
    p = subs.add_parser("annotate", help="Add a feature annotation")
    p.add_argument("file", help="GenBank file")
    p.add_argument("start", type=int, help="Start position (1-based)")
    p.add_argument("end", type=int, help="End position (1-based, inclusive)")
    p.add_argument("label", help="Feature label")
    p.add_argument("--type", default="misc_feature", help="Feature type (default: misc_feature)")
    p.add_argument("--color", default=None, help="SnapGene color as hex (e.g. #FF0000 for red)")
    p.add_argument("--note", default=None, help="Add a note/comment to the annotation")
    p.set_defaults(func=cmd_annotate)

    # remove
    p = subs.add_parser("remove", help="Remove annotation by label")
    p.add_argument("file", help="GenBank file")
    p.add_argument("label", help="Feature label to remove")
    p.set_defaults(func=cmd_remove)

    # orfs
    p = subs.add_parser("orfs", help="Find ORFs via ugenecl")
    p.add_argument("file", help="GenBank or FASTA file")
    p.add_argument("--min-length", type=int, default=100, help="Minimum ORF length in bp (default: 100)")
    p.set_defaults(func=cmd_orfs)

    # sites
    p = subs.add_parser("sites", help="Find restriction sites")
    p.add_argument("file", help="GenBank file")
    p.add_argument("--enzymes", help="Comma-separated enzyme names (default: common set)")
    p.set_defaults(func=cmd_sites)

    # revcomp
    p = subs.add_parser("revcomp", help="Reverse complement a sequence")
    p.add_argument("file", help="GenBank or FASTA file")
    p.add_argument("--out", help="Output file path (default: <name>_rc.gb)")
    p.set_defaults(func=cmd_revcomp)

    # diff
    p = subs.add_parser("diff", help="Compare two constructs")
    p.add_argument("file1", help="First GenBank file")
    p.add_argument("file2", help="Second GenBank file")
    p.set_defaults(func=cmd_diff)

    # search
    p = subs.add_parser("search", help="Search for a DNA sequence motif (both strands by default)")
    p.add_argument("file", help="GenBank or FASTA file")
    p.add_argument("sequence", help="DNA sequence to search for (case-insensitive, exact match)")
    p.add_argument("--fwd", action="store_true", help="Search forward strand only")
    p.add_argument("--rc", action="store_true", help="Search reverse complement strand only")
    p.set_defaults(func=cmd_search)

    # annotate-seq
    p = subs.add_parser("annotate-seq", help="Find a sequence and annotate its position")
    p.add_argument("file", help="GenBank file")
    p.add_argument("sequence", help="DNA sequence to search for (case-insensitive, exact match)")
    p.add_argument("label", help="Feature label")
    p.add_argument("--antisense", action="store_true", help="Reverse complement the input sequence before searching")
    p.add_argument("--type", default="misc_feature", help="Feature type (default: misc_feature)")
    p.add_argument("--color", default=None, help="SnapGene color as hex (e.g. #FF0000 for red)")
    p.add_argument("--note", default=None, help="Add a note/comment to the annotation")
    p.set_defaults(func=cmd_annotate_seq)

    # insert
    p = subs.add_parser("insert", help="Insert a DNA sequence at a given position")
    p.add_argument("file", help="GenBank file")
    p.add_argument("position", type=int, help="Insert position (1-based, inserts BEFORE this position)")
    p.add_argument("sequence", help="DNA sequence to insert (only ACGT, case-insensitive)")
    p.add_argument("--label", help="Label for the inserted region annotation (optional)")
    p.add_argument("--type", default="misc_feature", help="Feature type for annotation (default: misc_feature)")
    p.add_argument("--no-mark", action="store_true", help="Skip auto-annotation marking the edit")
    p.set_defaults(func=cmd_insert)

    # delete
    p = subs.add_parser("delete", help="Delete a DNA region (1-based, inclusive), shift downstream features, write back")
    p.add_argument("file", help="GenBank file")
    p.add_argument("start", type=int, help="Start position (1-based, inclusive)")
    p.add_argument("end", type=int, help="End position (1-based, inclusive)")
    p.add_argument("--no-mark", action="store_true", help="Skip auto-annotation marking the edit")
    p.set_defaults(func=cmd_delete)

    # replace
    p = subs.add_parser("replace", help="Replace a DNA region with a new sequence (1-based, inclusive), adjust features, write back")
    p.add_argument("file", help="GenBank file")
    p.add_argument("start", type=int, help="Start position (1-based, inclusive)")
    p.add_argument("end", type=int, help="End position (1-based, inclusive)")
    p.add_argument("sequence", help="Replacement DNA sequence (only ACGT, case-insensitive)")
    p.add_argument("--label", help="Label for the new region annotation (optional)")
    p.add_argument("--type", default="misc_feature", help="Feature type for annotation (default: misc_feature)")
    p.add_argument("--no-mark", action="store_true", help="Skip auto-annotation marking the edit")
    p.set_defaults(func=cmd_replace)

    # export
    p = subs.add_parser("export", help="Convert between sequence file formats (fasta, genbank/gb, tab)")
    p.add_argument("file", help="Input GenBank or FASTA file")
    p.add_argument("format", help="Output format: fasta, genbank (or gb), tab")
    p.add_argument("--out", help="Output file path (default: auto-generated from input filename)")
    p.set_defaults(func=cmd_export)

    # open
    p = subs.add_parser("open", help="Open file in default viewer")
    p.add_argument("file", help="File to open")
    p.set_defaults(func=cmd_open)

    # blast
    p = subs.add_parser("blast", help="Run a remote NCBI BLAST search")
    p.add_argument("file", help="GenBank or FASTA file")
    p.add_argument("--db", default="nt", help="BLAST database (default: nt). Options: nt, nr, refseq_rna, refseq_select_rna")
    p.add_argument("--program", default="blastn", help="BLAST program (default: blastn). Options: blastn, blastp, blastx, tblastn, tblastx")
    p.add_argument("--evalue", type=float, default=0.01, help="E-value threshold (default: 0.01)")
    p.add_argument("--max-hits", type=int, default=10, dest="max_hits", help="Max number of hits to return (default: 10)")
    p.add_argument("--region", type=int, nargs=2, metavar=("START", "END"), help="Blast only this subsequence (1-based, inclusive)")
    p.set_defaults(func=cmd_blast)

    # exons
    p = subs.add_parser("exons", help="Find exon boundaries by aligning an mRNA transcript against a genomic sequence")
    p.add_argument("file", help="Genomic sequence file (.gb, .fasta, .dna)")
    p.add_argument("--transcript", "-t", required=True, help="NCBI accession for the mRNA transcript (e.g. NM_004992.4)")
    p.add_argument("--annotate", "-a", action="store_true", help="Add exon annotations to the file")
    p.add_argument("--email", default="user@example.com", help="Email for NCBI Entrez (default: user@example.com)")
    p.set_defaults(func=cmd_exons)

    p = subs.add_parser("stitch", help="Extract regions by annotation, stitch together, optionally translate")
    p.add_argument("file", help="Sequence file (.gb, .fasta, .dna)")
    p.add_argument("labels", nargs="*", help="Feature labels to extract (in genomic order)")
    p.add_argument("--exons", "-e", action="store_true", help="Auto-find all exon-type annotations")
    p.add_argument("--translate", "-t", action="store_true", help="Translate stitched sequence to protein")
    p.set_defaults(func=cmd_stitch)

    # gibson
    p = subs.add_parser("gibson", help="Design Gibson assembly eBlocks")
    p.add_argument("file", help="Backbone GenBank file")
    p.add_argument("--enzymes", required=True, help="Comma-separated enzyme pair")
    p.add_argument("--insert", help="Raw insert DNA sequence")
    p.add_argument("--insert-from", help="GenBank file to extract features from")
    p.add_argument("--features", help="Comma-separated feature labels to stitch")
    p.add_argument("--preserve-sites", action="store_true", help="Include restriction sites in overlaps")
    p.add_argument("--tm-target", type=float, default=60.0, help="Target Tm for overlaps")
    p.add_argument("--label", default="insert", help="Name for output files")
    p.add_argument("--out-dir", help="Output directory (default: same as backbone)")
    p.set_defaults(func=cmd_gibson)

    # check
    p = subs.add_parser("check", help="Preflight validation of constructs")
    p.add_argument("file", help="GenBank file to check")
    p.add_argument("--backbone", help="Backbone file for Gibson-specific checks")
    p.add_argument("--enzymes", help="Comma-separated enzymes for Gibson checks")
    p.set_defaults(func=cmd_check)

    # branchpoint
    p = subs.add_parser("branchpoint", help="Predict branch point locations in introns using BPP and SVM-BPfinder")
    p.add_argument("file", nargs="?", help="GenBank or FASTA file containing the intron")
    p.add_argument("--intron", help="Feature label of the intron annotation")
    p.add_argument("--region", type=int, nargs=2, metavar=("START", "END"), help="Intron coordinates (1-based, inclusive)")
    p.add_argument("--seq", help="Raw intron DNA sequence (instead of file)")
    p.set_defaults(func=cmd_branchpoint)

    # splice-signals
    p = subs.add_parser("splice-signals", help="Annotate splice signals (5'SS, 3'SS, BPS, PPT) on an intron")
    p.add_argument("file", help="GenBank file")
    p.add_argument("--intron", type=int, nargs=2, required=False, metavar=("START", "END"),
                   help="Intron coordinates (1-based, inclusive). If not provided, introns are auto-detected from exon annotations.")
    p.set_defaults(func=cmd_splice_signals)

    # splicemap
    p_splicemap = subs.add_parser("splicemap",
        help="Map all splicing regulatory elements (5'SS, 3'SS, BPS, PPT, ESE, ESS)")
    p_splicemap.add_argument("file", help="GenBank file")
    p_splicemap.add_argument("-t", "--transcript",
        help="mRNA accession for exon discovery (e.g. NM_004992.4)")
    p_splicemap.add_argument("--no-ese", action="store_true",
        help="Skip ESE/ESS scanning (faster for large files)")
    p_splicemap.add_argument("--clean", action="store_true",
        help="Remove existing splicemap annotations only")
    p_splicemap.set_defaults(func=cmd_splicemap)

    # varmap
    p = subs.add_parser("varmap", help="Visualize variant positions mapped onto a DNA sequence")
    p.add_argument("file", help="Sequence file (.gb, .fasta, .dna)")
    p.add_argument("variants_csv", help="CSV file with variant positions (must have 'Reverted Positions' column)")
    p.add_argument("--region", type=int, nargs=2, metavar=("START", "END"), help="Display region (1-based, inclusive). Auto-detected if omitted.")
    p.add_argument("--intron", help="Intron feature label — positions in CSV are intron-relative")
    p.set_defaults(func=cmd_varmap)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
