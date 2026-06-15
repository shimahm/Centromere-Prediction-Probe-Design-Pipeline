# Centromere Prediction & Probe Design Pipeline
**Species:** *Brassica rapa* (BrapaRo18V23) and *Sinapis alba* (S.alba_v1)  
**Author:** Mahmoudi Lab  
**Date:** June 2026  
**Tools:** EDTA → CentIER → Probe Design (Primer3)

---

## Overview

```
Genome FASTA
     │
     ▼
[Step 1] EDTA — Transposable Element Annotation
     │
     ├──► TE library (.TElib.fa)
     ├──► Intact TEs (.intact.fa / .intact.gff3)
     └──► Full annotation (.TEanno.gff3)
          │
          ▼
[Step 2] CentIER — Centromere Location Prediction
     │
     └──► Centromere coordinates per chromosome
          │
          ▼
[Step 3] Probe Design — Centromere-Specific FISH Probes
     │
     └──► Primer pairs for FISH labeling
```

---

## Environment Setup

```bash
# Create and activate the CentIER conda environment
conda create -n centier -c bioconda -c conda-forge ltr_retriever=3.0.0
conda activate centier

# Install Python dependencies
pip install pyfastx==2.1.0
pip install pandas matplotlib scipy numpy

# Install LTR_Finder
conda install -c bioconda ltr_finder -y

# Compile GenomeTools 1.6.5 (without cairo to avoid X11 dependency)
wget https://github.com/genometools/genometools/archive/refs/tags/v1.6.5.tar.gz
tar -zxf v1.6.5.tar.gz
cd genometools-1.6.5
make clean
make -j4 cairo=no
export PATH=$PATH:$(pwd)/bin
cd ..

# Symlink gt into conda env so LTR_retriever can find it
GT_BIN=~/path/to/genometools-1.6.5/bin
ln -sf $GT_BIN/gt /home/$USER/.conda/envs/centier/bin/gt
ln -sf $GT_BIN/gt /home/$USER/.conda/envs/centier/bin/ltr_harvest

# Set GTDATA so the symlinked gt binary finds its data directory
export GTDATA=~/path/to/genometools-1.6.5/gtdata

# Clone CentIER
git clone https://github.com/simon19891216/CentIER.git
export PATH=$PATH:$(pwd)/CentIER

# Verify all dependencies
which gt && gt --version
which ltr_harvest
which ltr_finder
which LTR_retriever
python -c "import pyfastx; print('pyfastx ok')"
python CentIER/centIER.py -h
```

> **Note:** Add the `export GTDATA` and `export PATH` lines to your `~/.bashrc` or
> to `/home/$USER/.conda/envs/centier/etc/conda/activate.d/env_vars.sh`
> so they persist across sessions.

---

## Step 1 — EDTA: Transposable Element Annotation

### 1.1 Run EDTA

```bash
#!/bin/bash
source /opt/anaconda3/etc/profile.d/conda.sh
conda activate /home/mahmoudi/.conda/envs/EDTA_clean

WORKDIR=/tmp/$USER/edta_run
GENOME_DIR=/home/mahmoudi/Bigdata/computing/Shima/hi_fi_seq
GENOME=BrapaRo18V23_sequences.fa

rm -rf $WORKDIR
mkdir -p $WORKDIR/tmp
export TMPDIR=$WORKDIR/tmp
cp $GENOME_DIR/$GENOME $WORKDIR/
cd $WORKDIR

EDTA.pl \
  --genome $GENOME \
  --species others \
  --step all \
  --sensitive 1 \
  --anno 1 \
  --threads 30 \
  --overwrite 1
```

### 1.2 Key Output Files

| File | Purpose |
|---|---|
| `*.mod.EDTA.TElib.fa` | All TE consensus sequences |
| `*.mod.EDTA.intact.fa` | Full-length intact TEs (best for probes) |
| `*.mod.EDTA.intact.gff3` | Positions of intact elements |
| `*.mod.EDTA.TEanno.gff3` | All TE annotations genome-wide |
| `*.mod.EDTA.TEanno.sum` | TE summary statistics |

### 1.3 Survey TE Composition

```bash
# Check TE superfamily breakdown
grep ">" *.mod.EDTA.TElib.fa | \
    awk -F'#' '{print $2}' | sort | uniq -c | sort -rn | head -30
```

---

## Step 2 — CentIER: Centromere Location Prediction

### 2.1 Important: Check Chromosome ID Format

CentIER requires chromosome IDs in the format `>Chr1`, `>Chr2`, etc.
Check your genome headers first:

```bash
grep ">" BrapaRo18V23_sequences.fa | head -20
```

**If IDs are in a different format** (e.g. `>A01 dna:primary_assembly ...`), rename them:

```bash
# Example: rename >A01, >A02 ... to >Chr1, >Chr2 ...
# Also strips long descriptions after the space
awk '/^>/ {
    split($1, a, "A");
    num = sprintf("%d", a[2]);
    print ">Chr" num
} !/^>/ {print}' BrapaRo18V23_sequences.fa > BrapaRo18V23_ChrFixed.fa

# Verify
grep ">" BrapaRo18V23_ChrFixed.fa | head -10
grep -c ">" BrapaRo18V23_sequences.fa
grep -c ">" BrapaRo18V23_ChrFixed.fa   # must match original count
```

> **Note for BrapaRo18V23:** The gene GFF also uses `A01`-`A10` IDs,
> so the IDs already match between genome and GFF — no renaming needed
> if you pass the GFF. Only rename if running without a GFF.

### 2.2 Run CentIER

```bash
conda activate /home/mahmoudi/.conda/envs/centier

# Set required environment variables every session
export GTDATA=~/Bigdata/computing/Shima/hi_fi_seq/cent_detect/genometools-1.6.5/gtdata
export PATH=~/Bigdata/computing/Shima/hi_fi_seq/cent_detect/genometools-1.6.5/bin:$PATH

# Go to genome directory (CentIER writes temp files relative to working directory)
cd /home/mahmoudi/Bigdata/computing/Shima/hi_fi_seq/cent_detect

GENOME_DIR="/home/mahmoudi/Bigdata/computing/Shima/hi_fi_seq"
GENOME="BrapaRo18V23_sequences.fa"
GFF="/home/mahmoudi/Bigdata/computing/Shima/hi_fi_seq/2.Centromere_position_prediction/BrapaRo18V23_genes.gff"

# Run with nohup (job takes 1-3 hours on a 300 MB genome)
nohup python ~/Bigdata/computing/Shima/hi_fi_seq/cent_detect/CentIER/centIER.py \
  $GENOME_DIR/$GENOME \
  --gff $GFF \
  --step_len 200000 \
  -o CentIER_output_v2/ \
  > centier_run_final.log 2>&1 &

echo "PID: $!"
tail -f centier_run_final.log
```

### 2.3 CentIER Parameters Explained

| Parameter | Value used | Description |
|---|---|---|
| `genome` | `BrapaRo18V23_sequences.fa` | Input genome FASTA (required) |
| `--gff` | `BrapaRo18V23_genes.gff` | Gene annotation — helps identify gene-poor centromeric regions |
| `--step_len` | `200000` | Window step size (200 kb, same as used by colleague on Chiifu_V41) |
| `--mul_cents` | *(not used)* | Only add if species has metapolycentric chromosomes |
| `--SIGNAL_THRESHOLD` | *(default 0.7)* | Lower to 0.5 if centromeres are missed |
| `--matrix1/2` | *(not used)* | Optional Hi-C matrices for higher accuracy |

### 2.4 Expected Output Files

```
CentIER_output_v2/
├── *.fasta_centromere_range.txt       ← centromere coordinates per chromosome
├── *.fasta_all_centromere_seq.txt     ← centromere sequences
├── *.fasta_monomer_seq.txt            ← tandem repeat monomer sequences
├── *.fasta_monomer_in_centromere.txt  ← monomers within centromeres
├── *.fasta_ltr_position.txt           ← LTR positions
├── *.fasta_LTR_statistics.txt         ← LTR statistics
├── *.fasta_draw_cen.svg               ← visualization plot
└── *.fastaLTR-hmmresult.txt           ← HMM results
```

### 2.5 Temp Files Created During Run (Safe to Delete After)

```bash
# Remove temp files after successful run
rm -f BrapaRo18V23_sequences.fa.2.5.7.80.10.50.2000.dat
rm -f BrapaRo18V23_sequences.fa.2.5.7.80.10.50.2000.mask
rm -f BrapaRo18V23_sequences.fa_LTR_database.*
```

---

## Step 3 — Centromere-Specific Probe Design (S. alba example)

This step uses EDTA output to design FISH probes targeting centromeric
transposable elements. Demonstrated on *S. alba* scaffold-level assembly.

### 3.1 Survey and Extract Gypsy Elements

```bash
cd /path/to/edta_S.alba_v1_genomic

# Survey TE composition
grep ">" S.alba_v1_genomic.fna.mod.EDTA.TElib.fa | \
    awk -F'#' '{print $2}' | sort | uniq -c | sort -rn | head -30

# Extract Gypsy, Copia, and unknown LTR families
python3 << 'EOF'
from Bio import SeqIO

gypsy, copia, unk = [], [], []
for rec in SeqIO.parse("S.alba_v1_genomic.fna.mod.EDTA.TElib.fa", "fasta"):
    classification = rec.id.split("#")[-1] if "#" in rec.id else ""
    if "LTR/Gypsy" in classification:
        gypsy.append(rec)
    elif "LTR/Copia" in classification:
        copia.append(rec)
    elif "LTR/unknown" in classification:
        unk.append(rec)

SeqIO.write(gypsy, "Gypsy_lib.fa", "fasta")
SeqIO.write(unk,   "LTR_unknown_lib.fa", "fasta")
print(f"Gypsy:   {len(gypsy)}")
print(f"LTR/unk: {len(unk)}")
print(f"Copia:   {len(copia)}")

sizes = sorted([len(r.seq) for r in gypsy], reverse=True)
print(f"\nGypsy size range: {min(sizes)} - {max(sizes)} bp")
print(f"Gypsy >4000 bp (likely full-length): {sum(1 for s in sizes if s > 4000)}")
EOF
```

### 3.2 Find Centromeric Scaffolds (Gypsy-Enriched)

```bash
# Extract Gypsy annotations from intact GFF3
grep "LTR/Gypsy\|LTR_retrotransposon.*Gypsy\|gypsy\|Gypsy" \
    S.alba_v1_genomic.fna.mod.EDTA.intact.gff3 > Gypsy_intact.gff3

# Count Gypsy hits per scaffold (intact elements)
awk '{print $1}' Gypsy_intact.gff3 | sort | uniq -c | sort -rn | head -20

# Also count from full annotation (all copies including fragments)
grep -i "gypsy\|LTR/Gypsy" S.alba_v1_genomic.fna.mod.EDTA.TEanno.gff3 | \
    awk '{print $1}' | sort | uniq -c | sort -rn | head -20
```

### 3.3 Extract Intact Gypsy from Top Centromeric Scaffolds

The intact.fa header format is: `TE_00000012|MU105553.1:247582..254509#LTR/Gypsy`

```bash
python3 << 'EOF'
from Bio import SeqIO
from collections import Counter

# Top centromeric scaffolds identified from Gypsy density analysis
target = {
    "MU106082.1",   # 741 Gypsy hits — TOP
    "MU105712.1",   # 691 hits
    "MU105568.1",   # 616 hits
    "MU105933.1",   # 606 hits
    "MU105768.1",   # 559 hits
    "MU105667.1",   # 550 hits
    "MU105889.1",   # 505 hits
    "MU105651.1",   # 414 hits
    "MU105744.1",   # 411 hits
    "MU105578.1",   # 396 hits
}

intact_gypsy = []
all_gypsy    = []
skipped      = 0

for rec in SeqIO.parse("S.alba_v1_genomic.fna.mod.EDTA.intact.fa", "fasta"):
    rid = rec.id
    if "|" not in rid or "#" not in rid:
        skipped += 1
        continue
    try:
        scaffold       = rid.split("|")[1].split(":")[0]
        coords_class   = rid.split("|")[1].split(":")[1]
        coords         = coords_class.split("#")[0]
        classification = coords_class.split("#")[1]
        start, end     = coords.split("..")
        length         = abs(int(end) - int(start)) + 1
    except (IndexError, ValueError):
        skipped += 1
        continue

    if "Gypsy" in classification:
        all_gypsy.append((rec, scaffold, int(start), int(end), length))
        if scaffold in target:
            intact_gypsy.append((rec, scaffold, int(start), int(end), length))

print(f"Skipped malformed headers:             {skipped}")
print(f"Total intact Gypsy in genome:          {len(all_gypsy)}")
print(f"Intact Gypsy on centromeric scaffolds: {len(intact_gypsy)}")

SeqIO.write([r[0] for r in intact_gypsy], "Gypsy_centromeric_intact.fa",  "fasta")
SeqIO.write([r[0] for r in all_gypsy],    "Gypsy_all_intact.fa", "fasta")

if intact_gypsy:
    sizes = [x[4] for x in intact_gypsy]
    print(f"\nSize range:             {min(sizes):,} - {max(sizes):,} bp")
    print(f"Full-length (>5000 bp): {sum(1 for s in sizes if s > 5000)}")

    print("\n=== Intact Gypsy per centromeric scaffold ===")
    scaf_counts = Counter(x[1] for x in intact_gypsy)
    for scaf, cnt in scaf_counts.most_common():
        avg = sum(x[4] for x in intact_gypsy if x[1]==scaf) // cnt
        print(f"  {scaf:<15}  {cnt:>3} elements   avg size {avg:,} bp")
EOF
```

### 3.4 Validate Centromeric Specificity by BLAST

```bash
# Build genome BLAST database
makeblastdb -in S.alba_v1_genomic.fna -dbtype nucl -out genome_db -parse_seqids

# BLAST all intact Gypsy against genome
blastn \
    -query Gypsy_centromeric_intact.fa \
    -db genome_db \
    -perc_identity 75 \
    -outfmt "6 qseqid sseqid pident length qlen evalue" \
    -num_threads 8 \
    -out probe_vs_genome.blast

# Count hits per scaffold
python3 << 'EOF'
from collections import defaultdict

centromeric = {
    "MU106082.1","MU105712.1","MU105568.1","MU105933.1",
    "MU105768.1","MU105667.1","MU105889.1","MU105651.1",
    "MU105744.1","MU105578.1"
}

hits = defaultdict(set)
with open("probe_vs_genome.blast") as f:
    for line in f:
        cols = line.strip().split("\t")
        te   = cols[0].split("|")[0]
        scaf = cols[1]
        hits[te].add(scaf)

print(f"{'TE Family':<20} {'Total_scaffolds':>16} {'Centromeric':>12} {'Specificity%':>13}")
print("-" * 65)
for te, scaffolds in sorted(hits.items(), key=lambda x: len(x[1]), reverse=True):
    cen   = len(scaffolds & centromeric)
    total = len(scaffolds)
    spec  = cen / total * 100
    flag  = " ✓ SPECIFIC" if spec >= 70 else ""
    print(f"{te:<20} {total:>16} {cen:>12}    {spec:>6.1f}%{flag}")
EOF
```

### 3.5 Extract LTR Regions from Best Candidate (TE_00001488)

TE_00001488 was identified as the best probe candidate: highest centromeric
enrichment (16.1% of 22,681 hits on 10 centromeric scaffolds, with those
scaffolds occupying 6 of the top 8 positions by hit count).

```bash
# Extract LTR coordinates for TE_00001488 from intact GFF3
python3 << 'EOF'
target_te   = "TE_00001488"
ltr_entries = []

with open("S.alba_v1_genomic.fna.mod.EDTA.intact.gff3") as f:
    for line in f:
        if line.startswith("#"):
            continue
        cols = line.strip().split("\t")
        if len(cols) < 9:
            continue
        if cols[2] == "long_terminal_repeat" and target_te in cols[8]:
            scaffold = cols[0]
            start    = int(cols[3])
            end      = int(cols[4])
            length   = end - start + 1
            ltr_entries.append((scaffold, start, end, length))
            print(f"{scaffold}  {start:>10}-{end:>10}  len={length:>5}")

print(f"\nTotal LTR regions for {target_te}: {len(ltr_entries)}")

with open("TE_00001488_LTRs.bed", "w") as out:
    for i, (scaf, s, e, l) in enumerate(ltr_entries):
        out.write(f"{scaf}\t{s-1}\t{e}\tLTR_{i+1}\t.\t+\n")
print("Written: TE_00001488_LTRs.bed")
EOF

# Extract LTR sequences with bedtools
bedtools getfasta \
    -fi S.alba_v1_genomic.fna \
    -bed TE_00001488_LTRs.bed \
    -name \
    -fo TE_00001488_LTRs.fa

echo "LTR sequences extracted: $(grep -c '>' TE_00001488_LTRs.fa)"
```

### 3.6 Multiple Alignment and Primer Design

```bash
# Align all LTR sequences
mafft --auto --thread 8 TE_00001488_LTRs.fa > TE_00001488_LTRs_aligned.fa

# Build consensus and run Primer3
python3 << 'PYEOF'
from Bio import AlignIO, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import subprocess

# Load alignment
aln     = AlignIO.read("TE_00001488_LTRs_aligned.fa", "fasta")
aln_len = aln.get_alignment_length()
n_seq   = len(aln)
print(f"Alignment: {n_seq} seqs x {aln_len} columns")

# Per-position conservation
conserved = []
for i in range(aln_len):
    col = [str(aln[j].seq[i]).upper() for j in range(n_seq)]
    non_gap = [b for b in col if b != '-']
    if len(non_gap) < n_seq * 0.75:
        conserved.append(0.0)
    else:
        top = max(set(non_gap), key=non_gap.count)
        conserved.append(non_gap.count(top) / len(non_gap))

# Find best conserved window
for threshold in [0.85, 0.80, 0.75, 0.70]:
    windows = []
    for start in range(0, aln_len - 200, 10):
        chunk = conserved[start:start+200]
        avg   = sum(chunk) / len(chunk)
        if avg >= threshold:
            windows.append((start, start+200, avg))
    if windows:
        print(f"Threshold {threshold:.0%}: {len(windows)} windows found")
        best_start, best_end, best_avg = sorted(windows, key=lambda x: x[2], reverse=True)[0]
        print(f"Best window: pos {best_start}-{best_end}  avg={best_avg:.1%}")
        break

# Build full LTR consensus for Primer3
full_consensus = []
for i in range(aln_len):
    col = [str(aln[j].seq[i]).upper() for j in range(n_seq)
           if str(aln[j].seq[i]) != '-']
    if col:
        full_consensus.append(max(set(col), key=col.count))
full_seq = ''.join(full_consensus)

with open("LTR_consensus_full.fa", "w") as f:
    f.write(f">LTR_consensus_full_TE00001488\n{full_seq}\n")
print(f"Full LTR consensus: {len(full_seq)} bp  →  Written: LTR_consensus_full.fa")

# Primer3 input (relaxed settings for repetitive templates)
p3_input = f"""SEQUENCE_ID=TE_00001488_LTR_centromere
SEQUENCE_TEMPLATE={full_seq}
PRIMER_TASK=generic
PRIMER_PICK_LEFT_PRIMER=1
PRIMER_PICK_RIGHT_PRIMER=1
PRIMER_OPT_SIZE=22
PRIMER_MIN_SIZE=18
PRIMER_MAX_SIZE=28
PRIMER_OPT_TM=58.0
PRIMER_MIN_TM=54.0
PRIMER_MAX_TM=65.0
PRIMER_MIN_GC=35.0
PRIMER_MAX_GC=70.0
PRIMER_PRODUCT_SIZE_RANGE=150-600
PRIMER_NUM_RETURN=5
PRIMER_MAX_SELF_ANY_TH=47.0
PRIMER_MAX_SELF_END_TH=40.0
PRIMER_MAX_HAIRPIN_TH=30.0
PRIMER_MAX_POLY_X=5
PRIMER_EXPLAIN_FLAG=1
=
"""
with open("primer3_input.txt", "w") as f:
    f.write(p3_input)

result = subprocess.run(
    ["primer3_core", "primer3_input.txt"],
    capture_output=True, text=True
)
with open("primer3_output.txt", "w") as f:
    f.write(result.stdout)

# Parse results
data = {}
for line in result.stdout.split('\n'):
    if '=' in line:
        k, _, v = line.partition('=')
        data[k.strip()] = v.strip()

n_pairs = int(data.get('PRIMER_PAIR_NUM_RETURNED', 0))
print(f"\nPrimer pairs returned: {n_pairs}")

rows = []
for i in range(n_pairs):
    lf   = data.get(f'PRIMER_LEFT_{i}_SEQUENCE',    'N/A')
    rf   = data.get(f'PRIMER_RIGHT_{i}_SEQUENCE',   'N/A')
    tm_l = data.get(f'PRIMER_LEFT_{i}_TM',          'N/A')
    tm_r = data.get(f'PRIMER_RIGHT_{i}_TM',         'N/A')
    gc_l = data.get(f'PRIMER_LEFT_{i}_GC_PERCENT',  'N/A')
    gc_r = data.get(f'PRIMER_RIGHT_{i}_GC_PERCENT', 'N/A')
    size = data.get(f'PRIMER_PAIR_{i}_PRODUCT_SIZE', 'N/A')
    lpos = data.get(f'PRIMER_LEFT_{i}',              'N/A')
    rpos = data.get(f'PRIMER_RIGHT_{i}',             'N/A')
    print(f"\nPair {i+1}:")
    print(f"  F: 5'-{lf}-3'  pos={lpos}  Tm={float(tm_l):.1f}°C  GC={float(gc_l):.1f}%")
    print(f"  R: 5'-{rf}-3'  pos={rpos}  Tm={float(tm_r):.1f}°C  GC={float(gc_r):.1f}%")
    print(f"  Product: {size} bp")
    rows.append([str(i+1), lf, rf, size,
                 f"{float(tm_l):.1f}", f"{float(tm_r):.1f}",
                 f"{float(gc_l):.1f}", f"{float(gc_r):.1f}"])

with open("centromere_probe_primers_TE00001488.tsv", "w") as out:
    out.write("Pair\tForward\tReverse\tProduct_bp\tTm_F\tTm_R\tGC_F\tGC_R\n")
    for row in rows:
        out.write("\t".join(row) + "\n")
print("\nSaved: centromere_probe_primers_TE00001488.tsv")
PYEOF
```

---

## Final Output Summary

### CentIER outputs (Step 2)
| File | Contents |
|---|---|
| `*_centromere_range.txt` | Predicted centromere start/end per chromosome |
| `*_monomer_seq.txt` | Tandem repeat monomer sequences |
| `*_draw_cen.svg` | Visualization of k-mer signal + centromere calls |
| `*_LTR_statistics.txt` | LTR retrotransposon stats within centromeres |

### Probe design outputs (Step 3)
| File | Contents |
|---|---|
| `Gypsy_centromeric_intact.fa` | Intact Gypsy elements from centromeric scaffolds |
| `TE_00001488_LTRs.fa` | LTR sequences of best probe candidate |
| `LTR_consensus_full.fa` | Full-length LTR consensus for primer design |
| `centromere_probe_primers_TE00001488.tsv` | **Final primer pairs for FISH** |
| `probe_vs_genome.blast` | Specificity validation |

---

## Probe Selection Guidelines

| TE Family | Genome hits | % on centromeric scaffolds | Recommendation |
|---|---|---|---|
| TE_00001488 | 22,681 | 15.7% (top 6 of top 8) | ✅ PRIMARY PROBE |
| TE_00000525 | 8,426 | 9.2% | ✅ Secondary validation |
| TE_00002172 | 3,519 | 9.9% | ✅ Secondary validation |
| TE_00001972 | 3,410 | 6.2% | Moderate |
| TE_00001898 | 161 | 16.1% | ⚠ Too low copy for FISH |

> **Interpretation:** True centromeric probes in plants are pericentromeric
> (enriched but not exclusive). The pattern of top centromeric scaffolds
> dominating the hit list is the expected CRM-Gypsy signature.

---

## Experimental Follow-Up

Two recommended approaches for wet-lab probe production:

| Approach | Method | Signal type |
|---|---|---|
| PCR + nick translation | Amplify Pair 1 from genomic DNA, label with Dig or Biotin | Indirect FISH |
| Direct oligo | Synthesize the amplicon region as ATTO-488 or Cy3 labeled oligo | Direct FISH |

Expected FISH result: bright pericentromeric signals on most/all chromosomes,
consistent with Gypsy-CRM pericentromeric repeat patterns seen across Brassicaceae.

---

## Common Errors and Fixes

| Error | Cause | Fix |
|---|---|---|
| `cairo.h: No such file` | GenomeTools needs cairo | Compile with `make cairo=no` |
| `could not find gtdata/` | Symlinked `gt` can't find data dir | `export GTDATA=/path/to/genometools-1.6.5/gtdata` |
| `finder.scn does not exist` | `ltr_finder` not installed | `conda install -c bioconda ltr_finder` |
| `harvest.scn` is 0 KB | `gt ltrharvest` failing | Fix `GTDATA` export, recompile GenomeTools |
| `RuntimeError: / is not fasta` | Variables `$GENOME_DIR/$GENOME` not set | Re-export variables in current shell session |
| Output folder empty | Chr IDs don't match required format | Check IDs with `grep ">" genome.fa \| head` |
| `ModuleNotFoundError: pandas` | Python deps missing | `pip install pandas matplotlib scipy numpy` |
