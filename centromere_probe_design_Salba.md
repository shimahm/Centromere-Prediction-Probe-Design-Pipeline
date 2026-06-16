# Centromere-Specific Probe Design from EDTA Output
### *Sinapis alba* v1 genome — LTR/Gypsy (TE_00001488) FISH probe pipeline

---

## Overview

This pipeline identifies centromeric Gypsy LTR retrotransposons from EDTA repeat
annotation output and designs PCR primers for centromere-specific FISH probes.

**Input:** EDTA output directory from a scaffold-level genome assembly  
**Output:** PCR primer pairs targeting a high-copy centromeric Gypsy LTR family  
**Species:** *Sinapis alba* v1 (scaffold assembly, 12,410 sequences)  
**Environment:** Linux HPC, conda environment `centier`

---

## Dependencies

```bash
conda activate centier
conda install -c bioconda mafft primer3 bedtools cd-hit -y
pip install biopython
```

---

## Directory Structure

```
edta_S.alba_v1_genomic/
├── S.alba_v1_genomic.fna                          # Reference genome
├── S.alba_v1_genomic.fna.mod.EDTA.TElib.fa        # TE consensus library
├── S.alba_v1_genomic.fna.mod.EDTA.intact.fa       # Intact TE sequences
├── S.alba_v1_genomic.fna.mod.EDTA.intact.gff3     # Intact TE annotations
└── S.alba_v1_genomic.fna.mod.EDTA.TEanno.gff3    # All TE annotations
```

---

## Step 1 — Survey TE Composition

Examine the TE superfamily breakdown to identify centromere-relevant families.
In plants, **LTR/Gypsy** elements (especially CRM-clade) are the primary
centromeric repeat marker.

```bash
grep ">" S.alba_v1_genomic.fna.mod.EDTA.TElib.fa | \
    awk -F'#' '{print $2}' | sort | uniq -c | sort -rn | head -30
```

**Results for *S. alba*:**
```
1186 LTR/Copia
 579 LTR/Gypsy     ← target
 364 unknown
 351 DNA/DTH
 168 LTR/unknown
 ...
```

---

## Step 2 — Extract Gypsy Sequences from TElib

Parse the TE library by classification. The FASTA header format is:
`TE_XXXXX#LTR/Gypsy`

```bash
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

print(f"Gypsy:      {len(gypsy)}")
print(f"LTR/unk:    {len(unk)}")
print(f"Copia:      {len(copia)}")

sizes = sorted([len(r.seq) for r in gypsy], reverse=True)
print(f"\nGypsy size range: {min(sizes)} - {max(sizes)} bp")
print(f"Gypsy >4000 bp (likely full-length): {sum(1 for s in sizes if s > 4000)}")
EOF
```

**Results:** 579 Gypsy sequences, size range 80–11,590 bp, 170 full-length (>4 kb)

---

## Step 3 — Identify Centromeric Scaffolds by Gypsy Density

Centromeric scaffolds in scaffold-level assemblies are identified by high Gypsy
element density (copies per scaffold). Two complementary counts are used.

```bash
# Count intact full-length Gypsy per scaffold (from intact GFF3)
grep "LTR/Gypsy\|LTR_retrotransposon.*Gypsy\|gypsy\|Gypsy" \
    S.alba_v1_genomic.fna.mod.EDTA.intact.gff3 > Gypsy_intact.gff3

awk '{print $1}' Gypsy_intact.gff3 | sort | uniq -c | sort -rn | head -20

# Count all Gypsy hits per scaffold (including fragments, from TEanno)
grep -i "gypsy\|LTR/Gypsy" S.alba_v1_genomic.fna.mod.EDTA.TEanno.gff3 | \
    awk '{print $1}' | sort | uniq -c | sort -rn | head -20
```

**Top centromeric scaffolds identified (by TEanno Gypsy hit count):**

| Scaffold | Gypsy hits | Intact Gypsy | Avg element size |
|---|---|---|---|
| MU106082.1 | 741 | 12 | 8,596 bp |
| MU105712.1 | 691 | 9 | 8,205 bp |
| MU105568.1 | 616 | 4 | 7,570 bp |
| MU105933.1 | 606 | 8 | 7,365 bp |
| MU105768.1 | 559 | 4 | 9,855 bp |
| MU105667.1 | 550 | 5 | 11,120 bp |
| MU105889.1 | 505 | 5 | 8,471 bp |
| MU105651.1 | 414 | 8 | 7,090 bp |
| MU105744.1 | 411 | 4 | 7,161 bp |
| MU105578.1 | 396 | 7 | 6,724 bp |

---

## Step 4 — Extract Intact Gypsy Elements from Centromeric Scaffolds

The intact.fa header format is: `TE_XXXXX|scaffold:start..end#LTR/Gypsy`

```bash
python3 << 'EOF'
from Bio import SeqIO
from collections import Counter

target = {
    "MU106082.1", "MU105712.1", "MU105568.1", "MU105933.1",
    "MU105768.1", "MU105667.1", "MU105889.1", "MU105651.1",
    "MU105744.1", "MU105578.1",
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
        # Format: TE_00000012|MU105553.1:247582..254509#LTR/Gypsy
        right          = rid.split("|")[1]
        scaffold       = right.split(":")[0]
        coords_class   = right.split(":")[1]
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

print(f"Skipped malformed headers:              {skipped}")
print(f"Total intact Gypsy in genome:           {len(all_gypsy)}")
print(f"Intact Gypsy on centromeric scaffolds:  {len(intact_gypsy)}")

SeqIO.write([r[0] for r in intact_gypsy], "Gypsy_centromeric_intact.fa", "fasta")
SeqIO.write([r[0] for r in all_gypsy],    "Gypsy_all_intact.fa", "fasta")

sizes = [x[4] for x in intact_gypsy]
print(f"\nSize range:              {min(sizes):,} - {max(sizes):,} bp")
print(f"Full-length (>5000 bp):  {sum(1 for s in sizes if s > 5000)}")
print(f"Partial (2000-5000 bp):  {sum(1 for s in sizes if 2000 <= s <= 5000)}")
print(f"Fragment (<2000 bp):     {sum(1 for s in sizes if s < 2000)}")

print("\n=== Intact Gypsy per centromeric scaffold ===")
scaf_counts = Counter(x[1] for x in intact_gypsy)
for scaf, cnt in scaf_counts.most_common():
    avg = sum(x[4] for x in intact_gypsy if x[1]==scaf) // cnt
    print(f"  {scaf:<15}  {cnt:>3} elements   avg size {avg:,} bp")
EOF
```

**Results:** 66 intact Gypsy elements on 10 centromeric scaffolds (56 full-length >5 kb)

---

## Step 5 — Select Priority Families (Multi-Scaffold Elements)

Gypsy families present on multiple centromeric scaffolds are the most reliable
centromeric markers. From inspection of the element table:

| TE Family | Copies on cen-scaffolds | Scaffolds | Priority |
|---|---|---|---|
| TE_00000525 | 5 | MU105667, MU105712×2, MU105889×2 | ⭐⭐⭐ |
| TE_00001972 | 3 | MU105933, MU106082, MU105651 | ⭐⭐⭐ |
| TE_00002172 | 3 | MU105744, MU105889, MU105578 | ⭐⭐⭐ |
| TE_00001488 | 3 | MU106082, MU105889, MU105651 | ⭐⭐⭐ |
| TE_00000769 | 2 | MU105933, MU105712 | ⭐⭐ |
| TE_00001898 | 2 | MU106082×2 | ⭐⭐ |

```bash
python3 << 'EOF'
from Bio import SeqIO

priority = {
    "TE_00000525", "TE_00001488",
    "TE_00001972", "TE_00002172",
    "TE_00000769", "TE_00001898"
}

probe_candidates = []
for rec in SeqIO.parse("Gypsy_centromeric_intact.fa", "fasta"):
    te_id = rec.id.split("|")[0]
    if te_id in priority:
        probe_candidates.append(rec)

SeqIO.write(probe_candidates, "Gypsy_probe_candidates.fa", "fasta")
print(f"Probe candidate sequences saved: {len(probe_candidates)}")
for rec in probe_candidates:
    print(f"  {rec.id.split('|')[0]:<20}  {len(rec.seq):,} bp")
EOF
```

---

## Step 6 — Genome-wide BLAST to Assess Copy Number

High copy number = strong FISH signal. Build genome database and BLAST all
candidate families simultaneously.

```bash
makeblastdb -in S.alba_v1_genomic.fna \
            -dbtype nucl \
            -out genome_db \
            -parse_seqids

blastn \
    -query Gypsy_probe_candidates.fa \
    -db genome_db \
    -perc_identity 75 \
    -outfmt "6 qseqid sseqid pident length qlen evalue" \
    -num_threads 8 \
    -out probe_vs_genome.blast

echo "=== Genome-wide hit count per probe family ==="
awk '{print $1}' probe_vs_genome.blast | \
    sed 's/|.*//' | \
    sort | uniq -c | sort -rn
```

**Results:**

| TE Family | Genome-wide BLAST hits | Selection |
|---|---|---|
| TE_00001488 | 22,681 | ✅ **Selected** — highest copy, strong pericentromeric enrichment |
| TE_00000525 | 8,426 | Good |
| TE_00002172 | 3,519 | Moderate |
| TE_00001972 | 3,410 | Moderate |
| TE_00000769 | 2,756 | Moderate |
| TE_00001898 | 161 | Too low for reliable FISH signal |

> **Note:** These Gypsy elements are pericentromeric (distributed across the genome
> but highly enriched on centromeric scaffolds), consistent with CRM-clade Gypsy
> biology in Brassicaceae. Density analysis on centromeric scaffolds confirms
> enrichment even though absolute specificity is not 100%.

---

## Step 7 — Extract LTR Regions of TE_00001488

LTRs (~1,300 bp) are the most conserved part of the element — ideal for
designing primers that amplify across all copies.

```bash
# Parse intact GFF3 for LTR coordinates of TE_00001488
python3 << 'EOF'
target_te = "TE_00001488"
ltr_entries = []

with open("S.alba_v1_genomic.fna.mod.EDTA.intact.gff3") as f:
    for line in f:
        if line.startswith("#"): continue
        cols = line.strip().split("\t")
        if len(cols) < 9: continue
        if cols[2] == "long_terminal_repeat" and target_te in cols[8]:
            scaffold = cols[0]
            start    = int(cols[3])
            end      = int(cols[4])
            length   = end - start + 1
            ltr_entries.append((scaffold, start, end, length, cols[8]))

print(f"Total LTR regions for {target_te}: {len(ltr_entries)}")

with open("TE_00001488_LTRs.bed", "w") as out:
    for i, (scaf, s, e, l, attr) in enumerate(ltr_entries):
        out.write(f"{scaf}\t{s-1}\t{e}\tLTR_{i+1}\t.\t+\n")
print("BED file written: TE_00001488_LTRs.bed")
EOF

# Extract LTR sequences from genome
bedtools getfasta \
    -fi S.alba_v1_genomic.fna \
    -bed TE_00001488_LTRs.bed \
    -name \
    -fo TE_00001488_LTRs.fa

echo "LTR sequences extracted: $(grep -c '>' TE_00001488_LTRs.fa)"
```

**Results:** 16 LTR sequences extracted, length range 1,285–1,345 bp
(8 pairs from 8 intact elements across 7 scaffolds)

---

## Step 8 — Multiple Sequence Alignment of LTRs

```bash
mafft --auto \
      --thread 8 \
      TE_00001488_LTRs.fa > TE_00001488_LTRs_aligned.fa

echo "Alignment done: $(grep -c '>' TE_00001488_LTRs_aligned.fa) sequences"
```

---

## Step 9 — Find Conserved Regions and Build Consensus

Calculate per-position conservation across the alignment to identify the most
conserved region for primer anchoring.

```bash
python3 << 'PYEOF'
from Bio import AlignIO, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# Load alignment
aln     = AlignIO.read("TE_00001488_LTRs_aligned.fa", "fasta")
aln_len = aln.get_alignment_length()
n_seq   = len(aln)
print(f"Alignment: {n_seq} seqs x {aln_len} columns")

# Per-position conservation score
conserved = []
for i in range(aln_len):
    col     = [str(aln[j].seq[i]).upper() for j in range(n_seq)]
    non_gap = [b for b in col if b != '-']
    if len(non_gap) < n_seq * 0.75:
        conserved.append(0.0)
    else:
        top = max(set(non_gap), key=non_gap.count)
        conserved.append(non_gap.count(top) / len(non_gap))

# Find best conserved window (200 bp sliding window)
for threshold in [0.85, 0.80, 0.75, 0.70]:
    windows = []
    for start in range(0, aln_len - 200, 10):
        avg = sum(conserved[start:start+200]) / 200
        if avg >= threshold:
            windows.append((start, start+200, avg))
    if windows:
        print(f"Threshold {threshold:.0%}: {len(windows)} windows found")
        best_start, best_end, best_avg = sorted(windows, key=lambda x: x[2], reverse=True)[0]
        print(f"Best window: pos {best_start}-{best_end}  avg={best_avg:.1%}")
        break

# Build full LTR consensus (used as Primer3 template)
full_consensus = []
for i in range(aln_len):
    col = [str(aln[j].seq[i]).upper() for j in range(n_seq)
           if str(aln[j].seq[i]) != '-']
    if col:
        full_consensus.append(max(set(col), key=col.count))

full_seq = ''.join(full_consensus)
print(f"Full LTR consensus: {len(full_seq)} bp")

with open("LTR_consensus_full.fa", "w") as f:
    f.write(f">LTR_consensus_full_TE00001488\n{full_seq}\n")
print("Written: LTR_consensus_full.fa")
PYEOF
```

**Result:** Full LTR consensus of ~1,350 bp built from 16 copies (75–80% average conservation)

---

## Step 10 — Design Primers with Primer3

```bash
python3 << 'PYEOF'
from Bio import SeqIO
import subprocess

seq      = next(SeqIO.parse("LTR_consensus_full.fa", "fasta"))
template = str(seq.seq).upper()

p3_input = f"""SEQUENCE_ID=TE_00001488_LTR_centromere
SEQUENCE_TEMPLATE={template}
PRIMER_TASK=generic
PRIMER_PICK_LEFT_PRIMER=1
PRIMER_PICK_RIGHT_PRIMER=1
PRIMER_OPT_SIZE=22
PRIMER_MIN_SIZE=18
PRIMER_MAX_SIZE=26
PRIMER_OPT_TM=60.0
PRIMER_MIN_TM=57.0
PRIMER_MAX_TM=63.0
PRIMER_MIN_GC=40.0
PRIMER_MAX_GC=65.0
PRIMER_PRODUCT_SIZE_RANGE=200-500
PRIMER_NUM_RETURN=5
PRIMER_MAX_SELF_ANY_TH=45.0
PRIMER_MAX_SELF_END_TH=35.0
PRIMER_MAX_HAIRPIN_TH=24.0
PRIMER_MAX_POLY_X=4
PRIMER_EXPLAIN_FLAG=1
=
"""

with open("primer3_input.txt", "w") as f:
    f.write(p3_input)

result = subprocess.run(["primer3_core", "primer3_input.txt"],
                        capture_output=True, text=True)

with open("primer3_output.txt", "w") as f:
    f.write(result.stdout)

data = {}
for line in result.stdout.split('\n'):
    if '=' in line:
        k, _, v = line.partition('=')
        data[k.strip()] = v.strip()

n_pairs = int(data.get('PRIMER_PAIR_NUM_RETURNED', 0))
print(f"Primer pairs returned: {n_pairs}")

rows = []
for i in range(n_pairs):
    lf   = data.get(f'PRIMER_LEFT_{i}_SEQUENCE',    'N/A')
    rf   = data.get(f'PRIMER_RIGHT_{i}_SEQUENCE',   'N/A')
    tm_l = data.get(f'PRIMER_LEFT_{i}_TM',          'N/A')
    tm_r = data.get(f'PRIMER_RIGHT_{i}_TM',         'N/A')
    gc_l = data.get(f'PRIMER_LEFT_{i}_GC_PERCENT',  'N/A')
    gc_r = data.get(f'PRIMER_RIGHT_{i}_GC_PERCENT', 'N/A')
    size = data.get(f'PRIMER_PAIR_{i}_PRODUCT_SIZE', 'N/A')
    print(f"\nPair {i+1}:")
    print(f"  F: 5'-{lf}-3'  Tm={float(tm_l):.1f}°C  GC={float(gc_l):.1f}%")
    print(f"  R: 5'-{rf}-3'  Tm={float(tm_r):.1f}°C  GC={float(gc_r):.1f}%")
    print(f"  Product: {size} bp")
    rows.append([str(i+1), lf, rf, size,
                 f"{float(tm_l):.1f}", f"{float(tm_r):.1f}",
                 f"{float(gc_l):.1f}", f"{float(gc_r):.1f}"])

with open("centromere_probe_primers_TE00001488.tsv", "w") as out:
    out.write("Pair\tForward\tReverse\tProduct_bp\tTm_F\tTm_R\tGC_F\tGC_R\n")
    for row in rows:
        out.write("\t".join(row) + "\n")

print(f"\nSaved: centromere_probe_primers_TE00001488.tsv")
PYEOF
```

---

## Final Results — Primer Pairs

| Pair | Forward (5'→3') | Reverse (5'→3') | Product | Tm F/R | GC F/R | Recommendation |
|---|---|---|---|---|---|---|
| **1** | `ACGCCAAGTGTATTCCATCGAT` | `TACCAGTTCTAGGCCAGGGTAA` | 212 bp | 60.2/60.0 | 45.5/50.0 | ✅ **Primary** |
| 2 | `ACGCCAAGTGTATTCCATCGAT` | `TTACCAGTTCTAGGCCAGGGTA` | 213 bp | 60.2/60.0 | 45.5/50.0 | (duplicate of Pair 1) |
| **3** | `ATTTCTAGGGAGTCGACGAAGC` | `GTCTCTCCTTCCGTCTCTCCTA` | 306 bp | 59.9/59.8 | 50.0/54.5 | ✅ **Independent validation** |
| 4 | `ACATTTCCGTGCTTCCTATCCA` | `TACCAGTTCTAGGCCAGGGTAA` | 488 bp | 59.8/60.0 | 45.5/50.0 | Backup |
| 5 | `ACATTTCCGTGCTTCCTATCCA` | `TTACCAGTTCTAGGCCAGGGTA` | 489 bp | 59.8/60.0 | 45.5/50.0 | Backup |

---

## Output Files Summary

| File | Description |
|---|---|
| `Gypsy_lib.fa` | All 579 Gypsy consensus sequences from TElib |
| `Gypsy_intact.gff3` | GFF3 of intact Gypsy elements genome-wide |
| `Gypsy_centromeric_intact.fa` | 66 intact Gypsy elements on 10 centromeric scaffolds |
| `Gypsy_all_intact.fa` | All 631 intact Gypsy elements genome-wide |
| `Gypsy_probe_candidates.fa` | 18 sequences from 6 priority TE families |
| `probe_vs_genome.blast` | BLAST results: probe candidates vs whole genome |
| `TE_00001488_LTRs.bed` | BED coordinates of 16 LTR regions |
| `TE_00001488_LTRs.fa` | 16 extracted LTR sequences (~1,300 bp each) |
| `TE_00001488_LTRs_aligned.fa` | MAFFT alignment of 16 LTRs |
| `LTR_consensus_full.fa` | Full LTR consensus sequence (Primer3 template) |
| `primer3_input.txt` | Primer3 input file |
| `primer3_output.txt` | Full Primer3 output |
| `centromere_probe_primers_TE00001488.tsv` | **Final primer table (TSV)** |

---

## Probe Production for FISH

### Option A — PCR + Nick Translation (indirect FISH)
1. Amplify with Pair 1 or Pair 3 from *S. alba* genomic DNA
2. Label PCR product by nick translation (biotin-16-dUTP or digoxigenin-11-dUTP)
3. Use as FISH probe with anti-Dig or streptavidin-fluorochrome detection

### Option B — Direct labeled oligo (direct FISH)
Order a 40–50 bp oligo spanning the primer binding site with 5' ATTO-488 or Cy3 label.

### Expected FISH signal
- TE_00001488 has **22,681 genome-wide BLAST hits** (≥75% identity)
- Enriched on centromeric scaffolds → expected **pericentromeric dot signals**
  on most/all chromosomes of *S. alba* (2n=24)
- Signal pattern consistent with CRM-clade Gypsy pericentromeric repeats
  described in other Brassicaceae species

---

## Citation / Method Reference

> Repeat annotation was performed with EDTA v2 (Ou et al. 2019).
> Centromeric scaffold identification was based on LTR/Gypsy element density.
> Multiple sequence alignment was performed with MAFFT v7.525 (Katoh & Standley 2013).
> PCR primers were designed with Primer3 v2.6.1 (Untergasser et al. 2012).

---

*Pipeline developed for S. alba v1 genomic assembly (12,410 scaffolds)*  
*Date: June 2026*
