# Prime Editing pegRNA Library Construction

This repository contains a computational pipeline that designs prime-editing guide RNAs (pegRNAs) for comprehensive mutational scanning of a given coding sequence. The pipeline generates pegRNAs capable of mutating each amino acid in the coding sequence into all other codons (exhaustive codon-level enumeration), producing a candidate set that can be scored with PRIDICT/PRIDICT2.0 (https://www.pridict.it/) to select pegRNAs with high predicted editing efficiency. After running PRIDICT on all designed pegRNAs, the best pegRNA is chosen for each of the 20 possible amino acid substitutions at every position. The pipeline has been run and tested on exon 4 of the human low-density lipoprotein receptor (LDLR).

## Design pegRNAs (design_pegrna.py)

Design strategy:

**Exhaustive codon enumeration.** For each amino-acid (AA) position, the code iterates over all pairs of the existing AA and desired mutant codons. For mutant codons, all synonymous, missense, and nonsense variants that differ from the existing codon are considered—up to 63 codon alternatives per AA.

**Nearby PAM discovery.** For each AA position and desired codon change, the code searches for **NGG** PAMs (“GG” in practice) in the vicinity of the codon: up to **17 bp upstream** on the sense strand or **17 bp downstream** on the antisense strand. Each valid PAM defines a candidate spacer for a pegRNA that installs the desired edit.

**PAM/seed disruption and Hamming distance.**
- **Prefer PAM disruption via a synonymous change.** If possible, introduce a synonymous mutation that disrupts at least one G in the NGG PAM to avoid repeated editing of the same locus after installation.
- **If PAM disruption is not possible, disrupt the seed.** Search for synonymous changes within the seed sequence (the three nucleotides between the PAM and the nick site) to prevent re-targeting.
- **Enforce edit identifiability.** To distinguish genuine edits from potential sequencing errors, require a **Hamming distance ≥ 2** between the original and mutated sequences. When the combination of the AA-changing mutation and PAM/seed-disrupting mutation already yields ≥2 nucleotide differences, no further synonymous changes are added.
- **Special case when PAM lies within the edited codon.** If the PAM overlaps the AA being mutated and the AA-changing edit itself disrupts the PAM (e.g., GGC, Glycine → AGC, Glycine), the Hamming distance may be 1. In these cases, add an additional synonymous mutation within 1–3 upstream AAs to reach Hamming distance ≥ 2.

**Multiplicity and amino-acid edge cases.**
- When multiple synonymous options can disrupt the PAM or seed and/or satisfy Hamming ≥ 2, all such options are retained so PRIDICT can later prioritize the most efficient design.
- **Met** and **Trp**, each encoded by a single codon, naturally have no synonymous alternatives; no synonymous mutations are added in those instances.
- If no suitable NGG PAM is found within ±17 bp, the search is extended to ±35 bp for that target.

## PRIDICT post-processing and Final Selection (load_pridict_output.py)

1. **Score all candidates.** All (variant, PAM) designs are formatted and scored with PRIDICT/PRIDICT2.0 (the code is adjusted to incorporate differences in the output between PRIDICT and PRIDICT2.0).
2. **Filter by spacer.** Because PRIDICT may construct pegRNAs using any nearby PAM (i.e., spacer is not predetermined), post-processing is run to find those PRIDICT sequences whose spacer sequence matches the input spacer defined at design time.
3. **Codon optimality.** For each designed mutation, keep nucleotide changes that lead to optimal (most frequent in humans) codons. Exception: if the original AA was encoded by the optimal codon, select the *second-best* synonymous codon instead.
4. **Resolve multiple valid options.** Multiple candidates may remain per AA change due to different PAM/seed-disrupting choices and/or synonymous additions for Hamming distance. Iterate across all remaining options and select the single pegRNA with the highest PRIDICT editing score as the final design.
5. **Final checks.** As part of loading PRIDICT outputs, the pipeline re-verifies that Hamming distance ≥ 2 and that the intended AA change matches the actual change implied by the pegRNA sequence and warns the user if this is not the case. We verified that no designs are excluded at this stage.

In the specific case of LDLR exon 4, all considered PAM and seed mutations were inside the coding sequences.

## Inputs

The code expects the following inputs (files or variables set at the top of the scripts):

- **Target coding sequence**
  - A FASTA file with the coding DNA sequence for the exon/CDS to edit.

- **Target set specification**
  - Either provide an explicit list of desired amino-acid changes (table/CSV), or let the script enumerate all possible AA substitutions at each position.

- **PAM search parameters**
  - Window: by default ±17 bp around the target codon (with an optional extension in case PAM is not found).

- **Codon preference table (optional)**
  - The post-processing step find the optimal for each edit. The script contains an internal mapping of codons to the frequency of codon usage in humans.

- **PRIDICT outputs**
  - CSV exports from PRIDICT or PRIDICT2.0 for the designed candidate set.
  - The loader expects standard PRIDICT columns such as: `Spacer`, `PAM`, `Strand orientation`, `PBSlength`, `RTlength`, `Genomic index`, `ORF index`, `Seed mutated`, `Total Hamming distance`, `Original_Sequence`, `Edited_Sequence`, and a prediction/score column.

- **(Optional) Genomic↔ORF correspondence file**
  - If you already have a CSV mapping genomic indices to ORF indices, the scripts can read it directly; otherwise, the helpers can generate it from your gDNA and ORF inputs.

