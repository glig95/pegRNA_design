# pegRNA_design

'pegRNA_design_for_tiling' is a computational pipeline that designs prime-editing guide RNAs (pegRNAs) for comprehensive mutational scanning of a given coding sequence. The code designs pegRNAs that can be used for mutating each amino-acid in the coding sequence into all other codons. Such an exhaustive set of pegRNAs can then be used as an input for PRIDICT, a tool that infers the efficiency of the prime editing based on the pegRNA, to select the pegRNAs that would lead to high-efficiency genome editing. After running the PRIDICT for all designed pegRNAs, the best pegRNA is chosen for each of the 21 possible mutations for all amino-acid mutations in the given sequence. The code was run and tested on human low-density lipoproteinn receptor (LDLR) exon 4.

The code iterates over all pairs of existing amino-acids (AAs) and desired mutant codons. For mutant codons, we consider all synonymous, missense, and nonsense variants that are different from the existing codon (in total 63 for each AA). Then, for a given AA position in the exon, we search for the nearby PAM sequences (“NGG”, in practice “GG”) that are located up to 17 bp upstream of the desired AA change on the sense strand, or 17 bp downstream, on the antisense strand. When a potential PAM sequence is found, we proceed with the following steps:
1.	If possible, find a synonymous mutation that disrupts at least one G nucleotide in the PAM sequence to be introduced with pegRNA, to avoid repeated editing of the genome. 
2.	If there are no such synonymous mutations, search for possible synonymous mutations inside the seed sequence, i.e. the 3 nucleotides between the PAM and the nick-site. 

To ensure that potential sequencing errors can be easily distinguished from the successful genome edits, we imposed a requirement that the original and mutated exon 4 region must differ in at least 2 nucleotides (i.e. have Hamming distance bigger than 1). Introducing at least one nucleotide mutation in the AA of interest and one in the PAM region or the seed region already satisfied this requirement. Only in the special cases where the PAM sequence was inside the AA that is being mutated, and the introduced mutation disrupted the PAM sequence too (e.g. glycine GGC into AGC), the Hamming distance was 1. In these situations, we introduced an additional synonymous mutation in the region that is 1 to 3 AAs upstream of the AA of interest. 

For some AAs, there were multiple synonymous mutations that could disrupt PAM, or seed, or ensure that the Hamming distance is bigger than 1. In such cases, we listed all of them to give more opportunities of finding a pegRNA with high editing efficiency by PRIDICT. This resulted in an average number of approximately 15 pegRNAs for a given amino-acid position and the desired mutation. On the other hand, in certain cases, some AAs are encoded by only one triplet of nucleotides, and we did not include any synonymous mutations for them in those cases. Lastly, for certain amino-acids we could not find any PAMs that are 17 bp away from the mutation of interest or closer. In these cases only, we extended the search for PAM sequences to up to 30 bp away from the mutation.

After PRIDICT was run on such sequences, we loaded the PRIDICT output using 'Loading_PRIDICT_output' notebook. Briefly, for each designed mutation we took the best pegRNA as redesigned by PRIDICT. We kept only the nucleotide changes that lead to optimal (most frequent in humans) codons. The exceptions to this were the synonymous mutations (e.g. change of L to L) that already had the optimal codon in the original sequence - for them, we took the second-best codon. For methionine and tryptophan codons, which are always coded by the same triplet, we did not introduce synonymous mutations. After keeping only the best codons possible, we still had multiple options for each mutation. This was due to several possibilities for nearby PAMS and several possibilities for introducing synonymous mutations (that were used to satisfy Hamming distance bigger than 1 requirement). Then we iterated over these designs for each desired amino-acid change at a specific position and kept the one with the highest editing score calculated based on PRIDICT as the final pegRNA. As a part of 'Loading_PRIDICT_output' we also performed additional verifications for ensuring that the Hamming distance is bigger than 1 and that the intended amino-acid change matches the actual one based on pegRNA sequence.

In the specific case of LDLR exon 4, all considered PAM and seed mutations were inside the coding sequences. 
