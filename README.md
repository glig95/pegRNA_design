# pegRNA_design

This is a computational pipeline that designs pegRNAs for comprehensive mutational scanning of a given coding sequence. The pipeline designs all possible pegRNAs that could be used for a desired edit and such exhaustive set of pegRNAs can then be used as an input for PRIDICT, to select the pegRNAs that would lead to high-efficiency genome editing.
The code iterates over all pairs of existing amino-acids (AAs) and desired mutant codons. For mutant codons, we consider all synonymous, missense and nonsense variants that are different from the existing codon (in total 63 for each AA). Then, for a given AA position in the exon, we search for the nearby PAM sequences (“NGG”, in practice “GG”) that are located up to 17 bp upstream of the desired AA change on the sense strand, or 17 bp downstream, on the antisense strand. When a potential PAM sequence is found, we proceed with the following steps:
1.	If possible, find a synonymous mutation that disrupts at least one G nucleotide in the PAM sequence to be introduced with pegRNA, to avoid repeated editing of the genome. 
2.	If there are no such synonymous mutations, search for possible synonymous mutations inside the seed sequence, i.e. the 3 nucleotides between the PAM and the nick-site. 
In the specific case of LDLR exo 4, all considered PAM and seed mutations were inside the coding sequences. In only two cases (D70M and D70W), out of 83*63 = 5229, neither requirement 1. nor 2. could be satisfied for any proximal PAM region.
To ensure that potential sequencing errors can be easily distinguished from the successful genome edits, we imposed a requirement that the original and mutated exon 4 region must differ in at least 2 nucleotides (have Hamming distance bigger than 1). Introducing at least one nucleotide mutation in the AA of interest and one in the PAM region or the seed region already satisfied this requirement. Only in the special cases where PAM sequence is inside the AA that is being mutated, and the introduced mutation disrupts the PAM sequence too (e.g. glycine GGC into AGC), the Hamming distance was 1. In these situations, we introduced an additional synonymous mutation in the region that is 1 to 3 AAs upstream of the AA of interest. 
For some AAs, there were multiple synonymous mutations that could disrupt PAM, or seed, or ensure that the Hamming distance is bigger than 1. In such cases, we listed all of them to give more opportunities of finding a pegRNA with high editing efficiency as predicted by PRIDICT. This resulted in average number of approximately 15 pegRNAs for a given amino-acid position and the desired mutation. On the other hand, in certain cases some AAs are encoded by only one triplet of nucleotides, and we did not include any synonymous mutations for then in those cases. 
