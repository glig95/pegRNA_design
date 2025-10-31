import csv
import pandas as pd
import numpy as np
import time
from helpers import *

def hamming_distance(seq1, seq2): #TESTED
    """
    Calculates the Hamming distance between two nucleotide sequences.

    Arguments:
    seq1 -- a string representing the first nucleotide sequence
    seq2 -- a string representing the second nucleotide sequence

    Returns:
    The Hamming distance between the two sequences, i.e., the number of positions
    at which the corresponding nucleotides are different.
    """
    if len(seq1) != len(seq2):
        raise ValueError("Sequence lengths must be equal")
    
    # Check if sequences contain only valid nucleotides
    valid_nucleotides = {'A', 'T', 'G', 'C'}
    if set(seq1) - valid_nucleotides or set(seq2) - valid_nucleotides:
        raise ValueError("Sequences must contain only A, T, G, or C nucleotides")
    
    # Initialize the Hamming distance
    hamming_dist = 0
    
    # Iterate over the characters in the sequences
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            hamming_dist += 1
    
    return hamming_dist


def read_fasta_file(filename): #TESTED
    """
    Reads a FASTA file and returns the sequence as a string.

    Arguments:
    filename -- the name of the FASTA file to read

    Returns:
    A string representing the sequence in the FASTA file.
    """
    sequence = ""
    with open(filename, "r") as f:
        for line in f:
            if not line.startswith(">"):
                sequence += line.strip()
    return sequence.upper() #return the sequence in upper case to ensure compatibility with the downstream code (e.g. genetic table search)


def reverse_complement(sequence): #TESTED
    """
    Returns the reverse complement of a nucleotide sequence.

    Arguments:
    sequence -- a string representing the nucleotide sequence

    Returns:
    The reverse complement of the input sequence, as a string.
    """
    complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
    rev_comp = ""
    for base in reversed(sequence):
        if base not in complement:
            raise ValueError("Invalid nucleotide: {}".format(base))
        rev_comp += complement[base]
    return rev_comp


def search_sequence(sequence, position, target_sequence, n): #TESTED.  NOT USED
    """
    Searches for the first N occurrences of a given sequence in the negative direction
    starting from a given position.

    Arguments:
    sequence -- a string representing the sequence to search
    position -- an integer representing the starting position
    target_sequence -- a string representing the sequence to find
    n -- an integer representing the number of occurrences to find

    Returns:
    A list of integers representing the positions of the found sequences,
    or -1 for each sequence that was not found.
    """
    positions = []
    i = position - 1  # convert 1-based index to 0-based index
    while n > 0 and i >= 0:
        if sequence[i:i+len(target_sequence)] == target_sequence:
            positions.append(i)
            n -= 1
        i -= 1
    while n > 0:
        positions.append(-1)
        n -= 1
    return positions  # reverse the order of the positions


def translate(sequence): #TESTED
    """
    Translates a nucleotide sequence into an amino acid sequence.

    Arguments:
    sequence -- a string representing the nucleotide sequence

    Returns:
    A string representing the corresponding amino acid sequence
    """
    codon_table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }
    protein = ''
    if len(sequence) % 3 != 0:
        print("Warning: Sequence length is not a multiple of 3.")
    for i in range(0, len(sequence), 3):
        codon = sequence[i:i+3]
        if codon in codon_table:
            protein += codon_table[codon]
        else:
            protein += 'X'  # unknown amino acid
    return protein


def find_n_minima(array1, array2, n): #NOT TESTED, NOT USED, KEPT IF NEEDED LATER
    """
    Finds the first N minima coming from any of the input arrays and
    outputs from which array they came.

    Arguments:
    array1 -- a list representing the first input array
    array2 -- a list representing the second input array
    n -- an integer representing the number of minima to find

    Returns:
    A list of tuples, where each tuple contains the following elements:
    - the minimum value
    - the index of the minimum value in the combined array
    - an integer indicating from which array the minimum value came:
      - 0 for array1
      - 1 for array2
    """
    combined_array = array1 + array2
    indices = range(len(combined_array))
    minima = sorted(zip(combined_array, indices))
    result = []
    for value, index in minima:
        if len(result) == n:
            break
        if index < len(array1):
            result.append((value, index, 0))
        else:
            result.append((value, index - len(array1), 1))
    return result

def find_synonymous_coding_sequences(aa_sequence): #TESTED
    """
    Finds all nucleotide sequences that give the same amino acid sequence as the given
    amino acid sequence.

    Arguments:
    aa_sequence -- a string representing the amino acid sequence

    Returns:
    A list of nucleotide sequences that give the same amino acid sequence as the
    given amino acid sequence
    """
    # Define the reverse genetic code table
    genetic_code = {
        'A': ['GCT', 'GCC', 'GCA', 'GCG'],
        'C': ['TGT', 'TGC'],
        'D': ['GAT', 'GAC'],
        'E': ['GAA', 'GAG'],
        'F': ['TTT', 'TTC'],
        'G': ['GGT', 'GGC', 'GGA', 'GGG'],
        'H': ['CAT', 'CAC'],
        'I': ['ATT', 'ATC', 'ATA'],
        'K': ['AAA', 'AAG'],
        'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
        'M': ['ATG'],
        'N': ['AAT', 'AAC'],
        'P': ['CCT', 'CCC', 'CCA', 'CCG'],
        'Q': ['CAA', 'CAG'],
        'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
        'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
        'T': ['ACT', 'ACC', 'ACA', 'ACG'],
        'V': ['GTT', 'GTC', 'GTA', 'GTG'],
        'W': ['TGG'],
        'Y': ['TAT', 'TAC'],
        '_': ['TAA', 'TAG', 'TGA']
    }


    # Find all nucleotide sequences that give the same amino acid sequence
    nucleotide_seqs = ['']
    for aa in aa_sequence:
        new_nucleotide_seqs = []
        for seq in nucleotide_seqs:
            codons = genetic_code[aa]
            for codon in codons:
                if translate(codon) == aa:
                    new_nucleotide_seqs.append(seq + codon)
        nucleotide_seqs = new_nucleotide_seqs

    return nucleotide_seqs


def saturate_mutagenesis(codon): #TESTED
    """
    Performs saturate mutagenesis on the given codon, which involves generating
    all possible triple-nucleotide substitutions in the codon except for the given codon.

    Arguments:
    codon -- a string representing a codon of nucleotides

    Returns:
    A dictionary with keys as nucleotide triplets and values as amino acids
    translated from those triplets, for all possible triple-nucleotide substitutions
    in the given codon except the given input codon
    """
    substitutions = ['A', 'C', 'G', 'T']
    amino_acids = {}

    for sub1 in substitutions:
        for sub2 in substitutions:
            for sub3 in substitutions:
                if sub1+sub2+sub3 == codon: #If the current codon is same as the input codon, skip
                        continue
                new_codon = sub1 + sub2 + sub3 
                amino_acids[new_codon] = translate(new_codon)

    return amino_acids


def generate_stop_codons(codon): #TESTED
    """
    Generates stop codons.

    Arguments:
    codon -- a string representing a codon of nucleotides #Used only in case the input codon is a stop codon

    Returns:
    A dictionary with keys as nucleotide triplets and values as amino acids
    translated from those triplets, where the last nucleotide has been replaced
    to create stop codons.
    """
    stop_codons = ['TAA', 'TAG', 'TGA']
    amino_acids = {}

    for stop_codon in stop_codons:
        if stop_codon == codon:  # If the current codon is the same as the input codon, skip. Used only in case of stop codon as a the input codon
            continue
        amino_acids[stop_codon] = translate(stop_codon)

    return amino_acids


def generate_synonymous_codons(codon): #TESTED
    """
    Generates all possible synonymous codons of a given codon.

    Arguments:
    codon -- a string representing a codon of nucleotides

    Returns:
    A dictionary with keys as nucleotide triplets and values as amino acids
    translated from those triplets, where individual nucleotides have been replaced
    to create synonymous codons. It does not include the input codon.
    """
    nucleotides = ['A', 'C', 'G', 'T']
    amino_acids = {}

    for i in range(3):
        for nt in nucleotides:
            if nt == codon[i]:  # If the current nucleotide is the same as in the input codon, skip
                continue
            new_codon = codon[:i] + nt + codon[i+1:]
            if(translate(new_codon) == translate(codon)):
                amino_acids[new_codon] = translate(new_codon)

    return amino_acids


def add_row_to_csv(file_path, row): #TESTED
    """
    Adds a row to an existing CSV file.

    Args:
        file_path (str): The path to the CSV file.
        row (list): The row to add to the CSV file.

    Returns:
        None
    """
    with open(file_path, "a", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(row)


def find_occurrences(sequence, search_string, start_position, max_distance): #TESTED
    """
    Finds occurrences of a string in a given sequence of characters, searching in the reverse direction
    from a given starting position, up to a maximum distance.

    Args:
        sequence (str): The sequence of characters to search in.
        search_string (str): The string to search for.
        start_position (int): The starting position of the search.
        max_distance (int): The maximum distance to search from the starting position.

    Returns:
        A list of indices where the search string was found, or an empty list if it was not found.
    """
    occurrences = []
    search_length = len(search_string)

    # Iterate backwards from the start position up to the maximum distance
    for i in range(start_position, start_position- max_distance, -1):
        # Check if the search string matches the substring starting at the current position
        if sequence[i-search_length:i] == search_string:
            occurrences.append(i-search_length)

    return occurrences



def find_string_differences(str1, str2): #TESTED
    """
    Finds the region between the first and last differences between two strings of equal length.

    Args:
        str1 (str): The first string.
        str2 (str): The second string.

    Returns:
        A string with a structure 100bp of unmutated_region(original_sequence/mutated_sequence) 100 bp of unmutated_region
    """
    # Check that the strings are the same length
    if len(str1) != len(str2):
        raise ValueError("Input strings must be the same length")

    # Find the first difference
    first_diff_index = None
    for i in range(len(str1)):
        if str1[i] != str2[i]:
            first_diff_index = i
            break

    # If there are no differences, return empty strings
    if first_diff_index is None:
        return "", ""

    # Find the last difference
    last_diff_index = None
    for i in range(len(str1)-1, first_diff_index-1, -1):
        if str1[i] != str2[i]:
            last_diff_index = i
            break

    # Get the regions between the first and last differences
    region1 = str1[first_diff_index:last_diff_index+1]
    region2 = str2[first_diff_index:last_diff_index+1]
    
    before_region1 = str1[first_diff_index-100:first_diff_index]
    after_region2 = str1[last_diff_index+1:last_diff_index+101]
    
    output = before_region1 + '(' + region1 + '/' + region2 + ')' + after_region2
    
    return output


def find_difference_index(str1, str2):
    """Returns -1 if the two strings are not thee same lenght, 
    or if they are the same length but contain no differences or if they contain multiple differences
    If they contain exactly one difference, it returns the index of that difference"""
    
    if len(str1) != len(str2):
        return -1  # Strings are not of equal length, cannot find difference

    difference_index = None  # Index of the differing character, if any

    for i in range(len(str1)):
        if str1[i] != str2[i]:
            if difference_index is not None:
                return -1  # More than one difference found, return -1
            difference_index = i

    if difference_index is None:
        return -1  # No difference found, return -1

    return difference_index


def find_spacer(string, orientation, pam_pos): #TESTED
    """This function finds spacer based on PAM position
    It changes the first nucleotide of the spacer into G, same as PRIDICT"""
    
    pam_pos = pam_pos-1 #When load from excel where indexing is 1-based to Python where it's 0-based
    if orientation == '+':
        start = pam_pos - 21
        end = pam_pos - 1
        substring = string[start:end]
    elif orientation == '-':
        start = pam_pos + 3
        end = pam_pos + 23
        substring = string[start:end]
        substring = reverse_complement(substring)
    else:
        raise ValueError("Invalid orientation parameter: must be '+' or '-'")

    return 'G' + substring[1:]  # Change first character to 'G' (this is already performed by PRIDICT output)


def extract_data_for_custom_mutagenesis(excel_file):
    """
    Extract data from 'aa_pos', 'aa_ref_codon', and 'possible_alt_codons' columns in an Excel file.
    
    Args:
        excel_file (str): The path to the Excel file.

    Returns:
        tuple: A tuple containing lists for 'aa_pos', 'aa_ref_codon', and 'possible_alt_codons'.
    """
    # Load the Excel file into a DataFrame
    df = pd.read_excel(excel_file)

    # Extract the data from the specified columns
    aa_pos_data = df['aapos'].tolist()
    aa_ref_codon_data = df['ref_codon'].tolist()
    #possible_alt_codons_data = df['possible_alt_codons'].tolist()
    selected_alt_codon_data = df['selected_alt_codon'].tolist()
    return aa_pos_data, aa_ref_codon_data, selected_alt_codon_data


def find_genomic_index(orf_index, csv_file):
    """
    Find the genomic index for a given ORF index in a CSV file.

    Parameters:
    - orf_index: The ORF index to search for.
    - csv_file: The name of the CSV file containing 'Amino-acid', 'Genomic index', and 'ORF index' columns.

    Returns:
    - The genomic index corresponding to the given ORF index.
    - -1 in case there's no corresponding genomic index for the given ORF index

    Raises:
    - ValueError if there are multiple matches for the given ORF index.
    - FileNotFoundError if the specified CSV file is not found.
    """

    # Read the CSV file and initialize variables
    try:
        with open(csv_file, 'r') as csvfile:
            reader = csv.DictReader(csvfile)
            matches = []

            # Search for matches in the 'ORF index' column
            for row in reader:
#                 print(row)
#                 print(row['ORF index'])
#                 print(int(row['ORF index']))
                if row['ORF index'] == str(int(float(orf_index))): #There's an error if we try direct conversion from string to int
                    matches.append(row['Genomic index'])

            # Check for multiple matches
            if len(matches) == 0:
                return -1
            elif len(matches) > 1:
                raise ValueError(f"Multiple matches found for ORF index {orf_index} in {csv_file}.")

            # Return the genomic index
            return int(matches[0])

    except FileNotFoundError:
        raise FileNotFoundError(f"CSV file '{csv_file}' not found.")


def generate_pegRNA_reverse_strand(seq, aa_pos, pam_pos, mut_aa_nts, first_nt, print_search = 1): #TESTED
    """This is the similar function as generate pegRNA except it searches for pegRNA in the antisense direction.
        Positions in the function are still given with respect to sense strand.
        There are essentially 2 differences compared to the generate_pegRNA function:
        1. It searches for CC as PAM instead of GG
        2. If it doesn't find a PAM mutating synonymous codon, it proceeds to downstream codons (in the sense strand), instead of upstream in the generate_pegRNA function
        
        
    Args:
        seq(str) - the string sequence of the whole exon4 and the surrounding
        aa_pos (int) - position of the first nucleotide of the AA to be mutated (0-based indexing)
        pam_pos (int) - position of the first G nucleotide of the PAM sequence (0-based indexing), note that we neglect N in NGG
        mut_aa_nts (str) - nucleotide tripled coding for the desired AA mutation
        first_nt (int) - first nucleotide of the first amino-acid in the fasta file to be mutated (needed for defining reference ORF)
        print_search (bool) - whether to print the search process or not
        
    Returns:
        output - a list of tupples containing mutated sequence (str), and three flags: PAM_mutated (bool), seed_mutated (bool) and additional_mutations (to ensure Hamming distance > 1)
    """
    
    output = []
    pam_mutated = 0
    seed_mutated = 0
    additional_mutations = 0
    cannot_mutate_pam = 0
    #Check if PAM is already mutated by the AA mutation
    mut_seq = seq[0:aa_pos] + mut_aa_nts + seq[aa_pos+3:] #whole sequence
    
    if(print_search == 1):
        print('PAM position')
        print(pam_pos)
        print('Original AA:')
        print(seq[aa_pos:aa_pos+3])
        print('Mutated AA')
        print(mut_aa_nts)
    
    if(mut_seq[pam_pos:pam_pos+2] != 'CC'):
        pam_mutated = 1
        if(print_search == 1):
            print('PAM mutated by AA mutation')
        hamming = hamming_distance(seq, mut_seq)
        
        if(hamming > 1):
            output.append((mut_seq, pam_mutated, seed_mutated, additional_mutations))
            
        if(hamming == 1): #In rare cases where PAM is disrupted by the AA mutation that has one nucleotide change, we have to mutate one AA more to ensure Hamming > 1
            max_aa_cnt = 3 #We go up to max_aa amino-acids upstream of PAM to search for synonymous mutations
            cnt = 0
            #Identify index of AA of the second G in PAM.
            aa2_index = (pam_pos+1-first_nt)//3*3+first_nt
            
            while (cnt < max_aa_cnt):
                additional_mutations = 0
                potential_aa_nts = mut_seq[aa2_index+3+cnt*3:aa2_index+6+cnt*3] 
                potential_syn_aas = find_synonymous_coding_sequences(translate(potential_aa_nts))
                if(print_search == 1):
                    print(f'Currently looking {cnt+1} AAs upstream of PAM')
                if(len(potential_syn_aas) > 0):
                    if(print_search == 1):
                        print('Accepted seed mutations')
                    for pot_syn_aa in potential_syn_aas: #Search through all potential synonymous mutations
                        if(pot_syn_aa != potential_aa_nts): #Function that searches for synonymous mutation also outputs the original codon so we check here if they are different   
                            mut_seq = mut_seq[0:aa2_index+3+cnt*3] + pot_syn_aa + mut_seq[aa2_index+6+cnt*3:]
                            additional_mutated = 1
                            if(print_search == 1):
                                print(pot_syn_aa)
                            output.append((mut_seq, pam_mutated, seed_mutated, additional_mutations))
                cnt = cnt+1
                
    elif(pam_pos >= aa_pos and pam_pos+2 <= aa_pos+3): #In this case PAM is inside the AA that we want to mutate and we canot mutate it
        cannot_mutate_pam = 1
        if(print_search == 1):
            print('Cannot mutate PAM because it is inside AA that we are mutating')

    additional_mutations = 0
    #If PAM is not disrupted by this AA change, mutate it by inserting a synonymous mutation
    if(pam_mutated == 0):
        aa_starting_pos = [] #starting nucleotides of GG-containing AAs. Could be one or two 
        #We assume that all PAMS are inside coding sequences as it is generally the case

        #Identify index of AA of the first G in PAM
        aa1_index = (pam_pos-first_nt)//3*3+first_nt
        aa_starting_pos.append(aa1_index)

        #Identify index of AA of the second G in PAM. If not the same as in first G, store it for potential synonymous exchange
        aa2_index = (pam_pos+1-first_nt)//3*3+first_nt
        if(aa2_index != aa1_index):
            aa_starting_pos.append(aa2_index)

        #List all possible mutations that mutate PAM containing AAs into synonymous AAs
        original_pam = ""
        for ind in aa_starting_pos:
            original_pam += seq[ind:ind+3]
        
        if(cannot_mutate_pam == 0):
            potential_pam_mutations = find_synonymous_coding_sequences(translate(original_pam))

            #In rare cases, the first G is inside the AA we want to mutate, so that positon should not be touched
            if(len(aa_starting_pos) == 2 and aa1_index == aa_pos):
                new_pams = []
                for pot_pam in potential_pam_mutations:
                    new_pams.append(mut_aa_nts + pot_pam[3:6])
                potential_pam_mutations = list(set(new_pams)) #Keep only unique ones

            if(print_search == 1):
                print('Original PAM containing codon')
                print(original_pam)
                print('Potential PAM mutating codons')
                print(potential_pam_mutations)

            #Among those, keep only mutations that do not contain 'CC' from PAM any more
            good_pam_mutations = []
            relative_pam_index = pam_pos - aa1_index #PAM position relative to the start of the first PAM-containing AA
            for pot_pam in potential_pam_mutations:
                if(pot_pam[relative_pam_index:relative_pam_index+2] != 'CC'): #One C is OK, CC is not
                    #Check whether by introducing this PAM mutation we do not mess up the desired AA mutation
                    #(relevant only in cases PAM is inside AA of interest)
                    pot_mut_seq = mut_seq[0:aa_starting_pos[0]] + pot_pam + mut_seq[aa_starting_pos[-1]+3:]
                    if(pot_mut_seq[aa_pos:aa_pos+3] == mut_seq[aa_pos:aa_pos+3]):
                        good_pam_mutations.append(pot_pam)
                    

            if(len(good_pam_mutations) > 0):
                pam_mutated = 1
                for good_pam_mut in good_pam_mutations:
                    mut_seq = mut_seq[0:aa_starting_pos[0]] + good_pam_mut + mut_seq[aa_starting_pos[-1]+3:]
                    output.append((mut_seq, pam_mutated, seed_mutated, additional_mutations))
            elif(print_search == 1):
                print('Cannot find synonymous mutations that disrupt PAM')
    
    if(pam_mutated == 0):
        #Not possible to mutate PAM; try to mutate seed
        
        #Seed is defined as 3 nts upstream of PAM
        seed_pos = pam_pos + 3
        
        #Identify index of AA of the seed nucleoutides
        aa1_seed_index = (seed_pos-first_nt)//3*3+first_nt
        aa2_seed_index = (seed_pos+1-first_nt)//3*3+first_nt
        aa3_seed_index = (seed_pos+2-first_nt)//3*3+first_nt
        aa_seed_starting_pos = sorted(list(set([aa1_seed_index, aa2_seed_index, aa3_seed_index])))
        
        original_seed = ""
        for ind in aa_seed_starting_pos:
            original_seed += mut_seq[ind:ind+3]
        if(print_search == 1):   
            print('Original seed containing codons')
            print(original_seed)
        potential_seed_mutations = find_synonymous_coding_sequences(translate(original_seed))
        if(print_search == 1): 
            print('Potential seed disrupting mutations')
            print(potential_seed_mutations)
            
        accepted_seed_mutations = []
        for pot_mut in potential_seed_mutations:
            if(pot_mut != original_seed):
                seed_mut_seq = mut_seq[0:aa_seed_starting_pos[0]] + pot_mut + mut_seq[aa_seed_starting_pos[-1]+3:]
                
                #Check whether there is actual seed mutation
                if(mut_seq[seed_pos:seed_pos+3]!=seed_mut_seq[seed_pos:seed_pos+3]):
                    #Check whether the seed mutation does not mess up intented aa mutation
                    if(seed_mut_seq[aa_pos:aa_pos+3] == mut_seq[aa_pos:aa_pos+3]):
                        seed_mutated = 1
                        accepted_seed_mutations.append(pot_mut)
                        output.append((seed_mut_seq, pam_mutated, seed_mutated, additional_mutations))
                    
        if(print_search == 1): 
            print('Accepted_seed_mutations')
            print(accepted_seed_mutations)
        
    if(pam_mutated == 0 and seed_mutated == 0 and print_search == 1):
        print('Impossible to find mutations that mutate the given PAM or seed into a synonymous mutation')
        
    if(print_search == 1):
        print('################################################') #End of search for a given mutation
    return output



#Main function for tiling mutagenesis
def generate_all_pegRNA(first_nt, last_nt, limit, seq_name, output_csv_file):
    # Assuming you have already defined the functions read_fasta_file, reverse_complement,
    # add_row_to_csv, saturate_mutagenesis, find_occurrences, generate_pegRNA,
    # generate_pegRNA_reverse_strand, translate, find_string_differences, and hamming_distance.

    seq = read_fasta_file(seq_name)  # Reading the fasta file (has to be in .fa format and in the same folder as the code)
    rev_seq = reverse_complement(seq)

    # Prepare Excel header file
    add_row_to_csv(output_csv_file, ['Position of the first nt in the aa', 'Original_aa', 'Original_nts', 'Mutated_aa',
                              'Mutated_nts', 'Position of the first G in the PAM (sense direction)', 'Strand orientation',
                              'PAM mutated', 'Seed mutated', 'Additional synonymous mutations to satisfy Hamming > 1',
                              'Total Hamming distance', 'Difference in sequences'])

    for i in range(first_nt, last_nt, 3):
        aa_pos = i
        aa_nts = seq[i:i + 3]
        
        for mut_nts, mut_aa in saturate_mutagenesis(aa_nts).items():
            found_pam = 0  # A flag that sets whether we found a nearby PAM or not

            pam_positions = find_occurrences(seq, 'GG', aa_pos + 3, limit)  # We search for all GG occurrences within 'limit'-1 nucleotides from AA, including AA itself

            if pam_positions:
                found_pam = 1

            for pam_pos in pam_positions:
                output = generate_pegRNA(seq, aa_pos, pam_pos, mut_nts, first_nt, 0)
                for pegRNA in output:
                    seq_diff = find_string_differences(seq, pegRNA[0])
                    hamming = hamming_distance(seq, pegRNA[0])
                    if(hamming < 4): #We limit Hamming distance to 4. This is a sensible reduction of data to run PRIDICT on sinne efficiency of editing drops substantially with number of mutations
                        add_row_to_csv(output_csv_file, [aa_pos + 1, translate(aa_nts), aa_nts, mut_aa, mut_nts, pam_pos + 1,
                                              '+', pegRNA[1], pegRNA[2], pegRNA[3], hamming, seq_diff])
                if not output:
                    add_row_to_csv(output_csv_file, [aa_pos + 1, translate(aa_nts), aa_nts, mut_aa, mut_nts, pam_pos+1, '+',
                                              'COULD NOT FIND SYNONYMOUS PAM OR SEED DISRUPTING MUTATIONS'])

            # Searching for PAMs in the antisense direction
            rev_ind = len(seq) - i - 3
            pam_positions = find_occurrences(rev_seq, 'GG', rev_ind + 3, limit)
            if pam_positions:
                found_pam = 1

            for pam_pos in pam_positions:
                pam_pos = len(seq) - pam_pos - 3 + 1
                output = generate_pegRNA_reverse_strand(seq, aa_pos, pam_pos, mut_nts, first_nt, 0)
                for pegRNA in output:
                    seq_diff = find_string_differences(seq, pegRNA[0])
                    hamming = hamming_distance(seq, pegRNA[0])
                    if(hamming < 4): #We limit Hamming distance to 4. This is a sensible reduction of data to run PRIDICT on sinne efficiency of editing drops substantially with number of mutations
                        add_row_to_csv(output_csv_file, [aa_pos + 1, translate(aa_nts), aa_nts, mut_aa, mut_nts, pam_pos + 1,
                                              '-', pegRNA[1], pegRNA[2], pegRNA[3], hamming, seq_diff])
                if not output:
                    add_row_to_csv(output_csv_file, [aa_pos + 1, translate(aa_nts), aa_nts, mut_aa, mut_nts, pam_pos + 1, '-',
                                              'COULD NOT FIND SYNONYMOUS PAM OR SEED DISRUPTING MUTATIONS'])

            if found_pam == 0:
                add_row_to_csv(output_csv_file, [aa_pos + 1, translate(aa_nts), aa_nts, mut_aa, mut_nts,
                                          'NO NEARBY PAM SEQUENCES'])


#Main function for stop-gain variants
def generate_stop_pegRNA(first_nt, last_nt, limit, seq_name, output_csv_file):
    # Assuming you have already defined the functions read_fasta_file, reverse_complement,
    # add_row_to_csv, saturate_mutagenesis, find_occurrences, generate_pegRNA,
    # generate_pegRNA_reverse_strand, translate, find_string_differences, and hamming_distance.

    seq = read_fasta_file(seq_name)  # Reading the fasta file (has to be in .fa format and in the same folder as the code)
    rev_seq = reverse_complement(seq)

    # Prepare Excel header file
    add_row_to_csv(output_csv_file, ['Position of the first nt in the aa', 'Original_aa', 'Original_nts', 'Mutated_aa',
                              'Mutated_nts', 'Position of the first G in the PAM (sense direction)', 'Strand orientation',
                              'PAM mutated', 'Seed mutated', 'Additional synonymous mutations to satisfy Hamming > 1',
                              'Total Hamming distance', 'Spacer', 'Difference in sequences'])

    for i in range(first_nt, last_nt, 3):
        aa_pos = i
        aa_nts = seq[i:i + 3]
        
        for mut_nts, mut_aa in generate_stop_codons(aa_nts).items():
            found_pam = 0  # A flag that sets whether we found a nearby PAM or not

            pam_positions = find_occurrences(seq, 'GG', aa_pos + 4, limit)  # We search for all GG occurrences within 'limit'-1 nucleotides from AA, including AA itself

            if pam_positions:
                found_pam = 1

            for pam_pos in pam_positions:
                output = generate_pegRNA(seq, aa_pos, pam_pos, mut_nts, first_nt, 0)
                for pegRNA in output:
                    seq_diff = find_string_differences(seq, pegRNA[0])
                    hamming = hamming_distance(seq, pegRNA[0])
                    if(hamming <= 5): #We limit Hamming distance to 4. This is a sensible reduction of data to run PRIDICT on sinne efficiency of editing drops substantially with number of mutations
                        add_row_to_csv(output_csv_file, [aa_pos + 1, translate(aa_nts), aa_nts, mut_aa, mut_nts, pam_pos + 1,
                                              '+', pegRNA[1], pegRNA[2], pegRNA[3], hamming, find_spacer(seq, '+', pam_pos+1), seq_diff])
                if not output:
                    add_row_to_csv(output_csv_file, [aa_pos + 1, translate(aa_nts), aa_nts, mut_aa, mut_nts, pam_pos+1, '+',
                                              'COULD NOT FIND SYNONYMOUS PAM OR SEED DISRUPTING MUTATIONS'])

            # Searching for PAMs in the antisense direction
            rev_ind = len(seq) - i - 3
            pam_positions = find_occurrences(rev_seq, 'GG', rev_ind + 4, limit)
            if pam_positions:
                found_pam = 1

            for pam_pos in pam_positions:
                pam_pos = len(seq) - pam_pos - 3 + 1
                output = generate_pegRNA_reverse_strand(seq, aa_pos, pam_pos, mut_nts, first_nt, 0)
                for pegRNA in output:
                    seq_diff = find_string_differences(seq, pegRNA[0])
                    hamming = hamming_distance(seq, pegRNA[0])
                    if(hamming <= 5): #We limit Hamming distance to 4. This is a sensible reduction of data to run PRIDICT on sinne efficiency of editing drops substantially with number of mutations
                        add_row_to_csv(output_csv_file, [aa_pos + 1, translate(aa_nts), aa_nts, mut_aa, mut_nts, pam_pos + 1,
                                              '-', pegRNA[1], pegRNA[2], pegRNA[3], hamming, find_spacer(seq, '-', pam_pos+1), seq_diff])
                if not output:
                    add_row_to_csv(output_csv_file, [aa_pos + 1, translate(aa_nts), aa_nts, mut_aa, mut_nts, pam_pos + 1, '-',
                                              'COULD NOT FIND SYNONYMOUS PAM OR SEED DISRUPTING MUTATIONS'])

            if found_pam == 0:
                add_row_to_csv(output_csv_file, [aa_pos + 1, translate(aa_nts), aa_nts, mut_aa, mut_nts,
                                          'NO NEARBY PAM SEQUENCES'])


#Main function for synonymous mutagenesis
def generate_synonymous_pegRNA(first_nt, last_nt, limit, seq_name, output_csv_file):
    # Assuming you have already defined the functions read_fasta_file, reverse_complement,
    # add_row_to_csv, saturate_mutagenesis, find_occurrences, generate_pegRNA,
    # generate_pegRNA_reverse_strand, translate, find_string_differences, and hamming_distance.

    seq = read_fasta_file(seq_name)  # Reading the fasta file (has to be in .fa format and in the same folder as the code)
    rev_seq = reverse_complement(seq)
    
    # Prepare Excel header file
    add_row_to_csv(output_csv_file, ['Position of the first nt in the aa', 'Original_aa', 'Original_nts', 'Mutated_aa',
                              'Mutated_nts', 'Position of the first G in the PAM (sense direction)', 'Strand orientation',
                              'PAM mutated', 'Seed mutated', 'Additional synonymous mutations to satisfy Hamming > 1',
                              'Total Hamming distance', 'Spacer', 'Difference in sequences'])

    for i in range(first_nt, last_nt, 3):
        aa_pos = i
        aa_nts = seq[i:i + 3]
        
        for mut_nts, mut_aa in generate_synonymous_codons(aa_nts).items():
            found_pam = 0  # A flag that sets whether we found a nearby PAM or not

            pam_positions = find_occurrences(seq, 'GG', aa_pos + 4, limit)  # We search for all GG occurrences within 'limit'-1 nucleotides from AA, including AA itself

            if pam_positions:
                found_pam = 1

            for pam_pos in pam_positions:
                output = generate_pegRNA(seq, aa_pos, pam_pos, mut_nts, first_nt, 0)
                for pegRNA in output:
                    seq_diff = find_string_differences(seq, pegRNA[0])
                    hamming = hamming_distance(seq, pegRNA[0])
                    if(hamming <= 5): #We limit Hamming distance to 4. This is a sensible reduction of data to run PRIDICT on sinne efficiency of editing drops substantially with number of mutations
                        add_row_to_csv(output_csv_file, [aa_pos + 1, translate(aa_nts), aa_nts, mut_aa, mut_nts, pam_pos + 1,
                                              '+', pegRNA[1], pegRNA[2], pegRNA[3], hamming, find_spacer(seq, '+', pam_pos+1), seq_diff])
                if not output:
                    add_row_to_csv(output_csv_file, [aa_pos + 1, translate(aa_nts), aa_nts, mut_aa, mut_nts, pam_pos+1, '+',
                                              'COULD NOT FIND SYNONYMOUS PAM OR SEED DISRUPTING MUTATIONS'])

            # Searching for PAMs in the antisense direction
            rev_ind = len(seq) - i - 3
            pam_positions = find_occurrences(rev_seq, 'GG', rev_ind + 4, limit)
            if pam_positions:
                found_pam = 1

            for pam_pos in pam_positions:
                pam_pos = len(seq) - pam_pos - 3 + 1
                output = generate_pegRNA_reverse_strand(seq, aa_pos, pam_pos, mut_nts, first_nt, 0)
                for pegRNA in output:
                    seq_diff = find_string_differences(seq, pegRNA[0])
                    hamming = hamming_distance(seq, pegRNA[0])
                    if(hamming <= 5): #We limit Hamming distance to 4. This is a sensible reduction of data to run PRIDICT on sinne efficiency of editing drops substantially with number of mutations
                        add_row_to_csv(output_csv_file, [aa_pos + 1, translate(aa_nts), aa_nts, mut_aa, mut_nts, pam_pos + 1,
                                              '-', pegRNA[1], pegRNA[2], pegRNA[3], hamming, find_spacer(seq, '-', pam_pos+1), seq_diff])
                if not output:
                    add_row_to_csv(output_csv_file, [aa_pos + 1, translate(aa_nts), aa_nts, mut_aa, mut_nts, pam_pos + 1, '-',
                                              'COULD NOT FIND SYNONYMOUS PAM OR SEED DISRUPTING MUTATIONS'])

            if found_pam == 0:
                add_row_to_csv(output_csv_file, [aa_pos + 1, translate(aa_nts), aa_nts, mut_aa, mut_nts,
                                          'NO NEARBY PAM SEQUENCES'])


#Main function for custom mutagenesis
def generate_custom_pegRNA(first_nt, last_nt, limit, seq_name, output_csv_file, genome_orf_correspondence_file, custom_mutations_file):
    # Assuming you have already defined the functions read_fasta_file, reverse_complement,
    # add_row_to_csv, saturate_mutagenesis, find_occurrences, generate_pegRNA,
    # generate_pegRNA_reverse_strand, translate, find_string_differences, and hamming_distance.

    seq = read_fasta_file(seq_name)  # Reading the fasta file of the genomic loci (has to be in .fa format and in the same folder as the code)
    rev_seq = reverse_complement(seq)
    (aa_positions_list, ref_codons_list, alternative_codons_list) = extract_data_for_custom_mutagenesis(custom_mutations_file)
    
    # Prepare Excel header file
    add_row_to_csv(output_csv_file, ['Position of the first nt in the aa', 'Original_aa', 'Original_nts', 'Mutated_aa',
                              'Mutated_nts', 'Position of the first G in the PAM (sense direction)', 'Strand orientation',
                              'PAM mutated', 'Seed mutated', 'Additional synonymous mutations to satisfy Hamming > 1',
                              'Total Hamming distance', 'Spacer', 'Difference in sequences'])

    for i, ref_codon, alternative_codons in zip(np.array(aa_positions_list)*3, ref_codons_list, alternative_codons_list):
                                                #i is the nucleotide number, so we multiply by 3
        #i is the index of the first nt in the aa in the orf
        #we have to convert that to index of the first nt in the genomic file

        aa_pos = find_genomic_index(i/3, genome_orf_correspondence_file)
        if aa_pos in range(first_nt,  last_nt):
            if(aa_pos != -1):
                #If we found the corresponding aa_pos in the genomic_loci
                #The reason why it sometimes cannot be found is because an aa is split between exons
                #We just skip these aas
                aa_nts = seq[aa_pos:aa_pos + 3]
                #Perform check whether the reference codon in the custom mutations file matches the one in the sequence
                if aa_nts != ref_codon:
                    print('AA_POS, AA_NTS, REF_CODON')
                    print(aa_pos)
                    print(aa_nts)
                    print(ref_codon)
                    raise ValueError("The two variables are not the same")

                for mut_nts in alternative_codons.split(','):
                    mut_aa = translate(mut_nts)
                    found_pam = 0  # A flag that sets whether we found a nearby PAM or not

                    pam_positions = find_occurrences(seq, 'GG', aa_pos + 6, limit)  # We search for all GG occurrences within 'limit'-1 nucleotides upstream from AA, including AA itself and three nt downstream
                    #print('PAM pos fwd:')
                    #print(pam_positions)
                    if pam_positions:
                        found_pam = 1
                    for pam_pos in pam_positions:
                        output = generate_pegRNA(seq, aa_pos, pam_pos, mut_nts, first_nt, 0)
                        for pegRNA in output:
                            seq_diff = find_string_differences(seq, pegRNA[0])
                            hamming = hamming_distance(seq, pegRNA[0])
                            if(hamming < 250): #We limit Hamming distance to 4. This is a sensible reduction of data to run PRIDICT on sinne efficiency of editing drops substantially with number of mutations
                                add_row_to_csv(output_csv_file, [aa_pos + 1, translate(aa_nts), aa_nts, mut_aa, mut_nts, pam_pos + 1,
                                                      '+', pegRNA[1], pegRNA[2], pegRNA[3], hamming,  find_spacer(seq, '+', pam_pos+1), seq_diff])
                        if not output:
                            add_row_to_csv(output_csv_file, [aa_pos + 1, translate(aa_nts), aa_nts, mut_aa, mut_nts, pam_pos+1, '+',
                                                      'COULD NOT FIND SYNONYMOUS PAM OR SEED DISRUPTING MUTATIONS'])

                    # Searching for PAMs in the antisense direction
                    rev_ind = len(seq) - aa_pos - 3
                    pam_positions = find_occurrences(rev_seq, 'GG', rev_ind + 6, limit)
                    #print('PAM pos rev:')
                    #print(pam_positions)
                    if pam_positions:
                        found_pam = 1

                    for pam_pos in pam_positions:
                        pam_pos = len(seq) - pam_pos - 3 + 1
                        output = generate_pegRNA_reverse_strand(seq, aa_pos, pam_pos, mut_nts, first_nt, 0)
                        for pegRNA in output:
                            seq_diff = find_string_differences(seq, pegRNA[0])
                            hamming = hamming_distance(seq, pegRNA[0])
                            if(hamming <250): #We limit Hamming distance to 4. This is a sensible reduction of data to run PRIDICT on sinne efficiency of editing drops substantially with number of mutations
                                add_row_to_csv(output_csv_file, [aa_pos + 1, translate(aa_nts), aa_nts, mut_aa, mut_nts, pam_pos + 1,
                                                      '-', pegRNA[1], pegRNA[2], pegRNA[3], hamming,  find_spacer(seq, '-', pam_pos+1), seq_diff])
                        if not output:
                            add_row_to_csv(output_csv_file, [aa_pos + 1, translate(aa_nts), aa_nts, mut_aa, mut_nts, pam_pos + 1, '-',
                                                      'COULD NOT FIND SYNONYMOUS PAM OR SEED DISRUPTING MUTATIONS'])

                    if found_pam == 0:
                        add_row_to_csv(output_csv_file, [aa_pos + 1, translate(aa_nts), aa_nts, mut_aa, mut_nts,
                                                         'NO NEARBY PAM SEQUENCES'])


#Since we want to tile according to genomic locus sequence but by taking into account only ORF information (exons),
#we have to create correspondences between genomic indices and ORF indices 
#(custom mutations are given with respect to ORF indices, so for each ORF index we have to find the locus coordinate)
#we tile start +1 : end -1

def create_genomic_ORF_correspondences(first_nt_list, last_nt_list, input_gdna, input_orf, output_csv):
    """Creates correspondences between genomic DNA indices and the ORF
    Takes a list of first nucleotides of amino-acids inside each exon, 
    and a list of last nucleotides of amino-acids. 
    As well as two file names: input_gdna and input_orf
    
    Outputs an excel sheet with both input_gdna index and input_orf index for each amino-acid"""

    seq_gdna = read_fasta_file(input_gdna)
    seq_orf = read_fasta_file(input_orf)

    add_row_to_csv(output_csv, ['Amino-acid', 'ORF index', 'Genomic index'])
    
    
    current_exon = 0 #index of the starting exon
    current_nt = first_nt_list[0] #index in the gDNA
    
    i = 0 #starting index in the ORF; we're not skipping the first methionine for completeness
    
    while i<len(seq_orf):

        #Check if we're entering a new exon
        if(current_nt >= last_nt_list[current_exon]): 
            current_exon = current_exon + 1
            current_nt = first_nt_list[current_exon]
            
            #Check whether we skipped some AA due to going into a new exon
            if(seq_orf[i:i+9]!=seq_gdna[current_nt:current_nt+9]): #TODO: this check is dumb, I am just assuming that the three following AA would have been different if aa was not shared between exons 
                print(["Skipping amino acid ", translate(seq_orf[i:i+3]), " at nt orf position ", str(i)])
                i = i+3

            
        print("ORF is ", seq_orf[i:i+3])
        print("gDNA is ", seq_gdna[current_nt:current_nt+3])
        print(current_nt)
        print("ORF translation is ", translate(seq_orf[i:i+3]))
        print('############')
        
        if(seq_gdna[current_nt:current_nt + 3] == seq_orf[i:i+3]):
            #It's a match, write it down
            print(i)
            add_row_to_csv(output_csv, [translate(seq_orf[i:i+3]), i//3 + 1, current_nt])
        else:
            raise ValueError("No correspondence found")
            

        
        current_nt = current_nt + 3
        i = i + 3
        

#Calls for LDLR
first_nt_list = [187, 10863, 13304, 15860, 17205, 18032, 21292, 22154, 23918, 24174, 26733, 27499, 30730, 31010, 33814, 38648, 40153, 41919]
last_nt_list =  [253, 10982, 13423, 16237, 17324, 18151, 21408, 22276, 24085, 24398, 26849, 27636, 30870, 31159, 33981, 38722, 40308, 41954]
input_gdna = 'Input_sequences/LDLR gDNA.fa'
input_orf = 'Input_sequences/ORFs/LDLR ORF_for_custom.fa'
output_csv = 'Correspondences_LDLR.csv'
create_genomic_ORF_correspondences(first_nt_list, last_nt_list, input_gdna, input_orf, output_csv)


#function that goes over different regions (exons) generates pegRNAs. Eeach exon will have its own output
def generate_all_pegRNA_for_all_exons(first_nt_list, last_nt_list, limit_list, input_seq_list, output_file_list):
    for first_nt, last_nt, limit, input_seq, output_file in zip(first_nt_list, last_nt_list, limit_list, input_seq_list, output_file_list):
        generate_all_pegRNA(first_nt, last_nt, limit, input_seq, output_file)
        
#function that goes over different regions (exons) generates stop pegRNAs. Eeach exon will have its own output       
def generate_stop_pegRNA_for_all_exons(first_nt_list, last_nt_list, limit_list, input_seq_list, output_file_list):
    for first_nt, last_nt, limit, input_seq, output_file in zip(first_nt_list, last_nt_list, limit_list, input_seq_list, output_file_list):
        print(first_nt)
        generate_stop_pegRNA(first_nt, last_nt, limit, input_seq, output_file)
        
#function that goes over different regions (exons) generates synonymous pegRNAs. Eeach exon will have its own output       
def generate_synonymous_pegRNA_for_all_exons(first_nt_list, last_nt_list, limit_list, input_seq_list, output_file_list):
    for first_nt, last_nt, limit, input_seq, output_file in zip(first_nt_list, last_nt_list, limit_list, input_seq_list, output_file_list):
        print(first_nt)
        generate_synonymous_pegRNA(first_nt, last_nt, limit, input_seq, output_file)
        
#function that goes over different regions (exons) generates custom pegRNAs. Eeach exon will have its own output       
def generate_custom_pegRNA_for_all_exons(first_nt_list, last_nt_list, limit_list, input_seq_list, output_file_list, genome_ORF_correspondences_file_list, custom_mutations_file_list):
    for first_nt, last_nt, limit, input_seq, output_file, genome_ORF_correspondences_file, custom_mutations_file in zip(first_nt_list, last_nt_list, limit_list, input_seq_list, output_file_list, genome_ORF_correspondences_file_list, custom_mutations_file_list):
        #print(first_nt)
        generate_custom_pegRNA(first_nt, last_nt, limit, input_seq, output_file, genome_ORF_correspondences_file, custom_mutations_file)


def generate_pegRNA(seq, aa_pos, pam_pos, mut_aa_nts, first_nt, print_search = 1): #TESTED
    """This is the main function that generates pegRNA for a given AA position and  PAM position inside a given sequence
    
    Args:
        seq(str) - the string sequence of the whole exon4 and the surrounding
        aa_pos (int) - position of the first nucleotide of the AA to be mutated (0-based indexing)
        pam_pos (int) - position of the first G nucleotide of the PAM sequence (0-based indexing), note that we neglect N in NGG
        mut_aa_nts (str) - nucleotide tripled coding for the desired AA mutation
        first_nt (int) - first nucleotide of the first amino-acid in the fasta file to be mutated (needed for defining reference ORF)
        print_search (bool) - whether to print the search process or not
        
    Returns:
        output - a list of tupples containing mutated sequence (str), and three flags: PAM_mutated (bool), seed_mutated (bool) and additional_mutations (to ensure Hamming distance > 1)
    """
    
    output = []
    pam_mutated = 0
    seed_mutated = 0
    additional_mutations = 0
    cannot_mutate_pam = 0
    #Check if PAM is already mutated by the AA mutation
    mut_seq = seq[0:aa_pos] + mut_aa_nts + seq[aa_pos+3:] #whole sequence
    
    if(print_search == 1):
        print('PAM position')
        print(pam_pos)
        print('Original AA:')
        print(seq[aa_pos:aa_pos+3])
        print('Mutated AA')
        print(mut_aa_nts)
    
    if(mut_seq[pam_pos:pam_pos+2] != 'GG'):
        pam_mutated = 1
        if(print_search == 1):
            print('PAM mutated by AA mutation')
            
        hamming = hamming_distance(seq, mut_seq) #In rare cases where PAM is disrupted by the AA mutation that has one nucleotide change, we have to mutate one AA more to ensure Hamming > 1
        
        if(hamming > 1):
            output.append((mut_seq, pam_mutated, seed_mutated, additional_mutations))
            
        if(hamming == 1):
            max_aa_cnt = 3 #We go up to max_aa amino-acids upstream of PAM to search for synonymous mutations 
            cnt = 0
            
            #Identify index of AA of the first G in PAM
            aa1_index = (pam_pos-first_nt)//3*3+first_nt
            
            while (cnt < max_aa_cnt):
                additional_mutations = 0
                potential_aa_nts = mut_seq[aa1_index-3-cnt*3:aa1_index-cnt*3]
                potential_syn_aas = find_synonymous_coding_sequences(translate(potential_aa_nts))
                if(print_search == 1):
                    print(f'Currently looking {cnt+1} AAs upstream of PAM')
                if(len(potential_syn_aas) > 0):
                    if(print_search == 1):
                        print('Accepted additional mutations')
                    for pot_syn_aa in potential_syn_aas:
                        if(pot_syn_aa != potential_aa_nts): #Function that searches for synonymous mutation also outputs the original codon so we check here whether they are different   
                            mut_seq = mut_seq[0:aa1_index-3-cnt*3] + pot_syn_aa + mut_seq[aa1_index-cnt*3:]
                            additional_mutations = 1
                            if(print_search == 1):
                                print(pot_syn_aa)
                            output.append((mut_seq, pam_mutated, seed_mutated, additional_mutations))
                cnt = cnt+1
    elif(pam_pos >= aa_pos and pam_pos+2 <= aa_pos+3): #In this case PAM is inside the AA that we want to mutate and we canot mutate it
        cannot_mutate_pam = 1
        if(print_search == 1):
            print('Cannot mutate PAM because it is inside AA that we are mutating')
                    
    additional_mutations = 0
    #If PAM is not disrupted by this AA change, mutate it by inserting a synonymous mutation
    if(pam_mutated == 0):
        aa_starting_pos = [] #starting nucleotides of GG-containing AAs. Could be one or two 
        #We assume that all PAMS are inside coding sequences as it is generally the case

        #Identify index of AA of the first G in PAM
        aa1_index = (pam_pos-first_nt)//3*3+first_nt
        aa_starting_pos.append(aa1_index)

        #Identify index of AA of the second G in PAM. If not the same as in first G, store it for potential synonymous exchange
        aa2_index = (pam_pos+1-first_nt)//3*3+first_nt
        if(aa2_index != aa1_index):
            aa_starting_pos.append(aa2_index)

        #List all possible mutations that mutate PAM containing AAs into synonymous AAs
        original_pam = ""
        for ind in aa_starting_pos:
            original_pam += seq[ind:ind+3]
            
        if(cannot_mutate_pam == 0):
            potential_pam_mutations = find_synonymous_coding_sequences(translate(original_pam))

            #In rare cases, the second G is inside the AA we want to mutate, so that positon should not be touched
            if(len(aa_starting_pos) == 2 and aa2_index == aa_pos):
                new_pams = []
                for pot_pam in potential_pam_mutations:
                    new_pams.append(pot_pam[0:3]+mut_aa_nts)
                potential_pam_mutations = list(set(new_pams)) #Keep only unique ones

            if(print_search == 1):
                print('Original PAM containing codon')
                print(original_pam)
                print('Potential PAM mutating codons')
                print(potential_pam_mutations)

            #Among those, keep only mutations that do not contain 'GG' from PAM any more 
            good_pam_mutations = []
            relative_pam_index = pam_pos - aa1_index #PAM position relative to the start of the first PAM-containing AA
            for pot_pam in potential_pam_mutations:
                if(pot_pam[relative_pam_index:relative_pam_index+2] != 'GG'): #One G is OK, GG is not
                    #Check whether by introducing this PAM mutation we do not mess up the desired AA mutation
                    #(relevant only in cases PAM is inside AA of interest)
                    pot_mut_seq = mut_seq[0:aa_starting_pos[0]] + pot_pam + mut_seq[aa_starting_pos[-1]+3:]
                    if(pot_mut_seq[aa_pos:aa_pos+3] == mut_seq[aa_pos:aa_pos+3]):
                        good_pam_mutations.append(pot_pam)
                    
            
            if(len(good_pam_mutations) > 0):
                pam_mutated = 1
                for good_pam_mut in good_pam_mutations:
                    mut_seq = mut_seq[0:aa_starting_pos[0]] + good_pam_mut + mut_seq[aa_starting_pos[-1]+3:]
                    output.append((mut_seq, pam_mutated, seed_mutated, additional_mutations))

            elif(print_search == 1):
                print('Cannot find synonymous mutations that disrupt PAM')
                
    if(pam_mutated == 0):
        #Not possible to mutate PAM; try to mutate seed
        
        #Seed is defined as 3 nts upstream of PAM - be careful that PAM defined in this code does not include N from NGG
        seed_pos = pam_pos - 4
        
        #Identify index of AA of the seed nucleoutides
        aa1_seed_index = (seed_pos-first_nt)//3*3+first_nt
        aa2_seed_index = (seed_pos+1-first_nt)//3*3+first_nt
        aa3_seed_index = (seed_pos+2-first_nt)//3*3+first_nt #TODO BUT IT CAN BE ONLY TWO AAs NO? SAME FOR REVERSE pegRNA
        aa_seed_starting_pos = sorted(list(set([aa1_seed_index, aa2_seed_index, aa3_seed_index])))
        
        original_seed = ""
        for ind in aa_seed_starting_pos:
            original_seed += mut_seq[ind:ind+3]
            
        if(print_search == 1):     
            print('Original seed containing codons')
            print(original_seed)
        potential_seed_mutations = find_synonymous_coding_sequences(translate(original_seed))
        
        if(print_search == 1): 
            print('Potential seed disrupting mutations')
            print(potential_seed_mutations)
            
        #Find acceptable seed-disrupting mutations    
        accepted_seed_mutations = []
        for pot_mut in potential_seed_mutations:
            if(pot_mut != original_seed):
                seed_mut_seq = mut_seq[0:aa_seed_starting_pos[0]] + pot_mut + mut_seq[aa_seed_starting_pos[-1]+3:]
                
                #Check whether there is actual seed mutation
                if(mut_seq[seed_pos:seed_pos+3]!=seed_mut_seq[seed_pos:seed_pos+3]):
                    #Check whether the seed mutation does not mess up intented aa mutation
                    if(seed_mut_seq[aa_pos:aa_pos+3] == mut_seq[aa_pos:aa_pos+3]):
                        seed_mutated = 1
                        accepted_seed_mutations.append(pot_mut)
                        output.append((seed_mut_seq, pam_mutated, seed_mutated, additional_mutations))
        
        if(print_search == 1):
            print('Accepted_seed_mutations')
            print(accepted_seed_mutations)
        
    if(pam_mutated == 0 and seed_mutated == 0 and print_search == 1):
        print('Impossible to find mutations that mutate the given PAM or seed into a synonymous mutation')
        
    if(print_search == 1):
        print('################################################') #End of search for a given mutation
    return output


#Calls for LDLR
first_nt_list = [187, 10863, 13304, 15860, 17205, 18032, 21292, 22154, 23918, 24174, 26733, 27499, 30730, 31010, 33814, 38648, 40153, 41919]
last_nt_list =  [253, 10982, 13423, 16237, 17324, 18151, 21408, 22276, 24085, 24398, 26849, 27636, 30870, 31159, 33981, 38722, 40308, 41954]
last_nt_list = np.array(last_nt_list) + 1
limit_list = [20]*18 #effective limit is 20-3=17 !!!!!!!!!!!!!!
input_seq_list = ['Input_sequences/LDLR gDNA.fa']*18
output_file_list = ['LDLR_synonymous_output_correct_seed.csv']*18
generate_synonymous_pegRNA_for_all_exons(first_nt_list, last_nt_list, limit_list, input_seq_list, output_file_list)


