import csv
import os
import glob
from natsort import natsorted
import pandas as pd
import gc
import numpy as np
import random
random.seed(42)
import matplotlib.pyplot as plt
from helpers import *


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


#Codon frequencies based on Human table from GenScript https://www.genscript.com/tools/codon-frequency-table

#These have been double checked for accuracy 
most_frequent_codon = {"F": "TTC", "L": "CTG", "Y": "TAC", "_":"TGA", "H":"CAC", "Q":"CAG", "I":"ATC", "M":"ATG", 
                      "N":"AAC", "K":"AAG", "V":"GTG", "D":"GAC", "E":"GAG", "S":"AGC", "C":"TGC", "W":"TGG",
                       "P":"CCC", "R":"CGG", "T": "ACC", "A":"GCC", "G":"GGC"}

#None means there is only one codon
second_most_frequent_codon = {"F": "TTT", "L": "CTC", "Y": "TAT", "_":"TAA", "H":"CAT", "Q":"CAA", "I":"ATT", "M":None, 
                      "N":"AAT", "K":"AAA", "V":"GTC", "D":"GAT", "E":"GAA", "S":"TCC", "C":"TGT", "W":None,
                       "P":"CCT", "R":"AGA", "T": "ACA", "A":"GCT", "G":"GGG"}

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


def load_csv_to_list_of_dict(filename, read_index_from_column=0, header=None):
    """Loading the file that has been used as an input for PRIDICT as a list of dictionaries"""
    result = []
    
    # Add header if provided
    if header:
        result.append(header)
    
    with open(filename, 'r', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        for i, row in enumerate(reader):
            if read_index_from_column == 1 and 'index_pegRNA' in row:
                row['id'] = row['index_pegRNA']
            else:
                row['id'] = i + 1
            result.append(row)
    return result


def append_csv_data(data, folder_path, seq, spacer_given=0, old=1, long=0, efficiency_key='PRIDICT_editing_Score_deep'):
    """Function to append PRIDICT analysis to each designed pegRNA (we're taking the first row with matching spacer)
    seq is the full FASTA sequence of the region that is modified (including the flanking sites)
    We perform checks whether the spacer from PRIDICT matches the intended SPACER
    
    We use the flag 'old' to make the code compatible with the old output of PRIDICT
    in the old output, the name of the Spacer column was 'Spacer-Sequence'
    in the new output, the name of the Spacer column is 'Protospacer-Sequence'
    
    We use the flag 'long' to specify the form of the PRIDICT output file
    In case we used the long PAM distance, long should be 1
    In case we used the short PAM distance, long should be 0
    
    We use 'spacer_given' to say whether the spacer sequence has to be reconstructed based on seq or it's already given as a column in the input csv file
    If spacer_given is 0, find_spacer function is used to calculate expected_spacer.
    If spacer_given is 1, expected_spacer is directly assigned from the 'Spacer' column in the data.
    """
    for d in data:
        id_value = d['id']
        if long == 0:
            filename = folder_path + '/' + folder_path + str(id_value) + '_pegRNA_Pridict_full.csv'
        elif long == 1:
            filename = folder_path + '/' + folder_path + str(id_value) + '_long_pegRNA_Pridict_full.csv'

        if d['Strand orientation'] == '+' or d['Strand orientation'] == '-':
            try:
                if spacer_given == 0:
                    expected_spacer = find_spacer(seq, d['Strand orientation'], int(d['Position of the first G in the PAM (sense direction)']))
                elif spacer_given == 1:
                    expected_spacer = d['Spacer']

                with open(filename, 'r') as f:
                    reader = csv.DictReader(f)
                    correct_spacer = 0
                    while correct_spacer == 0:
                        new_data = next(reader)
                        if old == 1:
                            pridict_spacer = new_data['Spacer-Sequence']
                        elif old == 0:
                            pridict_spacer = new_data['Protospacer-Sequence']

                        if pridict_spacer == expected_spacer:
                            correct_spacer = 1

                    d.update(new_data)
            except Exception as e:
                print(f"Failed to find a matching spacer for {id_value}: {e}")
                d.update({efficiency_key: -1.0})
        else:
            print(f"Invalid strand orientation for {id_value}")
            d.update({efficiency_key: -1.0})

    return data

def group_dict_by_keys(list_of_dicts, keep_only_frequent_codons = 0, key_ids = [0,1,2]):
    """
    Takes a list of dictionaries and groups them by the first three keys.
    Returns a new list of dictionaries with the grouped dictionaries nested
    within their respective parent dictionaries.
    
    if keep_only_frequent_codons is 1, we keep only the otpimal codon for a mutation
    (the rest is not present in the output). For designing a synonymous change,
    in case the original codon is already optimal, we take the second most optimal.
    In case of W and M, we don't output anything (there are no possible mutations)
    """
    grouped_dict_list = []
    i = 0
    for d in list_of_dicts:
        #values = list(d.values())
        #key = (values[key_ids[0]], values[key_ids[1]], values[key_ids[2]]) #original aa, pos, mutated aa
        key = (d['Original_aa'], d['Position of the first nt in the aa'], d['Mutated_aa'])
        matching_dict = None
        
        for group in grouped_dict_list:
            group_key = (group['original_aa'], group['pos'], group['mutated_aa'])
            if group_key == key:
                matching_dict = group
                break
                
        if matching_dict is None:
            matching_dict = dict(zip(('original_aa', 'pos', 'mutated_aa'), key))
            grouped_dict_list.append(matching_dict)
            
        if(keep_only_frequent_codons == 0):
            matching_dict.setdefault('items', []).append(d)
        
        if(keep_only_frequent_codons == 1):
            try:
                if(key[0] != key[2]): #Non-synonynous change
                    if(d['Mutated_nts'] == most_frequent_codon[d['Mutated_aa']]): #if the codon we're looking at is optimal
                        matching_dict.setdefault('items', []).append(d)         

                elif(key[0] == key[2]): #Synonymous change
                    if(d['Original_nts'] != most_frequent_codon[d['Original_aa']]): #if existing codon is not optimal
                        if(d['Mutated_nts'] == most_frequent_codon[d['Mutated_aa']]): #if we're looking at the optimal codon
                            matching_dict.setdefault('items', []).append(d)  

                    if(d['Original_nts'] == most_frequent_codon[d['Original_aa']]): #if existing codon is already optimal
                        if(d['Mutated_nts'] == second_most_frequent_codon[d['Mutated_aa']]): #if we're looking at the second best codon
                            matching_dict.setdefault('items', []).append(d)
            except:
                #This is a dumb bug - when I skip an exon, I reprint a header in the csv file, which needs to be skipped here
                print(f"Skipping design for {d['id']} due to error in mutated aa value")
        if(i%1000 == 0):
            print(i)
        i += 1
    return grouped_dict_list  

def merge_grouped_dicts(list1, list2): #TESTED
    """Merges two grouped dictionaries
    Should be used before finding the best score
    We do not perform checks to see if we're adding duplicates - in any case only one will be chosen at the end by max score"""
    merged_list = []
    
    for dict1 in list1:
        merged_dict = dict1.copy()
        
        for dict2 in list2:
            if dict1['original_aa'] == dict2['original_aa'] and dict1['pos'] == dict2['pos'] and dict1['mutated_aa'] == dict2['mutated_aa']:
                merged_dict['items'] = dict1.get('items', []) + dict2.get('items', [])
        merged_list.append(merged_dict)
    
    for dict2 in list2:
        is_unique = True
        
        for dict1 in list1:
            if dict1['original_aa'] == dict2['original_aa'] and dict1['pos'] == dict2['pos'] and dict1['mutated_aa'] == dict2['mutated_aa']:
                is_unique = False
                break
        
        if is_unique:
            merged_list.append(dict2)
    
    return merged_list

    
def find_highest_score(scores, efficiency_key = 'PRIDICT_editing_Score_deep'): #TESTED
    """Finds the dictionary with the highest 'PRIDICT_editing_Score_deep' value' among all dictionaries with the same (original_aa, position, mutated_aa)"""
    highest_score = None
    
    for score in scores:
        if highest_score is None or float(score[efficiency_key]) > float(highest_score[efficiency_key]):
            highest_score = score
            
    return highest_score

def print_highest_score_values(data, output_file, efficiency_key = 'PRIDICT_editing_Score_deep'): #TESTED
    # Extract the header keys from the first nested dictionary
    items_keys = list(data[0]['items'][0].keys())
    max_score_list = [] #List of the pegRNAs with the highest scores
    # Open the output file and create a CSV writer
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)

        # Write the header row
        writer.writerow(items_keys)

        # Loop over the dictionaries in the list
        for i,d in enumerate(data):
            # Find the dictionary with the highest 'score' value
            max_score_pegRNA = find_highest_score(d['items'], efficiency_key)
            
            # Write its values to the CSV file
            if(max_score_pegRNA != None):
                writer.writerow(max_score_pegRNA.values())
                max_score_list.append(max_score_pegRNA)
                
    return max_score_list

def segregate_by_format(input_csv, column_name):
    """
    Segregates entries in a CSV file based on the format of a specified column.
    
    For each unique format in the specified column, a new CSV file is created,
    containing only entries from the original CSV file that match that format.
    
    The entries in the new CSV files appear in the same order as in the input file.

    Parameters:
    - input_csv (str): The path to the input CSV file.
    - column_name (str): The name of the column to segregate based on.

    Returns:
    - format_dataframes (dict): A dictionary where keys are unique formats and values are corresponding DataFrames.
    """

    # Read the input CSV file into a DataFrame
    df = pd.read_csv(input_csv)

    # Create a dictionary to store DataFrames for each unique format
    format_dataframes = {}
    unique_formats = set()
    for column_value in df[column_name].tolist():
        unique_formats.add(column_value[:column_value.rfind("_")])
    
    unique_formats = list(unique_formats)
    print(unique_formats)
    # Iterate through unique formats and create DataFrames
    for unique_format in unique_formats:
        # Filter DataFrame based on the unique format
        filtered_df = df[df[column_name].str.startswith(unique_format + '_')]

        # Add the DataFrame to the dictionary
        format_dataframes[unique_format] = filtered_df

        # Save the DataFrame to a new CSV file
        output_csv = f"{input_csv}_{unique_format}.csv"
        filtered_df.to_csv(output_csv, index=False)

    return format_dataframes


def find_spacer(string, orientation, pam_pos): #TESTED
    """This function finds spacer based on PAM position
    It changes the first nucleotide of the spacer into G, same as PRIDICT"""
    
    pam_pos = pam_pos-1 #We load from excel where indexing is 1-based to Python where it's 0-based
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


def reconstruct_whole_mut_seq(seq, original_subseq, mutated_subseq, strand_orientation):
    """Finds whole mutate sequence (same length as input fasta) 
    seq: original input fasta
    original_subseq: original seq fed to predict
    mutated_seq: mutated_seq fed to pridict
    strand orientation is needed because in case we use Rv strand, PRIDICT will use reverse complement for original_subseq
    and mutated_subseq while seq will be read from Fw strand"""
    
    #In newer versions of PRIDICT, the output contains small letters in the edited region
    #This is not compatible with the rest of the pipeline so we have convert it to upper
    mutated_subseq = mutated_subseq.upper() 
    
    if(strand_orientation == '-'):
        original_subseq = reverse_complement(original_subseq)
        mutated_subseq = reverse_complement(mutated_subseq)

    index = seq.find(original_subseq)
    L = len(mutated_subseq)
    
    seq1 = seq[:index]
    seq1 += mutated_subseq
    seq1 += seq[index+L:]
    
    return seq1

def mutated_seed(original_seq, mutated_seq, pam_pos, orientation, printing = 0): #TESTED
    """Verifies that the seed is mutated"""
    pam_pos = pam_pos - 1 #We load from excel where indexing is 1-based to Python where it's 0-based
    #Inputs should have the same lengths as the input FASTA so that pam_pos is correct
    #Note that pam_pos is the position of the first G in the case of + orientation
    #Or the position of first C in the + strand in the case of - orientation
    if(orientation == "+"):
        seed_pos = pam_pos - 4
    if(orientation == "-"):
        seed_pos = pam_pos + 3
    
    seed_original = original_seq[seed_pos:seed_pos + 3]

    seed_mutated = mutated_seq[seed_pos:seed_pos+3]
    
    if(printing == 1):
        print('seed original')
        print(seed_original)
        print('seed mutated')
        print(seed_mutated)
    if(seed_original != seed_mutated):
        return 1
    else:
        return 0


def find_mutated_codon(mutated_seq, aa_pos):
    return translate(mutated_seq[aa_pos:aa_pos+3])


def find_orf_index(correspondence_file, genomic_target_index):
    """
    Finds the corresponding 'ORF_index' for a given 'genomic_target_index' from a correspondence file.
    If no exact match is found, it returns the 'ORF_index' of the closest 'Genomic_index'.

    Parameters:
    - correspondence_file (str): The path to the correspondence csv file.
    - genomic_target_index (int): The target genomic index.

    Returns:
    - orf_index (int): The corresponding 'ORF_index'.
    """

    # Read the correspondence file into a pandas DataFrame
    df = pd.read_csv(correspondence_file)

    # Check if the exact 'Genomic_index' exists
    exact_match = df[df['Genomic index'] == genomic_target_index]

    if not exact_match.empty:
        # If an exact match is found, return the corresponding 'ORF_index'
        orf_index = exact_match['ORF index'].values[0]
    else:
        # If no exact match is found, find the closest 'ORF_index'
        closest_index = df['Genomic index'].sub(genomic_target_index).abs().idxmin()
        orf_index = df.loc[closest_index, 'ORF index']

    return orf_index


def generate_test_data(num_instances, L):
    """
    Generates synthetic pegRNA data for testing purposes.

    Parameters:
    - num_instances (int): Number of synthetic pegRNA instances to generate.

    Returns:
    - test_data (list): List containing synthetic pegRNA data.
    """
    test_data = []

    for i in range(num_instances):
        entry = {
            'editing_efficiency': random.uniform(0, 100),  # Random value between 0 and 1
            'Position of the first G in the PAM': random.uniform(0, L)       # Random value between 0 and 100 (adjust range as needed)
        }

        test_data.append(entry)

    return test_data

def create_pam_position_histogram(data):
    """
    Creates a histogram for the distribution of 'pam_position' values in the input dictionary.

    Parameters:
    - data (dict): A dictionary containing pegRNA data.

    Returns:
    - None: Displays the histogram plot.
    """

    # Extract 'pam_position' values from the data dictionary
    pam_positions = [entry['Position of the first G in the PAM (sense direction)'] for entry in data]

    # Create a histogram
    plt.hist(pam_positions, bins='auto', edgecolor='black')
    
    # Set labels and title
    plt.xlabel('pam_position')
    plt.ylabel('Frequency')
    plt.title('Distribution of pam_position Values')

    # Display the plot
    plt.show()
# Example usage:
# Generate synthetic data with 100 instances for testing
# L = 200
# test_data = generate_test_data(1000, L)
# Nfh = 40
# Nt = 50
# p = 60
# selected_instances = select_random_pegRNA(test_data, Nfh, L, p, Nt, efficiency_key = 'editing_efficiency')
# create_pam_position_histogram(test_data)
# create_pam_position_histogram(selected_instances)

def remove_bad_pegRNAs(grouped_dict_list, seq, efficiency_key = 'PRIDICT_editing_Score_deep'):
    """Removes the dictionaries that are considered bad by fullfilling one of the following criteria:
    1. The expected spacer seq does not match the spacer sequence assumed by PRIDICT
    2. The seed is not actually mutated
    3. (planned to be added later) There are unnecessary mutations in the PAM-containing AAs that do not disrupt PAM
    4. pegRNAs that have length loner than 300 bp, i.e. PBS + 2*RT >= 86"""
    new_grouped_dict_list = []
    i = 0
    for d in grouped_dict_list:
        #setting the header of the dict
        new_dict = {}
        new_dict['original_aa'] = d['original_aa']
        new_dict['pos'] = d['pos']
        new_dict['mutated_aa'] = d['mutated_aa']
        new_dict['items'] = []
        
        #if(d['original_aa']=='C' and d['pos'] == '270' and d['mutated_aa'] == 'K'): #Was used for testing
        #if(d['original_aa']=='S' and d['pos'] == '219' and d['mutated_aa'] == 'S'):
        
        #iterating over items to filter out bad designs
        if 'items' in d:
            for item in d['items']:
                if(item[efficiency_key] != -1): #This happens in some cases, we need to avoid them so that it runs through
                    if(int(item['Total Hamming distance']) > 1):
                        expected_mut_aa = item['Mutated_aa']
                        whole_mut_seq = reconstruct_whole_mut_seq(seq, item['Original_Sequence'], item['Edited_Sequence'], item['Strand orientation'])
                        aa_pos = int(item['Position of the first nt in the aa'])-1 #Converting to Python 0-based indexing
                        actual_mut_aa = find_mutated_codon(whole_mut_seq, aa_pos)

                        if(expected_mut_aa == actual_mut_aa):
                            expected_spacer = find_spacer(seq, item['Strand orientation'], int(item['Position of the first G in the PAM (sense direction)']))
                            pridict_spacer = item['Spacer-Sequence']
                            #print('=====')
                            #print(item['Position of the first G in the PAM (sense direction)'])
                            #print(expected_spacer)
                            #print(pridict_spacer)

                            #condition 1 and 4
                            if(expected_spacer == pridict_spacer and (int(item['PBSlength']) + 2*int(item['RTlength'])) < 86):
                                #condition 2 if seed is expected to be mutated
                                if(item['Seed mutated'] == '1' and mutated_seed(seq, whole_mut_seq, int(item['Position of the first G in the PAM (sense direction)']), item['Strand orientation'])):
                                    #Keep only these
                                    new_dict['items'].append(item)  
                                elif(item['Seed mutated'] == '0'): #if seed is not expected to be mutated, we just print out
                                    new_dict['items'].append(item)
                        else:
                            print('Hamming > 1 but AA mut not well designed')
                            print('[ID: ',item['id'])
                            print('Mutated nts as in item:', item['Mutated_nts'])
                            #print(seq[aa_pos-6:aa_pos+9])
                            print("Mutated nts in original seq", seq[aa_pos:aa_pos+3])
                            print(translate(seq[aa_pos:aa_pos+3]))
                            print('Reconstructed mut seq')
                            #print(whole_mut_seq[aa_pos-6:aa_pos+9])
                            print("Mutated nts in mutated seq", whole_mut_seq[aa_pos:aa_pos+3])
                            print(translate(whole_mut_seq[aa_pos:aa_pos+3]))
                            print('###############')
            new_grouped_dict_list.append(new_dict)
            i = i+1
    return new_grouped_dict_list


def sort_list_of_dicts(my_list, *sorting_keys):
    def custom_sort(item):
        return tuple(item[key] for key in sorting_keys)

    return natsorted(my_list, key=custom_sort) #Note the natural sort


def select_random_pegRNA(data, index_correspondences_file, Nfh, L, p, Nt, efficiency_key, output_file):
    """
    Selects Nt random instances of pegRNA data based on specified criteria.

    Parameters:
    - data (list): A list containing pegRNA data.
    - Nfh (int): "N first half", Number of instances with 'pam_position' < L/2 among the Nt random selections.
    - L (float): Length of the sequence.
    - p (float): The percentile value [0, 100], used to calculate the threshold 'editing_efficiency'.
    - Nt (int): "N total", Total number of random pegRNA instances to be selected.
    - efficiency_key (str): The key in the dictionary that corresponds to 'editing_efficiency' values.

    Returns:
    - selected_instances (list): A list of Nt random pegRNA instances meeting the specified criteria.
    """
    data = data.copy()
    # Extract 'editing_efficiency' values from the data list
    editing_efficiency_values = [float(entry[efficiency_key]) for entry in data]
    
    # Calculate the threshold corresponding to percentile p
    threshold = np.percentile(editing_efficiency_values, p)
    print("Using threshold {}".format(threshold))
    # Filter instances with 'editing_efficiency' above the threshold
    filtered_data = [entry for entry in data if float(entry[efficiency_key]) > threshold]

    # Initialize variables for the while loop
    selected_instances = []
    lower_half_count = 0
    upper_half_count = 0

    # While loop to select Nt random instances meeting the criteria
    while len(selected_instances) < Nt:
        # Randomly select an index from the filtered data
        random_index = random.randrange(len(filtered_data))
        random_instance = filtered_data.pop(random_index)
        #position of the first nucleotide in the aa that we chose to mutate (with respect to ORF, not genome)
        #we use this to know the relative position of the aa in the gene
        aa_pos = 3*float(find_orf_index(index_correspondences_file, int(random_instance['Position of the first nt in the aa'])))
        # Check if the 'pam_position' criteria is met
        if (aa_pos <= L/2 and lower_half_count < Nfh) or \
           (aa_pos > L/2 and upper_half_count < (Nt - Nfh)):
            selected_instances.append(random_instance)

            # Update counts based on 'pam_position'
            if aa_pos <= L/2:
                lower_half_count += 1
            else:
                upper_half_count += 1

    with open(output_file, 'w', newline='') as f:
        items_keys = selected_instances[0].keys()
        writer = csv.writer(f)
        # Write the header row
        writer.writerow(items_keys)
        # Loop over the dictionaries in the list
        for i,d in enumerate(selected_instances):
            writer.writerow(d.values())

    return selected_instances

def find_and_print_best_pegRNA(input_csv, PRIDICT_output_folder, genomic_file, ORF_file, index_correspondences_file, output_csv_file_all, output_csv_file_random, efficiency_key='PRIDICT2_0_editing_Score_deep_K562', select_random_stop=0, select_random_syn=0):
    # This is Main()
    # Iterate over each pair of corresponding elements in input_csv and PRIDICT_output_folder
    
    genomic_seq = read_fasta_file(genomic_file)
    merged_dicts = []
    
    for input_csv_single, PRIDICT_output_folder_single in zip(input_csv, PRIDICT_output_folder):
        print(PRIDICT_output_folder_single)
        data = load_csv_to_list_of_dict(input_csv_single, 1)
        data = append_csv_data(data, PRIDICT_output_folder_single, genomic_seq, spacer_given=1, old=1, long=0, efficiency_key=efficiency_key)
        grouped_dict_list = group_dict_by_keys(data, keep_only_frequent_codons=0, key_ids=[2, 1, 4])
        merged_dicts = merge_grouped_dicts(merged_dicts, grouped_dict_list)
        print("Len of merged dictionary is {}".format(len(merged_dicts)))
    #sorting the dicts
    merged_dicts = sort_list_of_dicts(merged_dicts, 'pos', 'original_aa', 'mutated_aa')
    merged_dicts_good = remove_bad_pegRNAs(merged_dicts, genomic_seq, efficiency_key = efficiency_key)
    list_of_best_pegRNA = print_highest_score_values(merged_dicts_good, output_csv_file_all, efficiency_key=efficiency_key)
    ORF_seq = read_fasta_file(ORF_file)
    #For random stop and random syn, print only some pegRNAs
    if select_random_stop == 1:
        # 40 pegRNAs randomly distributed in the first half of the ORF and 10 in the second half of the ORF
        random_data = select_random_pegRNA(list_of_best_pegRNA, index_correspondences_file, Nfh=40, L=len(ORF_seq), p=80, Nt=50,
                                           efficiency_key=efficiency_key, output_file=output_csv_file_random)
    if select_random_syn == 1:
        # 50 randomly distributed pegRNA
        random_data = select_random_pegRNA(list_of_best_pegRNA, index_correspondences_file, Nfh=25, L=len(ORF_seq), p=80, Nt=50,
                                           efficiency_key=efficiency_key, output_file=output_csv_file_random)
    return merged_dicts


def read_csv_to_dict(file_path):
    with open(file_path, 'r') as csv_file:
        reader = csv.DictReader(csv_file)
        return list(reader)

def write_dict_to_csv(file_path, data, fieldnames):
    with open(file_path, 'w', newline='') as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(data)

def merge_and_sort_csv(file1, file2, merge_keys, sort_keys, output_file, max_value_key):
    # Read CSV files into dictionaries
    dict1 = read_csv_to_dict(file1)
    dict2 = read_csv_to_dict(file2)

    # Merge dictionaries
    merged_dict = dict1 + dict2

    # Sort the merged dictionary based on specified keys
    sorted_merged_dict = sort_list_of_dicts(merged_dict, *sort_keys)

    # Identify entries with the same values for all given keys
    unique_entries = []
    seen_values = set()

    for entry in sorted_merged_dict:
        values = tuple(entry[key] for key in merge_keys)
        if values not in seen_values:
            seen_values.add(values)
            unique_entries.append(entry)

    # Find the entry with the highest value for a specific key among duplicates
    final_dict = []
    seen_max_values = {}

    for entry in unique_entries:
        key_value = float(entry[max_value_key])
        if key_value not in seen_max_values or key_value > float(seen_max_values[key_value][max_value_key]):
            seen_max_values[key_value] = entry

    # Keep only the entries with the highest values for the specific key
    final_dict = list(seen_max_values.values())

    # Write out the sorted and filtered data to a new CSV file
    fieldnames = final_dict[0].keys()
    write_dict_to_csv(output_file, final_dict, fieldnames)

    # Print the final list of dictionaries
    #for entry in final_dict:
        #print(entry)


def count_mutated_aa_combinations(data, title):
    # Create a dictionary to store counts for each unique combination of 'original_aa' and 'pos'
    counts = {}

    # Iterate through the list of dictionaries
    for entry in data:
        # Extract values for 'original_aa', 'pos', and 'mutated_aa'
        original_aa = entry.get('original_aa')
        pos = entry.get('pos')
        mutated_aa = entry.get('mutated_aa')

        # Create a tuple representing the unique combination of 'original_aa' and 'pos'
        key = (original_aa, pos)

        # If the combination is not in the counts dictionary, add it with an empty set as the value
        if key not in counts:
            counts[key] = set()

        # Add the 'mutated_aa' value to the set for the current combination
        counts[key].add(mutated_aa)

    # Create a new dictionary to store the counts
    result = {}

    # Iterate through the counts dictionary and count the unique 'mutated_aa' values for each combination
    for key, value_set in counts.items():
        result[key] = len(value_set)

    # Plotting the histogram
    counts_values = list(result.values())
    plt.hist(counts_values, bins=max(counts_values)-min(counts_values)+1, edgecolor='black', density=True)
    plt.xlabel('Number of Unique mutated_aa Values')
    plt.ylabel('Frequency')
    plt.title(title)
    plt.xlim([1,7])
    plt.show()

    return result

def load_excel_and_rename_keys(file_path):
    # Load the Excel sheet into a DataFrame
    df = pd.read_excel(file_path)

    # Rename the columns as per the specified mapping
    df = df.rename(columns={'aaref': 'original_aa', 'aapos': 'pos', 'selected_alt_codon': 'mutated_aa'})

    # Convert the DataFrame to a list of dictionaries
    data = df.to_dict(orient='records')

    return data


# example of calls
# whole LDLR receptor, run on 29.09.2024
# custom variants
data = load_csv_to_list_of_dict('LDLR_custom_output.csv')
seq = read_fasta_file('LDLR gDNA.fa') #Reading the fasta file (has to be in .fa format and in the same folder as the code)
data = append_csv_data(data, 'LDLR_custom_', seq, old = 1, long = 0, efficiency_key = 'PRIDICT2_0_editing_Score_deep_K562')
grouped_dict_list = group_dict_by_keys(data, keep_only_frequent_codons = 0)
data_additional_seed = load_csv_to_list_of_dict('LDLR_custom_update.csv')
seq = read_fasta_file('LDLR gDNA.fa') #Reading the fasta file (has to be in .fa format and in the same folder as the code)
data_additional_seed = append_csv_data(data_additional_seed, 'LDLR_custom_update_', seq, old = 1, long = 0, efficiency_key = 'PRIDICT2_0_editing_Score_deep_K562')
grouped_dict_list2 = group_dict_by_keys(data_additional_seed, keep_only_frequent_codons = 0)
merged_dicts = merge_grouped_dicts(grouped_dict_list, grouped_dict_list2)
merged_dicts = sort_list_of_dicts(merged_dicts, 'pos', 'original_aa', 'mutated_aa')
grouped_dict_list_clean = remove_bad_pegRNAs(merged_dicts, seq, efficiency_key = 'PRIDICT2_0_editing_Score_deep_K562')
print_highest_score_values(grouped_dict_list_clean, 'Library_custom_merged.csv', efficiency_key = 'PRIDICT2_0_editing_Score_deep_K562')


# whole LDLR receptor, run on 29.09.2024
# stop variants
data = load_csv_to_list_of_dict('LDLR_stop_output.csv')
seq = read_fasta_file('LDLR gDNA.fa') #Reading the fasta file (has to be in .fa format and in the same folder as the code)
data = append_csv_data(data, 'LDLR_stop_', seq, old = 1, long = 0, efficiency_key = 'PRIDICT2_0_editing_Score_deep_K562')
grouped_dict_list = group_dict_by_keys(data, keep_only_frequent_codons = 1)
data_additional_seed = load_csv_to_list_of_dict('LDLR_stop_update.csv')
data_additional_seed = append_csv_data(data_additional_seed, 'LDLR_stop_update_', seq, old = 1, long = 0, efficiency_key = 'PRIDICT2_0_editing_Score_deep_K562')
grouped_dict_list2 = group_dict_by_keys(data_additional_seed, keep_only_frequent_codons = 1)
merged_dicts = merge_grouped_dicts(grouped_dict_list, grouped_dict_list2)
sorted_dicts = sort_list_of_dicts(merged_dicts, 'pos', 'original_aa', 'mutated_aa')
grouped_dict_list_clean = remove_bad_pegRNAs(sorted_dicts, seq, efficiency_key = 'PRIDICT2_0_editing_Score_deep_K562')
print_highest_score_values(grouped_dict_list_clean, 'Library_stop_merged.csv', efficiency_key = 'PRIDICT2_0_editing_Score_deep_K562')


# whole LDLR receptor, run on 29.09.2024
# synonymous variants
data = load_csv_to_list_of_dict('LDLR_synonymous_output.csv')
seq = read_fasta_file('LDLR gDNA.fa') #Reading the fasta file (has to be in .fa format and in the same folder as the code)
data = append_csv_data(data, 'LDLR_syn_', seq, old = 1, long = 0, efficiency_key = 'PRIDICT2_0_editing_Score_deep_K562')
grouped_dict_list = group_dict_by_keys(data, keep_only_frequent_codons = 1)
data_additional_seed = load_csv_to_list_of_dict('LDLR_synonymous_update.csv')
data_additional_seed = append_csv_data(data_additional_seed, 'LDLR_synonymous_update_', seq, old = 1, long = 0, efficiency_key = 'PRIDICT2_0_editing_Score_deep_K562')
grouped_dict_list2 = group_dict_by_keys(data_additional_seed, keep_only_frequent_codons = 1)
merged_dicts = merge_grouped_dicts(grouped_dict_list, grouped_dict_list2)
merged_dicts = sort_list_of_dicts(merged_dicts, 'pos', 'original_aa', 'mutated_aa')
grouped_dict_list_clean = remove_bad_pegRNAs(merged_dicts, seq, efficiency_key = 'PRIDICT2_0_editing_Score_deep_K562')
print_highest_score_values(grouped_dict_list_clean, 'Library_synonymous_merged.csv', efficiency_key = 'PRIDICT2_0_editing_Score_deep_K562')

