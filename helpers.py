import pandas as pd
import csv
from natsort import natsorted

# Provided sorting function
def sort_list_of_dicts(my_list, *sorting_keys):
    def custom_sort(item):
        return tuple(item[key] for key in sorting_keys)

    return natsorted(my_list, key=custom_sort)

# Function to load CSV into a list of dictionaries
def load_csv_into_dict(file_path):
    with open(file_path, mode='r', newline='') as file:
        reader = csv.DictReader(file)
        data = list(reader)
    return data

# Function to save sorted CSV data to a new file
def save_sorted_csv(input_file_path, output_file_path):
    # Load the data from the input CSV file
    data = load_csv_into_dict(input_file_path)

    # Define sorting keys
    sorting_keys = ("Position of the first nt in the aa", "Original_aa", "Mutated_aa")

    # Sort the list of dictionaries
    sorted_data = sort_list_of_dicts(data, *sorting_keys)

    # Write the sorted data into a new CSV file
    with open(output_file_path, mode='w', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=sorted_data[0].keys())
        writer.writeheader()
        writer.writerows(sorted_data)



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
    
    i = 3 #starting index in the ORF; skipping first methionine
    
    while i<len(seq_orf):
        
        print("ORF should be ", seq_orf[i:i+3])
        print("gDNA is ", seq_gdna[current_nt:current_nt+3])
        print("current nt is ", current_nt)
        print("AA in ORF is", translate(seq_orf[i:i+3]))
        
        #Check if we're entering a new exon
        if(current_nt + 3 >= last_nt_list[current_exon]):
            #Taking the potential leftover from the current exon
            leftover = seq_gdna[current_nt:last_nt_list[current_exon]]
            current_exon = current_exon + 1
            current_nt = first_nt_list[current_exon]
            print("leftover is ", leftover)
            #Take the AA based on the leftover from the previous exon and beginnning of the new exon
            if(leftover!=""):
                beginning = seq_gdna[current_nt:current_nt+3-len(leftover)]
                print("beginning is ",beginning)
                split_aa = leftover + beginning
                if(split_aa == seq_orf[i:i+3]):
                    print("Found a match for split AA")
                    
                    add_row_to_csv(output_csv, [translate(seq_orf[i:i+3]), i//3 + 1, current_nt])
                    current_nt = first_nt_list[current_exon] + len(beginning)
                    print('############')
                else:
                    raise ValueError("No correspondence found")
                

            
            
        elif(seq_gdna[current_nt:current_nt + 3] == seq_orf[i:i+3]):
            #It's a match, write it down
            add_row_to_csv(output_csv, [translate(seq_orf[i:i+3]), i//3 + 1, current_nt])
            print('############')
            current_nt = current_nt + 3
        else:
            raise ValueError("No correspondence found")
            

        

        i = i + 3

