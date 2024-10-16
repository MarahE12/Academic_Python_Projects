##This code was generated by Marah Elhariry (Team leader role)


#import relevant modules
from Bio import SeqIO
import re

#### Step 1: open input.txt file,& store input data based on the prefixes ('1', '2', or '3') in a dictionary ####

def read_input_file (file_path):
    input_data = {'1': [], '2': [], '3': []} #dictionaries for the three input types
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip() #strip leading or trailing whitespace
            if line.startswith('1:'):
                input_data['1'].append(line[2:]) #if input type 1, append to the resepective dictionary
            elif line.startswith('2:'):
                input_data['2'].append(line[2:]) #if input type 2, append to the respective dictionary
            elif line.startswith('3:'):
                input_data['3'].append(line[2:]) #in input type 3, append to the respective dictionary
            else:
                #if input does not match the format or not exist, error message
                print('Unknown input type or incorrect format:', line)
                pass

    return input_data #return populated dictionary



#### Step 2: process and organise the different inputs ####

def process_inputs(input_data):
    records = {} #dictionary to store organised inputs

    if '1' in input_data:
        #if input type 1, store as 'protein_ID'
        records['Protein_ID'] = input_data['1']

    if '2' in input_data:
        #check if input type 2
        species_gene_pairs = []  #make a list for the pairs
        for pair in input_data['2']:
            #print("Contents of input_data['2']: ", input_data['2'])
            if pair.strip().startswith('OS=') and 'GN=' in pair:
                # Split species and gene names to make pairs
                species, gene = pair.split('OS=')[1].split('GN=')
                gene = gene.split()[0] #split species and gene names into a tuple
                species_gene_pairs.append((species, gene))  # Append a tuple of (species, gene)
        records['Species_Gene_Pair'] = species_gene_pairs #store as 'species_gene_pairs' in dictionary

    if '3' in input_data:
        # if input type 3, store as 'protein_motif'
        records['Protein_Motif'] = input_data['3']

    return records

##### Step 3: extracting from fasta based on processed inputs ####

def read_fasta (fasta_file, processed_inputs, output_file):
    #empty lists to store different types of records
    protein_id_record = []
    species_gene_record = []
    protein_motif_record = []

    print("Starting read_fasta() function...")

    for input_type, description in processed_inputs.items(): #refer and collect from dictionary

        ##input 1 (protein id)
        if input_type == 'Protein_ID': #from dictionary
            print("Processing Protein_ID...")
            identifiers = description
        # open fasta file, compare identifier with description and append matching records
            with open(fasta_file, 'r') as file:
                for record in SeqIO.parse(file, 'fasta'): #open fasta file
                    for identifier in identifiers:
                        # Remove extra spaces from the identifier before comparison
                        if identifier.strip() in record.id:
                            #print(f"Match found for identifier '{identifier.strip()}' in record description: '{record.description}'")
                            protein_id_record.append(record)

        #input 2 (species and gene)
        elif input_type == 'Species_Gene_Pair':
            print("Processing Species_Gene_Pair...")
            species_gene = description
        # open fasta file, compare species and gene names with description and append matching records
            with open(fasta_file, 'r') as file:
                for record in SeqIO.parse(file, 'fasta'): #open fasta file
                    for species, gene in species_gene: #account for both species and gene separately
                        if species in record.description and gene in record.description:
                            species_gene_record.append(record)


        #input 3 (protein motif)
        elif input_type == 'Protein_Motif':
            print("Processing Protein_Motif...")
            motifs = description
        # open fasta file, compare protein motif with sequences and append matching records
            with open(fasta_file, 'r') as file:
                for record in SeqIO.parse(file, 'fasta'): #open fasta file
                    sequence = str(record.seq)
                    #print(f"Sequence: {sequence}")
                    for motif in motifs:
                        #print(f"Motif: {motif}")
                        matches = re.finditer(motif.strip(), sequence.strip())
                        #print(f"Checking motif '{motif}' in record ID: '{record.id}' and description: '{record.description}'")
                        for match in matches:
                            protein_motif_record.append(record)
                            #print(f"Match found for motif '{motif}' in record ID: '{record.id}' and description: '{record.description}'")

    all_records = protein_id_record + species_gene_record + protein_motif_record

    #print('Fasta Record List Content: ', all_records)
    #create temporary FASTA file
    print("Creating temporary FASTA file...")
    with open(output_file, 'w') as temp_file:
        # write all output using the extracted fasta info from fasta file in fasta format
        SeqIO.write(all_records, temp_file, 'fasta')

    print("Exiting read_fasta() function...")
    return all_records #return all records based on processed inputs in fasta file


########################################################
## function calling ##
input_data = read_input_file('input.txt')
# #print(input_data)
processed_inputs = process_inputs(input_data)
# #print(('processed inputs:', processed_inputs))
fasta_file = 'uniprot-apicomplexa.fasta'
read_fasta(fasta_file, processed_inputs, 'temp_fasta_output.fasta')
