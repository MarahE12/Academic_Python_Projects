## In this assignment, I was tasked to write functions that can read in a DNA sequence, 
## translate a DNA sequence into protein sequence and handle different starting frames and strands
## I will then use these functions to process a chromosome FASTA file in "1_Identify_ORFs.py'

#produce a function to translate DNA into protein
def translate_dna(dna_seq):
    if isinstance(dna_seq, list):
        dna_seq = [seq.upper() for seq in dna_seq] #ensure DNA sequence in list is in upper case
    elif isinstance(dna_seq, str):
        dna_seq = dna_seq.upper() #ensure DNA sequence is in upper case
    else:
        print('Error: DNA sequence must be string or in a list.')
    #create a dictionary of the DNA codons and the amino acids they code for
    #the three stop codons 'TAA, TGA and TAG' will return a '*'
    codon_table = {'TTT': 'F', 'TCT': 'S', 'TAT': 'Y', 'TGT': 'C', 'TTC': 'F', 'TCC': 'S', 'TAC': 'Y', 'TGC': 'C',
           'TTA': 'L', 'TCA': 'S',  'TTG': 'L', 'TCG': 'S',  'TGG': 'W', 'CTT': 'L', 'CCT': 'P', 'CAT': 'H',
              'CGT': 'R', 'CTC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R', 'CTA': 'L', 'CCA': 'P', 'CAA': 'Q',
              'CGA': 'R', 'CTG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R', 'ATT': 'I', 'ACT': 'T', 'AAT': 'N',
           'AGT': 'S', 'ATC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S', 'ATA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',
           'ATG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R', 'GTT': 'V', 'GCT': 'A', 'GAT': 'D', 'GGT': 'G', 'GTC': 'V',
           'GCC': 'A', 'GAC': 'D', 'GGC': 'G', 'GTA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G', 'GTG': 'V', 'GCG': 'A',
           'GAG': 'E', 'GGG': 'G', 'TAA': '*', 'TGA': '*','TAG': '*'}

    prot_seq = ''  #create a string for the protein sequence

    for i in range(0, len(dna_seq), 3): #read the DNA sequence 3 bases at a time
            codon = dna_seq[i:i+3] #define a codon
            amino_acid = codon_table.get(codon, 'X')
            #using the .get function generate amino acid chain from codons using codon table.
            #If codon not present, return 'X'

            prot_seq += amino_acid #add the translated amino acids into the protein sequence string

    return prot_seq #return the translated protein sequence

#produce function to read in FASTA format files, distinguish between description and sequence and count to 1000 lines
def read_fasta(filename):
    record_description = ''  # create string refering to the description of FASTA file
    dna_seq = ''  # create string refering to the DNA seq in the FASTA file

    with open(filename, 'r') as fasta: #with statement to open the file and automatically close after
        count_lines = 0 #line-counter to count the lines in the file, starting from value '0'
        for line in fasta:
            line = line.strip() #remove any empty spaces, tabs etc.
            count_lines += 1 #count first line, add 1

 ###########if count_lines.startswith('>'):##### try and use this
            if count_lines == 1:
                record_description = line[1:] #identify line 1 as description and extract excluding the '>'
            else:
                dna_seq += line #append rest of lines into dna_seq

            if count_lines == 1000:
                break #stop counting lines after reaching 1000

    return record_description, dna_seq #return the 1000 lines of the fasta file


#produce function for reverse complement
def reverse_complement(seq):
        DNA = qr()  # convert DNA into upper case
        rev_DNA = DNA[::-1]  # reverse the DNA sequence
        reverse_seq = ''  # create an empty string for the reverse complement sequence
        table = rev_DNA.maketrans('ATCG', 'TAGC')  # create a translation table for the complement bases

        for base in rev_DNA:
            if base in 'ATCG': #check if base in the sequence is one of the four DNA bases
                reverse_seq += base.translate(table) #append translated complement base using the maketrans table
            else:
                reverse_seq += base #if base is not a DNA base, append as is
        return reverse_seq

#produce function to find open reading frames (ORFs) within protein sequence where ORFs sequence is longer than 50
def orf_finder(prot_seq):
    found_orfs = [] #create a list to store found ORFs

    #using .split function, split input at '*' character to create list of possible ORFs
    possible_orfs = prot_seq.split('*')

    for orf in possible_orfs: #loop through the possible ORFs
        if len(orf) > 50: #check if the length of the split possible ORF is greater than 50 amino acids
            found_orfs.append(orf) #if longer than 50 amino acids, append orf into found_orfs list

    return '\n'.join(found_orfs) #return all identified possible ORFs in protein sequence. One ORF per line





