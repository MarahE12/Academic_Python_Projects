#Declaration:
#I, Marah Elhariry, declare that I have employed Chat-GPT-3.5, to assist in the creation of this .py
# of script. Specifically, I used it to proactively identify and handle potential exceptions and errors that
#may arise during execution.


#Step 1: import toolkit functions
from toolkit import read_fasta
from toolkit import translate_dna
from toolkit import reverse_complement
from toolkit import orf_finder

#Step 2: open the fasta file
record_description, dna_seq = read_fasta('Oryza_sativa.IRGSP-1.0.dna.toplevel.fa1.fasta')
print(f'Description: {record_description}'
      f'\nDNA: {dna_seq}')

#Step 3: translate the DNA sequence from the fasta file into a protein sequence and print
protein_sequence = translate_dna(dna_seq)
print(f'Protein Sequence: {protein_sequence}')

#Step 4: read the first 3 frames of the protein sequence
protein_frame1 = protein_sequence[0:]
protein_frame2 = protein_sequence[1:]
protein_frame3 = protein_sequence[2:]
#print the first 3 frames
print(f'\nFrame 1: {protein_frame1} \n'
      f'Frame 2: {protein_frame2} \n'
      f'Frame 3: {protein_frame3}')

#Step 5: reverse complement the DNA and read the 3 frames of the protein translated reverse_complement
reverse_dna = reverse_complement(dna_seq)

protein_frame4 = translate_dna(reverse_dna[0:])
protein_frame5 = translate_dna(reverse_dna[1:])
protein_frame6 = translate_dna(reverse_dna[2:])

#print the 3 protein frames of the reverse complement
print(f'Frame 4: {protein_frame4} \n'
      f'Frame 5: {protein_frame5} \n'
      f'Frame 6: {protein_frame6}')

#Step 6: identify and print the possible ORFs of the 6 frames, one per line
print(f'\nFrame 1 ORFs: \n{orf_finder(protein_frame1)}')
print(f'\nFrame 2 ORFs: \n{orf_finder(protein_frame2)}')
print(f'\nFrame 3 ORFs: \n{orf_finder(protein_frame3)}')
print(f'\nFrame 4 ORFs: \n{orf_finder(protein_frame4)}')
print(f'\nFrame 5 ORFs: \n{orf_finder(protein_frame5)}')
print(f'\nFrame 6 ORFs: \n{orf_finder(protein_frame6)}')

#Step 7: write the final found ORFs output to a .txt file
with open('Rice_Chr1_ORFs.txt','w') as output_file:
      output_file.write(f'Frame 1 ORFs:\n----- \n {orf_finder(protein_frame1)}'
                        f'\n\nFrame 2 ORFs:\n----- \n{orf_finder(protein_frame2)}'
                        f'\n\nFrame 3 ORFs:\n----- \n{orf_finder(protein_frame3)}'
                        f'\n\nFrame 4 ORFs:\n----- \n{orf_finder(protein_frame4)}'
                        f'\n\nFrame 5 ORFs:\n----- \n{orf_finder(protein_frame5)}'
                        f'\n\nFrame 6 ORFs:\n----- \n{orf_finder(protein_frame6)}')
