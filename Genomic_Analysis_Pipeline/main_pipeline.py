#import functions and modules
from step1_fasta_search import read_input_file, process_inputs, read_fasta
from step2_BLAST import first10seq, run_blast_pipeline, process_blast_hits, process_blast_output
from step3_MSA import process_alignment
from step4_PT import create_phylogenetic_tree_NJ, create_phylogenetic_tree_UPGMA, literature_search
from matplotlib.backends.backend_pdf import PdfPages
import ssl

pipeline_completed = False  # Flag to track pipeline completion

# Step 1: Read input file, process inputs, and read FASTA
if not pipeline_completed:
    input_data = read_input_file('input.txt')  # Read input data from 'input.txt'
    processed_inputs = process_inputs(input_data)  # Process the input data
    fasta_file = 'uniprot-apicomplexa.fasta'  # Define the FASTA file
    read_fasta(fasta_file, processed_inputs, 'temp_fasta_output.fasta')  # Read the FASTA file
    pipeline_completed = True  # Mark step completion

# Step 2: Perform BLAST, process hits, generate histograms, etc.
if not pipeline_completed:
    first10seq('temp_fasta_output.fasta', 'first_10_seq.fasta', n=10)
    run_blast_pipeline("first_10_seq.fasta", "blast_output.tsv")
    process_blast_hits("blast_output.tsv", "blast_hits")

    pdf_output = 'BLAST_Results.pdf'  # Define the output PDF file name for BLAST results
    with PdfPages(pdf_output) as pdf:
        tsv_file = "blast_output.tsv"
        # Call the modified process_blast_output function
        generated_histogram_png_files = process_blast_output(tsv_file, pdf)

    pipeline_completed = True


# Step 3: MSA
if not pipeline_completed:
    folder_path = 'blast_hits'  # Define the folder path containing FASTA files
    output_folder = 'aln_files'  # Define the output folder for alignment files
    pdf_output = 'MSA_Results.pdf'  # Define the output PDF file name

    # Create a PDF document to store the MSA plots
    with PdfPages(pdf_output) as pdf:
        # Call the function to process alignments and retrieve generated MSA PNGs
        generated_msa_png_files = process_alignment(folder_path, output_folder, pdf)
        # Add the relevant PNGs to the PDF

    pipeline_completed = True


# Step 4: Phylogenetic Tree
if not pipeline_completed:
    ssl._create_default_https_context = ssl._create_unverified_context
    pdf_output = 'Phylogenetic_Results.pdf'

    with PdfPages(pdf_output) as pdf:
        create_phylogenetic_tree_NJ(pdf)
        create_phylogenetic_tree_UPGMA(pdf)
        input_file = "input.txt"
        output_file = "Pubmed_Report.txt"
        literature_search(input_file, output_file)

    pipeline_completed = True




