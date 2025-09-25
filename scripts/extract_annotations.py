#!/usr/bin/env python3

# Extract annotations from genbank files
# A function to create a table containing all annotation positions
# The script to extract intergenic regions is still not perfect and needs more attention...
# This is a modified (simplified) version. It does not have a function to write annotation stats of the genbank files.


import os
import argparse
import shutil
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature


# Constants
def output_dir(feature_lookup, extract_intergenic):
    OUTPUT_DIR = GENE_ANNOTATIONS_DIR = feature_lookup
    INTERGENIC_ANNOTATIONS_DIR = os.path.join(OUTPUT_DIR,"spacers")
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    os.makedirs(GENE_ANNOTATIONS_DIR, exist_ok=True)
    if extract_intergenic:
        os.makedirs(INTERGENIC_ANNOTATIONS_DIR, exist_ok=True)
    return GENE_ANNOTATIONS_DIR, INTERGENIC_ANNOTATIONS_DIR


def fix_genbank_file(input_genbank, output_genbank):
    # Some annotation tools write start and stop positions of the gene genbank file that are 
    with open(input_genbank, 'r') as infile, open(output_genbank, 'w') as outfile:
        for record in SeqIO.parse(infile, 'genbank'):           
            for feature in record.features:
                for location in feature.location.parts:
                    if location.start < 0:
                        # Adjust negative start position
                        new_start = max(1, location.start + 1)
                        location.start = new_start
                    if location.end < 0:
                        # Adjust negative end position
                        new_end = max(1, location.end + 1)
                        location.end = new_end
            
            SeqIO.write(record, outfile, "genbank")
        

def extract_annotations(input_dir, feature_lookup, extract_intergenic, min_length):
    gene_sequences = {}
    intergenic_sequences = {}
    illegible_files = []

    # Create output directories
    GENE_ANNOTATIONS_DIR, INTERGENIC_ANNOTATIONS_DIR = output_dir(feature_lookup, extract_intergenic)
    TEMP_DIR = "_temp"
    os.makedirs(TEMP_DIR, exist_ok=True)
    
    for genbank_file in os.listdir(input_dir):    
        try:
            if genbank_file.endswith(".gb"):
                sample_name = genbank_file.rsplit(".gb", 1)[0]
                input_file = os.path.join(input_dir, genbank_file)
                temp_genbank = os.path.join(input_dir, TEMP_DIR, f"_temp_{genbank_file}")
                fix_genbank_file(input_file, temp_genbank) # Fix start/stop positions in the input file first
                
                for genbank in SeqIO.parse(temp_genbank, "genbank"):
                    gene_names = []
                    gene_first_positions = []
                    gene_last_positions = []

                    # Here, we locate, record and extract all annotated genes
                    for feature in genbank.features:
                        if feature.type == feature_lookup:
                            gene_name = feature.qualifiers.get("gene", [""])[0]
                            # sometimes gene name is annotated as locus_tag
                            locus_name = feature.qualifiers.get(
                                "locus_tag", [""])[0]
                            standard_name = feature.qualifiers.get("standard_name", [""])[
                                0]  
                            # or as standard_name when using Geneious
                            product_name = feature.qualifiers.get("product", [""])[0]
                            # Sometimes the name is recorded as note
                            note_name = feature.qualifiers.get("note", [""])[0]
                            # Get gene name
                            if gene_name:
                                annotation_name = gene_name.replace(" ","_")
                            elif locus_name:
                                annotation_name = locus_name.replace(" ","_")
                            elif standard_name:
                                annotation_name = standard_name.replace(" ","_")
                            elif product_name:
                                annotation_name = product_name.replace(" ","_")
                            elif note_name:
                                annotation_name = note_name.replace(" ","_")
                            else:
                                print(f"Warning:\tCannot extract the qualifier '{feature_lookup}'. The annotations of the input {genbank_file} are not included in the final gene files.")
                                continue

                            gene_names.append(annotation_name)

                            # Get nucleotide position (convert it to a regular number not python number)
                            gene_first_positions.append(feature.location.start+1)
                            gene_last_positions.append(feature.location.end) 
                            gene_sequence = feature.extract(genbank.seq)
                            
                            # Store the annotation name in the gene_sequences dictionary
                            if annotation_name in gene_sequences:
                                if sample_name not in gene_sequences[annotation_name]:
                                    gene_sequences[annotation_name][sample_name] = gene_sequence
                                elif len(gene_sequences[annotation_name][sample_name]) < len(gene_sequence):
                                    gene_sequences[annotation_name][sample_name] = gene_sequence
                                    print(f"Warning:\tThere are multiple annotations of {annotation_name} in the sample {sample_name}. The one with the longest lenght is kept here. Please check your original genbank file to make sure if you want to extract this longest annotation.")
                            else:
                                gene_sequences[annotation_name] = {}
                                gene_sequences[annotation_name][sample_name] = gene_sequence

                        else:
                            continue
                        
                        # Check whether any input genbank file contains the type of annotation we are looking at at all. If the gene_names list is empty, then it will stop working on the current genbank file and run the next file.
                        if not gene_names:
                            print(
                                f"Error:\tThe input file: {genbank_file} does not contain {feature_lookup} annotation type.")
                            exit()
                        else:
                            # After getting all above information stored in lists, we can now combind them into a table for future usages
                            table = pd.DataFrame({"Gene": gene_names, "First_position": gene_first_positions,
                                                "Last_position": gene_last_positions, "Source": "Original"})
                            table = table.sort_values(
                                by="First_position").reset_index(drop=True)

                            # Now, we'll do the same with all intergenic non-coding regions (if users specify "True" in the function)

                            if extract_intergenic:
                                intergenic_regions = []
                                intergenic_first_positions = []
                                intergenic_last_positions = []

                                if genbank.annotations["topology"] == "linear":
                                    intergenic_regions.append(
                                        table.loc[0, 0]+"_non_coding_head")
                                    intergenic_first_positions.append(0)
                                    intergenic_last_positions.append(
                                        table.loc[0, 1]-1)

                                for n in range(len(table)):
                                    previous_gene = table.iloc[n, 0]
                                    next_gene = []
                                    intergenic_first_position = table.iloc[n, 2]+1
                                    try :
                                        if n == 0:
                                            if table.iloc[n, 2] < table.iloc[n+1, 1]:
                                                next_gene.append(table.iloc[n+1, 0])
                                                intergenic_last_position = table.iloc[n+1, 1]-1
                                            elif table.iloc[n,2] > table.iloc[n+1, 2] and table.iloc[n,2] > table.iloc[n+2, 1]:
                                                next_gene.append(table.iloc[n+2, 0])
                                                intergenic_last_position = table.iloc[n+2, 1]-1
                                        else:
                                            if table.iloc[n, 2] < table.iloc[n+1, 1] and table.iloc[n, 2] > table.iloc[n-1, 2]:
                                                next_gene.append(table.iloc[n+1, 0])
                                                intergenic_last_position = table.iloc[n+1, 1]-1
                                            elif table.iloc[n, 2] > table.iloc[n+1, 2] and table.iloc[n, 2] < table.iloc[n+2, 1]:
                                                next_gene.append(table.iloc[n+2, 0])
                                                intergenic_last_position = table.iloc[n+2, 1]-1
                                    except IndexError:
                                        if genbank.annotations["topology"] == "circular":
                                            next_gene.append(table.iloc[0,0])
                                            intergenic_last_position = table.iloc[0,1]-1
                                        else:
                                            next_gene.append("_non_coding_tail")
                                            intergenic_last_position = len(genbank.seq)

                                    if next_gene:
                                        if intergenic_last_position - intergenic_first_position > int(min_length):
                                            intergenic_name = previous_gene + "_" + next_gene[0]
                                            intergenic_regions.append(intergenic_name)
                                            intergenic_first_positions.append(
                                                intergenic_first_position)
                                            intergenic_last_positions.append(
                                                intergenic_last_position)

                                intergenic_table = pd.DataFrame(
                                    {"Intergenic_regions": intergenic_regions, "First_position": intergenic_first_positions, "Last_position": intergenic_last_positions})

                                for n, intergenic_region in enumerate(intergenic_regions):
                                    if intergenic_last_positions[n] > intergenic_first_positions[n]:
                                        intergenic_sequence = str(
                                            genbank.seq[intergenic_first_positions[n]-1:intergenic_last_positions[n]])
                                        if intergenic_region not in intergenic_sequences:
                                            intergenic_sequences[intergenic_region] = {
                                            }
                                        intergenic_sequences[intergenic_region
                                                            ][sample_name] = intergenic_sequence

                                    else:
                                        intergenic_sequence = str(genbank.seq[intergenic_first_positions[n]-1:])+str(
                                            genbank.seq[:intergenic_last_positions[n]])
                                        if intergenic_region not in intergenic_sequences:
                                            intergenic_sequences[intergenic_region] = {
                                            }
                                        intergenic_sequences[intergenic_region
                                                            ][sample_name] = intergenic_sequence
                                        
        except Exception as e:
            illegible_files.append(f"The input file {genbank_file} contains illigible format as {e}.")
            continue
    
    # Remove the temp directory we just created
    try:    
        shutil.rmtree(TEMP_DIR, )
    except Exception as e:
        print(f"Error:\tCannot remove the temporary directory: {e}")
        
    # Now, our gene and intergenic sequences dictionaries should contain every region from every sample with their nucleotide sequences
    # We can export those sequences into fasta files, one file per region, each sequence start with the original sample name

    for gene, samples in gene_sequences.items():
        output_file = ""
        if len(samples) == 1:
            SINGLE_SAMPLE = os.path.join(GENE_ANNOTATIONS_DIR, "single_sample")
            os.makedirs(SINGLE_SAMPLE, exist_ok=True)
            output_file = f"{SINGLE_SAMPLE}/{gene}.FNA"
        else:
            output_file = f"{GENE_ANNOTATIONS_DIR}/{gene}.FNA"
            
        with open(output_file, "w") as output_handle:
            for sample_name, sequences in samples.items():
                output_handle.write(f">{sample_name} {gene}\n{sequences}\n")

    if extract_intergenic:
        for intergenic, samples in intergenic_sequences.items():
            output_file = ""
            if len(samples) == 1:
                SINGLE_SAMPLE = os.path.join(INTERGENIC_ANNOTATIONS_DIR, "single_sample")
                os.makedirs(SINGLE_SAMPLE, exist_ok=True)
                output_file = f"{SINGLE_SAMPLE}/{intergenic}.FNA"
            else:
                output_file = f"{INTERGENIC_ANNOTATIONS_DIR}/{intergenic}.FNA"
            
            with open(output_file, "w") as output_handle:
                for sample_name, sequences in samples.items():
                    output_handle.write(f">{sample_name} {intergenic}\n{sequences}\n")
        print(
            f"Info:\tAll annotations defined as {feature_lookup} as well as intergenic spacers have been successfully extracted.")
    else:
        print(
            f"Info:\tAll annotations defined as {feature_lookup} have been successfully extracted.")
    

def main():
    CURRENT_DIR = os.getcwd()
    parser = argparse.ArgumentParser(description="Extract annotations from GenBank files into fasta files.")
    parser.add_argument("--input_dir", "-i", 
                        dest="input_dir", 
                        default=CURRENT_DIR, 
                        required=False, 
                        help="Directory containing GenBank files (.gb). Default = Current working directory.")
    
    parser.add_argument("--feature_lookup", "-f", 
                        dest="feature_lookup", 
                        default="gene", 
                        required=False, 
                        help="Type of annotation to extract (e.g., gene or CDS). Default = gene.")
    
    parser.add_argument("--extract_intergenic", "-s", 
                        dest="extract_intergenic", 
                        action = "store_true",
                        help= "If provided, the script will extract intergenic spacers.")
    
    parser.add_argument("--min_length", "-l", 
                        dest="min_length", 
                        default = 20, 
                        required=False, 
                        help="Minimum length of intergenic spacer to extract. Default = 21bp.")
    
    args = parser.parse_args()

    input_dir = args.input_dir
    feature_lookup = args.feature_lookup
    extract_intergenic = args.extract_intergenic
    min_length = args.min_length

    # Run the function
    extract_annotations(input_dir, feature_lookup, extract_intergenic, min_length)
    
    print(f"Info:\tFinished extracting the annotations.")


if __name__ == "__main__":
    main()