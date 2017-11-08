# -*- coding: utf-8 -*-
"""
Created on Mon Aug 28 13:20:19 2017

@author: chuck
"""
from Bio import SeqIO
from sys import argv


script, ref, error_locs, out = argv

genome = list(SeqIO.parse(ref, "fasta"))

error_locations = open(error_locs, 'r')

#Uses the "error locations" output file from the transcription error finder script as an input


error_types = ["A_to_C", "A_to_G", "A_to_T", "C_to_A", "C_to_G", "C_to_T", "G_to_C", "G_to_A", "G_to_T", "T_to_C", "T_to_G", "T_to_A"]
nucs = ["A", "C", "G", "T"]
error_dict = {}


for error in error_types:
    error_dict[error] = []

out_write = open(out, 'w')
    
for gene in genome:
    error_locations = open(error_locs, 'r')

    correct_orientation = True
    leading = True
    #Determine if gene is complemented relative to the reference. Only makes a difference when grabbing gene start and stop information from NCBI genes file
    if "complement" in gene.description:
        correct_orientation = False     
        gene_start = int(gene.description.split("complement(")[1].split("..")[0].strip("]").strip(">").strip("<"))
        gene_end =  int(gene.description.split("complement(")[1].split("..")[1].strip(")]").strip(">").strip("<"))
    else:
        gene_start = int(gene.description.split("location=")[1].split("..")[0].strip("]").strip(">").strip("<"))
        gene_end =  int(gene.description.split("location=")[1].split("..")[1].strip(")]").strip(">").strip("<"))   

    for error in error_locations:
        error_split = error.strip().split("\t")
        location = int(error_split[0])
        if location >= gene_start and location <= gene_end:
            
            if correct_orientation == True:
                loc_in_gene = location - gene_start + 1
                seq = gene.seq[loc_in_gene-9:loc_in_gene]
            else:
                loc_in_gene = gene_end - location + 1
                seq = gene.seq[loc_in_gene-9:loc_in_gene]
            if len(seq) != 0:
                error_dict[error_split[1]].append(seq)                
      
        else: 
            pass

    
### Write output of RNA:DNA hybrids
for error in error_types:
    temp_dict = {}
    for i in range(9):
        temp_dict[i] = {}
        temp_dict[i]["A"], temp_dict[i]["C"], temp_dict[i]["G"], temp_dict[i]["T"] = 0, 0, 0, 0
    if len(error_dict[error]) != 0:
        for sequence in error_dict[error]:

            for i in range(9):
                temp_dict[i][sequence[i]] += 1
            
            
    out_write.write("%s\n" % (error))
    
    for nuc in nucs:
        
        out_write.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (nuc, temp_dict[0][nuc], temp_dict[1][nuc], temp_dict[2][nuc],
                                                                             temp_dict[3][nuc], temp_dict[4][nuc], temp_dict[5][nuc],
                                                                             temp_dict[6][nuc], temp_dict[7][nuc], temp_dict[8][nuc]))
        
        
        
        
        
        
        
        
        
        
    
    