# -*- coding: utf-8 -*-
"""
Created on Thu Apr 28 15:00:00 2016

@author: chuck
"""

from sys import argv
from Bio import SeqIO

#####run from command line as: python script_name raw_errors_file coding_sequencing_file name_of_outfile
script, threshold_file, reference_file, outfile = argv

error_file = list(open(threshold_file, "r"))
reference = list(SeqIO.parse(reference_file, "fasta"))
genome = list(SeqIO.parse(reference_file, "fasta"))

polarity = list(open('/home/chuck/Documents/transcription_error_2.5/polarity2.txt', 'r'))

#coverage counters
total_cov = 0
total_A_cov, total_C_cov, total_G_cov, total_T_cov = 0, 0, 0, 0
#error counter
total_errors = 0

master_dict = {}

master_dict["A_to_C"], master_dict["A_to_G"], master_dict["A_to_T"] = 0, 0, 0
master_dict["C_to_A"], master_dict["C_to_G"], master_dict["C_to_T"] = 0, 0, 0
master_dict["G_to_C"], master_dict["G_to_A"], master_dict["G_to_T"] = 0, 0, 0
master_dict["T_to_C"], master_dict["T_to_G"], master_dict["T_to_A"] = 0, 0, 0


error_list = list(open(threshold_file, "r"))

error_locations = open(outfile.strip(".txt") + "error_locations.txt", "w")

#####Triplet counters with error in the center nucleotide#####
#list of triplets
trinucleotides = ['AAA', 'AAT', 'AAG', 'AAC', 'ATA', 'ATT', 'ATG', 'ATC', 'AGA', 'AGT', 'AGG', 'AGC', 'ACA', 'ACT', 'ACG', 'ACC', 
                  'TAA', 'TAT', 'TAG', 'TAC', 'TTA', 'TTT', 'TTG', 'TTC', 'TGA', 'TGT', 'TGG', 'TGC', 'TCA', 'TCT', 'TCG', 'TCC', 
                  'GAA', 'GAT', 'GAG', 'GAC', 'GTA', 'GTT', 'GTG', 'GTC', 'GGA', 'GGT', 'GGG', 'GGC', 'GCA', 'GCT', 'GCG', 'GCC', 
                  'CAA', 'CAT', 'CAG', 'CAC', 'CTA', 'CTT', 'CTG', 'CTC', 'CGA', 'CGT', 'CGG', 'CGC', 'CCA', 'CCT', 'CCG', 'CCC']
                  
error_types = ["A_to_C", "A_to_G", "A_to_T", "C_to_A", "C_to_G", "C_to_T", "G_to_C", "G_to_A", "G_to_T", "T_to_C", "T_to_G", "T_to_A"]
                  
trip_dict = {} #Define a dictionary to track the counts for each triplet
                  
for tri in trinucleotides:
    trip_dict[tri] = {}
    for error in error_types:
        trip_dict[tri].update({error:0})              

polarity_dict = {}

for gene in polarity:
    gensplit = gene.split("\t")
    polarity_dict[int(gensplit[3])] = gensplit[1]

leading_A, leading_C, leading_G, leading_T = 0, 0, 0, 0
lagging_A, lagging_C, lagging_G, lagging_T = 0, 0, 0, 0

################################################################################
##########Finds the errors in the protein coding regions of the genome##########
################################################################################

per_gene_list = []
per_gen_dict = {}
for gene in genome:
    
    correct_orientation = True
    leading = True

    
    if "complement" in gene.description:
        correct_orientation = False     
        gene_start = int(gene.description.split("complement(")[1].split("..")[0].strip("]").strip(">").strip("<"))
        gene_end =  int(gene.description.split("complement(")[1].split("..")[1].strip(")]").strip(">").strip("<"))
    else:
        gene_start = int(gene.description.split("location=")[1].split("..")[0].strip("]").strip(">").strip("<"))
        gene_end =  int(gene.description.split("location=")[1].split("..")[1].strip(")]").strip(">").strip("<"))
    
    if polarity_dict[gene_end] == "lead":
        leading = True
    else:
        leading = False
        
    location = gene_start - 1
    
    gene_coverage, gene_errors = 0, 0
    
    per_gen_dict["A_to_C"], per_gen_dict["A_to_G"], per_gen_dict["A_to_T"] = 0, 0, 0
    per_gen_dict["C_to_A"], per_gen_dict["C_to_G"], per_gen_dict["C_to_T"] = 0, 0, 0
    per_gen_dict["G_to_C"], per_gen_dict["G_to_A"], per_gen_dict["G_to_T"] = 0, 0, 0
    per_gen_dict["T_to_C"], per_gen_dict["T_to_G"], per_gen_dict["T_to_A"] = 0, 0, 0
    
    #Loop through each position in the current gene
    while location <= gene_end - 1:
        
        current_line = error_file[location].strip("\n").split("\t") #Split current line into each column
        base_position = int(current_line[0]) #Current position in reference
        coverage = int(current_line[2]) + int(current_line[3]) + int(current_line[4]) + int(current_line[5]) #calculate total coverage for position
        current_base = current_line[1] #Identity of the current base
        
        if coverage == 0: #continue to next base if there is no coverage
            location += 1
        
        elif coverage > 0: #If there is coverage, determine if there is an error
            
            #Coverage at each base
            A_cov, C_cov, G_cov, T_cov = int(current_line[2]), int(current_line[3]), int(current_line[4]), int(current_line[5])
            
            #Coverage at other bases
            A_error_check = C_cov + G_cov + T_cov
            C_error_check = A_cov + G_cov + T_cov
            G_error_check = C_cov + A_cov + T_cov
            T_error_check = C_cov + G_cov + A_cov
             
            total_cov += coverage
            #####Checks the current base and then determines if there is an error
            if current_base == "A":
                if correct_orientation == True:
                    total_A_cov += A_cov
                    if leading == True:
                        leading_A += A_cov
                    else:
                        lagging_A += A_cov
                        
                else:
                    total_T_cov += A_cov
                    if leading == True:
                        leading_T += A_cov
                    else:
                        lagging_T += A_cov                        
                
                if A_error_check == 0: ## Move on if there is no coverage at non-reference bases
                    pass
            
                elif A_error_check > 0: #If there is non-consensus coverage, identify and count the error
                    
                    gene_errors += 1                    
                    total_errors += 1                   

                    if correct_orientation == True:
                        if C_cov > 0:
                            master_dict["A_to_C"] += 1
                            per_gen_dict["A_to_C"] += 1
                            error_locations.write(str(base_position) + "\tA_to_C\t" + polarity_dict[gene_end] + "\n")
                        elif G_cov > 0:
                            master_dict["A_to_G"] += 1
                            per_gen_dict["A_to_G"] += 1
                            error_locations.write(str(base_position) + "\tA_to_G\t" + polarity_dict[gene_end] + "\n")
                        elif T_cov > 0:
                            master_dict["A_to_T"] += 1
                            per_gen_dict["A_to_T"] += 1
                            error_locations.write(str(base_position) + "\tA_to_T\t" + polarity_dict[gene_end] + "\n")
                        
                    else:
                        if C_cov > 0:
                            master_dict["T_to_G"] += 1
                            per_gen_dict["T_to_G"] += 1
                            error_locations.write(str(base_position) + "\tT_to_G\t" + polarity_dict[gene_end] + "\n")
                        elif G_cov > 0:
                            master_dict["T_to_C"] += 1
                            per_gen_dict["T_to_C"] += 1
                            error_locations.write(str(base_position) + "\tT_to_C\t" + polarity_dict[gene_end] + "\n")
                        elif T_cov > 0:
                            master_dict["T_to_A"] += 1
                            per_gen_dict["T_to_A"] += 1
                            error_locations.write(str(base_position) + "\tT_to_A\t" + polarity_dict[gene_end] + "\n")
                            
            elif current_base == "C":
                
                if correct_orientation == True:
                    total_C_cov += C_cov
                    if leading == True:
                        leading_C += C_cov
                    else:
                        lagging_C += C_cov
                else:
                    total_G_cov += C_cov
                    if leading == True:
                        leading_G += C_cov
                    else:
                        lagging_G += C_cov
                        
                if C_error_check == 0: ## Move on if there is no coverage at non-reference bases
                    pass
            
                elif C_error_check > 0: #If there is non-consensus coverage, identify and count the error
                    
                    gene_errors += 1
                    total_errors += 1                    
                        
                    if correct_orientation == True:
                        if A_cov > 0:
                            master_dict["C_to_A"] += 1
                            per_gen_dict["C_to_A"] += 1
                            error_locations.write(str(base_position) + "\tC_to_A\t" + polarity_dict[gene_end] + "\n")
                        elif G_cov > 0:
                            master_dict["C_to_G"] += 1
                            per_gen_dict["C_to_G"] += 1
                            error_locations.write(str(base_position) + "\tC_to_G\t" + polarity_dict[gene_end] + "\n")
                        elif T_cov > 0:
                            master_dict["C_to_T"] += 1
                            per_gen_dict["C_to_T"] += 1
                            error_locations.write(str(base_position) + "\tC_to_T\t" + polarity_dict[gene_end] + "\n")
                                                   
                    else:
                        if A_cov > 0:
                            master_dict["G_to_T"] += 1
                            per_gen_dict["G_to_T"] += 1
                            error_locations.write(str(base_position) + "\tG_to_T\t" + polarity_dict[gene_end] + "\n")
                        elif G_cov > 0:
                            master_dict["G_to_C"] += 1
                            per_gen_dict["G_to_C"] += 1
                            error_locations.write(str(base_position) + "\tG_to_C\t" + polarity_dict[gene_end] + "\n")
                        elif T_cov > 0:
                            master_dict["G_to_A"] += 1
                            per_gen_dict["G_to_A"] += 1
                            error_locations.write(str(base_position) + "\tG_to_A\t" + polarity_dict[gene_end] + "\n")
                            
                     
            elif current_base == "G":
                
                if correct_orientation == True:
                    total_G_cov += G_cov
                    if leading == True:
                        leading_G += G_cov
                    else:
                        lagging_G += G_cov                    
                else:
                    total_C_cov += G_cov
                    if leading == True:
                        leading_C += G_cov
                    else:
                        lagging_C += G_cov
                        
                if G_error_check == 0: ## Move on if there is no coverage at non-reference bases
                    pass
            
                elif G_error_check > 0: #If there is non-consensus coverage, identify and count the error
                    
                    gene_errors += 1                    
                    total_errors += 1
                        
                    if correct_orientation == True:
                        if A_cov > 0:
                            master_dict["G_to_A"] += 1
                            per_gen_dict["G_to_A"] += 1
                            error_locations.write(str(base_position) + "\tG_to_A\t" + polarity_dict[gene_end] + "\n")
                        elif C_cov > 0:
                            master_dict["G_to_C"] += 1
                            per_gen_dict["G_to_C"] += 1
                            error_locations.write(str(base_position) + "\tG_to_C\t" + polarity_dict[gene_end] + "\n")
                        elif T_cov > 0:
                            master_dict["G_to_T"] += 1
                            per_gen_dict["G_to_T"] += 1
                            error_locations.write(str(base_position) + "\tG_to_T\t" + polarity_dict[gene_end] + "\n")
                              
                    else:
                        
                        if A_cov > 0:
                            master_dict["C_to_T"] += 1
                            per_gen_dict["C_to_T"] += 1
                            error_locations.write(str(base_position) + "\tC_to_T\t" + polarity_dict[gene_end] + "\n")
                        elif C_cov > 0:
                            master_dict["C_to_G"] += 1
                            per_gen_dict["C_to_G"] += 1
                            error_locations.write(str(base_position) + "\tC_to_G\t" + polarity_dict[gene_end] + "\n")
                        elif T_cov > 0:
                            master_dict["C_to_A"] += 1
                            per_gen_dict["C_to_A"] += 1
                            error_locations.write(str(base_position) + "\tC_to_A\t" + polarity_dict[gene_end] + "\n")
                            

            elif current_base == "T":
                
                if correct_orientation == True:
                    total_T_cov += T_cov
                    if leading == True:
                        leading_T += T_cov
                    else:
                        lagging_T += T_cov
                else:
                    total_A_cov += T_cov
                    if leading == True:
                        leading_A += T_cov
                    else:
                        lagging_A += T_cov
                    
                if T_error_check == 0: ## Move on if there is no coverage at non-reference bases
                    pass
            
                elif T_error_check > 0: #If there is non-consensus coverage, identify and count the error
                    
                    gene_errors += 1                    
                    total_errors += 1
  
                    if correct_orientation == True:
                        if A_cov > 0:
                            master_dict["T_to_A"] += 1
                            per_gen_dict["T_to_A"] += 1
                            error_locations.write(str(base_position) + "\tT_to_A\t" + polarity_dict[gene_end] + "\n")
                        elif C_cov > 0:
                            master_dict["T_to_C"] += 1
                            per_gen_dict["T_to_C"] += 1
                            error_locations.write(str(base_position) + "\tT_to_C\t" + polarity_dict[gene_end] + "\n")
                        elif G_cov > 0:
                            master_dict["T_to_G"] += 1
                            per_gen_dict["T_to_G"] += 1
                            error_locations.write(str(base_position) + "\tT_to_G\t" + polarity_dict[gene_end] + "\n")
                            
                    else:
                        if A_cov > 0:
                            master_dict["A_to_T"] += 1
                            per_gen_dict["A_to_T"] += 1
                            error_locations.write(str(base_position) + "\tA_to_T\t" + polarity_dict[gene_end] + "\n")
                        elif C_cov > 0:
                            master_dict["A_to_G"] += 1
                            per_gen_dict["A_to_G"] += 1
                            error_locations.write(str(base_position) + "\tA_to_G\t" + polarity_dict[gene_end] + "\n")
                        elif G_cov > 0:
                            master_dict["A_to_C"] += 1
                            per_gen_dict["A_to_C"] += 1
                            error_locations.write(str(base_position) + "\tA_to_C\t" + polarity_dict[gene_end] + "\n")
                            
                        
            location += 1
            
            gene_coverage += coverage  
    per_gene_list.append("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (str(gene_start), str(polarity_dict[gene_end]), str(gene_coverage), str(gene_end-gene_start),
                                                                                               str(per_gen_dict["A_to_C"]), str(per_gen_dict["A_to_G"]), str(per_gen_dict["A_to_T"]),
                                                                                               str(per_gen_dict["C_to_A"]), str(per_gen_dict["C_to_G"]), str(per_gen_dict["C_to_T"]),
                                                                                               str(per_gen_dict["G_to_A"]), str(per_gen_dict["G_to_C"]), str(per_gen_dict["G_to_T"]),
                                                                                               str(per_gen_dict["T_to_A"]), str(per_gen_dict["T_to_C"]), str(per_gen_dict["T_to_G"])))
    print gene_start, total_errors, total_cov #prints progress gene by gene
    
    
    
out = open(outfile, "w")
out.write("A_to_C\tA_to_G\tA_to_T\tC_to_A\tC_to_G\tC_to_T\tG_to_A\tG_to_C\tG_to_T\tT_to_A\tT_to_C\tT_to_G\tCoding_Errors\tTotal_coding_bases\tA_Coverage\tC_Coverage\tG_Coverage\tT_Coverage\n")

out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (str(master_dict["A_to_C"]), str(master_dict["A_to_G"]), str(master_dict["A_to_T"]), 
                                                                                        str(master_dict["C_to_A"]), str(master_dict["C_to_G"]), str(master_dict["C_to_T"]),
                                                                                        str(master_dict["G_to_A"]), str(master_dict["G_to_C"]), str(master_dict["G_to_T"]), 
                                                                                        str(master_dict["T_to_A"]), str(master_dict["T_to_C"]), str(master_dict["T_to_G"]), 
                                                                                        str(total_errors), str(total_cov), str(total_A_cov), str(total_C_cov), str(total_G_cov), str(total_T_cov)))
    
gene_errors = open(outfile.strip(".txt") + "_geneerrors.txt", "w")

for gene in per_gene_list:
    gene_errors.write(gene)
    
error_locations.close()
gene_errors.close()
    
        

    
    
    
    
    
    
    