# def concatenate(handle, fo):
#     """Scans the num_species_way fasta format alignment file downloaded 
#     from UCSC genome browser,and concatenate the protein sequences.
#     parameters:
#     handle: opened FASTA file
#     fo: opened output file"""

#!/usr/bin/python
import sys
from Bio import SeqIO
from Bio.Alphabet import Gapped, IUPAC


isFirstBlock = True
alignment = {}
species_list = []
current_geneID = ''
i = 0

handle = open(sys.argv[1], 'rU')
fo = open(sys.argv[2], 'w')

for record in SeqIO.parse(handle, 'fasta'):
    fasta_descript, fasta_seq = record.id, record.seq
    fasta_geneID = '_'.join(fasta_descript.split()[0].split('_')[:-3])
    fasta_species = fasta_descript.split()[0].split('_')[-3]
    cds_index = int(fasta_descript.split()[0].split('_')[-2])
    cds_total = int(fasta_descript.split()[0].split('_')[-1])   


    if isFirstBlock:
        if fasta_species not in species_list:
            species_list.append(fasta_species)
        else:
            isFirstBlock = False

    if fasta_geneID != current_geneID:       
        if not isFirstBlock:
            if i % len(species_list) == 0:
                print str(i % len(species_list))
                # print("Error: {0} block size is not same the species number {1}".format(current_geneID, str(len(species_list))))
            for species in species_list:
                fo.write('"{0}","{1}","{2}"\n'.format(current_geneID, species, alignment[species]))
        current_geneID = fasta_geneID
        alignment = {}
        alignment[fasta_species] = fasta_seq
    else:
        if fasta_species not in alignment:
            alignment[fasta_species] = fasta_seq
        else:
            alignment[fasta_species] += fasta_seq

    i += 1

#write out the last gene's alignment
for species in species_list:
    fo.write('"{0}","{1}","{2}"\n'.format(current_geneID, species, alignment[species]))




    # i = 0
    # protein_list = []
    
    # for record in SeqIO.parse(handle, 'fasta'):
    #     fasta_descript, fasta_seq = record.id, record.seq
    #     #get the accessionID_speciesID for each record
    #     ucsc_id_species = '_'.join(fasta_descript.split('_')[:-2])

    #     #for error check, get the cds_index and total number.
    #     cds_index = int(fasta_descript.split('_')[-2])
    #     cds_total = int(fasta_descript.split('_')[-1])        
        
    #     #set up j to track hg19 line, where j==0
    #     j = i % num_species
    #     #initialize the protein_list as going through the first block 
    #     #of alignment.
    #     if len(protein_list) != num_species:
    #         #count cds fragments
    #         cds_counter = 1
    #         protein_list.append([ucsc_id_species, fasta_seq, cds_index, \
    #                             cds_total, cds_counter])
    #     #update the list element after the first block of alignment
    #     else:
    #         # write out access_id, ucsc_species, concatenated sequence, and
    #         # renew the list when new ucsc_id_species occurrs, and restart the concatenation.
    #         if protein_list[j][0] != ucsc_id_species:
    #             ucsc_id = '_'.join(protein_list[j][0].split('_')[:-1])
    #             species = protein_list[j][0].split('_')[-1]
    #             fo.write('"{0}","{1}","{2}"\n'.format(ucsc_id, species, protein_list[j][1]))
    #             #error check, to see the cds total number adds up the same.
    #             if protein_list[j][4] != protein_list[j][3]:
    #                 print("Error code:1; Record id:{0}; UCSC Species ID:{1}; \
    #                         Fragment counter:{2}; Fragment number:{3}\n".format(str(i), \
    #                         ucsc_id_species, str(protein_list[j][4]), cds_total))
    #             #reset cds_counter to 1
    #             cds_counter = 1
    #             protein_list[j] = ([ucsc_id_species, fasta_seq, cds_index, \
    #                                 cds_total, cds_counter])
    #         #If the ucsc_id_species keeps the same, then keep the concatenation.
    #         else:
    #             protein_list[j][1] += fasta_seq
    #             #increment the cds counter
    #             protein_list[j][4] += 1
    #             #error check, to see the fragements indces are ordered.
    #             if protein_list[j][4] != cds_index:
    #                 print("Error code:2; Record id:{0}; UCSC Species ID:{1}; \
    #                         Fragment counter:{2}; Fragment index:{3}\n".format(str(i), \
    #                         ucsc_id_species, str(protein_list[j][4]), cds_index))
    #                 #correct the tracker according to the cds_index in 
    #                 #order to avoid the accumulating error report.
    #                 protein_list[j][4] = cds_index
    #     i += 1
    # #write out the last access_id's concatenated alignment 
    # for j in range(num_species):
    #     ucsc_id = '_'.join(protein_list[j][0].split('_')[:-1])
    #     species = protein_list[j][0].split('_')[-1]
    #     fo.write('"{0}","{1}","{2}"\n'.format(ucsc_id, species, protein_list[j][1]))