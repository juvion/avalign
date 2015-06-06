#!/usr/bin/python
__author__ = "Xiaoju (Ju) Zhang"


from Bio import SeqIO
from Bio.Alphabet import Gapped, IUPAC


def meta_data(handle, fo):
    """Scans the fasta format alignment file downloaded from UCSC genome browser,
    and extracts UCSC ID, species ID, fragment ID, and aligned sequence from each entry."""

    header = "nucleotide_id\tsp_abbr\tfragment_id\tseq_in_align"
    fo.write("{0}\n".format(header))
    for record in SeqIO.parse(handle, 'fasta'):
        fasta_descript = record.id
        ucsc_id = fasta_descript.split('_')[0]
        species_id = fasta_descript.split('_')[1]
        fragment_id = fasta_descript.split('_')[2] + '_' + fasta_descript.split('_')[3]

        fo.write("{0}\t{1}\t{2}\t{3}\n".format(ucsc_id, species_id, fragment_id, record.seq))


def unique_ids(handle, num_species, fo):
    """Scans the num_species_way fasta format alignment file downloaded from UCSC genome browser,
    and extracts unique UCSC IDs."""

    current_id = ''
    i = 0
    for line in handle:
        #process only the description line of fasta
        if line[0] == '>':
            #check the ucsc id every number of species lines.
            if i % (num_species) == 0:
                ucsc_id = line[1:].split('_')[0]
                if current_id != ucsc_id:
                    fo.write("{0}\n".format(ucsc_id))
                    current_id = ucsc_id
            i += 1


def concatenate(handle, num_species, fo):
    """Scans the num_species_way fasta format alignment file downloaded from UCSC genome browser,
    and concatenate the protein fragment sequences."""

    i = 0
    protein_list = []
    for record in SeqIO.parse(handle, 'fasta'):
        fasta_descript, fasta_seq = record.id, record.seq
        ucsc_species_id = '_'.join(fasta_descript.split('_')[0:2])
        # species_id = fasta_descript.split('_')[1]
        fragment_index = int(fasta_descript.split('_')[2])
        fragment_num = int(fasta_descript.split('_')[3])        
        
        j = i % num_species
        #initialize the protein_list as going through the first block of alignment.
        if len(protein_list) != num_species:
            fragment_counter = 1
            protein_list.append([ucsc_species_id, fasta_seq, fragment_index, fragment_num, fragment_counter])
        #update the list element after the first block of alignment
        else:
            #renew the list when new ucsc_species_id occurrs, and restart the concatenation.
            if protein_list[j][0] != ucsc_species_id:
                # print protein_list[j][0], protein_list[j][1]
                fo.write("{0}\t{1}\t{2}\n".format(protein_list[j][0].split('_')[0], protein_list[j][0].split('_')[1], protein_list[j][1]))
                #error check, to see the fragments total number adds up the same.
                if protein_list[j][4] != protein_list[j][3]:
                    print("Error code:1; Record id:{0}; UCSC Species ID:{1}; Fragment counter:{2}; Fragment number:{3}\n".format(str(i), ucsc_species_id, str(protein_list[j][4]), fragment_num))
                fragment_counter = 1
                protein_list[j] = ([ucsc_species_id, fasta_seq, fragment_index, fragment_num, fragment_counter])
            #If the ucsc_species_id keeps the same, then keep the concatenation.
            else:
                protein_list[j][1] += fasta_seq
                #increment the fragment counter
                protein_list[j][4] += 1
                #error check, to see the fragements indces are ordered.
                if protein_list[j][4] != fragment_index:
                    print("Error code:2; Record id:{0}; UCSC Species ID:{1}; Fragment counter:{2}; Fragment index:{3}\n".format(str(i), ucsc_species_id, str(protein_list[j][4]), fragment_index))
        i += 1



if __name__ == "__main__":
    import sys
    # alignment_fasta = "/media/expand/avalign/ref/data/multiz100way/alignments/knownGene.exonAA.100way.fa"
    alignment_fasta = sys.argv[1]
    fo = open(sys.argv[2], 'w')
    handle = open(alignment_fasta, 'rU')
    concatenate(handle, 100, fo)
    fo.close()
