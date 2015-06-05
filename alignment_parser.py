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

    # #initialize a 2*num_species alignment_list filled with 0.
    # alignment_list = [[0] * 2 for i in range(num_species)]
    
    #set a counter to track record index 
    i = 0
    for record in SeqIO.parse(handle, 'fasta'):
        fasta_descript, fasta_seq = record.id, record.seq
        ucsc_species_id = '_'.join(fasta_descript.split('_')[0:2])
        # species_id = fasta_descript.split('_')[1]
        fragment_index = fasta_descript.split('_')[2]
        fragement_num = fasta_descript.split('_')[3]

        previous_fragment_index = 0
        #the alignment is in size of num_species repeatedly. Every 
        #num_species entries with be the same speices.
        j = i % num_species
        #concatenate the sequence if the current ucsc_species_id is identical
        #with the one in list positioned at jth position.
        if alignment_list[j][0] == ucsc_species_id:
            if fragment_index < fragement_num:
                alignment_list[j][1] += fasta_seq
            if fragment_index == fragement_num:
                alignment_list[j][1] += fasta_seq
                print alignment_list[j][0], alignment_list[j][1]
        #If the current entry is not associate with the jth in the list, then
        #starts with a new one.
        else:
            alignment_list[j] = [ucsc_species_id, fasta_seq]
        i += 1


if __name__ == "__main__":
    import sys
    # alignment_fasta = "/media/expand/avalign/ref/data/multiz100way/alignments/knownGene.exonAA.100way.fa"
    alignment_fasta = sys.argv[1]
    fo = open(sys.argv[2], 'w')
    handle = open(alignment_fasta, 'rU')
    concatenate(handle, 100, fo)
    fo.close()
