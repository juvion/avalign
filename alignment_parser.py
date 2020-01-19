#!/usr/bin/python
__author__ = "Xiaoju (Ju) Zhang"


from Bio import SeqIO
from Bio.Alphabet import Gapped, IUPAC


def meta_data(handle, fo):
    """Scans the fasta format alignment file downloaded from UCSC genome browser,
    and extracts UCSC ID, species ID, cds ID, and aligned sequence from each entry.
    parameters:
    handle: opened FASTA file
    fo: opened output file"""

    header = "nucleotide_id\tsp_abbr\tcds_id\tseq_in_align"
    fo.write("{0}\n".format(header))
    for record in SeqIO.parse(handle, 'fasta'):
        fasta_descript = record.id
        access_id = fasta_descript.split('_')[0]
        species_id = fasta_descript.split('_')[1]
        cds_id = fasta_descript.split('_')[2] + '_' + fasta_descript.split('_')[3]

        fo.write("{0}\t{1}\t{2}\t{3}\n".format(access_id, species_id, cds_id, record.seq))


def unique_ids(handle, num_species, fo, accession_type = 'ucsc'):
    """Scans the num_species_way fasta format alignment file downloaded from UCSC genome browser,
    and extracts unique IDs and the cds range information.
    parameters:
    accession_type: ucsc or refseq
    handle: opened FASTA file
    num_species: number of species in the alignment
    fo: opened output file"""

    #setting different index n for fetching the accesion id from FASTA header
    if accession_type == 'ucsc':
        n = 1
    elif accession_type == 'refseq':
        n = 2
    else:
        sys.exit("Error! Please choose correct accession_type: 'ucsc' or 'refseq'.")


    i = 0
    cds_range_list = []
    current_access_id = ''
    header = "accession_id\tCDS_ranges"
    #write the header
    fo.write("{0}\n".format(header))

    for record in SeqIO.parse(handle, 'fasta'):
        if i % num_species == 0:
            fasta_descript = record.description
            #different index n for fetching the accesion id from FASTA header
            access_id = '_'.join(fasta_descript.split('_')[:n])
            species_id = fasta_descript.split('_')[n]
            #error check if hg19 is not in the each alignment block's first line.
            if species_id != 'hg19':
                print("Error1: hg19 is not at the first line of the \
                        alignment. ucsc id:{0}\n".format(access_id))
            else:
                cds_start = fasta_descript[:-1].split(':')[-1].split('-')[0]
                cds_end = fasta_descript[:-1].split(':')[-1].split('-')[1]
                #creat the cds_range in the format matches AVA SQL: no space after ',' in the list.
                cds_range = [int(cds_start), int(cds_end)]
                #initialize the current_access_id when it encounters the first ucsc id of hg19.
                if len(current_access_id) == 0:
                    current_access_id = access_id
                    cds_range_list.append(cds_range)
                else:
                    #renew the cds_range_list and current_access_id if 
                    #alignment starts with a new access_id.
                    if access_id != current_access_id:
                        #print the range in the format matches AVA SQL: 
                        #no space after ',' in the list.
                        fo.write("{0}\t{1}\n".format(current_access_id, str(cds_range_list)))
                        current_access_id = access_id
                        cds_range_list = [cds_range]                        
                    #if the access_id keeps the same, append the cds range to cds_range_list.    
                    else:
                        cds_range_list.append(cds_range)
        i += 1

    #write the last block's record
    fo.write("{0}\t{1}\n".format(current_access_id, str(cds_range_list)))


def concatenate2(handle, fo):
    """Scans the num_species_way fasta format alignment file downloaded 
    from UCSC genome browser,and concatenate the protein sequences.
    This function is independent on the number of speices input.
    parameters:
    handle: opened FASTA file
    fo: opened output file"""

    isFirstBlock = True
    alignment = {}
    species_list = []
    current_geneID = ''
    i = 0
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
                if i % len(species_list) != 0:
                    print("Error: {0} block size is not same the species number {1}".format(current_geneID, str(len(species_list))))
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


def concatenate(handle, num_species, fo):
    """Scans the num_species_way fasta format alignment file downloaded 
    from UCSC genome browser,and concatenate the protein sequences.
    parameters:
    handle: opened FASTA file
    num_species: number of species in the alignment
    fo: opened output file"""

    i = 0
    protein_list = []
    
    for record in SeqIO.parse(handle, 'fasta'):
        fasta_descript, fasta_seq = record.id, record.seq
        #get the accessionID_speciesID for each record
        ucsc_id_species = '_'.join(fasta_descript.split('_')[:-2])

        #for error check, get the cds_index and total number.
        cds_index = int(fasta_descript.split('_')[-2])
        cds_total = int(fasta_descript.split('_')[-1])        
        
        #set up j to track hg19 line, where j==0
        j = i % num_species
        #initialize the protein_list as going through the first block 
        #of alignment.
        if len(protein_list) != num_species:
            #count cds fragments
            cds_counter = 1
            protein_list.append([ucsc_id_species, fasta_seq, cds_index, \
                                cds_total, cds_counter])
        #update the list element after the first block of alignment
        else:
            # write out access_id, ucsc_species, concatenated sequence, and
            # renew the list when new ucsc_id_species occurrs, and restart the concatenation.
            if protein_list[j][0] != ucsc_id_species:
                ucsc_id = '_'.join(protein_list[j][0].split('_')[:-1])
                species = protein_list[j][0].split('_')[-1]
                fo.write('"{0}","{1}","{2}"\n'.format(ucsc_id, species, protein_list[j][1]))
                #error check, to see the cds total number adds up the same.
                if protein_list[j][4] != protein_list[j][3]:
                    print("Error code:1; Record id:{0}; UCSC Species ID:{1}; \
                            Fragment counter:{2}; Fragment number:{3}\n".format(str(i), \
                            ucsc_id_species, str(protein_list[j][4]), cds_total))
                #reset cds_counter to 1
                cds_counter = 1
                protein_list[j] = ([ucsc_id_species, fasta_seq, cds_index, \
                                    cds_total, cds_counter])
            #If the ucsc_id_species keeps the same, then keep the concatenation.
            else:
                protein_list[j][1] += fasta_seq
                #increment the cds counter
                protein_list[j][4] += 1
                #error check, to see the fragements indces are ordered.
                if protein_list[j][4] != cds_index:
                    print("Error code:2; Record id:{0}; UCSC Species ID:{1}; \
                            Fragment counter:{2}; Fragment index:{3}\n".format(str(i), \
                            ucsc_id_species, str(protein_list[j][4]), cds_index))
                    #correct the tracker according to the cds_index in 
                    #order to avoid the accumulating error report.
                    protein_list[j][4] = cds_index
        i += 1
    #write out the last access_id's concatenated alignment 
    for j in range(num_species):
        ucsc_id = '_'.join(protein_list[j][0].split('_')[:-1])
        species = protein_list[j][0].split('_')[-1]
        fo.write('"{0}","{1}","{2}"\n'.format(ucsc_id, species, protein_list[j][1]))


def match_cds_ranges(fi_ref, fi_query, fo):
    """Match the CDS alignment location between knownGene.exonAA.100way.fa with refGene.exonAA.100way.fa, 
    where the former one uses ucsc accession id, and the latter one uses RefSeq. The matched record will
    confirm which ucsc id should be used for duplicated RefSeq in the alignments."""

    header = "quiry_access_id\tref_access_id\tcds_range"
    fo.write("{}\n".format(header))

    dict_ref = {}
    for line in fi_ref:
        ref_access_id = line.strip().split('\t')[0]
        ref_cds_range = line.strip().split('\t')[1]
        dict_ref.update({ref_cds_range:ref_access_id})

    for line in fi_query:
        access_id = line.strip().split('\t')[0]
        cds_range = line.strip().split('\t')[1]
        if cds_range in dict_ref.keys():
            ref_access_id = dict_ref.pop(cds_range)
        else:
            ref_access_id = 'Null'
        fo.write("{0}\t{1}\t{2}\n".format(access_id, ref_access_id, cds_range))


def match_2merge(fi_ref, fi_query, fo):
    """Compare query file's cds range set and RefSeq ID with reference file and find the match.
    Record the match types: ID match, cds range match or Perfect match."""

    dict_ref = {}
    for line in fi_ref:
        ucsc_id = line.strip().split('\t')[0]
        refseq_cds_ranges = '.'.join(line.strip().split('\t')[-2:])
        if refseq_cds_ranges in dict_ref.keys():
            (dict_ref[refseq_cds_ranges]).append(ucsc_id)
        else:
            dict_ref[refseq_cds_ranges] = [ucsc_id] 

    for line in fi_query:
        refseq_cds_ranges = '.'.join(line.strip().split('\t'))
        if refseq_cds_ranges in dict_ref.keys():
            match_ucsc_id = dict_ref[refseq_cds_ranges]
            fo.write("{0}\t{1}\n".format(refseq_cds_ranges, str(match_ucsc_id)))
        else:
            fo.write("{}\tnull\n".format(refseq_cds_ranges))


def check_duplicate(fi_ref, fi_query, num_species):
    """Read RefSeq/Exon coordinates/UCSC ID table, take the matched UCSC ID with 2+ hits, and iterate
    the original knownGene_exonAA alignment, and ckeck those UCSC ID's alignment descriptions are same or not.
    The description includes assembly, size of the assembly, location, size, corrdinates. 
    input: work_1b/knownGene_match_ucsc_id-refseq-gene_symbl-exon_range.out, alignments/knownGene.exonAA.100way.fa""" 

    #read in list of ucsc_id with 2+ members.
    dup_ucsc_id_lists = []
    dup_ucsc_ids = []
    for line in fi_query:
        #convert the string "['uc001bwr.3', 'uc001bws.3', 'uc009vue.3']"" into a list
        ucsc_id_list = line.strip().split('\t')[1].translate(None, "[]''").split(', ')
        if len(ucsc_id_list) > 1:
            #for duplicated ucsc_id, pool all the ucsc_ids, and collect all the ucsc_id lists.
            dup_ucsc_ids += ucsc_id_list
            dup_ucsc_id_lists.append(ucsc_id_list)


    #process the alignment lines. Concatenate all lines for the same ucsc_id, and create dictionary.
    i = 0
    ucsc_id_descript_dict = {}
    descriptions = '' 
    current_access_id = 'nil'
    for record in SeqIO.parse(fi_ref, 'fasta'):
        j = i % num_species

        if j == 0:
            #for every hg19 record (when j==0), extract ucsc_id and check if it is in dup_ucsc_ids
            ucsc_id = record.description.split('_')[0]
            # print '_'.join(record.description.split('_')[1:])
            if ucsc_id in dup_ucsc_ids:
                #only start to process when ucsc_id with duplicate is found.
                #allow the following (num_species - 1) lines to be processed.
                ticket = [True] * (num_species - 1)                
                # print ucsc_id
                #trim off the gene id part for the fasta_descrip
                fasta_descript = '_'.join(record.description.split('_')[1:])
                #append new item for ucsc_id_descript_dict, and renew current_access_id, description 
                #when alignment gene changes. Fist item add will be {'':''}
                if current_access_id != ucsc_id:
                    ucsc_id_descript_dict.update({current_access_id:descriptions})
                    current_access_id = ucsc_id
                    descriptions = fasta_descript
                #concatenate the fasta_descript for the same ucsc_id
                else:
                    descriptions += fasta_descript
            else:
                #disallow the following (num_species - 1) lines to be processed.
                ticket = [False] * (num_species - 1)
        #for the lines are not hg19, check the ticket to see if it is allowed to process.
        else:
            if ticket.pop():
                #trim off the gene id part for the fasta_descrip
                ucsc_id = record.description.split('_')[0]
                fasta_descript = '_'.join(record.description.split('_')[1:])
                descriptions += fasta_descript

        i += 1
    #write in the last instance of current_access_id:descriptions
    ucsc_id_descript_dict.update({current_access_id:descriptions})

    #extract all the concatenated line (value of dictionary) into a list, and compare the duplicity 
    #for the itmes in the new list.    
    for ucsc_id_list in dup_ucsc_id_lists:
        description_list = [ucsc_id_descript_dict.pop(ucsc_id) for ucsc_id in ucsc_id_list]
        if description_list[1:] != description_list[:-1]:
            print("Mismatch is found in {}".format(str(ucsc_id_list)))
        else:
            print("Passed check for {}".format(str(ucsc_id_list)))


def get_range(csv_file, fo):
    """Read in AVA dumped data set, and get the gene symbol with the widest coverage of the gene
    from all the mRNA ranges for the same gene, by taking the lowest start and largest end."""
    import csv

    gene_range_dict = {}
    with open(csv_file, 'rU') as csvfi:
        csv_handle = csv.reader(csvfi, delimiter = ',', quotechar = '"')
        for record in csv_handle:
            mrna_ranges = record[4]
            #omit the genes with no defined mRNA range
            if mrna_ranges == 'NULL':
                pass
            else:
                #get the gene symbol and chromosome ID
                gene_symbl, chr_id = record[0], record[3]
                #extract the start and end coordinates.
                mrna_coord = mrna_ranges.split(',')
                mrna_start, mrna_end = int(mrna_coord[0].strip('[')), int(mrna_coord[-1].strip(']'))
                
                key = "{0}\tchr{1}".format(gene_symbl, chr_id)
                value = [mrna_start, mrna_end]
                #update the start and end for the most outward coordinates.
                if gene_symbl not in gene_range_dict.keys():
                    gene_range_dict.update({key:value})
                else:
                    if mrna_start < gene_range_dict[gene_symbl][0]:
                        gene_range_dict[gene_symbl][0] = mrna_start
                    if mrna_end > gene_range_dict[gene_symbl][1]:
                        gene_range_dict[gene_symbl][1] = mrna_end

    gene_range_flank_dict = {}
    for key in gene_range_dict:
        start = gene_range_dict[key][0] - 2000
        if start < 0:
            start = 0
        end = gene_range_dict[key][1] + 2000

        fo.write("{0}\t{1}\t{2}\n".format(key, str(start), str(end)))


    # gene_size_dict = {}
    # for key in gene_range_dict:
    #     gene_size_dict[key] = gene_range_dict[key][1] - gene_range_dict[key][0]

    # sorted_gene_size = sorted(gene_size_dict.items(), key=operator.itemgetter(1))
    # for item in sorted_gene_size:
    #     print item[0], item[1]

def maf_iterate(handle, gene_ranges, fo1, fo2, block_id=1):
    """Iterates over lines in a MAF file. Extract interested range of sequences, 
    with species, chromosome, block coordinates, and sign an incremented fragment
    id to the block."""


    print gene_ranges
    # gets the first gene range, and remove it from the gene_ranges.
    gene_range = gene_ranges.pop(0)
    #ucsc is 0 indexing basis, AVA is 1 indexing basis, convert AVA index to ucsc.
    gene_start, gene_end = gene_range[0] - 1, gene_range[1] - 1
    
    is_match = False
    in_a_block = False

    while True:
        # allows parsing of the last block
        try:
            line = handle.next()
        except StopIteration:
            line = ""

        if in_a_block:
            if line.startswith('s'):
            #process the sequence 's' line, where all informations are presented
                line_split = line.strip().split()
                if len(line_split) != 7:
                    raise ValueError("Error parsing alignment - 's' line must have 7 fields")

                species, chromosome = line_split[1].split('.')[0], line_split[1].split('.')[1]    
                strand, sequence = line_split[4], line_split[6] 
            
                #check block coordinates of hg19 matches gene range.
                if species == "hg19":
                    block_start, block_size = int(line_split[2]), int(line_split[3])           
                    block_end = block_start + block_size - 1
                    if block_start <= gene_end and block_end >= gene_start:
                        #remove the first entry of gene ranges
                        is_match = True
                        fo1.write("{0},{1},{2},{3}\n".format(str(block_id), chromosome.strip('chr'), str(block_start), str(block_end)))
                    else:
                        is_match = False
                        if block_end > gene_end:
                            # iteration passes the matched block, get next gene range by pop(0)
                            if len(gene_ranges) > 0:
                                gene_range = gene_ranges.pop(0)
                                #ucsc is 0 indexing basis, AVA is 1 indexing basis, convert AVA index to ucsc.
                                gene_start, gene_end = gene_range[0] - 1, gene_range[1] - 1
                            else:
                                break
                if is_match:
                    fo2.write("{0},{1},{2},{3}\n".format(str(block_id), species, strand, sequence))
            elif line.startswith("e") or \
                line.startswith("i") or \
                line.startswith("q"):
                # not implemented
                pass
            #empty line as break of block
            elif not line.strip():
                in_a_block = False
            else:
                raise ValueError("Error parsing alignment - unexpected line:\n%s" % (line,))
        elif line.startswith("a"):
            # start a bundle of records
            in_a_block = True
            if is_match:
                block_id += 1
            if len(line.strip().split()[1:]) != line.count("="):
                raise ValueError("Error parsing alignment - invalid key in 'a' line")
        elif line.startswith("#"):
            # ignore comments
            pass
        elif not line:
            break
    return block_id


def chr_iterate(bed_input, out_dir, block_id, start_chr, end_chr):
    """iterates over all the chromosomes from hg19. Read in bed_input file "chr# gend_start gene_end gene_id" 
    and gzip maf files; outputs block range file and block sequence file for each chromosome.
    block_id will be defined as the starting block_id, chr_start will be the first chromosome to be processed, it 
    will be converted to list index for all chromosomes, end_chr will be included into the processing."""
    import os
    import gzip

    #gets the directory of the script
    path = os.path.dirname(os.path.abspath(__file__))
    #gets directory that the chr*.maf.gz sit.
    chr_path = '/' + path.strip('/tools') + '/ref/multiz100way/maf/'


    bed_list = []
    for line in bed_input:
        if line.strip():
            line_split = line.strip().split()
            chr_index = line_split[0]
            gene_start = int(line_split[1])
            gene_end = int(line_split[2])
            bed_list.append((gene_start, gene_end, chr_index))
    

    #note that AVA call mitochondra chromosome as chrMT, ucsc calls it chrM.
    chrs = ['chr' + str(x) for x in range(1,23)] + ['chrX', 'chrY', 'chrMT']
    # chrs = ['chr' + str(x) for x in range(22,23)]
    #range of chromosomes, X will be 23, and Y will be 24, and MT will be 25.
    for chr in chrs[start_chr -1 : end_chr]:
        print chr
        #gets the gene ranges from bed_list for the iterated chr
        gene_ranges = [x[:2] for x in bed_list if x[2] == chr]
        #sorts the list from smallest start to largest.
        gene_ranges.sort(key=lambda x: x[0])

        fo1_file = 'multiz100way_genome_alignment_coord' + chr + '.csv'
        fo2_file = 'multiz100way_genome_alignment_seq' + chr + '.csv'
        fo1 = open(os.path.join(out_dir, fo1_file), 'w')
        fo2 = open(os.path.join(out_dir, fo2_file), 'w')
        
        maf_gz_file = chr + ".maf.gz"
        handle = gzip.open(os.path.join(chr_path, maf_gz_file), 'rb')
        #run the maf_iterate on the iterated chr, and get the ending block_id for next iteration.
        block_id = maf_iterate(handle, gene_ranges, fo1, fo2, block_id)
        fo1.close()
        fo2.close()


if __name__ == "__main__":
    import sys

    bed_input = open(sys.argv[1], 'rU')
    out_dir = str(sys.argv[2])
    block_id = int(sys.argv[3])
    start_chr = int(sys.argv[4])
    end_chr = int(sys.argv[5])
    chr_iterate(bed_input, out_dir, block_id, start_chr, end_chr)

# if __name__ == "__main__":
#     import sys
#     import csv
#     import operator 

#     csv_file = sys.argv[1]
#     fo = open(sys.argv[2], 'w')
#     get_range(csv_file, fo)
#     fo.close()