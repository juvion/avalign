#!/bin/python

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

