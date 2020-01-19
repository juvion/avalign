#!/usr/bin/python
__author__ = "Xiaoju (Ju) Zhang"

import sys
from Bio import SeqIO
from Bio.Alphabet import Gapped, IUPAC


qseq = sys.argv[1]
handle = open(sys.argv[2], 'rU')

for record in SeqIO.parse(handle, 'fasta'):
    fasta_descript, fasta_seq = record.id, record.seq
    fasta_rc_seq=fasta_seq[::-1]
    if qseq in fasta_seq:
        print("sequence is found in + strand {}".format(fasta_descript))
    if qseq in fasta_rc_seq:
        print("sequence is found in - strand {}".format(fasta_descript))
