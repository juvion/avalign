#!/usr/bin/python
__author__ = "Xiaoju (Ju) Zhang"

import alignment_parser
import sys


alignment_fasta = sys.argv[1]
fo = open(sys.argv[2], 'w')
handle = open(alignment_fasta, 'rU')
alignment_parser.unique_ids(handle, 100, fo)
fo.close()
