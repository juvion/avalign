#!/usr/bin/python
__author__ = "Xiaoju (Ju) Zhang"

import alignment_parser
import sys
# import operator 

csv_file = sys.argv[1]
fo = open(sys.argv[2], 'w')
alignment_parser.get_range(csv_file, fo)
fo.close()
