#!/usr/bin/python
__author__ = "Xiaoju (Ju) Zhang"

import alignment_parser
import sys


fi_ref = open(sys.argv[1], 'rU')
fi_query = open(sys.argv[2], 'rU')
fo = open(sys.argv[3], 'w')

alignment_parser.match_cds_ranges(fi_ref, fi_query, fo)
fi_ref.close()
fi_query.close()
fo.close()
