#!/usr/bin/python
__author__ = "Xiaoju (Ju) Zhang"

import alignment_parser
import sys


fi_ref = open(sys.argv[1], 'rU')
fi_query = open(sys.argv[2], 'rU')


alignment_parser.check_duplicate(fi_ref, fi_query, 100)
fi_ref.close()
fi_query.close()

