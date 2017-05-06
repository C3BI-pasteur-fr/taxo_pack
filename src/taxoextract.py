#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Corinne Maufrais
# Institut Pasteur, Centre d'informatique pour les biologistes
# corinne.maufrais@pasteur.fr
#

# version 2.1

import os
import sys
import argparse

try:
    LIB = os.environ['RANKOPTIMIZERLIB']
except:
    LIB = '/usr/local/bin'

if LIB not in sys.path:
    sys.path.append( LIB )

try:
    LIB2 = os.environ['BLAST2TC_LIB']
except:
    LIB2 = '/usr/local/bin'

if LIB2 not in sys.path:
    sys.path.append( LIB2 )



if __name__=='__main__':
    parser = argparse.ArgumentParser(prog='taxoextract.py',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="Extract a list of entries from a taxoptimizer output file")

    usage = "taxoextract [options] -i <file>/-j <file> -b <file>"

    general_options = parser.add_argument_group(title="Options", description=None)

    general_options.add_argument("-i", "--one_column_file", dest="offsetfh",
                                 metavar="File",
                                 type=file,
                                 required=True,
                                 help="List of taxoptomizer's lines to parse (offset numbers given by kronaexract).")
    
    general_options.add_argument("-b", "--taxofile", dest="taxofh",
                                 help="taxoptimizer output file (format: blast m8 + OC + taxonomy +/-  DE ).",
                                 metavar="FILE",
                                 type=file,
                                 required=True,)
    
    general_options.add_argument("-a", "--append", dest="ap_mode",
                                 help="append mode. Result will be append to the output file \
                                 instead of creating a new one (Default: write mode)",
                                 action="store_true", 
                                 default=False,
                                 )
    
    general_options.add_argument("-o", "--out_file", dest="outfile",
                                help="Output file", 
                                metavar="FILE",
                                required=True)
    
    
    args = parser.parse_args()    
    


    if args.ap_mode is True:
        extractfh = open( args.outfile, 'a' )
    else:
        extractfh = open( args.outfile, 'w' )
    
    offline = args.offsetfh.readline()
    while offline:
        off_field = offline.strip()
        try:
            position = int(off_field)
        except:
            print >>sys.stderr, "format error: one column file is expected (int)" 
            sys.exit(1)
        args.taxofh.seek( position )
        blLine = args.taxofh.readline()
        print >>extractfh, blLine[:-1]
        offline = args.offsetfh.readline()
