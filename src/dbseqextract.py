#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Corinne Maufrais
# Institut Pasteur, Centre d'informatique pour les biologistes
# corinne.maufrais@pasteur.fr
#
#

# version 1.0

import sys
import os
import argparse
import re

# ############### ENV


class GoldenError:
    def __init__(self, err):
        self.err = err

    def __repr__(self):
        return "[GoldenError] " + self.err

try:
    import Golden
except ImportError, err:
    print >>sys.stderr, GoldenError("%s\n Install the mandatory program golden (https://github.com/C3BI-pasteur-fr/golden)" % err)
    sys.exit(1)

try:
    GOLDENDATA = os.environ['GOLDENDATA']
except:
#    # GOLDENDATA = "/local/gensoft2/exe/golden/1.1a/share/golden/db/"
#    # GOLDENDATA = "/mount/banques/prod/index/golden/"
    print >>sys.stderr, GoldenError('Set the mandatory GOLDENDATA environment variable. Consult https://github.com/C3BI-pasteur-fr/golden.')
    sys.exit(1)

# ############### Golden

def doGoldenAndParse(db, ac):
    # ########################## db ref
    if db in['sp', 'sw', 'swissprot', 'tr', 'trembl']:
        db = 'uniprot'
    elif db in ['emb', 'dbj']:
        db = 'embl'
    elif db in ['gb']:
        db = 'genbank'
    elif db in ['gp']:
        db = 'genpept'
    elif db in ['ref']:
        db = 'refseq'
    elif db in ['rdpii']:
        db = 'rdpii'
    elif db[0:8] == 'embl_wgs':
        db = 'embl_wgs'
    elif db[0:8] == 'genbank_wgs':
        db = 'genbank_wgs'
    elif db in ['pir', 'pdb', 'tpg', 'tpe', 'tpd', 'prf']:
        return '', '', '', ''
    try:
        flatFile = Golden.access(db, ac)
    except IOError, err:
        print >>sys.stderr, err, db, ac
        sys.exit()
    sequence = ''
    if flatFile:
        sequence = parse(flatFile)  # orgName, taxId, taxoLight, description
        flatFile = ''  # buffer free
    return sequence
    
# ############### Parser

def parseUniprot(flatFile):
    """
    parse uniprot or embl like flat file
    """
    sequence = ''
    lineFld = flatFile.split('\n')
    vuSQ = False
    for line in lineFld:
        if line == '\n':
            continue
        tag = line[0:2].strip()
        if tag == 'SQ':
            vuSQ = True
        elif vuSQ:
            sequence += re.sub('\d', '', line).strip().replace(' ', '')
    return sequence


def parseGenbank(flatFile):
    """
    parse genbank like flat file
    """
    sequence = ''
    vuSQ = False
    lineFld = flatFile.split('\n')
    for line in lineFld:
        if line == '\n':
            continue
        tag = line[0:12].strip()
        if tag == 'ORIGIN':
            vuSQ = True
        elif vuSQ:
            sequence += re.sub('\d', '', line).strip().replace(' ', '')
    return sequence


def parse(flatFile):
    """
    parse db flat file (uniprot, embl, genbank) into a DBRecord object
    """
    if flatFile[:2] == 'ID':
        return parseUniprot(flatFile)
    elif flatFile[:5] == 'LOCUS':
        return parseGenbank(flatFile)
    return ''


# ############### 




def offset_extract(infh):
    """
    parse kronaextract file (Query ID    Offset)
    return a containers 
    """
    line = infh.readline()
    containers = {}
    while line:
        fld = line.strip().split()
        containers[fld[0]] = int(fld[1])
        
        line = infh.readline()
    return containers


def write_fasta(fhout, comment, sequence, nb_pb_in_line):
    """
    write fasta file
    """
    nb_car = len(sequence)

    print_f = '>%s\n' % (comment)

    i = 0
    while i < nb_car:
        if i < nb_pb_in_line:
            # print_f += 'C_albi     ' + ref_snps_reference_sequence[i:i + nb_pb_in_line] + '\n'
            print_f += sequence[i:i + nb_pb_in_line] + '\n'
        else:
            print_f += sequence[i:i + nb_pb_in_line] + '\n'
        i += nb_pb_in_line
    print >>fhout, print_f[:-1]



def deb_extract(infh, outfh, containers):
    """
    Extract a list of sequences form a database and rewrite a part DNA sequence (HSP start/stop) in file
    """
    for posline in containers.values():
        infh.seek(posline)
        blast_line = infh.readline()
        fld = blast_line.split('\t')
        # Fields: Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score
        
        query_id, subject_id = fld[0], fld[1]
        q_start, q_end = fld[6], fld[7]
        s_start, s_end = int(fld[8]), int(fld[9])
        
        acc, db = subject_id_analyser(subject_id)
        subject_sequence = doGoldenAndParse(db, acc)
        part_subject_sequence = subject_sequence[s_start -1, s_end]
        write_fasta(outfh, '%s [%s:%s] (%s [%s:%s])' % (subject_id, s_start, s_end, query_id, q_start, q_end), part_subject_sequence, 60)


def subject_id_analyser(subject_id):
    """
    Return accession number and database's name extracted from a blast's subject_id
    """
    fld = subject_id.split('|')
    acc = ''
    db = ''
    if len(fld) == 5:
        db = fld[2]
        acc = fld[3].split('.')[0]
    elif len(fld) == 2 or len(fld) == 3:
        db = fld[0]
        if len(fld) == 3 and fld[1] == '':
            acc = fld[2].split('.')[0]
        else:
            acc = fld[1].split('.')[0]
    return acc, db



if __name__ == '__main__':

    usage = "dbseqextract [options] -i <FILE> "
    epilog = """
    dbseqextract extract fasta sequences from databases.
    dbseqextract is part of the taxoptimizer suite. It reports the fasta sequences of a list of taxon given by the kronaextract program.
    
    """
    parser = argparse.ArgumentParser(prog='dbseqextract.py',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter, usage=usage, epilog=epilog)
    general_options = parser.add_argument_group(title="Options", description=None)
    general_options.add_argument("-i", "--in",
                                 dest="krona_extract_file",
                                 help="kronaextract ouptut file (2 column with database's ID and offset number)",
                                 metavar="file",
                                 type=file,
                                 required=True)
    general_options.add_argument("-j", "--taxo_in",
                                 dest="taxoptimizer_file",
                                 help="taxoptimizer or blast report where offset positions have been found (taxoptimizer in/output file)",
                                 metavar="file",
                                 type=file,
                                 required=True)
    general_options.add_argument("-o", "--out", dest="outfile",
                                 help="Output file (fasta format).",
                                 type=argparse.FileType('w'),
                                 metavar="file",
                                 required=True)
    args = parser.parse_args()

    containers =  offset_extract(args.krona_extract_file)
    deb_extract(args.taxoptimizer_file, args.outfile, containers)