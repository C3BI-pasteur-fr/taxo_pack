#!/usr/bin/env python
# -*- coding: utf-8 -*-


# Corinne Maufrais
# Institut Pasteur, Centre d'informatique pour les biologistes
# maufrais@pasteur.fr
#

# version 2.1


import os
import sys
import argparse
import pickle
import re

class MappingError:
    def __init__(self, err):
        self.err = err

    def __repr__(self):
        return "[mappinghsp4taxo]" + self.err
    
try:
    LIB = os.environ['RANKOPTIMIZERLIB']
except:
    print >>sys.stderr, MappingError("The mandatory RANKOPTIMIZERLIB environment variable is not defined")
    sys.exit(1)

# try:
#     LIB2 = os.environ['BLAST2TC_LIB']
# except:
#     print >>sys.stderr, MappingError("The mandatory BLAST2TC_LIB environment variable is not defined")
#     sys.exit(1)


import Golden

try:
    GOLDENDATA = os.environ['GOLDENDATA']
except:
    #GOLDENDATA = "/mount/banques/prod/index/golden"
    #os.environ['GOLDENDATA'] = GOLDENDATA
    print >>sys.stderr, MappingError('The mandatory GOLDENDATA environment variable is not defined. Consult https://github.com/C3BI-pasteur-fr/golden.')
    #sys.exit(1)


def parseUniprot(flatFile, DE):
    """
    parse uniprot or embl like flat file
    """
    description = ''
    taxId = ''
    taxoLight = ''
    orgName = ''
    lineFld = flatFile.split('\n')
    vuOS = False
    vuOC = False
    vuOCXX = False
    for line in lineFld:
        if line == '\n':
            continue
        tag = line[0:2].strip()
        value = line[5:]
        value = value.replace('\n', '')
        if DE and tag == 'DE':
            description += value
        elif tag == 'OS' and not vuOS:
            fldOS = value.split('(')
            orgName += fldOS[0].strip()
            vuOS = True
        elif tag == 'OC' and not vuOCXX:
            taxoLight += value
            vuOC = True
        elif tag == 'XX' and vuOC:
            vuOCXX = True
        elif tag == 'OX':
            taxId += value
        elif tag in ['RN', 'DR', 'CC', 'FH', 'SQ']:
            return orgName, taxId, taxoLight, description
    return orgName, taxId, taxoLight, description


def parseUniprotSeq(flatFile, file_format, acc):
    """
    parse uniprot or embl like flat file
    """
    ID = ''
    SEQ = ''
    DE = ''
    seqVu = False
    lineFld = flatFile.split('\n')
    for line in lineFld:
        if line == '\n':
            continue
        tag = line[0:2].strip()
        value = line[5:]
        value = value.replace('\n', '')
        if tag == 'ID':
            ID += value.split()[0]
        elif tag == 'DE':
            DE += value
        elif tag == 'SQ' and not seqVu:
            seqVu = True
        elif seqVu:
            SEQ += value.replace(' ', '')
    if file_format == 'fasta':
        sequence = '>%s %s\n' % (acc, ID)
        i = 0
        while i < len(SEQ):
            sequence += SEQ[i:i+80] + '\n'
            i += 80
        return sequence
    elif file_format == 'raw':
        return SEQ


def parseGenbank(flatFile, DE):
    """
    parse genbank like flat file
    """
    description = ''
    taxId = ''
    taxoLight = ''
    orgName = ''
    taxoL = False
    lineFld = flatFile.split('\n')
    for line in lineFld:
        if line == '\n':
            continue
        tag = line[0:12].strip()
        value = line[12:]
        value = value.replace('\n', '')
        if DE and tag == 'DEFINITION':
            description += value
        elif tag == 'ORGANISM':  # one line ##### pb qqfois plusieurs ligne
            orgName += value.strip()
            taxoL = True
        elif taxoL and not tag:  # multiple line without tag
            taxoLight += value
        elif tag == 'TaxID':  # gi number
            taxId += value
        elif tag in ['REFERENCE', 'COMMENT', 'FEATURES', 'ORIGIN']:
            return orgName, taxId, taxoLight, description
        elif not tag:  # fin taxonomy
            taxoL = False
    return orgName, taxId, taxoLight, description


def parseGenbankSeq(flatFile, file_format, acc):
    """
    parse genbank like flat file
    """
    ID = ''
    SEQ = ''
    DE = ''
    seqVu = False
    lineFld = flatFile.split('\n')
    for line in lineFld:
        if line == '\n':
            continue
        tag = line[0:12].strip()
        value = line[12:]
        value = value.replace('\n', '')
        if tag == 'LOCUS':
            ID += value.split()[0]
        elif tag == 'DEFINITION':
            DE += value
        elif tag == 'ORIGIN' and not seqVu:
            seqVu = True
        elif seqVu:
            SEQ += line[10:].replace(' ', '')
    if file_format == 'fasta':
        sequence = '>%s %s\n' % (acc, ID)
        i = 0
        while i < len(SEQ):
            sequence += SEQ[i:i+80] + '\n'
            i += 80
        return sequence[:-1]
    elif file_format == 'raw':
        return SEQ


def format_golden_out(flatFile, file_format, acc):
    """
    parse db flat file (uniprot, embl, genbank)
    """
    if flatFile[:2] == 'ID':
        return parseUniprotSeq(flatFile, file_format, acc)
    elif flatFile[:5] == 'LOCUS':
        return parseGenbankSeq(flatFile, file_format, acc)
    return ''

##############################################################################
#
#            Golden: Keep this for compatibility and performance testing.
#
##############################################################################

def doGolden(db, ac, file_format='fasta'):
    # ## db ref
    if db in ['sp', 'sw', 'swissprot', 'tr', 'trembl']:
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
    elif db in ['pir','pdb','tpg', 'tpe', 'tpd', 'prf']:
        return ''

    try:
        flatFile = Golden.access(db, ac)
    except IOError, err:
        print >>sys.stderr, err, "%s:%s" % (db, ac)
        sys.exit()

    if flatFile:
        return format_golden_out(flatFile, file_format, ac)  # orgName, taxId
        flatFile = ''  # buffer free
    else:
        return ''


##############################################################################


def parse_krona_extract(fhin, blast_column):
    DBInfo = {}
    maxi = 0
    subjct_with_max_query = []
    line = fhin.readline()
    while line:
        fld = line.strip().split('\t')
        
        query = fld[0].strip()
        try:
            sbjctInfo = fld[blast_column].split(args.separator)
        except:
            if blast_column < len(fld):
                print >>sys.stderr, MappingError("Parsing error: %s" % fld[blast_column])
            else:
                print >>sys.stderr, MappingError("Parsing error: column %s doesn't exist in line nb %s" % (blast_column+1, line))
            continue
        sbjctAcc = ''
        sbjctDB = ''
        if len(sbjctInfo) == 5:
            # PF8/Pathoquest
            sbjctDB = sbjctInfo[2]
            sbjctAcc = sbjctInfo[3].split('.')[0]
        elif len(sbjctInfo) == 2 or len(sbjctInfo) == 3:
            # Lionel extraction
            sbjctDB = sbjctInfo[0]
            sbjctAcc = sbjctInfo[1]
        else:
            sbjctDB = sbjctInfo[0]
            sbjctAcc = sbjctInfo[1]
        if not sbjctAcc or not sbjctDB:
            print >>sys.stderr, MappingError("Parsing error: %s" % fld[blast_column])
            continue

        if sbjctAcc in DBInfo:
            DBInfo[sbjctAcc]['nb'] += 1
            DBInfo[sbjctAcc]['blast'].append(fld)
            #print sbjctAcc, DBInfo[sbjctAcc]['nb']
        else:
            DBInfo[sbjctAcc] = {'nb': 1, 'db': sbjctDB, 'blast': [fld]}

        if DBInfo[sbjctAcc]['nb'] > maxi:
            maxi = DBInfo[sbjctAcc]['nb']
            subjct_with_max_query = [sbjctAcc]
        elif DBInfo[sbjctAcc]['nb'] == maxi:
            subjct_with_max_query.append(sbjctAcc)

        line = fhin.readline()
    return DBInfo, subjct_with_max_query


def parse_krona_extract_with_html(fhin, blast_column):
    DBInfo = {}
    maxi = 0
    subjct_with_max_query = []
    line = fhin.readline()
    while line:
        fld = line.strip().split('\t')
        query = fld[0].strip()
        try:
            hsp = fld[1]
            db_info = re.search('>.*<', hsp).group(0)
            sbjctInfo = db_info[1:-1].split(args.separator)
        except:
            print >>sys.stderr, MappingError("Parsing error: %s" % hsp)
            continue
        sbjctAcc = ''
        sbjctDB = ''
        if len(sbjctInfo) == 5:
            # PF8/Pathoquest
            sbjctDB = sbjctInfo[2]
            sbjctAcc = sbjctInfo[3].split('.')[0]
        elif len(sbjctInfo) == 2 or len(sbjctInfo) == 3:
            # Lionel extraction
            sbjctDB = sbjctInfo[0]
            if '_merge.png' in sbjctInfo[1]:
                sbjctAcc = sbjctInfo[1].replace('_merge.png', '')
            else:
                sbjctAcc = sbjctInfo[1]
        else:
            sbjctDB = sbjctInfo[0]
            sbjctAcc = sbjctInfo[1]
        if not sbjctAcc or not sbjctDB:
            print >>sys.stderr, MappingError("Parsing error: %s" % fld[blast_column])
            continue

        if sbjctAcc in DBInfo:
            DBInfo[sbjctAcc]['nb'] += 1
            DBInfo[sbjctAcc]['blast'].append(fld)
            #print sbjctAcc, DBInfo[sbjctAcc]['nb']
        else:
            DBInfo[sbjctAcc] = {'nb': 1, 'db': sbjctDB, 'blast': [fld]}

        if DBInfo[sbjctAcc]['nb'] > maxi:
            maxi = DBInfo[sbjctAcc]['nb']
            subjct_with_max_query = [sbjctAcc]
        elif DBInfo[sbjctAcc]['nb'] == maxi:
            subjct_with_max_query.append(sbjctAcc)

        line = fhin.readline()
    return DBInfo, subjct_with_max_query



def print_fasta_subjetcs(fhop, fhon, DBInfo):
    #for accMax in subjct_with_max_query:
    for acc in DBInfo.keys():
        sbjctSeq = doGolden(DBInfo[acc]['db'], acc, file_format='fasta')
        
        if not sbjctSeq:
            print >>sys.stderr, ('Golden error: %s %s' % (DBInfo[acc]['db'], acc))
            continue
        elif DBInfo[acc]['db'] in ['sp', 'sw', 'swissprot', 'tr', 'trembl', 'uniprot']:
            print >>fhop, sbjctSeq
        elif DBInfo[acc]['db'] in ['embl', 'emb', 'dbj', 'gb', 'genbank']:
            print >>fhon, sbjctSeq

def dump_seq_subjetcs(DBInfo):
    #for accMax in subjct_with_max_query:
    for acc in DBInfo.keys():
        sbjctSeq = doGolden(DBInfo[acc]['db'], acc, file_format='fasta')
        
        if not sbjctSeq:
            print >>sys.stderr, ('Golden error: %s %s' % (DBInfo[acc]['db'], acc))
            continue
        else:
            DBInfo[acc]['seq'] = sbjctSeq
    return DBInfo


def printANDdump_fasta_subjetcs(fhop, fhon, DBInfo):
    for acc in DBInfo.keys():
        sbjctSeq = doGolden(DBInfo[acc]['db'], acc, file_format='fasta')
        
        if not sbjctSeq:
            print >>sys.stderr, ('Golden error: %s %s' % (DBInfo[acc]['db'], acc))
            continue
        elif DBInfo[acc]['db'] in ['sp', 'sw', 'swissprot', 'tr', 'trembl', 'uniprot']:
            print >>fhop, sbjctSeq
            DBInfo[acc]['seq'] = sbjctSeq
        elif DBInfo[acc]['db'] in ['embl', 'emb', 'dbj', 'gb', 'genbank']:
            print >>fhon, sbjctSeq
            DBInfo[acc]['seq'] = sbjctSeq


if __name__ == '__main__':
    usage = ""

    parser = argparse.ArgumentParser(prog='extractHSPsubject.py',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter, usage=usage)
    general_options = parser.add_argument_group(title="Options", description=None)

    general_options.add_argument("-i", "--two_column_file", dest="krona_extract_fh",
                                 help="kronaextract output file: query and HSP in two columns", 
                                 type=file,
                                 required=True)
    general_options.add_argument("-o", "--outprefix",
                                 dest='outname',
                                 help='Output prefix; split fasta seq in one prot and one dna files',
                                 default=False)
    general_options.add_argument("-c", "--blast_column",
                                 dest="blast_column", metavar="int",
                                 help="Column number to parse in ptaxoptimizer file (default: 2)",
                                 default=1)
    general_options.add_argument("-s", "--separator",
                                 dest="separator", metavar="STRING",
                                 help="Separator character (default '|')",
                                 default='|')
    general_options.add_argument("-d", "--dump_seq",
                                 dest="dump_seq",
                                 action="store_true",
                                 help="dump_seq",
                                 default=False)
    general_options.add_argument("-m", "--html",
                                 dest="html_input",
                                 action="store_true",
                                 help="html_input",
                                 default=False)


    args = parser.parse_args()
    
    if args.html_input:
        DBInfo, subjct_with_max_query = parse_krona_extract_with_html(args.krona_extract_fh, args.blast_column)
    else:
        DBInfo, subjct_with_max_query = parse_krona_extract(args.krona_extract_fh, args.blast_column)

    if args.dump_seq:
        if not args.outname:
            DBInfo = dump_seq_subjetcs(DBInfo)
        else:
            fhop = open('%s.fap' % args.outname, 'w')
            fhon = open('%s.fan' % args.outname, 'w')
            printANDdump_fasta_subjetcs(fhop, fhon, DBInfo)
        pickle.dump((DBInfo, subjct_with_max_query), open('%s.dmp' % args.krona_extract_fh.name, 'w'))
    elif args.outname:
        fhop = open('%s.fap' % args.outname, 'w')
        fhon = open('%s.fan' % args.outname, 'w')
        print_fasta_subjetcs(fhop, fhon, DBInfo)

    

