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
import genoGd



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
    sys.exit(1)


######################## analyse krona utput and golden


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


def parseUniprotSeq(flatFile, file_format):
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
        sequence = '>%s %s\n' % (ID, DE)
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


def parseGenbankSeq(flatFile, file_format):
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
        sequence = '>%s %s\n' % (ID, DE)
        i = 0
        while i < len(SEQ):
            sequence += SEQ[i:i+80] + '\n'
            i += 80
        return sequence[:-1]
    elif file_format == 'raw':
        return SEQ


def format_golden_out(flatFile, file_format):
    """
    parse db flat file (uniprot, embl, genbank)
    """
    if flatFile[:2] == 'ID':
        return parseUniprotSeq(flatFile, file_format)
    elif flatFile[:5] == 'LOCUS':
        return parseGenbankSeq(flatFile, file_format)
    return ''


############################   plot data   ################################

def _blastScoreVsColor(val):
    if 0 <= val < 50:
        return 'blue'
    elif 50 <= val < 100:
        return 'grey'
    elif 100 <= val < 150:
        return 'orange'
    elif 150 <= val < 200:
        return 'green'
    elif 200 <= val:
        return 'red'
    else:
        return 'black'


def _blastScoreVsColorInv(val):
    if 0 <= val < 50:
        return 'lightblue'
    elif 50 <= val < 100:
        return 'lightgrey'
    elif 100 <= val < 150:
        return 'yellow'
    elif 150 <= val < 200:
        return 'lightgreen'
    elif 200 <= val:
        return 'lightred'
    else:
        return 'black'



def plot(pictureFileName, sbjctLen, sbjctAcc, blHits, sizeX=700, fontsize=4, typePlot='png'):
    """ picture output """
    # blHits: list of split (blast m8 line)
    # blHits = [qry, sbjct, % identity, aln length, mismatches, gap, q. start, q. end, s. start, s. end, e-value, bit score]
    gdPlot = genoGd.genoGD(fontsize=fontsize, sizeX=sizeX, nbQuery=len(blHits))
    gdPlot.plotSeq(seqlen=sbjctLen, start=0, stop=sizeX, header=sbjctAcc, color='red', up=True)
    gdPlot._nbligneP()
    rpt = float(sizeX) / sbjctLen
    nbhit = 0
    for blh in blHits:
        head = True
        # if len(blLine[0])<= 10:
        #    header = blLine[0]
        # else:
        #    header = blLine[0][0:10]+'..'
        header = blh[0][8:18]+'..'
        if int(blh[8]) < int(blh[9]):
            gdPlot.plotHit(int(int(blh[8]) * rpt), int(int(blh[9]) * rpt), header, _blastScoreVsColor(float(blh[10])), delta=2, head=head)
            nbhit += 1
            gdPlot._nbligneP()
        else:
            gdPlot.plotHit(int(int(blh[9]) * rpt), int(int(blh[8]) * rpt), header, _blastScoreVsColorInv(float(blh[10])), delta=2, head=head)
            nbhit += 1
            gdPlot._nbligneP()
        head = False

    gdPlot._nbligneP()
    gdPlot.plotSeq(seqlen=sbjctLen, start=0, stop=sizeX, header=sbjctAcc, color='red', up=False)
    gdPlot.legend()
    gdPlot.legendInv()
    if typePlot == 'png' and nbhit > 0:
        picturefh = open(pictureFileName, "w")
        gdPlot.picture.writePng(picturefh)
        picturefh.close()


def plotMerge(pictureFileName, sbjctLen, sbjctAcc, blHits, sizeX=700, fontsize=4, typePlot='png'):
    """ picture output """
    # blHits: list of split (blast m8 line)
    # blHits = [qry, sbjct, % identity, aln length, mismatches, gap, q. start, q. end, s. start, s. end, e-value, bit score]
    gdPlot = genoGd.genoGD(fontsize=fontsize, sizeX=sizeX, nbQuery=1)
    gdPlot.plotSeq(seqlen=sbjctLen, start=0, stop=sizeX, header=sbjctAcc, color='red', up=True)
    gdPlot._nbligneP()
    gdPlot._nbligneP()
    rpt = float(sizeX) / sbjctLen
    nbhit = 1
    for blh in blHits:
        head = False
        header = ''
        if int(blh[8]) < int(blh[9]):
            gdPlot.plotHit(int(int(blh[8]) * rpt), int(int(blh[9]) * rpt), header, 'lightred', delta=2, head=head)
        else:
            gdPlot.plotHit(int(int(blh[9]) * rpt), int(int(blh[8]) * rpt), header, 'blue', delta=2, head=head)
        head = False

    gdPlot._nbligneP()
    # gdPlot.plotSeq(seqlen=sbjctLen, start=0, stop=sizeX, header=sbjctAcc, color='red', up = False)
    # gdPlot._nbligneP()
    gdPlot.legendMerge()
    if typePlot == 'png' and nbhit > 0:
        picturefh = open(pictureFileName, "w")
        gdPlot.picture.writePng(picturefh)
        picturefh.close()


def plotMergeFewPlus(pictureFileName, sbjctLen, sbjctAcc, blHits, sizeX=700, fontsize=4, typePlot='png'):
    """ picture output """
    # blHits: list of split (blast m8 line)
    # blHits = [qry, sbjct, % identity, aln length, mismatches, gap, q. start, q. end, s. start, s. end, e-value, bit score]
    gdPlot = genoGd.genoGD(fontsize=fontsize, sizeX=sizeX, nbQuery=len(blHits)/60)
    gdPlot.plotSeq(seqlen=sbjctLen, start=0, stop=sizeX, header=sbjctAcc, color='red', up=True)
    gdPlot._nbligneP()
    gdPlot._nbligneP()
    rpt = float(sizeX) / sbjctLen
    nbhit = 1
    i = 0
    nmax = 0
    while i < len(blHits):
        head = False
        header = ''
        end = int(blHits[i][9])
        gdPlot.plotHit(int(int(blHits[i][8]) * rpt), int(int(blHits[i][9]) * rpt), '', _blastScoreVsColor(float(blHits[i][11])), delta=2, head=False)
        n = 0
        while (i < len(blHits)-1) and (int(blHits[i+1][8]) < end):
            gdPlot._nbligneP()
            gdPlot.plotHit(int(int(blHits[i+1][8]) * rpt), int(int(blHits[i+1][9]) * rpt), header, _blastScoreVsColor(float(blHits[i+1][11])), delta=2, head=head)
            i += 1
            n += 1
            if n > nmax:
                nmax = n
        gdPlot.nbligne -= n
        i += 1

    gdPlot.nbligne += nmax
    # gdPlot.plotSeq(seqlen=sbjctLen, start=0, stop=sizeX, header=sbjctAcc, color='red', up = False)
    # gdPlot._nbligneP()
    gdPlot.legendMerge()
    if typePlot == 'png' and nbhit > 0:
        picturefh = open(pictureFileName, "w")
        gdPlot.picture.writePng(picturefh)
        picturefh.close()


def plotMergeFewMinus(pictureFileName, sbjctLen, sbjctAcc, blHits, sizeX=700, fontsize=4, typePlot='png'):
    """ picture output """
    # blHits: list of split (blast m8 line)
    # blHits = [qry, sbjct, % identity, aln length, mismatches, gap, q. start, q. end, s. start, s. end, e-value, bit score]
    gdPlot = genoGd.genoGD(fontsize=fontsize, sizeX=sizeX, nbQuery=len(blHits)/60)
    gdPlot.plotSeq(seqlen=sbjctLen, start=0, stop=sizeX, header=sbjctAcc, color='red', up=True)
    gdPlot._nbligneP()
    gdPlot._nbligneP()
    rpt = float(sizeX) / sbjctLen
    nbhit = 1
    i = len(blHits)-1
    nmax = 0
    while i >= 0:
        head = False
        header = ''
        begin = int(blHits[i][8])
        gdPlot.plotHit(int(int(blHits[i][9]) * rpt), int(int(blHits[i][8]) * rpt), '', _blastScoreVsColorInv(float(blHits[i][11])), delta=2, head=False)
        n = 0
        while (i > 1) and (int(blHits[i-1][9]) < begin):
            gdPlot._nbligneP()
            gdPlot.plotHit(int(int(blHits[i-1][9]) * rpt), int(int(blHits[i-1][8]) * rpt), header, _blastScoreVsColorInv(float(blHits[i-1][11])), delta=2, head=head)
            i -= 1
            n += 1
            if n > nmax:
                nmax = n
        gdPlot.nbligne -= n
        i -= 1
    gdPlot.nbligne += nmax
    # gdPlot.plotSeq(seqlen=sbjctLen, start=0, stop=sizeX, header=sbjctAcc, color='red', up = False)
    # gdPlot._nbligneP()
    gdPlot.legendMerge()
    if typePlot == 'png' and nbhit > 0:
        picturefh = open(pictureFileName, "w")
        gdPlot.picture.writePng(picturefh)
        picturefh.close()


def plotMergeFew(pictureFileName, sbjctLen, sbjctAcc, blHits, sizeX=700, fontsize=4, typePlot='png'):
    """ picture output """
    # blHits: list of split (blast m8 line)
    # blHits = [qry, sbjct, % identity, aln length, mismatches, gap, q. start, q. end, s. start, s. end, e-value, bit score]
    gdPlot = genoGd.genoGD(fontsize=fontsize, sizeX=sizeX, nbQuery=len(blHits)/60)
    gdPlot.plotSeq(seqlen=sbjctLen, start=0, stop=sizeX, header=sbjctAcc, color='red', up=True)
    gdPlot._nbligneP()
    gdPlot._nbligneP()
    rpt = float(sizeX) / sbjctLen
    nbhit = 1
    i = 0
    nmax = 0
    while i < len(blHits):
        head = False
        header = ''
        if int(blHits[i][8]) < int(blHits[i][9]):
            end = int(blHits[i][9])
            gdPlot.plotHit(int(int(blHits[i][8]) * rpt), int(int(blHits[i][9]) * rpt), '', _blastScoreVsColor(float(blHits[i][11])), delta=2, head=False)
        else:
            end = int(blHits[i][8])
            gdPlot.plotHit(int(int(blHits[i][9]) * rpt), int(int(blHits[i][8]) * rpt), '', _blastScoreVsColorInv(float(blHits[i][11])), delta=2, head=False)
        n = 0
        while (i < len(blHits)-1) and (int(blHits[i+1][8]) < end or int(blHits[i+1][9]) < end):
            gdPlot._nbligneP()
            if int(blHits[i+1][8]) < int(blHits[i+1][9]):
                gdPlot.plotHit(int(int(blHits[i+1][8]) * rpt), int(int(blHits[i+1][9]) * rpt), header, _blastScoreVsColor(float(blHits[i+1][11])), delta=2, head=head)
            else:
                gdPlot.plotHit(int(int(blHits[i+1][9]) * rpt), int(int(blHits[i+1][8]) * rpt), header, _blastScoreVsColorInv(float(blHits[i+1][11])), delta=1, head=head)
            i += 1
            n += 1
            if n > nmax:
                nmax = n
        gdPlot.nbligne -= n
        i += 1

    gdPlot.nbligne += nmax
    # gdPlot.plotSeq(seqlen=sbjctLen, start=0, stop=sizeX, header=sbjctAcc, color='red', up = False)
    # gdPlot._nbligneP()
    gdPlot.legendMerge()
    if typePlot == 'png' and nbhit > 0:
        picturefh = open(pictureFileName, "w")
        gdPlot.picture.writePng(picturefh)
        picturefh.close()



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
        return format_golden_out(flatFile, file_format)  # orgName, taxId
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
        fld = line.split('\t')

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

        if not sbjctAcc or not sbjctDB:
            print >>sys.stderr, MappingError("Parsing error: %s" % fld[blast_column])
            continue

        if sbjctAcc in DBInfo:
            DBInfo[sbjctAcc]['nb'] += 1
            DBInfo[sbjctAcc]['blast'].append(fld[:12])
            #print sbjctAcc, DBInfo[sbjctAcc]['nb']
        else:
            DBInfo[sbjctAcc] = {'nb': 1, 'db': sbjctDB, 'blast': [fld[:12]]}

        if DBInfo[sbjctAcc]['nb'] > maxi:
            maxi = DBInfo[sbjctAcc]['nb']
            subjct_with_max_query = [sbjctAcc]
        elif DBInfo[sbjctAcc]['nb'] == maxi:
            subjct_with_max_query.append(sbjctAcc)

        line = fhin.readline()
    return DBInfo, subjct_with_max_query



def print_subjetcs(fho, DBInfo):
    #for accMax in subjct_with_max_query:
    for acc in DBInfo.keys():
        sbjctSeq = doGolden(DBInfo[acc]['db'], acc, file_format='fasta')
        
        if not sbjctSeq:
            print >>sys.stderr, ('Golden error: %s %s' % (DBInfo[acc]['db'], acc))
            continue
        else:
            print >>fho, sbjctSeq
            #print >>fho, '>%s' % acc
            bls = DBInfo[acc]['blast']
            #bls.sort(lambda x, y: cmp(int(x[8]), int(y[8])))
            #for bl in bls:
            #    print >>fho, bl

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



if __name__ == '__main__':
    usage = ""

    parser = argparse.ArgumentParser(prog='mappinghsp4taxo.py',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter, usage=usage)
    general_options = parser.add_argument_group(title="Options", description=None)

    general_options.add_argument("-i", "--two_column_file", dest="krona_extract_fh",
                                 help="kronaextract output file: query and HSP in two columns", 
                                 type=file,
                                 required=True)
    general_options.add_argument("-o", "--print_seq",
                                 action='store',
                                 dest='fhout',
                                 type=argparse.FileType('w'),
                                 help='Output ')
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

    general_options.add_argument("-p", "--picture_file", dest="picture_prefix",
                                 help="Draw picture(s) with prefix <str>. png format", 
                                 metavar="STRING", required=True)

    general_options.add_argument("-x", "--strand_orientation", dest='strand_orientation',
                                 help="Split mapping acccording to strand's orientation (minus and plus orientations)",
                                action="store_true", default=False)

    general_options.add_argument("-m", "--merge", dest='merge',
                                 help="Draw all hsp in one line",
                                 action="store_true", default=False)

    general_options.add_argument("-M", "--merge_few", dest='merge_few',
                                 help="Draw non overlapping hsp in few lines   --> in developpment",
                                action="store_true", default=False)


    args = parser.parse_args()

    DBInfo, subjct_with_max_query = parse_krona_extract(args.krona_extract_fh, args.blast_column)
    DBInfo = dump_seq_subjetcs(DBInfo)
    #if args.fhout:
    #    print_subjetcs(args.fhout, DBInfo)
    #pickle.dump((DBInfo, subjct_with_max_query), open('subject.dmp', 'w'))

    for subject in subjct_with_max_query:
        plot_name_prefix = args.picture_prefix + '_' + subject
        plus = []
        minus = []
        DBInfo[subject]['blast'].sort(lambda x, y: cmp(int(x[8]), int(y[8])))
        if args.strand_orientation or (args.merge or args.merge_few):
            for bl in DBInfo[subject]['blast']:
                if int(bl[8]) <= int(bl[9]):
                    plus.append(bl)
                else:
                    minus.append(bl)
            if (args.merge or args.merge_few):
                new = []
                i = 0
                j = 0
                minus.sort(lambda x, y: cmp(int(x[9]), int(y[9])))
                while j < len(minus):
                    while i < len(plus)-1 and j < len(minus):
                        if int(minus[j][9]) > int(plus[i][8]):
                            if int(minus[j][9]) > int(plus[i+1][8]):
                                new.append(plus[i])
                                new.append(minus[j])
                                i += 1
                                j += 1
                            else:
                                new.append(plus[i])
                                i += 1
                        else:
                            new.append(minus[j])
                            j += 1
            else:
                minus.reverse()

        sbjctSeq = DBInfo[subject]['seq']
        if not sbjctSeq:
            print >>sys.stderr, ('Golden error: %s %s' % (DBInfo[subject]['db'], subject))
            continue
        sbjctLen = len(sbjctSeq)
        print sbjctLen
        if args.picture_prefix and (not args.merge and not args.merge_few):
            if not args.strand_orientation:
                plot(plot_name_prefix + '.png', sbjctLen, subject, DBInfo[subject]['blast'], sizeX=700, fontsize=8, typePlot='png')
            else:
                plot(plot_name_prefix + '_plus.png', sbjctLen, subject, plus, sizeX=700, fontsize=8, typePlot='png')
                plot(plot_name_prefix + '_minus.png', sbjctLen, subject, minus, sizeX=700, fontsize=8, typePlot='png')

        if args.picture_prefix and args.merge:
            plot_name_prefix += '_merge'

            if not args.strand_orientation:
                plotMerge(plot_name_prefix + '.png', sbjctLen, subject, DBInfo[subject]['blast'], sizeX=700, fontsize=8, typePlot='png')
            else:
                plotMerge(plot_name_prefix + '_plus.png', sbjctLen, subject, plus, sizeX=700, fontsize=8, typePlot='png')
                plotMerge(plot_name_prefix + '_minus.png', sbjctLen, subject, minus, sizeX=700, fontsize=8, typePlot='png')

        if args.picture_prefix and args.merge_few:
            plot_name_prefix += '_noOverlap'
            if not args.strand_orientation:
                plotMergeFew(plot_name_prefix + '.png', sbjctLen, subject, new, sizeX=700, fontsize=8, typePlot='png')
            else:
                plotMergeFewPlus(plot_name_prefix + '_plus.png', sbjctLen, subject, plus, sizeX=700, fontsize=8, typePlot='png')
                plotMergeFewMinus(plot_name_prefix + '_minus.png', sbjctLen, subject, minus, sizeX=700, fontsize=8, typePlot='png')
