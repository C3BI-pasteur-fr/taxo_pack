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

try:
    LIB2 = os.environ['BLAST2TC_LIB']
except:
    print >>sys.stderr, MappingError("The mandatory BLAST2TC_LIB environment variable is not defined")
    sys.exit(1)


import Golden
# os.environ['GOLDENDATA'] = "/mount/banques/prod/index/golden"

try:
    GOLDENDATA = os.environ['GOLDENDATA']
except:
    #GOLDENDATA = "/mount/banques/prod/index/golden"
    #os.environ['GOLDENDATA'] = GOLDENDATA
    print >>sys.stderr, MappingError('The mandatory GOLDENDATA environment variable is not defined. Consult https://github.com/C3BI-pasteur-fr/golden.')
    sys.exit(1)


import genoGd





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
        return sequence
    elif file_format == 'raw':
        return SEQ


def parse(flatFile, file_format):
    """
    parse db flat file (uniprot, embl, genbank)
    """
    if flatFile[:2] == 'ID':
        return parseUniprotSeq(flatFile, file_format)
    elif flatFile[:5] == 'LOCUS':
        return parseGenbankSeq(flatFile, file_format)
    return ''


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
    # elif db in ['pir','pdb','tpg', 'tpe', 'tpd', 'prf']:
    #    return '','','',''

    try:
        flatFile = Golden.access(db, ac)
    except IOError:
        return ''

    if flatFile:
        return parse(flatFile, file_format)  # orgName, taxId, taxoLight
    else:
        return ''


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
            gdPlot.plotHit(int(int(blh[8]) * rpt), int(int(blh[9]) * rpt), header, _blastScoreVsColor(float(blh[11])), delta=2, head=head)
            nbhit += 1
            gdPlot._nbligneP()
        else:
            gdPlot.plotHit(int(int(blh[9]) * rpt), int(int(blh[8]) * rpt), header, _blastScoreVsColorInv(float(blh[11])), delta=2, head=head)
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


def plotMerge2plus(pictureFileName, sbjctLen, sbjctAcc, blHits, sizeX=700, fontsize=4, typePlot='png'):
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


def plotMerge2minus(pictureFileName, sbjctLen, sbjctAcc, blHits, sizeX=700, fontsize=4, typePlot='png'):
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


def plotMerge2(pictureFileName, sbjctLen, sbjctAcc, blHits, sizeX=700, fontsize=4, typePlot='png'):
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


if __name__ == '__main__':
    usage = "mappinghsp4taxo [options] -i <file>/-j <file> -b <file>"

    parser = argparse.ArgumentParser(prog='mappinghsp4taxo.py',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter, usage=usage)
    general_options = parser.add_argument_group(title="Options", description=None)

    general_options.add_argument("-i", "--one_column_file", dest="oneofflinefh",
                                 help="Offset's list of line to parse taxoptomizer output file (cf option -b). File contains taxoptimizer's line offset (integer).",
                                 metavar="file",
                                 required=True)
    general_options.add_argument("-j", "--two_column_file", dest="twoofflinefh",
                                 help="kronaextract output file. FIle contains two column; first: reads name and second: taxoptimizer's line offset (integer).",
                                 metavar="file",
                                 required=True)
    general_options.add_argument("-b", "--taxofile", dest="taxofh",
                                 help="[p]taxoptimizer output file (format: blast m8 + OC + taxonomy +/-  DE).",
                                 metavar="file",
                                 required=True)

    general_options.add_argument("-e", "--extract_file", dest="extractfh",
                                 help="Extract lines of interest from [p]taxoptimizer file and write them in <file>",
                                 type=argparse.FileType('w'),
                                 metavar="file")

    general_options.add_argument("-c", "--blast_column",
                                 dest="blast_column", metavar="int",
                                 help="Column number to parse in ptaxoptimizer file (default: 2)", default=1)
    general_options.add_argument("-s", "--separator",
                                 dest="separator", metavar="STRING",
                                 help="Separator character (default '|')",
                                 default='|')

    general_options.add_argument("-p", "--picture_file", dest="picture_prefix",
                                 help="Draw picture(s) with prefix <str>. png format", metavar="STRING", default='')

    general_options.add_argument("-x", "--strand_orientation", dest='strand_orientation',
                                 help="Split mapping acccording to strand's orientation (minus and plus orientations)",
                                 metavar="BOOLEAN", action="store_true", default=False)

    general_options.add_argument("-m", "--merge", dest='merge',
                                 help="Draw all hsp in one line",
                                 metavar="BOOLEAN", action="store_true", default=False)

    general_options.add_argument("-M", "--merge2", dest='merge2',
                                 help="Draw non overlapping hsp in few lines   --> in developpment",
                                 metavar="BOOLEAN", action="store_true", default=False)

    general_options.add_argument("-g", "--galaxy_restrict", dest='galaxy_restrict',
                                 help="Restrict output(s) to one reference sequence",
                                 metavar="BOOLEAN", action="store_true", default=False)

    args = parser.parse_args()

    if args.oneofflinefh:
        try:
            offlinefh = args.oneofflinefh
            offset_column = 0
        except IOError, err:
            print >>sys.stderr, err
            sys.exit()
    elif args.twoofflinefh:
        try:
            offlinefh = args.twoofflinefh
            offset_column = 1
        except IOError, err:
            print >>sys.stderr, err
            sys.exit()

    DBInfo = {}

    maxi = 0
    sbjctAccMax = []
    offline = offlinefh.readline()
    while offline:
        off_field = offline.split('\t')
        if len(off_field) != offset_column + 1:
            print >>sys.stderr, "format error: cf option error: -i/--one_column_file or -j/two_column_file \n"
            sys.exit(1)
        position = int(off_field[offset_column])
        args.taxofh.seek(position)
        blLine = args.taxofh.readline()
        if args.extract_file:
            print >>args.extractfh, blLine[:-1]
        blfld = blLine.split('\t')
        query = blfld[0]
        try:
            sbjctInfo = blfld[args.blast_column].split(args.separator)
        except:
            if args.blast_column < len(blfld):
                print >>sys.stderr, RankOptimizerError("Parsing error: %s" % blfld[args.blast_column])
            else:
                print >>sys.stderr, RankOptimizerError("Parsing error: column %s doesn't exist in line nb %s" % (args.blast_column+1, blLine))
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
            print >>sys.stderr, RankOptimizerError("Parsing error: %s" % blfld[args.blast_column])
            continue

        if sbjctAcc in DBInfo:
            DBInfo[sbjctAcc]['nb'] += 1
            DBInfo[sbjctAcc]['blast'].append(blfld[:12])
        else:
            DBInfo[sbjctAcc] = {'nb': 1, 'db': sbjctDB, 'blast': [blfld[:12]]}

        if DBInfo[sbjctAcc]['nb'] > maxi:
            maxi = DBInfo[sbjctAcc]['nb']
            sbjctAccMax = [sbjctAcc]
        elif DBInfo[sbjctAcc]['nb'] == maxi:
            sbjctAccMax.append(sbjctAcc)

        offline = offlinefh.readline()

    if args.galaxy_restrict:
        # fichiers multiples mal geres par galaxy
        sbjctAccMax = sbjctAccMax[0:1]

    for elemMax in sbjctAccMax:
        if not args.galaxy_restrict:
            plot_name_prefix = args.picture_prefix + '_' + elemMax
        else:
            plot_name_prefix = args.picture_prefix
        plus = []
        minus = []
        DBInfo[elemMax]['blast'].sort(lambda x, y: cmp(int(x[8]), int(y[8])))
        if args.strand_orientation or (args.merge or args.merge2):
            for bl in DBInfo[elemMax]['blast']:
                if int(bl[8]) <= int(bl[9]):
                    plus.append(bl)
                else:
                    minus.append(bl)
            if (args.merge or args.merge2):
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

        sbjctSeq = doGolden(DBInfo[elemMax]['db'], elemMax, file_format='raw')
        if not sbjctSeq:
            print >>sys.stderr, ('Golden error: %s %s' % (DBInfo[elemMax]['db'], elemMax))
            continue
        sbjctLen = len(sbjctSeq)
        # sbjctLen = 2900
        if args.picture_prefix and (not args.merge and not args.merge2):
            if not args.strand_orientation:
                plot(plot_name_prefix + '.png', sbjctLen, elemMax, DBInfo[elemMax]['blast'], sizeX=700, fontsize=8, typePlot='png')
            else:
                plot(plot_name_prefix + '_plus.png', sbjctLen, elemMax, plus, sizeX=700, fontsize=8, typePlot='png')
                plot(plot_name_prefix + '_minus.png', sbjctLen, elemMax, minus, sizeX=700, fontsize=8, typePlot='png')

        if args.picture_prefix and args.merge:
            # if not args.galaxy_restrict:
            plot_name_prefix += '_merge'

            if not args.strand_orientation:
                plotMerge(plot_name_prefix + '.png', sbjctLen, elemMax, DBInfo[elemMax]['blast'], sizeX=700, fontsize=8, typePlot='png')
            else:
                plotMerge(plot_name_prefix + '_plus.png', sbjctLen, elemMax, plus, sizeX=700, fontsize=8, typePlot='png')
                plotMerge(plot_name_prefix + '_minus.png', sbjctLen, elemMax, minus, sizeX=700, fontsize=8, typePlot='png')

        if args.picture_prefix and args.merge2:
            # if not args.galaxy_restrict:
            plot_name_prefix += '_noOverlap'
            if not args.strand_orientation:
                plotMerge2(plot_name_prefix + '.png', sbjctLen, elemMax, new, sizeX=700, fontsize=8, typePlot='png')
            else:
                plotMerge2plus(plot_name_prefix + '_plus.png', sbjctLen, elemMax, plus, sizeX=700, fontsize=8, typePlot='png')
                plotMerge2minus(plot_name_prefix + '_minus.png', sbjctLen, elemMax, minus, sizeX=700, fontsize=8, typePlot='png')
