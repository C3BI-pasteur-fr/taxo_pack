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
        return "[plotHSPsubject]" + self.err
    
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


##############################################################################

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



if __name__ == '__main__':
    usage = ""

    parser = argparse.ArgumentParser(prog='plotHSPsubject.py',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter, usage=usage)
    general_options = parser.add_argument_group(title="Options", description=None)
    
    general_options.add_argument("-i", "--dump_name", dest="dump_name",
                                 help="dump_name", 
                                 type=str,
                                 required=True)
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

    DBInfo, subjct_with_max_query = pickle.load(open('%s.dmp' % args.dump_name, 'r'))
    
    for subject in DBInfo.keys():
        plot_name_prefix = args.picture_prefix + '/' + subject
        plus = []
        minus = []
        DBInfo[subject]['blast'].sort(lambda x, y: cmp(int(x[8]), int(y[8])))
        if args.strand_orientation: # or (args.merge or args.merge_few):
            for bl in DBInfo[subject]['blast']:
                if int(bl[8]) <= int(bl[9]):
                    plus.append(bl)
                else:
                    minus.append(bl)
            if (plus and minus):
                if  (args.merge or args.merge_few):
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
                        print 'a verifier le j+=1 ajouter ici'
                        j += 1   
                else:
                    minus.reverse()
            else:
                print 'a faire'   
        
        sbjctSeq = DBInfo[subject]['seq']
        if not sbjctSeq:
            print >>sys.stderr, ('Golden error: %s %s' % (DBInfo[subject]['db'], subject))
            continue
        sbjctLen = len(sbjctSeq)
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
