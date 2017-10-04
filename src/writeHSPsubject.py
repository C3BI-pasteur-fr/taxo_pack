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
import BlastParser

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



##############################################################################

def parse_blast(fhin, DBInfo):
    blastRecordS, format = BlastParser.parser(fhin)
    ## Check supported database
    if format not in ['m0','m8']:
        print >>sys.stderr, BlastParser.BlastError("Not an accepted blast format" )
        sys.exit()
    for blastRecord in blastRecordS.blastS:
        blastRecord.database = blastRecordS.database
        for hit in blastRecord.hitS:
            if hit.dbInfo in DBInfo.keys():
                if 'blast_m0' not in DBInfo[hit.dbInfo]:
                    DBInfo[hit.dbInfo]['blast_m0'] = {}
                for el_hsp in DBInfo[hit.dbInfo]['blast']:
                    el_hsp_qry =  el_hsp[0].strip()
                    if el_hsp_qry == blastRecord.queryDef:
                        DBInfo[hit.dbInfo]['blast_m0'][el_hsp_qry] = {'all': []}
                        for hsp in hit.hspS:
                            DBInfo[hit.dbInfo]['blast_m0'][el_hsp_qry]['all'].append([hsp.hspNum, hsp.hspExpect, hsp.hspQueryFrom,
                                                                       hsp.hspQueryTo, hsp.hspHitFrom, hsp.hspHitTo,
                                                                       hsp.hspIdentity, hsp.hspPositive,  hsp.hspGap,
                                                                       hsp.hspQseq, hsp.hspHseq, hsp.hspMidline,
                                                                       hsp.hspAlignLen, hsp.hspQueryFrame, 
                                                                       hsp.hspHitFrame])
#                             DBInfo[hit.dbInfo]['blast_m0'][el_hsp_qry]['all']=[['hspNum','hspExpect','hspQueryFrom',
#                                                                        'hspQueryTo', 'hspHitFrom', 'hspHitTo', 
#                                                                        'hspIdentity', 'hspPositive','hspGap',
#                                                                        'hspQseq', 'hspHseq','hspMidline', 
#                                                                        'hspAlignLen', 'hspQueryFrame',
#                                                                        'hspHitFrame'], ...]

    return DBInfo
##############################################################################

def read_fasta(fhin):
    """
    Read a fasta file: return sequence and description
    """
    error = True
    comment = ''
    sequence = ''
    pos = fhin.tell()
    line = fhin.readline()
    flag = None
    while line != '':
        if line[0] == '>':
            # si ligne d'entete deja vu
            if flag is not None:
                # recule le pointeur de lecture au debut de la ligne
                fhin.seek(pos)
                flag = None
                break
            else:
                flag = 1
                comment = line[1:].split()[0]
                error = False
        elif not error:
            line = line.replace(' ', '')
            sequence = sequence + line[:-1]
        else:
            raise IOError("Sequence file '%s' has an invalid format." % fhin.name)
        # recuperation du curseur de lecture
        pos = fhin.tell()
        line = fhin.readline()
    return comment, sequence


def extract_fasta(query_seq_fh, query):
    pos0 = query_seq_fh.tell()
    comment, sequence = read_fasta(query_seq_fh)
    while comment and sequence:
        if comment == query:
            query_seq_fh.seek(pos0)
            return sequence
        comment, sequence = read_fasta(query_seq_fh)
    query_seq_fh.seek(pos0)

def extract_seq(sbjct_seq):
    fld = sbjct_seq.split('\n')
    return ''.join(fld[1:])

def extrat(seq, start, stop):
    return seq[start-1: stop]


def complement(seq):
    newseq = seq.upper()
    seq = newseq.replace('A', 't').replace('T', 'a').replace('G', 'c').replace('C', 'g')
    return seq.upper()



def reverse(seq):
    return seq[::-1]


def reverse_complement(seq):
    #return reverse(complement(seq))
    return complement(reverse(seq))


def translation(sequence, genetic_code):
    prot = ""
    for i in range(0,len(sequence),3):
        codon = sequence[i:i+3].upper()
        if len(codon)%3 == 0:
            try:
                prot = prot + genetic_code[codon]
            except:
                prot = prot + '.'
        else:
            prot = prot + len(codon)*'+'
    return prot

def split_hsp_qry(DBInfo, sbjct_acc):
    if 'blast_m0' in DBInfo[sbjct_acc]:
        hsp_queries = DBInfo[sbjct_acc]['blast_m0'].keys()
        for hsp_qry in hsp_queries:
            DBInfo[sbjct_acc]['blast_m0'][hsp_qry]['plus'] = []
            DBInfo[sbjct_acc]['blast_m0'][hsp_qry]['minus'] = []
            for hsp in DBInfo[sbjct_acc]['blast_m0'][hsp_qry]['all']:
                if hsp[4] < hsp[5]:
                    DBInfo[sbjct_acc]['blast_m0'][hsp_qry]['plus'].append(hsp)
                else:
                    DBInfo[sbjct_acc]['blast_m0'][hsp_qry]['minus'].append(hsp)
    return DBInfo


def all_hsp_qry(DBInfo, sbjct_acc):
    if 'blast_m0' in DBInfo[sbjct_acc]:
        hsp_queries = DBInfo[sbjct_acc]['blast_m0'].keys()
        for hsp_qry in hsp_queries:
            all_plus = DBInfo[sbjct_acc]['blast_m0'][hsp_qry]['plus']
            all_plus.sort(lambda x, y: cmp(x[4], y[4]))
            all_minus = DBInfo[sbjct_acc]['blast_m0'][hsp_qry]['minus']
            all_minus.sort(lambda x, y: cmp(x[5], y[5]))
            DBInfo[sbjct_acc]['blast_m0'][hsp_qry]['all_sort'] = []
            p, m  = 0, 0
            
            if all_plus and all_minus:
                while p < len(all_plus) and m < len(all_minus):
                    if all_plus[p][4] < all_minus[m][5]:
                        DBInfo[sbjct_acc]['blast_m0'][hsp_qry]['all_sort'].append(all_plus[p])
                        p += 1
                    else:
                        DBInfo[sbjct_acc]['blast_m0'][hsp_qry]['all_sort'].append(all_minus[m])
                        m += 1
            if p < len(all_plus):
                DBInfo[sbjct_acc]['blast_m0'][hsp_qry]['all_sort'].extend(all_plus[p:])
            if m < len(all_minus):
                DBInfo[sbjct_acc]['blast_m0'][hsp_qry]['all_sort'].extend(all_minus[m:])
    return DBInfo


def merge_hsp_plus(DBInfo, sbjct_acc):
    if 'blast_m0' in DBInfo[sbjct_acc]:
        hsp_queries = DBInfo[sbjct_acc]['blast_m0'].keys()
        for hsp_qry in hsp_queries:
            plus = DBInfo[sbjct_acc]['blast_m0'][hsp_qry]['plus']
            plus.sort(lambda x, y: cmp(x[4], y[4]))
            if plus:
                DBInfo[sbjct_acc]['blast_m0'][hsp_qry]['plus_merge'] = [plus[0]]
            else:
                DBInfo[sbjct_acc]['blast_m0'][hsp_qry]['plus_merge'] = []
            i = 1
            while i < len(plus):
                if plus[i-1][5] <= plus[i][4]:
                    "hsp(i-1) .. hsp(i) ===> cote a cote: garde"
                    #print 'plus garde', sbjct_acc, hsp_qry, plus[i-1][4], plus[i-1][5], plus[i][4], plus[i][5]
                    DBInfo[sbjct_acc]['blast_m0'][hsp_qry]['plus_merge'].append(plus[i])
                elif plus[i-1][5] < plus[i][5]:
                    "hsp(i-1) .."
                    ".. hsp(i)   ===>  quinconce garde"
                    DBInfo[sbjct_acc]['blast_m0'][hsp_qry]['plus_merge'].append(plus[i])
                    #print 'plus garde', hsp_qry, plus[i-1][4], plus[i-1][5], plus[i][4], plus[i][5]

                else:
                    if plus[i-1][9] != plus[i][9]:
                        DBInfo[sbjct_acc]['blast_m0'][hsp_qry]['plus_merge'].append(plus[i])
                    #print 'plus gete', hsp_qry, plus[i-1][4], plus[i-1][5], plus[i][4], plus[i][5]
                    pass
                i += 1
    return DBInfo


def merge_hsp_minus(DBInfo, sbjct_acc):
    if 'blast_m0' in DBInfo[sbjct_acc]:
        hsp_queries = DBInfo[sbjct_acc]['blast_m0'].keys()
        for hsp_qry in hsp_queries:
            minus = DBInfo[sbjct_acc]['blast_m0'][hsp_qry]['minus']
            minus.sort(lambda x, y: cmp(x[4], y[4]))
            if minus:
                DBInfo[sbjct_acc]['blast_m0'][hsp_qry]['minus_merge'] = [minus[0]]
            else:
                DBInfo[sbjct_acc]['blast_m0'][hsp_qry]['minus_merge'] = []
            i = 1
            while i < len(minus):
                
                if minus[i-1][5] >= minus[i][4]:
                    "hsp(i-1) .. hsp(i) ===> cote a cote: garde"
                    DBInfo[sbjct_acc]['blast_m0'][hsp_qry]['minus_merge'].append(minus[i])
                    #print 'minus garde', hsp_qry, minus[i-1][4], minus[i-1][5], minus[i][4], minus[i][5]
                    
                elif minus[i-1][4] < minus[i][4]:
                    "hsp(i-1) .."
                    ".. hsp(i)   ===>  quinconce garde"
                    DBInfo[sbjct_acc]['blast_m0'][hsp_qry]['minus_merge'].append(minus[i])
                    #print 'minus garde', hsp_qry, minus[i-1][4], minus[i-1][5], minus[i][4], minus[i][5]
                else:
                    if minus[i-1][9] != minus[i][9]:
                        DBInfo[sbjct_acc]['blast_m0'][hsp_qry]['minus_merge'].append(minus[i])
                    #print 'minus gete',  hsp_qry, minus[i-1][4], minus[i-1][5], minus[i][4], minus[i][5]
                i += 1
    return DBInfo


def merge_query_plus(DBInfo, sbjct_acc, new_DBInfo):
    new_DBInfo[sbjct_acc]['blast_m0']['all_plus_merge']=[]
    if 'blast_m0' in DBInfo[sbjct_acc]:
        hsp_queries = DBInfo[sbjct_acc]['blast_m0'].keys()
        all_plus = []
        for hsp_qry in hsp_queries:
            plus_merge = DBInfo[sbjct_acc]['blast_m0'][hsp_qry]['plus_merge']
            all_plus.extend(plus_merge)
        all_plus.sort(lambda x, y: cmp(x[4], y[4]))
        if all_plus:
            new_DBInfo[sbjct_acc]['blast_m0']['all_plus_merge'].append(all_plus[0])
        i = 1
        while i < len(all_plus):
            if all_plus[i-1][5] <= all_plus[i][4]:
                "hsp(i-1) .. hsp(i) ===> cote a cote: garde"
                new_DBInfo[sbjct_acc]['blast_m0']['all_plus_merge'].append(all_plus[i])
                #print >>sys.stderr, 'plus garde', sbjct_acc, all_plus[i-1][4], all_plus[i-1][5], all_plus[i][4], all_plus[i][5]

            elif all_plus[i-1][5] < all_plus[i][5]:
                "hsp(i-1) .."
                ".. hsp(i)   ===>  quinconce garde"
                new_DBInfo[sbjct_acc]['blast_m0']['all_plus_merge'].append(all_plus[i])
                #print >>sys.stderr, 'plus garde', sbjct_acc, all_plus[i-1][4], all_plus[i-1][5], all_plus[i][4], all_plus[i][5]

            else:
                if all_plus[i-1][9] != all_plus[i][9]:
                    new_DBInfo[sbjct_acc]['blast_m0']['all_plus_merge'].append(all_plus[i])
            i += 1
    return new_DBInfo


def merge_query_minus(DBInfo, sbjct_acc, new_DBInfo):
    new_DBInfo[sbjct_acc]['blast_m0']['all_minus_merge']  = []
    if 'blast_m0' in DBInfo[sbjct_acc]:
        hsp_queries = DBInfo[sbjct_acc]['blast_m0'].keys()
        all_minus = []
        for hsp_qry in hsp_queries:
            minus_merge = DBInfo[sbjct_acc]['blast_m0'][hsp_qry]['minus_merge']
            all_minus.extend(minus_merge)
        all_minus.sort(lambda x, y: cmp(x[4], y[4]))
        if all_minus:
            new_DBInfo[sbjct_acc]['blast_m0']['all_minus_merge'].append(all_minus[0])
        
        i = 1
        while i < len(all_minus):
            
            if all_minus[i-1][5] >= all_minus[i][4]:
                "hsp(i-1) .. hsp(i) ===> cote a cote: garde"
                new_DBInfo[sbjct_acc]['blast_m0']['all_minus_merge'].append(all_minus[i])
            elif all_minus[i-1][4] < all_minus[i][4]:
                "hsp(i-1) .."
                ".. hsp(i)   ===>  quinconce garde"
                new_DBInfo[sbjct_acc]['blast_m0']['all_minus_merge'].append(all_minus[i])
            else:
                if all_minus[i-1][9] != all_minus[i][9]:
                    new_DBInfo[sbjct_acc]['blast_m0']['all_minus_merge'].append(all_minus[i])
            i += 1
    return new_DBInfo


def merge_query_all(DBInfo, sbjct_acc):
    merge_DBInfo = {sbjct_acc: {'seq': DBInfo[sbjct_acc]['seq'], 'blast_m0': {}}}
    merge_DBInfo = merge_query_plus(DBInfo, sbjct_acc, merge_DBInfo)
    merge_DBInfo = merge_query_minus(DBInfo, sbjct_acc, merge_DBInfo)
    all_plus_merge = merge_DBInfo[sbjct_acc]['blast_m0']['all_plus_merge']
    all_plus_merge.sort(lambda x, y: cmp(x[4], y[4]))
    all_minus_merge = merge_DBInfo[sbjct_acc]['blast_m0']['all_minus_merge']
    all_minus_merge.sort(lambda x, y: cmp(x[5], y[5]))
    merge_DBInfo[sbjct_acc]['blast_m0']['all_merge'] = []
    p, m  = 0, 0
    if all_plus_merge and all_minus_merge:
        while p < len(all_plus_merge) and m < len(all_minus_merge):
            if all_plus_merge[p][4] < all_minus_merge[m][5]:
                merge_DBInfo[sbjct_acc]['blast_m0']['all_merge'].append(all_plus_merge[p])
                p += 1
            else:
                merge_DBInfo[sbjct_acc]['blast_m0']['all_merge'].append(all_minus_merge[m])
                m += 1
    if p < len(all_plus_merge):
        merge_DBInfo[sbjct_acc]['blast_m0']['all_merge'].extend(all_plus_merge[p:])
    if m < len(all_minus_merge):
        merge_DBInfo[sbjct_acc]['blast_m0']['all_merge'].extend(all_minus_merge[m:])
    return merge_DBInfo


def lfind(word, car = '-'):
    "Return indexes in s where caracter is found--> list"
    lid = []
    for i in range(len(word)):
        if word[i] == car:
            lid.append(i)
    return lid

def word_insert(word, lpos, car):
    delta = 0
    nw = word
    for pos in lpos:
        nw = nw[:pos+delta] + car + nw[pos+delta:]
        #delta += 1
    return nw


def word_replace(tag, word, start, stop, nword):
    return word[:start]+nword+ word[stop:]


def cons_max_pos(l):
    maxi = 0
    vu = ['.']
    value = '.'
    for e in l:
        if e in vu:
            pass
        else:
            vu.append(e)
            ctp = l.count(e)
            if maxi < ctp:
                maxi = ctp
                value = e
    return value 

def consensus_max(str_qry, other):
    cons = ''
    i = 0
    while i< len(str_qry):
        cons += cons_max_pos([e[i] for e in [ str_qry] + other])
        i +=1
    return cons


def cons_max_iupac(l, jok):
    IUPAC = {'A_G':'R', 'C_T':'Y', 
             'A_T':'W', 'C_G':'S', 
             'G_T':'K', 'A_C':'M',
             'C_G_T': 'B', 'A_G_T': 'D',
             'A_C_T': 'H',  'A_C_G': 'V'}
         

    vu = []
    for e in l:
        if e != '.' and e != '-':
            vu.append(e.upper())
    vu = sorted(vu,key=l.count, reverse=True)
    s = set(vu)
    if len(s) > 4:
        # prot
        return jok # vu[0] 
    elif len(s) == 4 :
        return 'N'
    elif len(s) == 1 :
        return vu[0]
    elif len(s) != 0: 
        new_l = list(s)
        new_l.sort()
        new_k = '_'.join(new_l)
        if new_k in IUPAC:
            return IUPAC[new_k]
        else:
            # prot
            return jok
    else:
        return '.'

def consensus_iupac(str_qry, other, jok):
    cons = ''
    i = 0
    while i< len(str_qry):
        cons += cons_max_iupac([e[i] for e in [ str_qry] + other], jok)
        i +=1
    return cons


def cons_max_joker(l, jok):
    vu = []
    for e in l:
        if e.upper() not in vu and e != '.':
            vu.append(e.upper())
    if len(vu) == 1 :
        return vu[0]
    if len(vu) == 0 :
        return '.'
    else:
        return jok


def consensus_joker(str_qry, other, jok):
    cons = ''
    i = 0
    while i< len(str_qry):
        cons += cons_max_joker([e[i] for e in [ str_qry] + other], jok)
        i +=1
    return cons


def anlz_blast(sbjct_acc, DBInfo, merge_DBInfo):
    qj, hj = 1, 1
    
    len_dbseq = len(extract_seq(DBInfo[sbjct_acc]['seq']))

    DBInfo = split_hsp_qry(DBInfo, sbjct_acc)
    DBInfo = all_hsp_qry(DBInfo, sbjct_acc)
    DBInfo = merge_hsp_plus(DBInfo, sbjct_acc)
    DBInfo = merge_hsp_minus(DBInfo, sbjct_acc)
    merge_DBInfo = merge_query_all(DBInfo, sbjct_acc)
    return DBInfo, merge_DBInfo

def align_queries_merge(DBInfo, sbjct_acc, hj=1): ## WARNING: only positive frame 
    str_qry = ''
    str_db = ''
    str_mid =''
    start = 0
    other= []
   
#       DBInfo[hit.dbInfo]['blast_m0'][el_hsp_qry]['all']=[['hspNum','hspExpect','hspQueryFrom',
#                                                          'hspQueryTo', 'hspHitFrom', 'hspHitTo', 
#                                                          'hspIdentity', 'hspPositive','hspGap',
#                                                          'hspQseq', 'hspHseq','hspMidline', 
#                                                          'hspAlignLen', 'hspQueryFrame',
#                                                          'hspHitFrame'], ...]
    if 'blast_m0' in DBInfo[sbjct_acc]:
        x = 0
        for hsp in DBInfo[sbjct_acc]['blast_m0']['all_merge']:
            minus_vu = False
            if hsp[4] > hsp[5]:
                hsp_start = hsp[5]
                hsp_end = hsp[4]
                minus_vu = True
            else:
                hsp_start = hsp[4]
                hsp_end = hsp[5]
            if start <= hsp_start/hj:
                if minus_vu:
                    str_db += extract_seq(DBInfo[sbjct_acc]['seq'])[start:(hsp_start-1)/hj] + reverse_complement(hsp[10]).lower()
                    str_mid += ((hsp_start -1)/hj - start)*' ' + reverse_complement(hsp[11])
                    str_qry += ((hsp_start  -1)/hj - start)*'.' + reverse_complement(hsp[9])
                else:
                    str_db += extract_seq(DBInfo[sbjct_acc]['seq'])[start:(hsp_start-1)/hj] + hsp[10]
                    str_mid += ((hsp_start -1)/hj - start)*' ' + hsp[11]
                    str_qry += ((hsp_start  -1)/hj - start)*'.' + hsp[9]
                start = hsp_end/hj
            else:
                x +=1 
                deb, toto = 0, hsp_start/hj
                for i in range(len(hsp[10])):
                    if hsp[10][i] != '-': 
                        toto += 1
                    deb += 1
                    if toto == start:
                        break
                if hj == 1:
                    deb +=1
                # there is a 'gap insertion' in old alignment if gap exists in new hsp
                if minus_vu:
                    l_pos_gap_db = lfind(hsp[10][:-deb], '-')
                    for i in range(len(l_pos_gap_db)):
                        l_pos_gap_db[i] = len(hsp[10]) - l_pos_gap_db[i] 
                    l_pos_gap_qry = lfind(str_db[:hsp_start/hj], '-')
                    for i in range(len(l_pos_gap_qry)):
                        l_pos_gap_qry[i] = len(hsp[10]) - l_pos_gap_qry[i] 
                else:
                    l_pos_gap_db = lfind(hsp[10][:deb], '-')
                    l_pos_gap_qry = lfind(str_db[:hsp_start/hj], '-')
                
                down = len(str_db)
                up = (hsp_start+str_db[:hsp_end].count('-')-1)/hj
                
                str_db =  word_replace( 'str_db:', str_db, up, down, word_insert(str_db[up:], l_pos_gap_db, '-') )
                str_mid = word_replace( 'str_mid:', str_mid, up, down, word_insert(str_mid[up:], l_pos_gap_db, '-') )
                str_qry = word_replace( 'str_qry:', str_qry, up, down, word_insert(str_qry[up:], l_pos_gap_db, '-') )
                k = 0
                while k < len(other):
                    other[k] = word_replace( 'other:', other[k] , up, down, word_insert(other[k][up:down], l_pos_gap_db, '-') )
                    k +=1
                # 'new hsp query in other'
                if minus_vu:
                    new_hsp_qry_add = reverse_complement(hsp[9])
                else:
                    new_hsp_qry_add = hsp[9]
                if hsp_end/hj < start:
                    new = (hsp_start/hj + len(l_pos_gap_qry) -1/hj)*'.'+ new_hsp_qry_add
                    other.append(new)
                else:
                    #new = (down - len(lid))*'.'+ hsp[9]
                    new = (hsp_start/hj + len(l_pos_gap_qry) -1/hj)*'.'+ new_hsp_qry_add
                    other.append(new)
                
                if minus_vu:
                    str_db +=  reverse_complement(hsp[10][:-deb])#.upper()
                    str_mid += reverse_complement(hsp[11][:-deb])
                    str_qry += len(hsp[9][:-deb])*'.'
                else:
                    str_db +=  hsp[10][deb:]
                    str_mid += hsp[11][deb:]
                    str_qry += len(hsp[9][deb:])*'.'
                if hsp_end/hj >= start:
                    start = hsp_end/hj

    str_db += extract_seq(DBInfo[sbjct_acc]['seq'])[start:]
    str_mid += (len(extract_seq(DBInfo[sbjct_acc]['seq'])) - start)*' ' 
    str_qry += (len(extract_seq(DBInfo[sbjct_acc]['seq'])) - start)*'.' 
    j = 0
    while j < len(other):
        other[j] = other[j] + (len(str_db) -len (other[j]))*'.' 
        j+=1

    return str_db, str_mid, str_qry, other


def align_all_hsps(DBInfo, sbjct_acc, hj=1): ## WARNING: only positive frame 
    """ alignment of hsps: query VS seqref: mode 1: hspblastMerge"""
    str_qry = ''
    str_db = ''
    str_mid =''
    start = 0
    other= []
   
#       DBInfo[hit.dbInfo]['blast_m0'][el_hsp_qry]['all']=[['hspNum','hspExpect','hspQueryFrom',
#                                                          'hspQueryTo', 'hspHitFrom', 'hspHitTo', 
#                                                          'hspIdentity', 'hspPositive','hspGap',
#                                                          'hspQseq', 'hspHseq','hspMidline', 
#                                                          'hspAlignLen', 'hspQueryFrame',
#                                                          'hspHitFrame'], ...]
    if 'blast_m0' in DBInfo[sbjct_acc]:
        hsp_queries = DBInfo[sbjct_acc]['blast_m0'].keys()
        all_plus = []
        all_minus = []
        for hsp_qry in hsp_queries:
            all_plus.extend(DBInfo[sbjct_acc]['blast_m0'][hsp_qry]['plus'])
            all_minus.extend(DBInfo[sbjct_acc]['blast_m0'][hsp_qry]['minus'])
        all_plus.sort(lambda x, y: cmp(x[4], y[4]))
        all_minus.sort(lambda x, y: cmp(x[5], y[5]))
        all_hsps = []
        p, m  = 0, 0
        if all_plus and all_minus:
            while p < len(all_plus) and m < len(all_minus):
                if all_plus[p][4] < all_minus[m][5]:
                    all_hsps.append(all_plus[p])
                    p += 1
                else:
                    all_hsps.append(all_minus[m])
                    m += 1
        if p < len(all_plus):
            all_hsps.extend(all_plus[p:])
        if m < len(all_minus):
            all_hsps.extend(all_minus[m:])
        
        x = 0
        for hsp in all_hsps:
            minus_vu = False
            if hsp[4] > hsp[5]:
                hsp_start = hsp[5]
                hsp_end = hsp[4]
                minus_vu = True
            else:
                hsp_start = hsp[4]
                hsp_end = hsp[5]
            if start <= hsp_start/hj:
                if minus_vu:
                    str_db += extract_seq(DBInfo[sbjct_acc]['seq'])[start:(hsp_start-1)/hj] + reverse_complement(hsp[10]).lower()
                    str_mid += ((hsp_start -1)/hj - start)*' ' + reverse_complement(hsp[11])
                    str_qry += ((hsp_start  -1)/hj - start)*'.' + reverse_complement(hsp[9])
                else:
                    str_db += extract_seq(DBInfo[sbjct_acc]['seq'])[start:(hsp_start-1)/hj] + hsp[10]
                    str_mid += ((hsp_start -1)/hj - start)*' ' + hsp[11]
                    str_qry += ((hsp_start  -1)/hj - start)*'.' + hsp[9]
                start = hsp_end/hj
            else:
                x +=1 
                deb, toto = 0, hsp_start/hj
                for i in range(len(hsp[10])):
                    if hsp[10][i] != '-': 
                        toto += 1
                    deb += 1
                    if toto == start:
                        break
                if hj == 1:
                    deb +=1
                # there is a 'gap insertion' in old alignment if gap exists in new hsp
                if minus_vu:
                    l_pos_gap_db = lfind(hsp[10][:-deb], '-')
                    for i in range(len(l_pos_gap_db)):
                        l_pos_gap_db[i] = len(hsp[10]) - l_pos_gap_db[i] 
                    l_pos_gap_qry = lfind(str_db[:hsp_start/hj], '-')
                    for i in range(len(l_pos_gap_qry)):
                        l_pos_gap_qry[i] = len(hsp[10]) - l_pos_gap_qry[i] 
                else:
                    l_pos_gap_db = lfind(hsp[10][:deb], '-')
                    l_pos_gap_qry = lfind(str_db[:hsp_start/hj], '-')
                
                down = len(str_db)
                up = (hsp_start+str_db[:hsp_end].count('-')-1)/hj
                
                str_db =  word_replace( 'str_db:', str_db, up, down, word_insert(str_db[up:], l_pos_gap_db, '-') )
                str_mid = word_replace( 'str_mid:', str_mid, up, down, word_insert(str_mid[up:], l_pos_gap_db, '-') )
                str_qry = word_replace( 'str_qry:', str_qry, up, down, word_insert(str_qry[up:], l_pos_gap_db, '-') )
                k = 0
                while k < len(other):
                    other[k] = word_replace( 'other:', other[k] , up, down, word_insert(other[k][up:down], l_pos_gap_db, '-') )
                    k +=1
                # 'new hsp query in other'
                if minus_vu:
                    new_hsp_qry_add = reverse_complement(hsp[9])
                else:
                    new_hsp_qry_add = hsp[9]
                if hsp_end/hj < start:
                    new = (hsp_start/hj + len(l_pos_gap_qry) -1/hj)*'.'+ new_hsp_qry_add
                    other.append(new)
                else:
                    #new = (down - len(lid))*'.'+ hsp[9]
                    new = (hsp_start/hj + len(l_pos_gap_qry) -1/hj)*'.'+ new_hsp_qry_add
                    other.append(new)
                
                if minus_vu:
                    str_db +=  reverse_complement(hsp[10][:-deb])#.upper()
                    str_mid += reverse_complement(hsp[11][:-deb])
                    str_qry += len(hsp[9][:-deb])*'.'
                else:
                    str_db +=  hsp[10][deb:]
                    str_mid += hsp[11][deb:]
                    str_qry += len(hsp[9][deb:])*'.'
                if hsp_end/hj >= start:
                    start = hsp_end/hj

    str_db += extract_seq(DBInfo[sbjct_acc]['seq'])[start:]
    str_mid += (len(extract_seq(DBInfo[sbjct_acc]['seq'])) - start)*' ' 
    str_qry += (len(extract_seq(DBInfo[sbjct_acc]['seq'])) - start)*'.' 
    j = 0
    while j < len(other):
        other[j] = other[j] + (len(str_db) -len (other[j]))*'.' 
        j+=1

    return str_db, str_mid, str_qry, other



def align_by_query(DBInfo, sbjct_acc, hsp_qry, hj=1): ## WARNING: only positive frame 
    """ alignment of hsps: query VS seqref: mode 1: hspblastMerge"""
    str_qry = ''
    str_db = ''
    str_mid =''
    start = 0
    other= []
   
#       DBInfo[hit.dbInfo]['blast_m0'][el_hsp_qry]['all']=[['hspNum','hspExpect','hspQueryFrom',
#                                                          'hspQueryTo', 'hspHitFrom', 'hspHitTo', 
#                                                          'hspIdentity', 'hspPositive','hspGap',
#                                                          'hspQseq', 'hspHseq','hspMidline', 
#                                                          'hspAlignLen', 'hspQueryFrame',
#                                                          'hspHitFrame'], ...]
    x = 0
    for hsp in DBInfo[sbjct_acc]['blast_m0'][hsp_qry]['all_sort']:
        minus_vu = False
        if hsp[4] > hsp[5]:
            hsp_start = hsp[5]
            hsp_end = hsp[4]
            minus_vu = True
        else:
            hsp_start = hsp[4]
            hsp_end = hsp[5]
        if start <= hsp_start/hj:
            if minus_vu:
                str_db += extract_seq(DBInfo[sbjct_acc]['seq'])[start:(hsp_start-1)/hj] + reverse_complement(hsp[10]).lower()
                str_mid += ((hsp_start -1)/hj - start)*' ' + reverse_complement(hsp[11])
                str_qry += ((hsp_start  -1)/hj - start)*'.' + reverse_complement(hsp[9])
            else:
                str_db += extract_seq(DBInfo[sbjct_acc]['seq'])[start:(hsp_start-1)/hj] + hsp[10]
                str_mid += ((hsp_start -1)/hj - start)*' ' + hsp[11]
                str_qry += ((hsp_start  -1)/hj - start)*'.' + hsp[9]
            start = hsp_end/hj
        else:
            x +=1 
            deb, toto = 0, hsp_start/hj
            for i in range(len(hsp[10])):
                if hsp[10][i] != '-': 
                    toto += 1
                deb += 1
                if toto == start:
                    break
            if hj == 1:
                deb +=1
            # there is a 'gap insertion' in old alignment if gap exists in new hsp
            if minus_vu:
                l_pos_gap_db = lfind(hsp[10][:-deb], '-')
                for i in range(len(l_pos_gap_db)):
                    l_pos_gap_db[i] = len(hsp[10]) - l_pos_gap_db[i] 
                l_pos_gap_qry = lfind(str_db[:hsp_start/hj], '-')
                for i in range(len(l_pos_gap_qry)):
                    l_pos_gap_qry[i] = len(hsp[10]) - l_pos_gap_qry[i] 
            else:
                l_pos_gap_db = lfind(hsp[10][:deb], '-')
                l_pos_gap_qry = lfind(str_db[:hsp_start/hj], '-')
            
            down = len(str_db)
            up = (hsp_start+str_db[:hsp_end].count('-')-1)/hj
            
            str_db =  word_replace( 'str_db:', str_db, up, down, word_insert(str_db[up:], l_pos_gap_db, '-') )
            str_mid = word_replace( 'str_mid:', str_mid, up, down, word_insert(str_mid[up:], l_pos_gap_db, '-') )
            str_qry = word_replace( 'str_qry:', str_qry, up, down, word_insert(str_qry[up:], l_pos_gap_db, '-') )
            k = 0
            while k < len(other):
                other[k] = word_replace( 'other:', other[k] , up, down, word_insert(other[k][up:down], l_pos_gap_db, '-') )
                k +=1
            # 'new hsp query in other'
            if minus_vu:
                new_hsp_qry_add = reverse_complement(hsp[9])
            else:
                new_hsp_qry_add = hsp[9]
            if hsp_end/hj < start:
                new = (hsp_start/hj + len(l_pos_gap_qry) -1/hj)*'.'+ new_hsp_qry_add
                other.append(new)
            else:
                #new = (down - len(lid))*'.'+ hsp[9]
                new = (hsp_start/hj + len(l_pos_gap_qry) -1/hj)*'.'+ new_hsp_qry_add
                other.append(new)
            
            if minus_vu:
                str_db +=  reverse_complement(hsp[10][:-deb])#.upper()
                str_mid += reverse_complement(hsp[11][:-deb])
                str_qry += len(hsp[9][:-deb])*'.'
            else:
                str_db +=  hsp[10][deb:]
                str_mid += hsp[11][deb:]
                str_qry += len(hsp[9][deb:])*'.'
            if hsp_end/hj >= start:
                start = hsp_end/hj

    str_db += extract_seq(DBInfo[sbjct_acc]['seq'])[start:]
    str_mid += (len(extract_seq(DBInfo[sbjct_acc]['seq'])) - start)*' ' 
    str_qry += (len(extract_seq(DBInfo[sbjct_acc]['seq'])) - start)*'.' 
    j = 0
    while j < len(other):
        other[j] = other[j] + (len(str_db) -len (other[j]))*'.' 
        j+=1

    return str_db, str_mid, str_qry, other



def print_by_query(alnfopfx, DBInfo, sbjct_acc, colwidth, hj):
    
    len_dbseq = len(extract_seq(DBInfo[sbjct_acc]['seq']))
    if 'blast_m0' in DBInfo[sbjct_acc]:
        hsp_queries = DBInfo[sbjct_acc]['blast_m0'].keys()
        for hsp_qry in hsp_queries:
            alnfh = open(alnfopfx + '%s_%s.aln' % (sbjct_acc,hsp_qry.replace('/','')), 'w')
            #fafh = open(alnfopfx + '%s_%s.fa' % (sbjct_acc,hsp_qry.replace('/','')), 'w')
            str_db, str_mid, str_qry, other = align_by_query(DBInfo, sbjct_acc, hsp_qry, hj)
            str_cons_max = consensus_max(str_qry, other)
            str_cons_joker = consensus_joker(str_qry, other, '*')
            str_cons_upac = consensus_iupac(str_qry, other, '*')
            i = 0
            hfr = 1
            n_gap = 0
            str_len_db = len(str(len_dbseq))
            print >>alnfh, "%s  %s  %s" % (20 *"=", hsp_qry, 20 *"=" )
            while i < len(str_qry) :
                x = 1
                ht = hfr + (colwidth - str_db[i:i+colwidth].count('-'))*hj -1
                if ht > len_dbseq:
                    ht = len_dbseq
                print >>alnfh, ('%%-8s:  %%-%ss %%s %%s' % str_len_db) % (sbjct_acc[0:8], hfr, str_db[i:i+colwidth].lower(), ht)
                hfr = ht+1
                print >>alnfh, ('%%-11s%%-%ss %%s %%s' % str_len_db) % ('',' ',  str_mid[i:i+colwidth], ' ')
                print >>alnfh, ('%%-11s%%-%ss %%s %%s' % str_len_db) % ('HSP_%s:' % x, ' ',  str_qry[i:i+colwidth].lower(), ' ')
                for str_qn in other:
                    x+=1
                    if str_qn[i:i+colwidth] != len(str_qn[i:i+colwidth])*'.':
                    #if 1:
                        print >>alnfh, ('%%-11s%%-%ss %%s %%s' % str_len_db) % ('HSP_%s:' % x, ' ',  str_qn[i:i+colwidth].lower(), ' ')
                    
                print >>alnfh, ('%%-11s%%-%ss %%s %%s' % str_len_db) % ('Cons_max:', ' ',  str_cons_max[i:i+colwidth].upper(), ' ')
                print >>alnfh, ('%%-11s%%-%ss %%s %%s' % str_len_db) % ('Cons_joker:', ' ',  str_cons_joker[i:i+colwidth].upper(), ' ')
                print >>alnfh, ('%%-11s%%-%ss %%s %%s' % str_len_db) % ('Cons_uapc:', ' ',  str_cons_upac[i:i+colwidth].upper(), ' ')
        
                print >>alnfh, ''
                i = i + colwidth
            alnfh.close()
            #fafh.close()
    else:
        alnfh = open(alnfopfx + '%s_%s.aln' % (sbjct_acc, 'no_blast_m0'), 'w')
        #fafh = open(alnfopfx + '%s_%s.fa' % (sbjct_acc, 'no_blast_m0'), 'w')
        print >>alnfh, ''
        alnfh.close()
        #fafh.close()


def print_all_hsps(alnfh, DBInfo, sbjct_acc, colwidth, hj):
    len_dbseq = len(extract_seq(merge_DBInfo[sbjct_acc]['seq']))
    str_db, str_mid, str_qry, other = align_all_hsps(DBInfo, sbjct_acc, hj)
    str_cons_max = consensus_max(str_qry, other)
    str_cons_joker = consensus_joker(str_qry, other, '*')
    str_cons_upac = consensus_iupac(str_qry, other, '*')
    i = 0
    hfr = 1
    n_gap = 0
    str_len_db = len(str(len_dbseq))
    while i < len(str_qry) :
        x = 1
        ht = hfr + (colwidth - str_db[i:i+colwidth].count('-'))*hj -1
        if ht > len_dbseq:
            ht = len_dbseq
        print >>alnfh, ('%%-8s:  %%-%ss %%s %%s' % str_len_db) % (sbjct_acc[0:8], hfr, str_db[i:i+colwidth].lower(), ht)
        hfr = ht+1
        print >>alnfh, ('%%-11s%%-%ss %%s %%s' % str_len_db) % ('',' ',  str_mid[i:i+colwidth], ' ')
        print >>alnfh, ('%%-11s%%-%ss %%s %%s' % str_len_db) % ('HSP_%s:' % x, ' ',  str_qry[i:i+colwidth].lower(), ' ')
        
        for str_qn in other:
            x+=1
            if str_qn[i:i+colwidth] != len(str_qn[i:i+colwidth])*'.':
            #if 1:
                print >>alnfh, ('%%-11s%%-%ss %%s %%s' % str_len_db) % ('HSP_%s:' % x, ' ',  str_qn[i:i+colwidth].lower(), ' ')
            
        print >>alnfh, ('%%-11s%%-%ss %%s %%s' % str_len_db) % ('Cons_max:', ' ',  str_cons_max[i:i+colwidth].upper(), ' ')
        print >>alnfh, ('%%-11s%%-%ss %%s %%s' % str_len_db) % ('Cons_joker:', ' ',  str_cons_joker[i:i+colwidth].upper(), ' ')
        print >>alnfh, ('%%-11s%%-%ss %%s %%s' % str_len_db) % ('Cons_uapc:', ' ',  str_cons_upac[i:i+colwidth].upper(), ' ')

        print >>alnfh, ''
        i = i + colwidth
    



def print_queries_merge(alnfh, merge_DBInfo, sbjct_acc, colwidth, hj):
    len_dbseq = len(extract_seq(merge_DBInfo[sbjct_acc]['seq']))
    str_db, str_mid, str_qry, other = align_queries_merge(merge_DBInfo, sbjct_acc, hj)
    str_cons_max = consensus_max(str_qry, other)
    str_cons_joker = consensus_joker(str_qry, other, '*')
    str_cons_upac = consensus_iupac(str_qry, other, '*')
    i = 0
    hfr = 1
    n_gap = 0
    str_len_db = len(str(len_dbseq))
    while i < len(str_qry) :
        x = 1
        ht = hfr + (colwidth - str_db[i:i+colwidth].count('-'))*hj -1
        if ht > len_dbseq:
            ht = len_dbseq
        print >>alnfh, ('%%-8s:  %%-%ss %%s %%s' % str_len_db) % (sbjct_acc[0:8], hfr, str_db[i:i+colwidth].lower(), ht)
        hfr = ht+1
        print >>alnfh, ('%%-11s%%-%ss %%s %%s' % str_len_db) % ('',' ',  str_mid[i:i+colwidth], ' ')
        print >>alnfh, ('%%-11s%%-%ss %%s %%s' % str_len_db) % ('HSP_%s:' % x, ' ',  str_qry[i:i+colwidth].lower(), ' ')
        
        for str_qn in other:
            x+=1
            if str_qn[i:i+colwidth] != len(str_qn[i:i+colwidth])*'.':
            #if 1:
                print >>alnfh, ('%%-11s%%-%ss %%s %%s' % str_len_db) % ('HSP_%s:' % x, ' ',  str_qn[i:i+colwidth].lower(), ' ')
            
        print >>alnfh, ('%%-11s%%-%ss %%s %%s' % str_len_db) % ('Cons_max:', ' ',  str_cons_max[i:i+colwidth].upper(), ' ')
        print >>alnfh, ('%%-11s%%-%ss %%s %%s' % str_len_db) % ('Cons_joker:', ' ',  str_cons_joker[i:i+colwidth].upper(), ' ')
        print >>alnfh, ('%%-11s%%-%ss %%s %%s' % str_len_db) % ('Cons_uapc:', ' ',  str_cons_upac[i:i+colwidth].upper(), ' ')

        print >>alnfh, ''
        i = i + colwidth

if __name__ == '__main__':
    usage = ""

    parser = argparse.ArgumentParser(prog='mappinghsp4taxo.py',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter, usage=usage)
    general_options = parser.add_argument_group(title="Options", description=None)

    general_options.add_argument("-p", "--outprefix", dest="outprefix",
                                 help="outprefix (ex: path)", 
                                 metavar="STRING",
                                 default='')
#     general_options.add_argument("-f", "--query_fasta", dest="query_fasta_fh",
#                                  help="Query in fasta file. Warning: all query are necessary", 
#                                  type=file,
#                                  required=True)
    general_options.add_argument("-b", "--blast_prot", dest="blast_prot_fh",
                                 help="Proteic Blast m1", 
                                 type=file,
                                 required=True)
    general_options.add_argument("-B", "--blast_nuc", dest="blast_nuc_fh",
                                 help="Nuleic Blast m1", 
                                 type=file,
                                 required=True)
    general_options.add_argument("-i", "--dump_name", dest="dump_name",
                                 help="dump_name", 
                                 type=str,
                                 required=True)

    args = parser.parse_args()

    DBInfo, subjct_with_max_query = pickle.load(open('%s.dmp' % args.dump_name, 'r'))
    
    adn2prot = {1:{'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S',
           'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','TGT':'C','TGC':'C','TGA':'*','TGG':'W',
           'CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
           'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
            'ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T',
            'AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
            'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
            'GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G'}
            }
    
    DBInfo = parse_blast(args.blast_prot_fh, DBInfo)
    DBInfo = parse_blast(args.blast_nuc_fh, DBInfo)
    merge_DBInfo = {}
    
    for sbjct_acc in DBInfo.keys():
        fh_aln_merge = open(args.outprefix + '%s.merge.aln' % sbjct_acc, 'w')
        fh_aln_all = open(args.outprefix + '%s.all_hsps.aln' % sbjct_acc, 'w')
        DBInfo, merge_DBInfo= anlz_blast(sbjct_acc, DBInfo, merge_DBInfo)
        print_queries_merge(fh_aln_merge, merge_DBInfo, sbjct_acc, colwidth=70, hj=1)
        print_by_query(args.outprefix, DBInfo, sbjct_acc, colwidth=70, hj=1)
        print_all_hsps(fh_aln_all, DBInfo, sbjct_acc, colwidth=70, hj=1)
        fh_aln_merge.close()
        fh_aln_all.close()