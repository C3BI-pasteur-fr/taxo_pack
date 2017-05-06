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
#from bsddb3 import db as bdb
from bsddb import db as bdb

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
    print >>sys.stderr, GoldenError('The mandatory GOLDENDATA environment variable is not defined. Consult https://github.com/C3BI-pasteur-fr/golden.')
    sys.exit(1)

# ###################


class TaxOptimizerError:
    def __init__(self, err):
        self.err = err

    def __repr__(self):
        return "[taxoptimizer] " + self.err
    
class ParserError:
    def __init__(self, err):
        self.err = err

    def __repr__(self):
        return "[ParserError] " + self.err


class ParserWarning:
    def __init__(self, err):
        self.err = err

    def __repr__(self):
        return "[ParserWarning] " + self.err


class InputLine:
    # vlegrand@pasteur.fr development
    def __init__(self, orig_line, acc, skip_db):
        self.orig_line = orig_line
        self.acc = acc
        self.skip_db = skip_db


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


def parse(flatFile, DE):
    """
    parse db flat file (uniprot, embl, genbank) into a DBRecord object
    """
    if flatFile[:2] == 'ID':
        return parseUniprot(flatFile, DE)
    elif flatFile[:5] == 'LOCUS':
        return parseGenbank(flatFile, DE)
    return '', '', '', ''


##############################################################################
#
#            Golden: Keep this for compatibility and performance testing.
#
##############################################################################


def doGolden(db, ac, DE):
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
    if flatFile:
        orgName, taxId, taxoLight, description = parse(flatFile, DE)  # orgName, taxId, taxoLight, description
        flatFile = ''  # buffer free
        return orgName, taxId, taxoLight, description
    else:
        return '', '', '', ''


# Builds query input string
def buildQueryStr(txt_line, db, acc, l_cards, cnt_cards, l_lines):
    # vlegrand@pasteur.fr development
    skip_db = False
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
        # return '', '', '', ''
        skip_db = True
    if not skip_db:
        l_cards += db
        l_cards += ":"
        l_cards += acc
        l_cards += "\n"
    cnt_cards += 1
    li = InputLine(txt_line, acc, skip_db)
    l_lines.append(li)
    return l_cards, cnt_cards, l_lines

##############################################################################
#
#            Golden Multi. Using new version of golden.
#            #vlegrand@pasteur.fr
#
#  l_input : Depending on input_flag, a list of bank:AC bank:locus separated by '\n'
#  so, bank:AC\nbank:locus...
##############################################################################


def doGoldenMulti(allTaxo, l_cards, DE, allTaxId, osVSoc_bdb):
    idx_res = 0
    lst_input = l_cards.split("\n")
    try:
        flatFile = Golden.access_new(l_cards)
        while (flatFile is not None):
            db_acc = lst_input[idx_res]
            l_db_acc = db_acc.split(":")
            acc = l_db_acc[1]
            if acc not in allTaxo:
                allTaxo[acc] = {'db': l_db_acc[0]}
                allTaxo[acc]['orgName'], allTaxo[acc]['taxId'], allTaxo[acc]['taxoLight'], allTaxo[acc]['DE'] = parse(flatFile, DE)  # orgName, taxId, taxoLight, description
                allTaxo, allTaxId = extractTaxoFrom_osVSocBDB_multi(acc, allTaxo, allTaxId, osVSoc_bdb)

            flatFile = Golden.access_new(l_cards)
            idx_res += 1
        return allTaxo
    except IOError, err:
        print >>sys.stderr, err, l_cards
        sys.exit()


# display results from allTaxo dictionnary.
def printResults(l_lines, allTaxo, outfh, notaxfhout, splitFile):
    for li in l_lines:
        taxonomy = ''
        if not li.skip_db:
            if 'taxoFull' in allTaxo[li.acc]:
                taxonomy = allTaxo[li.acc]['taxoFull']
            else:
                taxonomy = allTaxo[li.acc]['taxoLight']

        if taxonomy:
            print >>outfh, li.orig_line, "\t%s\t%s\t%s" % (allTaxo[li.acc]['orgName'], taxonomy, allTaxo[li.acc]['DE'])
        else:
            if notaxfhout:
                print >>notaxfhout, li.orig_line
            if not splitFile:
                print >>outfh, li.orig_line

##############################################################################
#
#            Taxonomy
#
##############################################################################


def extractTaxoFrom_osVSocBDB(acc, allTaxo, allTaxId, BDB):
    taxonomy = allTaxo[acc]['taxoLight']
    orgName = allTaxo[acc]['orgName']
    if orgName and orgName not in allTaxId:
        taxoFull = BDB.get(str(orgName))
        if taxoFull:
            allTaxo[acc]['taxoFull'] = taxoFull
            allTaxId[orgName] = taxoFull
            allTaxo[acc]['taxoLight'] = ''
            taxonomy = taxoFull
    elif orgName:
        allTaxo[acc]['taxoFull'] = allTaxId[orgName]
        allTaxo[acc]['taxoLight'] = ''
        taxonomy = allTaxId[orgName]
    return taxonomy, allTaxo, allTaxId


def extractTaxoFrom_osVSocBDB_multi(acc, allTaxo, allTaxId, BDB):
    orgName = allTaxo[acc]['orgName']
    if orgName and orgName not in allTaxId:
        taxoFull = BDB.get(str(orgName))
        if taxoFull:
            allTaxo[acc]['taxoFull'] = taxoFull
            allTaxId[orgName] = taxoFull
            allTaxo[acc]['taxoLight'] = ''
    elif orgName:
        allTaxo[acc]['taxoFull'] = allTaxId[orgName]
        allTaxo[acc]['taxoLight'] = ''
    return allTaxo, allTaxId


def extractTaxoFrom_accVSos_ocBDB(acc, allTaxo, BDB):
    # 'acc: os_@#$_oc'
    # allTaxo[acc]['orgName'], allTaxo[acc]['taxId'], allTaxo[acc]['taxoLight'], allTaxo[acc]['DE']
    allTaxo[acc] = {'taxId': '', 'taxoLight': '', 'DE': '', 'orgName': ''}
    os_oc = BDB.get(acc)
    if os_oc:
        os_oc_fld = os_oc.split('_@#$_')
        allTaxo[acc]['orgName'] = os_oc_fld[0]
        allTaxo[acc]['taxoFull'] = os_oc_fld[1]

        return allTaxo[acc]['taxoFull'], allTaxo
    else:
        return '', allTaxo


def column_analyser(fldcolumn, db):
    acc = ''
    if len(fldcolumn) == 5:
        if not db:
            db = fldcolumn[2]
        acc = fldcolumn[3].split('.')[0]
    elif len(fldcolumn) == 2 or len(fldcolumn) == 3:
        if not db:
            db = fldcolumn[0]
        if len(fldcolumn) == 3 and fldcolumn[1] == '':
            acc = fldcolumn[2].split('.')[0]
        else:
            acc = fldcolumn[1].split('.')[0]
    elif len(fldcolumn) == 1 and db:
        acc = fldcolumn[0]
    return acc, db


def main_ncbi(tabfh, outfh, osVSoc_bdb, column, separator, max_cards, notaxofh=None, db=None, splitfile=False, description=False):
    allTaxo = {}
    allTaxId = {}
    try:
        line = tabfh.readline()
        lineNb = 1
    except EOFError, err:
        print >>sys.stderr, err
        sys.exit()
    l_cards = ""
    cnt_cards = 0
    l_lines = []
    while line:
        fld = line.split()
        if line == '\n':
            line = tabfh.readline()
            continue
        try:
            fldcolumn = fld[column - 1].split(separator)
        except:
            print >>sys.stderr, TaxOptimizerError("Parsing: column error: couldn't parse line: \n%s ...\n --> %s" % (line[0:50], lineNb))
            sys.exit()

        acc, db = column_analyser(fldcolumn, db)

        if not acc or not db:
            if not splitfile:
                print >>outfh, line[:-1]
            if notaxofh:
                print >>notaxofh, line[:-1]
                print >>sys.stderr, TaxOptimizerError("Parsing: no acc and db in %s with separator=%s (line %s)" % (fld[column - 1], separator, lineNb))
            try:
                line = tabfh.readline()
                lineNb += 1
            except EOFError, err:
                print >>sys.stderr, err
                print >>sys.stderr, TaxOptimizerError("in line %s" % (lineNb))
                sys.exit()
            continue
        elif db not in ['silva', 'gg']:
            l_cards, cnt_cards, l_lines = buildQueryStr(line[:-1], db, acc, l_cards, cnt_cards, l_lines)
            if cnt_cards == max_cards:
                allTaxo = doGoldenMulti(allTaxo, l_cards, description, allTaxId, osVSoc_bdb)
                printResults(l_lines, allTaxo, outfh, notaxofh, splitfile)
                l_cards = ""
                l_lines = []
                cnt_cards = 0

        try:
            line = tabfh.readline()
            lineNb += 1
        except EOFError, err:
            print >>sys.stderr, err
            print >>sys.stderr, TaxOptimizerError("in line %s" % (lineNb))

            sys.exit()
        lineNb += 1

    if cnt_cards != 0:
        allTaxo = doGoldenMulti(allTaxo, l_cards, description, allTaxId, osVSoc_bdb)
        printResults(l_lines, allTaxo, outfh, notaxofh, splitfile)


def main_vipr(tabfh, outfh, osVSoc_bdb, column, separator, max_cards, notaxofh=None, db=None, splitfile=False, description=False, viprVSdb={}):
    allTaxo = {}
    allTaxId = {}
    try:
        line = tabfh.readline()
        lineNb = 1
    except EOFError, err:
        print >>sys.stderr, err
        sys.exit()
    l_cards = ""
    cnt_cards = 0
    l_lines = []
    while line:
        fld = line.split()
        if line == '\n':
            line = tabfh.readline()
            continue
        try:
            fldcolumn = fld[column - 1].split(separator)
        except:
            print >>sys.stderr, TaxOptimizerError("Parsing: column error: couldn't parse line: \n%s ...\n --> %s" % (line[0:50], lineNb))
            sys.exit()

        acc, db_not_used = column_analyser(fldcolumn, db)

        if not acc or not db:
            if not splitfile:
                print >>outfh, line[:-1]
            if notaxofh:
                print >>notaxofh, line[:-1]
                print >>sys.stderr, TaxOptimizerError("Parsing: no acc and db in %s with separator=%s (line %s)" % (fld[column - 1], separator, lineNb))
            try:
                line = tabfh.readline()
                lineNb += 1
            except EOFError, err:
                print >>sys.stderr, err
                print >>sys.stderr, TaxOptimizerError("in line %s" % (lineNb))
                sys.exit()
            continue
        else:
            # subtitute Vipr by gb
            new_acc = viprVSdb[acc]
            l_cards, cnt_cards, l_lines = buildQueryStr(line[:-1], db, new_acc, l_cards, cnt_cards, l_lines)
            if cnt_cards == max_cards:
                allTaxo = doGoldenMulti(allTaxo, l_cards, description, allTaxId, osVSoc_bdb)
                printResults(l_lines, allTaxo, outfh, notaxofh, splitfile)
                l_cards = ""
                l_lines = []
                cnt_cards = 0

        try:
            line = tabfh.readline()
            lineNb += 1
        except EOFError, err:
            print >>sys.stderr, err
            print >>sys.stderr, TaxOptimizerError("in line %s" % (lineNb))

            sys.exit()
        lineNb += 1

    if cnt_cards != 0:
        allTaxo = doGoldenMulti(allTaxo, l_cards, description, allTaxId, osVSoc_bdb)
        printResults(l_lines, allTaxo, outfh, notaxofh, splitfile)



def main_gg_silva(tabfh, outfh, accVosocBDB, column, separator, notaxofh=None, db=None, splitfile=False, description=False):
    allTaxo = {}
    try:
        line = tabfh.readline()
        lineNb = 1
    except EOFError, err:
        print >>sys.stderr, err
        sys.exit()
    DE = ''
    while line:
        fld = line.split()
        if line == '\n':
            line = tabfh.readline()
            continue
        try:
            fldcolumn = fld[column - 1].split(separator)
        except:
            print >>sys.stderr, TaxOptimizerError("Parsing: column error: couldn't parse line: \n%s ...\n --> %s" % (line[0:50], lineNb))
            sys.exit()

        acc, db = column_analyser(fldcolumn, db)
        if not acc or not db:
            if not splitfile:
                print >>outfh, line[:-1]
            if notaxofh:
                print >>notaxofh, line[:-1]
            print >>sys.stderr, TaxOptimizerError("Parsing: no acc and db in %s with separator=%s (line %s)" % (fld[column - 1], separator, lineNb))

            try:
                line = tabfh.readline()
                lineNb += 1
            except EOFError, err:
                # print >>sys.stderr, err
                print >>sys.stderr, TaxOptimizerError("in line %s" % (lineNb))
                sys.exit()
            continue
        else:
            taxonomy = ''
            if acc in allTaxo:
                if 'taxoFull' in allTaxo[acc]:
                    taxonomy = allTaxo[acc]['taxoFull']
                else:
                    taxonomy = allTaxo[acc]['taxoLight']
                if description:
                    DE = allTaxo[acc]['DE']
            else:
                taxonomy = ''
                allTaxo[acc] = {'db': db}
                taxonomy, allTaxo = extractTaxoFrom_accVSos_ocBDB(acc, allTaxo, accVosocBDB)

            if taxonomy:
                print >>outfh, line[:-1], "\t%s\t%s\t%s" % (allTaxo[acc]['orgName'], taxonomy, DE)
            else:
                if notaxofh:
                    print >>notaxofh, line[:-1]
                if not splitfile:
                    print >>outfh, line[:-1]

        try:
            line = tabfh.readline()
            lineNb += 1
        except EOFError, err:
            print >>sys.stderr, err
            print >>sys.stderr, TaxOptimizerError("in line %s" % (lineNb))

            sys.exit()
        lineNb += 1


def doViprVSgb(fhin):
    viprVSgb = {}
    line = fhin.readline()
    while line:
        fld = line.strip().split()
        viprVSgb[fld[0]] = fld[1]
        line = fhin.readline()
    return viprVSgb


##############################################################################
#
#            MAIN
#
##############################################################################


if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog='taxoptimizer.py',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="Parse a blast output report and add NCBI, SILVA or Greengenes Taxonomy database information in each HSP.")

    general_options = parser.add_argument_group(title="Options", description=None)

    general_options.add_argument("-i", "--in", dest="tabfh",
                                 help="Tabulated input file. (Recommended, Blast m8 file)",
                                 metavar="File",
                                 type=file,
                                 required=True)
    general_options.add_argument("-o", "--out",
                                 action='store',
                                 dest='outfh',
                                 metavar="File",
                                 type=argparse.FileType('w'),
                                 help='Output file',
                                 required=True)
    general_options.add_argument("-b", "--bdb",
                                 dest="bdbfile",
                                 help="Berleley DB file generated by taxodb_ncbi (https://github.com/C3BI-pasteur-fr/taxodb_ncbi) or \
                                 taxo_rrna programs (https://github.com/C3BI-pasteur-fr/taxo_rrna).\
                                 If  NCBITAXODB_BDB, GGTAXODB_BDB or/and SILVATAXODB_BDB unix environment variable(s) are set, this option is not mandatory.",
                                 metavar="File",
                                 )
    general_options.add_argument("-w", "--crossId",
                                 dest="viprID",
                                 help="tabulated file with two column vipr_id and gb_id ",
                                 type=file,
                                 metavar="File",
                                 )
    
    general_options.add_argument("-t", "--bdb_type",
                                 dest="bdbtype",
                                 help="Berleley DB file type. Link to -b option.",
                                 type=str,
                                 default='ncbi',
                                 choices=['ncbi', 'gg', 'silva', 'silva_ssu', 'silva_lsu', 'vipr_nt'],
                                 required=True
                                 )
    general_options.add_argument("-c", "--column",
                                 action='store',
                                 dest='column',
                                 type=int,
                                 help='Column\'s number with ID and database informations for all HSPs',
                                 default=2)
    general_options.add_argument("-s", "--separator",
                                 dest="separator", metavar="str", type=str,
                                 help="Separator in database AC",
                                 default='|')

    general_options.add_argument("-e", "--description",
                                 dest="description",
                                 help="Add database description (DE) in the output.",
                                 action='store_true',
                                 default=False,)
    general_options.add_argument("-x", "--splitfile",
                                 dest="splitfile",
                                 help="Only show lines with a taxonomy correspondance. Could be used with -m option.",
                                 action='store_true',
                                 default=False,)
    general_options.add_argument("-f", "--no_taxo_file",
                                 action='store',
                                 dest='notaxofh',
                                 metavar="File",
                                 type=argparse.FileType('w'),
                                 help='Only show lines without a taxonomy correspondance. Could be used with -x option.',)
    general_options.add_argument('-d', '--database', metavar='str',
                                 dest='database',
                                 type=str,
                                 help="Supposed that all blast HSPs match this database. Not used for Silva and Greengenes databases",
                                 default=None,
                                 )
    general_options.add_argument('-n', '--nbcpu', metavar='int',
                                 dest='nb_cpu',
                                 type=int,
                                 help="Reserved option for mtaxoptimizer mpi program",
                                 )
    golden_options = parser.add_argument_group(title="Golden options", description=None)
    golden_options.add_argument("-m", "--max_cards",
                                action='store',
                                dest='max_cards',
                                type=int,
                                help='Maximum cards number used by Golden2.0 in analyses',
                                default=500)

    args = parser.parse_args()

    # ===== Tabulated file parsing
    NCBITAXODB_BDB = 'taxodb.bdb'
    GGTAXODB_BDB = '%s_accVosoc.bdb' % args.bdbtype
    SILVATAXODB_BDB = '%s_accVosoc.bdb' % args.bdbtype

    if args.bdbfile:
        bdbfile = args.bdbfile
    else:
        if args.bdbtype == 'ncbi':
            try:
                TABLE = os.environ['TAXOBDBDATA']
                bdbfile = TABLE + NCBITAXODB_BDB
            except:
                print >>sys.stderr, TaxOptimizerError("NCBI TaxoDB database is mandatory (bekley db format) or set the TAXOBDBDATA environment variable") 
                sys.exit(1)
        elif args.bdbtype == 'gg':
            try:
                TABLE = os.environ['GGBDBDATA']
                bdbfile = TABLE + GGTAXODB_BDB
            except:
                print >>sys.stderr, TaxOptimizerError("Greengenes TaxoDB database is mandatory (bekley db format)  or set the GGBDBDATA environment variable")
                sys.exit(1)
        elif args.bdbtype == 'silva':
            try:
                TABLE = os.environ['SILVABDBDATA']
                bdbfile = TABLE + SILVATAXODB_BDB
            except:
                print >>sys.stderr, TaxOptimizerError("Silva TaxoDB database is mandatory (bekley db format)  or set the SILVABDBDATA environment variable")
                sys.exit(1)


    if args.bdbtype == 'ncbi':
        osVSoc_bdb = bdb.DB()
        try:
            osVSoc_bdb.open(bdbfile, None, bdb.DB_HASH, bdb.DB_RDONLY)
        except StandardError, err:
            print >>sys.stderr, TaxOptimizerError("NCBI TaxoDB database open error, %s" % err)
            sys.exit()

        main_ncbi(args.tabfh, args.outfh, osVSoc_bdb, args.column, args.separator, args.max_cards, args.notaxofh, args.database, args.splitfile, args.description)
        osVSoc_bdb.close()
    elif args.bdbtype in ['gg', 'silva']:
        accVosocBDB = bdb.DB()
        try:
            accVosocBDB.open(bdbfile, None, bdb.DB_HASH, bdb.DB_RDONLY)
        except StandardError, err:
            print >>sys.stderr, TaxOptimizerError("Taxonomy Berkeley database open error, %s" % err)
            sys.exit()
        main_gg_silva(args.tabfh, args.outfh, accVosocBDB, args.column, args.separator, args.notaxofh, args.database, args.splitfile, args.description)
    elif args.bdbtype in ['vipr_nt']:
        osVSoc_bdb = bdb.DB()
        try:
            osVSoc_bdb.open(bdbfile, None, bdb.DB_HASH, bdb.DB_RDONLY)
        except StandardError, err:
            print >>sys.stderr, TaxOptimizerError("NCBI TaxoDB database open error, %s" % err)
            sys.exit()
        viprVSgb = doViprVSgb(args.viprID)
        main_vipr(args.tabfh, args.outfh, osVSoc_bdb, args.column, args.separator, args.max_cards, args.notaxofh, 'gb', args.splitfile, args.description, viprVSgb)
        osVSoc_bdb.close()