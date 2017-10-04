#! /usr/local/bin/python

# Corinne Maufrais
# Institut Pasteur, Groupe Logisiels et Banques de sonnees
# maufrais@pasteur.fr
#

import xml.sax, sys
from xml.sax.handler import ContentHandler
import re

class BlastError:
    def __init__( self, err ):
        self.err = err

    def __repr__( self ):
        return "[BlastError] " + self.err


class BlastWarning:
    def __init__( self, err ):
        self.err = err

    def __repr__( self ):
        return "[BlastWarning] " + self.err
    
class ParserError:
    def __init__( self, err ):
        self.err = err

    def __repr__( self ):
        return "[ParserError] " + self.err

class _XMLHandler( ContentHandler ):
     """Generic SAX Parser
     Redefine the methods startElement, characters and endElement.
     """

     def __init__( self ):
          "initialisation"
          self._tag = []
          self._value = ''

     def _secure_name( self, name ):
          """Replace '-' with '_' in Xml tag names"""
          return name.replace('-', '_')

     def startDocument( self ):
          "fonction appelee lorsque le parser rencontre le premier element"
          pass

     def startElement( self, name, attrs ):
          """fonction appelee lorsque le parser rencontre une balise ouvrante:
          name -- tag name
          attr -- tag attributes --> not in blast
          """
          self._tag.append( name )
          # Try to call a method
          try:
               eval( self._secure_name('self._start_' + name) )()
          except AttributeError:
               # Doesn't exist (yet)
               pass

     def endElement( self, name ):
          """fonction appelee lorsque le parser rencontre une balise fermante:
          name -- tag name
          """
          # Strip character buffer
          self._value = self._value.strip()
          # Try to call a method (defined in subclasses)
          try:
               eval( self._secure_name('self._end_' + name) )()
          except AttributeError: # Method doesn't exist (yet ?)
               pass

          # Reset character buffer
          self._value = ''

     def characters( self, chrs ):
          """fonction appelee lorsque le parser rencontre des donnees dans
          un element:
          chrs -- characters read
          """
          self._value += chrs

     def endDocument( self ):
          "fonction appelee lorsque le parser rencontre le dernier element"
          pass

class XMLBlastParser( _XMLHandler ):
    """Parse XML BLAST data into a Record.Blast object

    Methods:
    parse           Parses BLAST XML data.

    All XML 'action' methods are private methods and may be:
    _start_TAG      called when the start tag is found
    _end_TAG        called when the end tag is found
    """
    def __init__( self ):
        """ Constructor """
        # Calling superclass method
        _XMLHandler.__init__( self )
        # Create a parser
        self._parser = xml.sax.make_parser()
        # Tell the parser to use my handler
        self._parser.setContentHandler( self )
        # To avoid ValueError: unknown url type: NCBI_BlastOutput.dtd
        self._parser.setFeature( xml.sax.handler.feature_validation, 0 )
        self._parser.setFeature( xml.sax.handler.feature_namespaces, 0 )
        self._parser.setFeature( xml.sax.handler.feature_external_pes, 0 )
        self._parser.setFeature( xml.sax.handler.feature_external_ges, 0 )

        self._allBlast = BlastRecordS()

    def parse( self, handler ):
        """Parses the XML data
        handler -- file handler or StringIO
        """
        # Parse the input
        self._parser.parse( handler )
        return self._allBlast

    # Database Header
    def _end_BlastOutput_program( self ):
        self._allBlast.program = self._value
    
    def _end_BlastOutput_version( self ):
        self._allBlast.version = self._value
        if self._allBlast.version[-1] == '+':
            self._allBlast.vplus = True
    
    def _end_BlastOutput_db( self ):
        """ the database(s) searched """
        self._allBlast.database = self._value
        self._allBlast.databaseAll = self._value

    # Query Header
    def _start_Iteration( self ):
        self._allBlast.blastS.append( BlastRecord() )
        self._blast = self._allBlast.blastS[-1]
    
    def _end_Iteration_iter_num( self ):
        self._blast.blastNum = int( self._value )
    
    def _end_Iteration_query_ID( self ):
        """ the length of the query """
        self._blast.queryId = self._value 
    
    def _end_Iteration_query_def( self ):
        self._blast.queryDef = self._value 
    
    def _end_Iteration_query_len( self ):
        """ the length of the query """
        self._blast.queryLength = int( self._value )

    # Parameters
    # Hits
    def _start_Hit( self ):
        self._blast.hitS.append( BlastHitRecord() )
        self._hit = self._blast.hitS[-1]

    def _end_Hit_num(self):
        self._hit.hitNum = int( self._value )

    def _end_Hit_id( self ):
        self._hit.hitId = self._value

    def _end_Hit_def( self ):
        """ definition line of the database sequence  """
        self._hit.hitDef = self._value
        if len( self._value.split( ' ' )[0].split('|') ) > 1:
            self._hit.dbInfo = self._value.split( ' ' )[0]
        self._hit.comment = ' '.join( self._value.split( ' ' )[1:] )

    def _end_Hit_accession( self ):
        self._hit.hitAcc = self._value

    def _end_Hit_len( self ):
        self._hit.hitLength = int( self._value )

    # HSPs
    def _start_Hsp( self ):
        self._hit.hspS.append( BlastHspRecord( self._blast.queryLength ) )
        self._hsp = self._hit.hspS[-1]

    def _end_Hsp_evalue( self ):
        """ expect value value of the HSP """
        self._hsp.hspExpect = float( self._value )

    def _end_Hsp_num( self ):
        self._hsp.hspNum = int( self._value )

    def _end_Hsp_bit_score( self ):
        self._hsp.hspScore = float( self._value )

    def _end_Hsp_query_from( self ):
        self._hsp.hspQueryFrom = int( self._value )

    def _end_Hsp_query_to( self ):
        self._hsp.hspQueryTo = int( self._value )

    def _end_Hsp_hit_from( self ):
        self._hsp.hspHitFrom = int( self._value )

    def _end_Hsp_hit_to( self ):
        self._hsp.hspHitTo = int( self._value )

    def _end_Hsp_hit_frame(self):
        self._hsp.hspHitFrame = int( self._value )
    
    def _end_Hsp_query_frame(self):
        self._hsp.hspQueryFrame = int( self._value )

    def _end_Hsp_identity( self ):
        self._hsp.hspIdentity = int( self._value )

    def _end_Hsp_positive( self ):
        self._hsp.hspPositive = int( self._value )

    def _end_Hsp_gaps( self ):
        self._hsp.hspGap = int( self._value )

    def _end_Hsp_qseq( self ):
        self._hsp.hspQseq = self._value

    def _end_Hsp_hseq( self ):
        self._hsp.hspHseq = self._value

    def _end_Hsp_midline( self ):
        self._hsp.hspMidline = self._value

    def _end_Hsp_align_len( self ):
        self._hsp.hspAlignLen = int( self._value )

class PairwiseBlastParser:
    def __init__( self ):
        """ Constructor """
        self._allBlast = BlastRecordS()
        
    def parse( self, handler ):
        """Parses the  pairwise data
        handler -- file handler or StringIO
        """
        # Parse the input

        #header
        line = handler.readline()
        fld = line.split()
        self._allBlast.program = fld[0]
        self._allBlast.version = fld[1]
        if '+' in self._allBlast.version:
            self._allBlast.vplus = True
        
        #corps
        line = self._parse( handler )
        if 'psiblast' in line:
            return line
        while line and ('BLAST' in line or 'Query=' in line or 'Searching' in line ):
            # 'multiple blast' or 'psiblast' and 'phiblast'  with multipass
            line = self._parse( handler )
        #footer
        return self._allBlast

    def _parse( self, handler ):
        self._allBlast.blastS.append( BlastRecord( database=self._allBlast.database, 
                                                   dbSeq=self._allBlast.dbSeq, 
                                                   dbLeter=self._allBlast.dbLeter, 
                                                   program=self._allBlast.program,
                                                   version=self._allBlast.version) )
        _blast = self._allBlast.blastS[-1]
        _blast.blastNum = len( self._allBlast.blastS )
        
        #header
        line = handler.readline()
        
        while line and (line[0:19] != 'Sequences producing' and 'No hits found' not in line):
            if 'Database:' in line:
                if 'All non-redundant GenBank' in line:
                    self._allBlast.database = 'nrprot'
                    self._allBlast.databaseAll = 'nrprot'
                    line = handler.readline()
                    while 'letters' not in line:
                        line = handler.readline()
                        line = handler.readline()
                else:
                    fld = line[0:-1].split(':')
                    self._allBlast.database = fld[1].split()[0].strip().lower()
                    self._allBlast.databaseAll = ' '.join( fld[1:] )
                    while 'letters' not in line:
                        line = handler.readline()
                fld2 = line.split(';')
                self._allBlast.dbSeq = fld2[0].split('sequences')[0].strip()
                self._allBlast.dbLeter = fld2[1].split('total letters')[0].strip()
                _blast.dbSeq = self._allBlast.dbSeq
                _blast.dbLeter = self._allBlast.dbLeter
            if 'Query=' in line:
                _blast.queryDef = line[6:-1].strip()
                line = handler.readline()
                while 'letters' not in line and 'Length' not in line:
                    _blast.queryDef += ' ' + line[:-1].strip()
                    line = handler.readline()
                if 'letters' in line:
                    reg = re.compile( r'\s*\((\d*)\s*')
                    m = re.match( reg , line )
                    _blast.queryLength = int( m.group(1) )
                elif 'Length' in line:
                    fld = line[0:-1].split('=')
                    _blast.queryLength = int( fld[1] )
            #psiblast pb
#            if 'Results from round' in line:
#                psiblast = True
#                return 'psiblast'
            line = handler.readline()
            
        
        if 'No hits found' in line:
            line = handler.readline()
            line = handler.readline()
            line = handler.readline()
            if 'Database:' in line: # multiple blast not well form: pblastall of N. Joly
                while line and 'BLAST'  not in line:
                    pos = handler.tell()
                    line = handler.readline()
                handler.seek( pos )
            return line
        
        nbHitComment = False
        if 'Sequences producing' in line:
            fld = line[0:-1].split()
            if fld[-1] == 'N':
                nbHitComment = True
        #comment (newline between comment and sequences )
        line = handler.readline()
        line = handler.readline()
        while line and line != '\n':
            fld = line[0:-1].split()
            if fld and not _blast.hasHit( fld[0] ):
                _blast.hitS.append( BlastHitRecord() )
                hit = _blast.hitS[-1]
                hit.hitNum = len( _blast.hitS )
                hit.dbInfo = fld[0].strip()
                hit.hitDef = ' '.join( fld[:-2] )
                hit.db = line.strip().split()
                hit.hspS.append( BlastHspRecord( _blast.queryLength ) )
                first_hsp = hit.hspS[-1]
                try:
                    if nbHitComment:
                        first_hsp.hspExpect = float( fld[-2].strip() )
                    else:
                        first_hsp.hspExpect = float( fld[-1].strip() )
                except:
                    # e-164  pas reconnu par python mais 1.e-164 oui
                    if nbHitComment:
                        hspExpect = '1' + fld[-2].strip()
                    else:
                        hspExpect = '1' + fld[-1].strip()
                    first_hsp.hspExpect = float( hspExpect )
                if nbHitComment:
                    first_hsp.hspScore = float( fld[-3].strip() )
                else:
                    first_hsp.hspScore = float( fld[-2].strip() )
            line = handler.readline()
        
        #pairwise
        pos = handler.tell()
        line = handler.readline()
        #flush
        hit = None
        hitParsed, hitCvu, vuHsp, vuQ, vuS, ancre = False, False, False, False, False, -1
        hspQueryFrame, hspHitFrame = 0, 0
        # (pairwise finish by 'Database' for single BLAST )
        # (pairwise finish by 'BLAST?' for multiple BLAST )
        # (pairwise finish by 'Searching' for 'psiblast' and 'phiblast' with multipass)
        while line and 'Database:' not in line and 'BLAST' not in line and 'Query=' not in line and 'Searching' not in line:
            if line[0] == '>' and not hitParsed:
                hspFrame = ''
                # first HSP for new Hit
                fld = line[1:-1].split()
                if _blast.hasHit( fld[0] ):
                    hitParsed = True
                    hit = _blast.getHit( fld[0].strip() )
                    hit.hitDef = line[1:-1].strip()
                    hit.comment = ' '.join( fld[1:] )
                    # erase hspS obtained in comment_line
                    hit.hspS = []
                    pos = handler.tell()
                    line = handler.readline()
                    while 'Length' not in line:
                        hit.hitDef += ' ' + line[1:-1].strip()
                        hit.comment += ' ' + line[:-1].strip()
                        pos = handler.tell()
                        line = handler.readline()
                    handler.seek( pos )
                #else:
                #    print sys.stderr, BlastWarning( "Hit '%s' in pairwise but not in comment"  % fld[0] )
                pos = handler.tell()
            elif line[0] == '>' and  hitParsed:
                # new Hit
                hitParsed, hitCvu, vuHsp = False, False, False
                handler.seek( pos )
                vuQ, vuS, ancre = False, False, -1
                self.hsp( hit, _blast.queryLength, hspScore, hspScore2, hspExpect, hspQueryFrom, hspQueryTo, hspHitFrom, hspHitTo, hspIdentity, hspPositive, hspGap, hspAlignLen, hspQseq, hspHseq, hspMidline, hspQueryFrame, hspHitFrame)
            elif hitParsed:
                # other HSP in Hit
                if 'Length' in line:
                    fldl = line.split('=')
                    hit.hitLength = int( fldl[1].strip() )
                elif 'Score' in line:
                    if not vuHsp:
                        vuHsp = True
                    else:
                        vuQ, vuS, ancre = False, False, -1
                        self.hsp( hit, _blast.queryLength, hspScore, hspScore2, hspExpect, hspQueryFrom, hspQueryTo, hspHitFrom, hspHitTo, hspIdentity, hspPositive, hspGap, hspAlignLen, hspQseq, hspHseq, hspMidline, hspQueryFrame, hspHitFrame)
                        
                    flde = line.split( ',' )
                    fldex = flde[1].split( '=' )
                    hspExpect = fldex[-1].strip()
                    fldsc = flde[0].split( '=' )
                    hspScore = fldsc[1].strip().split( 'bits' )[0][:-1]
                    hspScore2 = fldsc[1].strip().split( '(' )[1][:-1]
                elif 'Identities' in line and vuHsp:
                    fldi = line.split( ',' )
                    hspIdentity, hspPositive, hspGap = None, None, None
                    for elem in fldi:
                        if 'Identities' in elem:
                            hspIdentity = elem.split( '=' )[1].split( '/' )[0].strip() 
                        elif 'Positives' in elem:
                            hspPositive = elem.split( '=' )[1].split( '/' )[0].strip()
                        elif 'Gaps' in elem:  
                            hspGap = elem.split( '=' )[1].split( '/')[0].strip() 
                    hspAlignLen = fldi[0].split( '=' )[1].split( '/')[1].split(' ')[0]
                elif 'Frame' in line and vuHsp:
                    fldf = line.split('=')
                    fldf2 = fldf[1].split('/')
                    if self._allBlast.program in ['BLASTX', 'TBLASTX']:
                        hspQueryFrame = int(fldf2[0].strip())
                        if len(fldf2) > 1:
                            hspHitFrame = int(fldf2[1].strip())
                    elif self._allBlast.program == 'TBLASTN':
                        hspHitFrame = int(fldf2[0].strip())
                    
                elif 'Query' in line and not vuQ:
                    vuQ = True
                    fldq = line.split( ' ' )
                    if self._allBlast.vplus: ## new NCBI blast: 2.2.21+ .... New formating
                        hspQueryFrom = fldq[2]
                        hspQseq = fldq[-3].strip()
                    else:
                        hspQueryFrom = fldq[1]
                        hspQseq = fldq[-2].strip()
                    hspQueryTo = fldq[-1][0:-1]
                    if ancre == -1:
                        ancre = line.find( hspQseq )
                    line = handler.readline()
                    hspMidline = line[ancre:-1]
                elif 'Query' in line and vuQ:
                    fldq = line.split( ' ' )
                    if self._allBlast.vplus:
                        hspQseq += fldq[-3].strip()
                    else:
                        hspQseq += fldq[-2].strip()
                    hspQueryTo = fldq[-1][0:-1]
                    line = handler.readline()
                    hspMidline += line[ancre:-1]
                elif 'Sbjct' in line and not vuS:
                    vuS = True
                    flds = line.split( ' ' )
                    if self._allBlast.vplus:
                        hspHitFrom = fldq[2]
                        hspHseq = flds[-3].strip()
                    else:
                        hspHitFrom = flds[1]
                        hspHseq = flds[-2].strip()
                    hspHitTo = flds[-1]
                    
                elif 'Sbjct' in line and vuS:
                    flds = line.split( ' ' )
                    if self._allBlast.vplus:
                        hspHseq += flds[-3].strip()
                    else:
                        hspHseq += flds[-2].strip()
                    hspHitTo = flds[-1]
                pos = handler.tell()
            else:
                pos = handler.tell()
            line = handler.readline()
        if hit:
            self.hsp( hit, _blast.queryLength, hspScore, hspScore2, hspExpect, hspQueryFrom, hspQueryTo, hspHitFrom, hspHitTo, hspIdentity, hspPositive, hspGap, hspAlignLen, hspQseq, hspHseq, hspMidline, hspQueryFrame, hspHitFrame)
        if 'BLAST'  in line or 'Query=' in line:
            handler.seek( pos )
        elif 'Database:' in line: # multiple blast not well form: pblastall of N. Joly
            while line and 'BLAST'  not in line:
                pos = handler.tell()
                line = handler.readline()
            handler.seek( pos )
        return line

    def _parseUniq( self, handler, database='',  dbSeq='', dbLeter='', program='', version= '' ):
        _blast = BlastRecord( database=database, dbSeq=dbSeq, dbLeter=dbLeter, program=program,version=version)
        #_blast.blastNum = len( self._allBlast.blastS )
        
        #header
        line = handler.readline()
        
        while line and (line[0:19] != 'Sequences producing' and 'No hits found' not in line):
            if 'Database:' in line:
                if 'All non-redundant GenBank' in line:
                    _blast.database = 'nrprot'
                    _blast.databaseAll = 'nrprot'
                    line = handler.readline()
                    while 'letters' not in line:
                        line = handler.readline()
                        line = handler.readline()
                else:
                    fld = line[0:-1].split(':')
                    _blast.database = fld[1].split()[0].strip().lower()
                    _blast.databaseAll = ' '.join( fld[1:] )
                    line = handler.readline()
                fld2 = line.split(';')
                _blast.dbSeq = fld2[0].split('sequences')[0].strip()
                _blast.dbLeter = fld2[1].split('total letters')[0].strip()
            if 'Query=' in line:
                _blast.queryDef = line[6:-1].strip()
                line = handler.readline()
                while 'letters' not in line and 'Length' not in line:
                    _blast.queryDef += ' ' + line[:-1].strip()
                    line = handler.readline()
                if 'letters' in line:
                    reg = re.compile( r'\s*\((\d*)\s*')
                    m = re.match( reg , line )
                    _blast.queryLength = int( m.group(1) )
                elif 'Length' in line:
                    fld = line[0:-1].split('=')
                    _blast.queryLength = int( fld[1] )
            #psiblast pb
#            if 'Results from round' in line:
#                psiblast = True
#                return 'psiblast'
            line = handler.readline()
            
        
        if 'No hits found' in line:
            line = handler.readline()
            line = handler.readline()
            line = handler.readline()
            if 'Database:' in line: # multiple blast not well form: pblastall of N. Joly
                while line and 'BLAST'  not in line:
                    pos = handler.tell()
                    line = handler.readline()
                handler.seek( pos )
            return line, _blast
        
        nbHitComment = False
        if 'Sequences producing' in line:
            fld = line[0:-1].split()
            if fld[-1] == 'N':
                nbHitComment = True
        #comment (newline between comment and sequences )
        line = handler.readline()
        line = handler.readline()
        while line and line != '\n':
            fld = line[0:-1].split()
            if fld and not _blast.hasHit( fld[0] ):
                _blast.hitS.append( BlastHitRecord() )
                hit = _blast.hitS[-1]
                hit.hitNum = len( _blast.hitS )
                hit.dbInfo = fld[0].strip()
                hit.hitDef = ' '.join( fld[:-2] )
                hit.hspS.append( BlastHspRecord( _blast.queryLength ) )
                first_hsp = hit.hspS[-1]
                try:
                    if nbHitComment:
                        first_hsp.hspExpect = float( fld[-2].strip() )
                    else:
                        first_hsp.hspExpect = float( fld[-1].strip() )
                except:
                    # e-164  pas reconnu par python mais 1.e-164 oui
                    if nbHitComment:
                        hspExpect = '1' + fld[-2].strip()
                    else:
                        hspExpect = '1' + fld[-1].strip()
                    first_hsp.hspExpect = float( hspExpect )
                if nbHitComment:
                    first_hsp.hspScore = float( fld[-3].strip() )
                else:
                    first_hsp.hspScore = float( fld[-2].strip() )
            line = handler.readline()
        
        #pairwise
        pos = handler.tell()
        line = handler.readline()
        #flush
        hit = None
        hitParsed, hitCvu, vuHsp, vuQ, vuS, ancre = False, False, False, False, False, -1
        hspQueryFrame, hspHitFrame = 0, 0
        # (pairwise finish by 'Database' for single BLAST )
        # (pairwise finish by 'BLAST?' for multiple BLAST )
        # (pairwise finish by 'Searching' for 'psiblast' and 'phiblast' with multipass)
        while line and 'Database:' not in line and 'BLAST' not in line and 'Query=' not in line and 'Searching' not in line:
            if line[0] == '>' and not hitParsed:
                hspFrame = ''
                # first HSP for new Hit
                fld = line[1:-1].split()
                if _blast.hasHit( fld[0] ):
                    hitParsed = True
                    hit = _blast.getHit( fld[0].strip() )
                    hit.hitDef = line[1:-1].strip()
                    hit.comment = ' '.join( fld[1:] )
                    # erase hspS obtained in comment_line
                    hit.hspS = []
                    pos = handler.tell()
                    line = handler.readline()
                    while 'Length' not in line:
                        hit.hitDef += ' ' + line[1:-1].strip()
                        hit.comment += ' ' + line[:-1].strip()
                        pos = handler.tell()
                        line = handler.readline()
                    handler.seek( pos )
                #else:
                #    print sys.stderr, BlastWarning( "Hit '%s' in pairwise but not in comment"  % fld[0] )
                pos = handler.tell()
            elif line[0] == '>' and  hitParsed:
                # new Hit
                hitParsed, hitCvu, vuHsp = False, False, False
                handler.seek( pos )
                vuQ, vuS, ancre = False, False, -1
                self.hsp( hit, _blast.queryLength, hspScore, hspScore2, hspExpect, hspQueryFrom, hspQueryTo, hspHitFrom, hspHitTo, hspIdentity, hspPositive, hspGap, hspAlignLen, hspQseq, hspHseq, hspMidline, hspQueryFrame, hspHitFrame)
            elif hitParsed:
                # other HSP in Hit
                if 'Length' in line:
                    fldl = line.split('=')
                    hit.hitLength = int( fldl[1].strip() )
                elif 'Score' in line:
                    if not vuHsp:
                        vuHsp = True
                    else:
                        vuQ, vuS, ancre = False, False, -1
                        self.hsp( hit, _blast.queryLength, hspScore, hspScore2, hspExpect, hspQueryFrom, hspQueryTo, hspHitFrom, hspHitTo, hspIdentity, hspPositive, hspGap, hspAlignLen, hspQseq, hspHseq, hspMidline, hspQueryFrame, hspHitFrame)
                        
                    flde = line.split( ',' )
                    fldex = flde[1].split( '=' )
                    hspExpect = fldex[-1].strip()
                    fldsc = flde[0].split( '=' )
                    hspScore = fldsc[1].strip().split( 'bits' )[0][:-1]
                    hspScore2 = fldsc[1].strip().split( '(' )[1][:-1]
                elif 'Identities' in line and vuHsp:
                    fldi = line.split( ',' )
                    hspIdentity, hspPositive, hspGap = None, None, None
                    for elem in fldi:
                        if 'Identities' in elem:
                            hspIdentity = elem.split( '=' )[1].split( '/' )[0].strip() 
                        elif 'Positives' in elem:
                            hspPositive = elem.split( '=' )[1].split( '/' )[0].strip()
                        elif 'Gaps' in elem:  
                            hspGap = elem.split( '=' )[1].split( '/')[0].strip() 
                    hspAlignLen = fldi[0].split( '=' )[1].split( '/')[1].split(' ')[0]
                elif 'Frame' in line and vuHsp:
                    fldf = line.split('=')
                    fldf2 = fldf[1].split('/')
                    if _blast.program in ['BLASTX', 'TBLASTX']:
                        hspQueryFrame = int(fldf2[0].strip())
                        if len(fldf2) > 1:
                            hspHitFrame = int(fldf2[1].strip())
                    elif _blast.program == 'TBLASTN':
                        hspHitFrame = int(fldf2[0].strip())
                    
                elif 'Query' in line and not vuQ:
                    vuQ = True
                    fldq = line.split( ' ' )
                    if _blast.vplus: ## new NCBI blast: 2.2.21+ .... New formating
                        hspQueryFrom = fldq[2]
                        hspQseq = fldq[-3].strip()
                    else:
                        hspQueryFrom = fldq[1]
                        hspQseq = fldq[-2].strip()
                    hspQueryTo = fldq[-1][0:-1]
                    if ancre == -1:
                        ancre = line.find( hspQseq )
                    line = handler.readline()
                    hspMidline = line[ancre:-1]
                elif 'Query' in line and vuQ:
                    fldq = line.split( ' ' )
                    if _blast.vplus:
                        hspQseq += fldq[-3].strip()
                    else:
                        hspQseq += fldq[-2].strip()
                    hspQueryTo = fldq[-1][0:-1]
                    line = handler.readline()
                    hspMidline += line[ancre:-1]
                elif 'Sbjct' in line and not vuS:
                    vuS = True
                    flds = line.split( ' ' )
                    if _blast.vplus:
                        hspHitFrom = fldq[2]
                        hspHseq = flds[-3].strip()
                    else:
                        hspHitFrom = flds[1]
                        hspHseq = flds[-2].strip()
                    hspHitTo = flds[-1]
                    
                elif 'Sbjct' in line and vuS:
                    flds = line.split( ' ' )
                    if _blast.vplus:
                        hspHseq += flds[-3].strip()
                    else:
                        hspHseq += flds[-2].strip()
                    hspHitTo = flds[-1]
                pos = handler.tell()
            else:
                pos = handler.tell()
            line = handler.readline()
        if hit:
            self.hsp( hit, _blast.queryLength, hspScore, hspScore2, hspExpect, hspQueryFrom, hspQueryTo, hspHitFrom, hspHitTo, hspIdentity, hspPositive, hspGap, hspAlignLen, hspQseq, hspHseq, hspMidline, hspQueryFrame, hspHitFrame)
        if 'BLAST'  in line or 'Query=' in line:
            handler.seek( pos )
        elif 'Database:' in line: # multiple blast not well form: pblastall of N. Joly
            while line and 'BLAST'  not in line:
                pos = handler.tell()
                line = handler.readline()
            handler.seek( pos )
        return line, _blast


    # HSPs
    def hsp( self, hit, queryLength, hspScore, hspScore2, hspExpect, hspQueryFrom, hspQueryTo, hspHitFrom, hspHitTo, hspIdentity, hspPositive, hspGap, hspAlignLen, hspQseq, hspHseq, hspMidline, hspQueryFrame, hspHitFrame):
        hit.hspS.append( BlastHspRecord( queryLength ) )
        hsp = hit.hspS[-1]
        hsp.hspNum = len( hit.hspS )
        hsp.hspScore = float( hspScore )
        hsp.hspScore2 = int( hspScore2 )
        try:
            hsp.hspExpect = float( hspExpect )
        except:
            # e-164  pas reconnu par python mais 1.e-164 oui
            hspExpect = '1' + hspExpect.strip()
            hsp.hspExpect = float( hspExpect )
        hsp.hspQueryFrom = int( hspQueryFrom )
        hsp.hspQueryTo = int( hspQueryTo )
        hsp.hspHitFrom = int( hspHitFrom )
        hsp.hspHitTo = int( hspHitTo )
        hsp.hspIdentity = int( hspIdentity )
        if hspPositive:
            hsp.hspPositive = int( hspPositive )
        else:
            # pas de hspPositice dans blastn
            hsp.hspPositive = None
        if hspGap:
            hsp.hspGap = int( hspGap )
        hsp.hspAlignLen = int( hspAlignLen )
        hsp.hspQseq = hspQseq
        hsp.hspHseq = hspHseq
        hsp.hspMidline = hspMidline
        hsp.hspQueryFrame = hspQueryFrame
        hsp.hspHitFrame = hspHitFrame
        
        
class TabularBlastParser:
    def __init__( self ):
        """ Constructor : Return a list of BlastRecord object"""
        self._allBlast = BlastRecordS()

    def parse( self, handler ):
        """Parses the  tabular data handler 
        """
        # Parse the input
        
        self._allBlast.blastS.append( BlastRecord() )
        _blast = self._allBlast.blastS[-1]
        _blast.blastNum = len( self._allBlast.blastS )
        
        line = handler.readline()
        Id=''
        while line and line != '\n':
            fld = line[0:-1].split()
            if _blast.queryId == None:
                _blast.queryId = fld[0]
                _blast.queryDef = fld[0]
                Id = fld[0]
            elif _blast.queryId != fld[0] and Id != fld[0]:
                self._allBlast.blastS.append( BlastRecord() )
                _blast = self._allBlast.blastS[-1]
                _blast.blastNum = len( self._allBlast.blastS )
                _blast.queryId = fld[0]
                _blast.queryDef = fld[0]
                Id=fld[0]
            if fld[2] == '100.00':
                _blast.queryLength = int( fld[7] )
            if _blast.hasHit( fld[1] ):
                hit = _blast.getHit( fld[1] )
                self._hsp( hit, _blast, fld )
            else:
                _blast.hitS.append( BlastHitRecord() )
                hit = _blast.hitS[-1]
                hit.hitAcc = fld[1]
                hit.dbInfo = fld[1]
                hit.hitNum = len( _blast.hitS )
                self._hsp( hit, _blast, fld )
            line = handler.readline()
            
        return self._allBlast

    # HSPs
    def _hsp( self, hit, _blast, field ):
        hit.hspS.append( BlastHspRecord( _blast.queryLength ) )
        hsp = hit.hspS[-1]
        hsp.hspExpect = float( field[10] )
        hsp.hspNum = len( hit.hspS )
        hsp.hspScore = float( field[11] )
        hsp.hspQueryFrom = int( field[6] )
        hsp.hspQueryTo = int(field[7])
        hsp.hspHitFrom = int(field[8])
        hsp.hspHitTo = int( field[9] )
        hsp.hspGap = int( field[5] )
        #hsp.hspQseq = None
        #hsp.hspHseq = None
        #hsp.hspMidline = None
        hsp.hspAlignLen = int( field[3] )
        hsp.hspIdentity = hsp.hspAlignLen - int( field[4] ) - hsp.hspGap
        #hsp.hspPositive = None

class BlastRecordS:
    def __init__( self, database = '', dbSeq='',  dbLeter='', program='', version=''):
        self.database = database
        self.dbSeq = dbSeq
        self.dbLeter = dbLeter
        self.databaseAll = ''
        self.program = program
        self.version = version
        self.vplus = False ## new NCBI blast: 2.2.21+ .... New formating
        self.blastS=[]

    def hasDatabase( self ):
        return self.database != None
        
class BlastRecord( BlastRecordS ):
    def __init__( self, database = '', dbSeq='',  dbLeter='', program='', version=''):
        BlastRecordS.__init__(self, database=database, 
                                  dbSeq=dbSeq, 
                                  dbLeter=dbLeter, 
                                  program=program,
                                  version=version)
        self.blastNum = None
        self.queryId = None
        self.queryDef = ''
        self.queryLength = None
        self.hitS = []
        self.vplus = False ## new NCBI blast: 2.2.21+ .... New formating
        self.program = program
        self.version = version
        self.dbLeter = dbLeter
        self.database = database
        self.dbSeq= dbSeq
        
        
    def hasHitS( self ):
        return self.hitS != []

    def hasHit( self, dbInfo ):
        for hit in self.hitS:
            if hit.dbInfo == dbInfo:
                return True
        return False

    def getHit( self, dbInfo ):
        for hit in self.hitS:
            if hit.dbInfo == dbInfo:
                return hit

    def getHitS( self ):
        return self.hitS

    def getLength( self ):
        return self.queryLength

class BlastHitRecord:
    def __init__( self ):
        self.hitNum = None
        self.hitId = None
        self.hitDef = ''
        self.hitAcc = None
        self.hitLength = None
        self.dbInfo = ''
        self.comment = ''
        self.hspS=[]
        self.hspSplus = []
        self.hspSminus = []

    def getLength( self ):
        return self.hitLength

    def getHspS( self ):
        return self.hspS

    def getNum( self ):
        return self.hitNum

    def extractDBinfo( self ):
        """ return db, acc """
        fld = []
        if not self.dbInfo:
            if self.hitId[0:2] == 'gi':
                return self.hitId.split('|')[2], self.hitAcc
        else:
            field = self.dbInfo.split()
            fld = field[0].split( '|' )
        
        if len(fld) > 1:
            return fld[0], fld[1]
        elif len(fld) == 1:
            return None, fld[0]
        else:
            return None, None

    def extractHspsNotSupperposed( self ):
        import copy
        if len (self.hspS) == 1:
            return self.hspS
        work = copy.deepcopy(self.hspS)
        i = 0
        while i< len(work):
            j = i+1
            while j <  len( work ):
                if work[i].isSupperposed( work[j] ):
                    if work[i].isBetter( work[j]):
                        work.pop( j)
                        j-=1
                    else:
                        work.pop( i)
                        i=-1
                        break
               
                j+=1
            i+=1
        i = 0
        while i< len(work):
            j = i+1
            while j <  len( work ):
                if not work[i].isCompatible( work[j] ):
                    if work[i].isBetter( work[j]):
                        work.pop( j)
                        j-=1
                    else:
                        work.pop( i)
                        i=-1
                        break
               
                j+=1
            i+=1        
        return work

    def splitHspSbyHitFrame( self ):
        self.hspSplus = []
        self.hspSminus = []
        for hsp in self.hspS:
            if hsp.hspHitFrame < 0:
                self.hspSminus.append(hsp)
            else:
                self.hspSplus.append(hsp)
        
    def splitHspSbyQueryFrame( self ):
        self.hspSplus = []
        self.hspSminus = []
        for hsp in self.hspS:
            if hsp.hspQueryFrame < 0:
                self.hspSminus.append(hsp)
            else:
                self.hspSplus.append(hsp)
                
    def sortByHspHitFrom( self, hit = True ):
        self.splitHspSbyHitFrame()
        self.hspSplus.sort(lambda x, y: cmp(x.hspHitFrom, y.hspHitFrom))
        self.hspSminus.sort(lambda x, y: cmp(x.hspHitTo, y.hspHitTo))
        self.hspS = self.hspSminus + self.hspSplus    
        
    def sortByHspQueryFrom( self ):
        self.splitHspSbyQueryFrame()
        self.hspSplus.sort(lambda x, y: cmp(x.hspQueryFrom, y.hspQueryFrom))
        self.hspSminus.sort(lambda x, y: cmp(x.hspQueryTo, y.hspQueryTo))
        self.hspS = self.hspSminus + self.hspSplus
    
class BlastHspRecord( BlastRecord ):
    def __init__( self , queryLength = None ):
        self.hspExpect = None
        self.hspNum = None
        self.hspScore = None
        self.hspScore2 = ''
        self.hspQueryFrom = None
        self.hspQueryTo = None
        self.hspHitFrom = None
        self.hspHitTo = None
        self.hspIdentity = None
        self.hspPositive = None
        self.hspGap = 0
        self.hspQseq = ''
        self.hspHseq = None
        self.hspMidline = None
        self.hspAlignLen = None
        self.hspQueryFrame = 0
        self.hspHitFrame = 0
        self.queryLength = queryLength

    def isQuery( self ):
        return  self.hspHitFrom == 1 and self.hspHitTo == self.queryLength and  self.hspIdentity == self.queryLength

    def partOfQuery( self ):
        return ( self.hspQueryTo - self.hspQueryFrom +1 ) * 100 / float( self.queryLength )

    def isSupperposed( self, hsp ) :
        if self.hspQueryFrom <= hsp.hspQueryFrom <= self.hspQueryTo: #ok
            return True
        elif hsp.hspQueryFrom <= self.hspQueryFrom <= hsp.hspQueryTo: #ok
            return True
        elif self.hspHitFrom <= hsp.hspHitFrom <= self.hspHitTo: #ok
            return True
        elif hsp.hspHitFrom <= self.hspHitFrom <= hsp.hspHitTo: #ok
            return True
        else:
            return False

    def isCompatible( self, hsp ):
#        if (self.hspQueryFrame == 0 and hsp.hspHitFrame == 0) or ( abs(self.hspQueryFrame) > 0 and abs(self.hspHitFrame) > 0 ):
#            if self.isQueryCompatible( hsp ) and self.isDBCompatible( hsp ):
#                return True
#        elif (abs(self.hspQueryFrame) > 0 and hsp.hspHitFrame == 0):
#            if self.isQueryCompatible( hsp ):
#                return True
#        elif (abs(self.hspHitFrame) > 0 and hsp.hspQueryFrame == 0):
#            if self.isDBCompatible( hsp ):
#                return True
        if self.isQueryCompatible( hsp ) and self.isDBCompatible( hsp ):
            return True
        return False
    
    def isQueryCompatible( self, hsp ):      
        if self.hspQueryFrame >= 0 and  hsp.hspQueryFrame >= 0:
            if (self.hspQueryTo < hsp.hspQueryFrom) or (hsp.hspQueryTo < self.hspQueryFrom):
                return True
            else:
                return False
        elif self.hspQueryFrame < 0 and  hsp.hspQueryFrame < 0:
            if (self.hspQueryTo > hsp.hspQueryFrom) or (hsp.hspQueryTo > self.hspQueryFrom):
                return True
            else:
                return False
        return True
      
    def isDBCompatible( self, hsp ):
        if self.hspHitFrame >= 0 and  hsp.hspHitFrame >= 0:
            if (self.hspHitTo < hsp.hspHitFrom) or (hsp.hspHitTo < self.hspHitFrom):
                return True
            else:
                return False
        elif self.hspHitFrame < 0 and  hsp.hspHitFrame < 0:
            if (self.hspHitTo > hsp.hspHitFrom) or (hsp.hspHitTo > self.hspHitFrom):
                return True
            else:
                return False
        return True
      
    def isQueryOverlap( self, hsp, delta ):
        if self.hspQueryFrame >= 0 and hsp.hspQueryFrame >= 0:
            if ((self.hspQueryFrom < hsp.hspQueryFrom <= self.hspQueryTo < hsp.hspQueryTo) and ((hsp.hspQueryFrom + delta) >= self.hspQueryTo) ) or ((hsp.hspQueryFrom < self.hspQueryFrom<= hsp.hspQueryTo < self.hspQueryTo ) and ((self.hspQueryFrom + delta) >= hsp.hspQueryTo) ):
                return True
            else:
                return False
        elif self.hspQueryFrame < 0 and  hsp.hspQueryFrame < 0:
            if ((self.hspQueryFrom > hsp.hspQueryFrom >= self.hspQueryTo > hsp.hspQueryTo) and ((hsp.hspQueryFrom - delta) <= self.hspQueryTo) ) or ((hsp.hspQueryFrom > self.hspQueryFrom >= hsp.hspQueryTo > self.hspQueryTo ) and ((self.hspQueryFrom - delta) <= hsp.hspQueryTo) ):
                return True
            else:
                return False
        return False
        
        
    def isDBOverlap( self, hsp, delta ):
        if self.hspHitFrame >= 0 and hsp.hspHitFrame >= 0:
            if ((self.hspHitFrom < hsp.hspHitFrom <= self.hspHitTo < hsp.hspHitTo) and ((hsp.hspHitFrom + delta) >= self.hspHitTo) ) or ((hsp.hspHitFrom < self.hspHitFrom<= hsp.hspHitTo < self.hspHitTo ) and ((self.hspHitFrom + delta) >= hsp.hspHitTo) ):
                return True
            else:
                return False
        elif self.hspHitFrame < 0 and hsp.hspHitFrame < 0:
            if ((self.hspHitFrom > hsp.hspHitFrom >= self.hspHitTo > hsp.hspHitTo) and ((hsp.hspHitFrom - delta) <= self.hspHitTo) ) or ((hsp.hspHitFrom > self.hspHitFrom >= hsp.hspHitTo > self.hspHitTo ) and ((self.hspHitFrom - delta) <= hsp.hspHitTo) ):
                return True
            else:
                return False
        return False
    
    def isDBincluded( self, hsp):
        if self.hspHitFrame >= 0 and hsp.hspHitFrame >= 0:
            if ((self.hspHitFrom < hsp.hspHitFrom <= hsp.hspHitTo < self.hspHitTo)  ) or ((hsp.hspHitFrom < self.hspHitFrom<= self.hspHitTo < hsp.hspHitTo )):
                return True
            else:
                return False
        elif self.hspHitFrame < 0 and hsp.hspHitFrame < 0:
            if ((self.hspHitFrom > hsp.hspHitFrom >= hsp.hspHitTo > self.hspHitTo) ) or ((hsp.hspHitFrom > self.hspHitFrom >= self.hspHitTo > hsp.hspHitTo )):
                return True
            else:
                return False
        return False
    
    def dbNotReverse( self, hsp ):
        #mandatory: hsps are in same frame (+/-)
        if self.hspQueryFrame >=0:
            if self.hspHitFrame >=0:
                if self.hspQueryFrom <= hsp.hspQueryFrom and self.hspHitFrom <= hsp.hspHitFrom:
                    return True
            else:
                if self.hspQueryFrom <= hsp.hspQueryFrom and self.hspHitFrom >= hsp.hspHitFrom:
                    return True
        else:
            if self.hspHitFrame >=0: 
                if self.hspQueryFrom >= hsp.hspQueryFrom and self.hspHitFrom <= hsp.hspHitFrom:
                    return True
            else:
                if self.hspQueryFrom >= hsp.hspQueryFrom and self.hspHitFrom >= hsp.hspHitFrom:
                    return True
        return False
    
    def isBetter( self, hsp ):
        if self.hspExpect <= hsp.hspExpect:
            return True

###############################

def checkBlastFormat( handler ):
    #### Very, Very succinct
    deb = handler.tell()

    line = handler.readline()
    if line[0] == '<':
        if 'xml' in line.lower():
            handler.seek( deb )
            return 'm7'
    elif len( line.split() ) == 12: ## before 2.19 no gap inside line => 11 fields
        handler.seek( deb )
        return 'm8'
    elif line[0] == '#' or 'BLAST' not in line or '{' in line:
        return False
    else:
        format = False
        line = handler.readline()
        while line:
            if line[0] == '>' or 'No hits found' in line:
                # m1, m2, m3, m4, m5, m6 no pairwise alignment
                handler.seek( deb )
                return 'm0'
            line = handler.readline()
        return format
    return False

def parser( fhin ):

    """  Input blast parsing """
    

    # ===== Check blast output format

    # pairwise (0), xml (7) or tabular (8) format: very succinct cheching

    format = checkBlastFormat( fhin )
    if format == 'm0':
        pairwiseparser = PairwiseBlastParser()
        blastRecordS = pairwiseparser.parse( fhin )
        
        #if blastRecordS == 'psiblast':
        #    format = 'psiblast'
    elif format == 'm7':
        xmlparser = XMLBlastParser()
        blastRecordS = xmlparser.parse( fhin )
    elif format == 'm8':
        tabularparser = TabularBlastParser()
        blastRecordS = tabularparser.parse( fhin )
    elif format == 'noHit':
        print >>sys.stderr, BlastError( "No hits found in blast report" )
        sys.exit()
    else:
        print >>sys.stderr, BlastError( "Not a supported Blast output format: m0, m7, m8" )
        sys.exit()

    fhin.close()
    return blastRecordS, format

def parser2( fhin ):

    """  Input blast parsing """
    

    # ===== Check blast output format

    # pairwise (0), xml (7) or tabular (8) format: very succinct cheching

    format = checkBlastFormat( fhin )
    if format == 'm0':
        pairwiseparser = PairwiseBlastParser()
        blastRecordS = pairwiseparser.parse( fhin )
        
        #if blastRecordS == 'psiblast':
        #    format = 'psiblast'
    elif format == 'm7':
        xmlparser = XMLBlastParser()
        blastRecordS = xmlparser.parse( fhin )
    elif format == 'm8':
        tabularparser = TabularBlastParser()
        blastRecordS = tabularparser.parse( fhin )
    elif format == 'noHit':
        print >>sys.stderr, BlastError( "No hits found in blast report" )
        sys.exit()
    else:
        print >>sys.stderr, BlastError( "Not a supported Blast output format: m0, m7, m8" )
        sys.exit()

    fhin.close()
    return blastRecordS, format


def printLine( outfh, line, car=65 ):
    i = 0
    first = True
    while i < len(line):
        st = line[i:i+car]
        while (i+car)< len(line) and st[-1] != ' ':
            st = st[:-1]
            i -=1
        if first: 
            print >>outfh, '%s' % (st.strip())
            first = False
        else:
            print >>outfh, '           %s' % (st.strip())
        i += car

def printSeq( outfh, program, hitseq, hspHitFrom, hspHitTo, qryseq, hspQueryFrom, hspQueryTo, middle, car=60 ):
    qry3 = ['BLASTX', 'TBLASTX']
    hit3 = ['TBLASTX', 'TBLASTN']
    qj, hj = 1,1
    if program in qry3:
        qj = 3
    if program in hit3:
        hj = 3
    i = 0
    hfr = hspHitFrom
    qfr = hspQueryFrom
    while i < len(middle):
        hsq = hitseq[i:i+car]
        mid = middle[i:i+car]
        qsq = qryseq[i:i+car]
        if hspQueryFrom<hspQueryTo:
            if qfr + (car -qsq.count('-'))*qj -1 < hspQueryTo: 
                
                print >>outfh, ('Query: %%-%ss %%s %%s' %(max(len(str(hspHitFrom)), len(str(hspHitTo)), len(str(hspQueryFrom)), len(str(hspQueryTo))))) % (qfr,  qsq, qfr + (car -qsq.count('-'))*qj-1)
            else:
                print >>outfh, ('Query: %%-%ss %%s %%s' %(max(len(str(hspHitFrom)), len(str(hspHitTo)), len(str(hspQueryFrom)), len(str(hspQueryTo))))) % (qfr,  qsq, hspQueryTo)
        else:
            if qfr - (car-qsq.count('-'))*qj +1 > hspQueryTo: 
                print >>outfh, ('Query: %%-%ss %%s %%s' %(max(len(str(hspHitFrom)), len(str(hspHitTo)), len(str(hspQueryFrom)), len(str(hspQueryTo))))) % (qfr, qsq, qfr - (car -qsq.count('-'))*qj+1)
            else:
                print >>outfh, ('Query: %%-%ss %%s %%s' %(max(len(str(hspHitFrom)), len(str(hspHitTo)), len(str(hspQueryFrom)), len(str(hspQueryTo))))) % (qfr, qsq, hspQueryTo)
        
        print >>outfh, '       %s %s' % (max(len(str(hspHitFrom)), len(str(hspHitTo)), len(str(hspQueryFrom)), len(str(hspQueryTo)))*' ', mid)
        
        if hspHitFrom<hspHitTo:
            if hfr + (car-hsq.count('-'))*hj-1 < hspHitTo:
                print >>outfh, ('Sbjct: %%-%ss %%s %%s' %(max(len(str(hspHitFrom)), len(str(hspHitTo)), len(str(hspQueryFrom)), len(str(hspQueryTo))))) % (hfr,hsq , hfr + (car-hsq.count('-'))*hj-1)
            else:
                print >>outfh, ('Sbjct: %%-%ss %%s %%s' %(max(len(str(hspHitFrom)), len(str(hspHitTo)), len(str(hspQueryFrom)), len(str(hspQueryTo))))) % (hfr,hsq , hspHitTo)
        else:
            if hfr - (car-hsq.count('-'))*hj+1 > hspHitTo:
                print >>outfh, ('Sbjct: %%-%ss %%s %%s' %(max(len(str(hspHitFrom)), len(str(hspHitTo)), len(str(hspQueryFrom)), len(str(hspQueryTo))))) % (hfr,hsq , hfr - (car-hsq.count('-'))*hj +1)
            else:
                print >>outfh, ('Sbjct: %%-%ss %%s %%s' %(max(len(str(hspHitFrom)), len(str(hspHitTo)), len(str(hspQueryFrom)), len(str(hspQueryTo))))) % (hfr,hsq , hspHitTo)
        print >>outfh
        
        if hspQueryFrom<hspQueryTo:
            qfr = qfr + (car -qsq.count('-'))*qj
        else:
            qfr = qfr - (car - qsq.count('-'))*qj
        if hspHitFrom<hspHitTo:
            hfr = hfr + (car-hsq.count('-'))*hj
        else:
            hfr = hfr - (car-hsq.count('-'))*hj
        i += car


def PairwiseBlastPrinterS( blastRecordS, fhout=sys.stdout ):
    for blastRecord in  blastRecordS.blastS:
        print >>fhout, '%s %s' %  (blastRecordS.program, blastRecordS.version )
        print >>fhout
        print >>fhout
        print >>fhout, 'Reference: Altschul, Stephen F., Thomas L. Madden, Alejandro A. Schaffer,'
        print >>fhout, 'Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997),'
        print >>fhout, '"Gapped BLAST and PSI-BLAST: a new generation of protein database search'
        print >>fhout, 'programs",  Nucleic Acids Res. 25:3389-3402.'
        print >>fhout
        print >>fhout, 'Query= %s' % blastRecord.queryDef
        print >>fhout, '         (%s letters)' % blastRecord.queryLength
        print >>fhout, 'Database: %s version' % blastRecordS.database.capitalize()
        print >>fhout, '           %s sequences; %s total letters' %( blastRecordS.dbSeq, blastRecordS.dbLeter )
        print >>fhout
        print >>fhout, 'Searching..................................................done'
        print >>fhout
        print >>fhout
        print >>fhout
        print >>fhout, '                                                                 Score    E'
        print >>fhout, 'Sequences producing significant alignments:                      (bits) Value'
        print >>fhout
        for hit in blastRecord.hitS:
            print >>fhout, '%s...    %s   %s' % (hit.hitDef[0:64], int(hit.hspS[0].hspScore), hit.hspS[0].hspExpect)
        print >>fhout
        for hit in blastRecord.hitS:
            for hsp in hit.hspS:
                printLine( fhout, '>'+hit.hitDef,  car=65 )
                #print >>fhout, '>%s' % (hit.hitDef)
                print >>fhout, '          Length = %s' % (hit.hitLength)
                print >>fhout
                print >>fhout, ' Score = %s bits (%s), Expect = %s' % (hsp.hspScore,hsp.hspScore2, hsp.hspExpect )
                try:
                    identities = ' Identities = %s/%s (%s%%)' % (hsp.hspIdentity, hsp.hspAlignLen, int(hsp.hspIdentity*100/hsp.hspAlignLen))
                except:
                    identities = ''
                try:  
                    positives = ', Positives = %s/%s (%s%%)' % (hsp.hspPositive, hsp.hspAlignLen, int(hsp.hspPositive*100/hsp.hspAlignLen))
                except:
                    positives = ''
                try:
                    gaps = ', Gaps = %s/%s (%s%%)' % (hsp.hspGap, hsp.hspAlignLen, int(hsp.hspGap*100/hsp.hspAlignLen) )
                except:
                    gaps = ''
                print >>fhout, '%s%s%s' % (identities, positives, gaps)
                if hsp.hspQueryFrame and hsp.hspHitFrame:
                    print >>fhout, ' Frame = %s / %s' % (hsp.hspQueryFrame, hsp.hspHitFrame)
                elif hsp.hspQueryFrame:
                    print >>fhout, ' Frame = %s' % hsp.hspQueryFrame
                elif hsp.hspHitFrame:
                    print >>fhout, ' Frame = %s' % hsp.hspHitFrame
                print >>fhout
                printSeq( fhout, blastRecord.program, hsp.hspHseq, hsp.hspHitFrom, hsp.hspHitTo,hsp.hspQseq, hsp.hspQueryFrom, hsp.hspQueryTo,hsp.hspMidline, car=60 )
            
                print >>fhout
    print >>fhout, '  Database: %s version' % blastRecordS.database.capitalize()
    
def PairwiseBlastPrinter( blastRecord, fhout=sys.stdout ):
    print >>fhout, '%s %s' %  (blastRecord.program, blastRecord.version )
    print >>fhout
    print >>fhout
    print >>fhout, 'Reference: Altschul, Stephen F., Thomas L. Madden, Alejandro A. Schaffer,'
    print >>fhout, 'Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997),'
    print >>fhout, '"Gapped BLAST and PSI-BLAST: a new generation of protein database search'
    print >>fhout, 'programs",  Nucleic Acids Res. 25:3389-3402.'
    print >>fhout
    print >>fhout, 'Query= %s' % blastRecord.queryDef
    print >>fhout, '         (%s letters)' % blastRecord.queryLength
    print >>fhout, 'Database: %s version' % blastRecord.database.capitalize()
    print >>fhout, '           %s sequences; %s total letters' %( blastRecord.dbSeq, blastRecord.dbLeter )
    print >>fhout
    print >>fhout, 'Searching..................................................done'
    print >>fhout
    print >>fhout
    print >>fhout
    print >>fhout, '                                                                 Score    E'
    print >>fhout, 'Sequences producing significant alignments:                      (bits) Value'
    print >>fhout
    for hit in blastRecord.hitS:
        print >>fhout, '%s...    %s   %s' % (hit.hitDef[0:64], int(hit.hspS[0].hspScore), hit.hspS[0].hspExpect)
    print >>fhout
    for hit in blastRecord.hitS:
        header = False
        for hsp in hit.hspS:
            if not hsp.hspIdentity:
                continue
            if not header:
                printLine( fhout, '>'+hit.hitDef,  car=65 )
                print >>fhout, '          Length = %s' % (hit.hitLength)
                print >>fhout
                header = True
            #print >>fhout, '>%s' % (hit.hitDef)
            print >>fhout, ' Score = %s bits (%s), Expect = %s' % (hsp.hspScore,hsp.hspScore2, hsp.hspExpect )
            try:
                identities = ' Identities = %s/%s (%s%%)' % (hsp.hspIdentity, hsp.hspAlignLen, int(hsp.hspIdentity*100/hsp.hspAlignLen))
            except:
                identities = ''
            try:  
                positives = ', Positives = %s/%s (%s%%)' % (hsp.hspPositive, hsp.hspAlignLen, int(hsp.hspPositive*100/hsp.hspAlignLen))
            except:
                positives = ''
            try:
                if hsp.hspGap:
                    gaps = ', Gaps = %s/%s (%s%%)' % (hsp.hspGap, hsp.hspAlignLen, int(hsp.hspGap*100/hsp.hspAlignLen) )
                else:
                    gaps = ''
            except:
                gaps = ''
            print >>fhout, '%s%s%s' % (identities, positives, gaps)
            if hsp.hspQueryFrame and hsp.hspHitFrame:
                print >>fhout, ' Frame = %s / %s' % (hsp.hspQueryFrame, hsp.hspHitFrame)
            elif hsp.hspQueryFrame:
                print >>fhout, ' Frame = %s' % hsp.hspQueryFrame
            elif hsp.hspHitFrame:
                print >>fhout, ' Frame = %s' % hsp.hspHitFrame
            print >>fhout
            printSeq( fhout, blastRecord.program, hsp.hspHseq, hsp.hspHitFrom, hsp.hspHitTo,hsp.hspQseq, hsp.hspQueryFrom, hsp.hspQueryTo,hsp.hspMidline, car=60 )
            
            print >>fhout
            print >>fhout
    #print >>fhout, '  Database: %s version' % blastRecord.database.capitalize()
        
def XMLBlastPrinter( blastRecord, fhout=sys.stdout ):
    pass
        
def printer( blastRecordS, format ):
#    print '%s +blastRecordS database: ' % format, blastRecordS.database
#    print '                 dbSeq:  ' , blastRecordS.dbSeq
#    print '                 dbLetter:  ' , blastRecordS.dbLeter
#    print '                 program:  ', blastRecordS.program
#    print '                 version:  ', blastRecordS.version
#    print '                 vplus:    ', blastRecordS.vplus
#    for blastRecord in  blastRecordS.blastS:
#        print '%s    -blastRecord blastNum:    ' % format, blastRecord.blastNum
#        print '                   queryId:     ', blastRecord.queryId
#        print '                   queryDef:    ', blastRecord.queryDef
#        print '                   queryLength: ', blastRecord.queryLength
#        
#        for hit in blastRecord.hitS: # list of BlastHit
#            print '%s        -hit hitNum:    ' % format, hit.hitNum
#            print '               hitId:     ', hit.hitId
#            print '               hitDef:    ', hit.hitDef
#            print '               hitAcc:    ', hit.hitAcc
#            print '               hitLength: ', hit.hitLength
#            print '               dbInfo:    ', hit.dbInfo
#            print '               comment:   ', hit.comment
#            for hsp in hit.hspS:
#                print '%s           -hsp hspExpect:    ' % format, hsp.hspExpect
#                print '                  hspNum:       ', hsp.hspNum
#                print '                  hspScore:     ', hsp.hspScore
#                print '                  hspQueryFrom: ', hsp.hspQueryFrom
#                print '                  hspQueryTo:   ', hsp.hspQueryTo
#                print '                  hspHitFrom:   ', hsp.hspHitFrom
#                print '                  hspHitTo:     ', hsp.hspHitTo
#                print '                  hspIdentity:  ', hsp.hspIdentity
#                print '                  hspPositive:  ', hsp.hspPositive
#                print '                  hspGap:       ', hsp.hspGap
#                print '                  hspQseq:      ', hsp.hspQseq
#                print '                  hspHseq:      ', hsp.hspHseq
#                print '                  hspMidline:   ', hsp.hspMidline
#                print '                  hspAlignLen:  ', hsp.hspAlignLen
#                print '                  hspFrame:     ', hsp.hspQueryFrame, '/', hsp.hspHitFrame            
#                print '                  queryLength:  ', hsp.queryLength
#                
#        print
#        print
    PairwiseBlastPrinterS( blastRecordS )    

    
if __name__ == '__main__':

    import sys
    import os
    # ===== Blast parsing
    blastFile = sys.argv[1]
    try:
        fhin = open( blastFile )
        if os.path.getsize( blastFile ) == 0:
            print >>sys.stderr, ParserError( "Empty file: '%s'" % blastFile)  
            sys.exit( 0 )
    except IOError, err:
        print >>sys.stderr, err
        sys.exit( 0 )
    blastRecordS, format = parser( fhin )
    print format
    printer( blastRecordS, format )
