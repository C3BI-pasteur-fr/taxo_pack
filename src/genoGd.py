import gd, os, sys




fontlist = ['courier.ttf', 'gdMediumBoldFont', ]


for f in fontlist:
    try:
        FONT = os.path.basename(f)
        break
    except IOError, str:
        print >>sys.stderr, str

fontpath = '/Users/maufrais/Developpements/Blast2LAT/lib'        

try:
    GDFONTPATH = os.environ['RANKOPTIMIZERLIB']
except:
    print >>sys.stderr, RankOptimizerError("The mandatory RANKOPTIMIZERLIB environment variable is not defined")
    sys.exit(1)


if GDFONTPATH not in sys.path:
    sys.path.append( GDFONTPATH )

os.environ['GDFONTPATH'] = GDFONTPATH


try:
    FONT
except NameError:
    print "No fonts found"
    sys.exit(1)


class genoGD:
    def __init__( self, fontsize = 10, angle = 0, sizeX = 700, limit = 10 , nbQuery = None):
        self.fontsize = fontsize
        self.angle = angle
        self.limit = limit
        self.nbQuery = nbQuery

        self.recX = self.limit + (self.fontsize*10)
        self.pixY = 2*self.fontsize
        self.nbligne = 1

        self.sizeX = sizeX
        self.sizeY = self.estimatSizeY()

        self.picture = gd.image( ( self.sizeX+self.recX+2*self.limit, self.sizeY ) )

        self.color = {'white' : self.picture.colorAllocate( (255, 255, 255) ),
                      'grey': self.picture.colorAllocate( (125, 125, 125) ),
                      'lightgrey': self.picture.colorAllocate( (211, 211, 211) ),
                      'green' : self.picture.colorAllocate( (0,205,0) ),
                      'lightgreen' : self.picture.colorAllocate( (152,251,152) ),   
                      'orange' : self.picture.colorAllocate( (255,187,0) ),
                      'blue' : self.picture.colorAllocate( (0, 0, 255) ),
                      'lightblue' : self.picture.colorAllocate( (0, 100, 200) ),
                      'red' : self.picture.colorAllocate( (255, 0, 0) ),
                      'lightred' : self.picture.colorAllocate( (255, 128, 128) ),
                      'pink' : self.picture.colorAllocate( (255, 0, 255) ),
                      'yellow' : self.picture.colorAllocate( (255, 255, 0) ),
                      'black' : self.picture.colorAllocate( (0, 0, 0) )
                      }

        
        #self.picture.interlace(1)
        #self.picture.setThickness(2)
    
    
    def estimatSizeY( self ):
        if self.nbQuery is None:
            print "Error of plot's size: nothing"
            sys.exit(1)

        legendSize = self.pixY * 9
        #c = ((self.sizeX - (2*self.limit)) /10.)
        #b = nbQuery *c 
        a = self.nbQuery * self.pixY
        return int ( a ) + legendSize
    
    def plotSeq( self, seqlen = 700, start = 0, stop = 700, header = 'TUTUTUTU', color = 'red', up = True ):
        """
        Display sequence
        """
        self.picture.string_ttf(FONT, self.fontsize, self.angle, ( self.limit , (self.pixY*self.nbligne)+4), header, self.color['red'] )
        #self.picture.filledRectangle( ( self.recX, self.pixY-4) , (self.recX+self.sizeX, (self.pixY*self.nbligne)+4), self.color[color] )
        self.picture.filledRectangle( ( start+self.recX, self.pixY*self.nbligne-4) , (stop+self.recX, (self.pixY*self.nbligne)+4), self.color[color] )
        unit = float((stop)-(start))/seqlen
        delta = 1
        for i in range ( 10, seqlen+1, 10 ):
            tmp = seqlen / i
            if (tmp/5) >10:
                delta = i*5
            elif tmp > 10:
                delta = i
        if up:
            y1 = self.pixY*self.nbligne-2
            self.nbligne += 1
            y2 = self.pixY*self.nbligne-6
            yleg= y2 + self.pixY/2 +1
        else:
            y1 = self.pixY*self.nbligne+2
            self.nbligne -= 1
            y2 = self.pixY*self.nbligne+6
            yleg= y2 - self.pixY/4
        
        for i in range (seqlen/delta+1):
            pos = i * delta * unit
            x1 = x2 = start+self.recX + int(pos)
            self.picture.filledRectangle(( x1,y1),( x2,y2),self.color['red'])
            nm = len(str(i*delta))*self.fontsize
            if i%5 == 0:
                self.picture.string_ttf(FONT, self.fontsize, self.angle, ( x1-nm/3, yleg), str(i*delta), self.color['black'] )
                                     
    def _nbligneP(self):
        self.nbligne += 1
        
    def _nbligneM(self):
        self.nbligne -= 1
        
    def plotHit( self, start = 0, stop = 700, header='TOTO', color = 'blue', head = True, delta = 4 ):
        if head:
            self.picture.string_ttf(FONT, self.fontsize, self.angle, ( self.limit , (self.pixY)*self.nbligne), header, self.color['black'] )
        self.picture.filledRectangle( ( start+self.recX, self.pixY*self.nbligne-delta) , (stop+self.recX, (self.pixY*self.nbligne)+delta), self.color[color] )
        

    def legend( self, start = 0, stop = 700, legendScore = [('0','blue'),('50','grey'),('100','orange'),('150','green'),('200','red') ] ):
        """
        specific legend
        """
        self.nbligne += 3

        self.picture.string_ttf(FONT, self.fontsize, self.angle, ( self.limit , (self.pixY*self.nbligne)+4), "Scoring +", self.color['red'] )
        #self.nbligne += 1
        #self.picture.string_ttf(FONT, self.fontsize, self.angle, ( self.limit , (self.pixY*(self.nbligne+1))+4), "colors", self.color['red'] )
        nbcol = len( legendScore )
        lencol = (stop - start) / nbcol
        pos = start
        prec = 0
        for vcol in legendScore:
            self.picture.filledRectangle( ( pos+start+self.recX, self.pixY*self.nbligne-4) , (stop+self.recX, (self.pixY*self.nbligne)+4), self.color[vcol[1]] )
            #self.picture.string_ttf(FONT, self.fontsize, self.angle, ( pos+start+self.recX , (self.pixY*(self.nbligne+1))+4), "S>%s"%vcol[0], self.color['black'] )
            pos += lencol
    
    def legendInv( self, start = 0, stop = 700, legendScore = [('0','lightblue'),('50','lightgrey'),('100','yellow'),('150','lightgreen'),('200','lightred') ] ):
        """
        specific legend
        """
        self.nbligne += 1

        self.picture.string_ttf(FONT, self.fontsize, self.angle, ( self.limit , (self.pixY*self.nbligne)+4), "Scoring -", self.color['red'] )
        #self.nbligne += 1
        self.picture.string_ttf(FONT, self.fontsize, self.angle, ( self.limit , (self.pixY*(self.nbligne+1))+4), "colors", self.color['red'] )

        nbcol = len( legendScore )
        lencol = (stop - start) / nbcol
        pos = start
        prec = 0
        for vcol in legendScore:
            self.picture.filledRectangle( ( pos+start+self.recX, self.pixY*self.nbligne-4) , (stop+self.recX, (self.pixY*self.nbligne)+4), self.color[vcol[1]] )
            self.picture.string_ttf(FONT, self.fontsize, self.angle, ( pos+start+self.recX , (self.pixY*(self.nbligne+1))+4), "S>%s"%vcol[0], self.color['black'] )
            pos += lencol

    def legendMerge( self, start = 0, stop = 700, legendScore = [('minus','blue'),('plus','lightred')]):
        self.nbligne += 1

        self.picture.string_ttf(FONT, self.fontsize, self.angle, ( self.limit , (self.pixY*self.nbligne)+4), "Scoring", self.color['red'] )
        #self.nbligne += 1
        self.picture.string_ttf(FONT, self.fontsize, self.angle, ( self.limit , (self.pixY*(self.nbligne+1))+4), "colors", self.color['red'] )

        nbcol = len( legendScore )
        lencol = (stop - start) / nbcol
        pos = start
        prec = 0
        for vcol in legendScore:
            self.picture.filledRectangle( ( pos+start+self.recX, self.pixY*self.nbligne-4) , (stop+self.recX, (self.pixY*self.nbligne)+4), self.color[vcol[1]] )
            self.picture.string_ttf(FONT, self.fontsize, self.angle, ( pos+start+self.recX , (self.pixY*(self.nbligne+1))+4), "%s"%vcol[0], self.color['black'] )
            pos += lencol
        

if __name__ == '__main__':
    
    #start = 0, stop = 700, header='TOTO', color = 'blue'
    #'white','green' ,'blue','lithblue' ,'lithyellow' ,'yellow' ,'red' ,'black' 
    lenSeq = 3945
    Hits = [[(0,340,'toto1','grey'), (600,920,'toto1','blue'), (1000,1250,'toto1','yellow'), (2000,2500,'toto1', 'pink'),(2600, lenSeq,'toto1', 'red')],
            [(0,340,'toto2','green'), (600,920,'toto2','orange'), (1000,1250,'toto2','yellow'), (2000,2500,'toto2', 'pink'),(2600, lenSeq,'toto2', 'red')],
            [(0,340,'toto3','green'), (600,920,'toto3','blue'), (1000,1250,'toto3','yellow'), (2000,2500,'toto3', 'pink'),(2600, lenSeq,'toto3', 'red')]]
    
    
    expend = False
    #expend = True
    if expend:
        sizeX = 700 *2
        nbQuery = len( Hits )*5 +1
    else:
        sizeX = 700
        nbQuery = len( Hits )
    rpt = float(sizeX) / lenSeq
    
    gdPlot = genoGD( sizeX = sizeX, nbQuery = nbQuery )
    gdPlot.plotSeq()
    for hit in Hits:
        head = True
        for h in hit:
            gdPlot.plotHit(int(h[0]*rpt), int(h[1]*rpt), h[2], h[3], head=head)
            head = False
            if expend:
                gdPlot._nbligneP()
        gdPlot._nbligneP()
    gdPlot.legend()
    f=open( 'toto.png',"w")
    gdPlot.picture.writePng(f)
    f.close()
