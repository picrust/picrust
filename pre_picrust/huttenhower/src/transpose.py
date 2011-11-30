#!/usr/bin/env python

import sys

def get( a, i ):
	
	return ( a[i] if ( i < len( a ) ) else "" )

astrLines = map( lambda strLine: strLine.split( "\t" ),
	filter( lambda strLine: strLine, (strLine.rstrip( "\n\r" ) for strLine in sys.stdin) ) )
for iRow in range( len( astrLines[0] ) ):
	print( "\t".join( get( astrLines[iCol], iRow ) for iCol in range( len( astrLines ) ) ) )
