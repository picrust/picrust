#!/usr/bin/env python

from fir import CFir
import csv
import fir
import fir_support
import logging
import sys

if len( sys.argv ) != 2:
	raise Exception( "Usage: taxa2fir.py <tree.fir> < <taxa.pcl>" )
strTree = sys.argv[1]
logging.basicConfig( format = "%(asctime)s %(levelname)s:%(message)s", level = logging.WARNING )

pTree = CFir( )
pTree.open( open( strTree ) )
hashTree = {}
fir_support.funcSpecies( pTree.get_roots( ), hashTree, True )

astrTaxa = []
ahashFir = []
hashCache = {}
fFirst = True
for astrLine in csv.reader( sys.stdin, csv.excel_tab ):
	strID = astrLine[0]
	if fFirst:
		fFirst = False
		continue

	astrID = strID.split( "|" )
	for i in range( len( astrID ) - 2, -1, -1 ):
		hashNodes = hashCache.get( astrID[i] )
		if hashNodes != None:
			break
		apNodes = fir_support.funcFind( astrID[i], hashTree )
		if apNodes:
			strCache = astrID[i]
			break
	if not ( hashNodes or apNodes ):
		logging.warn( "Unknown node: %s" % strID )
		continue

	if hashNodes:
		logging.info( "Matched %s to cache" % strID )
	else:
		logging.info( "Matched %s to %s" % (strID, [( "%s:%s" % (p.get_id( ), p.get( fir_support.c_strName )) ) for p in apNodes]) )
		hashNodes = {}
		dWeight = 1.0 / len( apNodes )
		for pNode in apNodes:
			hashNodes[pNode] = dWeight + hashNodes.get( pNode, 0 )
		hashCache[strCache] = hashNodes
	astrTaxa.append( strID )
	ahashFir.append( hashNodes )
setpNodes = set()
for hashNodes in ahashFir:
	setpNodes |= set(hashNodes.keys( ))
apNodes = list(setpNodes)

csvw = csv.writer( sys.stdout, csv.excel_tab )
csvw.writerow( ["TID"] + [p.get_id( ) for p in apNodes] )
for iRow in range( len( astrTaxa ) ):
	if not ( iRow % 500 ):
		logging.warn( "%s/%s" % (iRow, len( astrTaxa )) )
	csvw.writerow( [astrTaxa[iRow]] + [ahashFir[iRow].get( p, "" ) for p in apNodes] )
