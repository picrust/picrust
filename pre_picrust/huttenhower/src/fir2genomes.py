#!/usr/bin/env python

from fir import CFir
import csv
import fir_support
import logging
import sys

if len( sys.argv ) != 2:
	raise Exception( "Usage: fir2genomes.py <tree.fir> < <genomes.tsv>" )
strTree = sys.argv[1]
logging.basicConfig( format = "%(asctime)s %(levelname)s:%(message)s", level = logging.INFO )

pTree = CFir( )
pTree.open( open( strTree ) )
hashTree = {}
fir_support.funcSpecies( pTree.get_roots( ), hashTree )

hashGenomes = {}
hashCounts = {}
astrGenes = None
for astrLine in csv.reader( sys.stdin, csv.excel_tab ):
	strID, astrData = astrLine[0], astrLine[1:]
	if not astrGenes:
		astrGenes = astrData
		continue
	import random
#	if random.random( ) > 0.005:
#		continue

	apNodes = fir_support.funcFind( strID, hashTree )
	if not apNodes:
		logging.warn( "Unknown node: %s" % strID )
		continue
	logging.info( "Matched %s to %s" % (strID, [( "%s:%s" % (p.get_id( ), p.get( fir_support.c_strName )) ) for p in apNodes]) )
	adData = [float(s) for s in astrData]
	for pNode in apNodes:
		dWeight = 1.0 / len( apNodes )
		while pNode:
			hashCounts[pNode] = hashCounts.get( pNode, 0 ) + dWeight
			adGenome = hashGenomes.setdefault( pNode, [0] * len( astrGenes ) )
			for i in range( len( adData ) ):
				if adData[i]:
					adGenome[i] += dWeight * adData[i]
			dWeight *= 2 ** -pNode.m_dLength
			pNode = pNode.m_pParent
csvw = csv.writer( sys.stdout, csv.excel_tab )
csvw.writerow( ["TID"] + astrGenes )
for pNode, dTotal in hashCounts.items( ):
	if ( not pNode.get_id( ) ) or ( dTotal == 0 ):
		continue
	csvw.writerow( [pNode.get_id( )] + [( d / dTotal ) for d in hashGenomes[pNode]] )
