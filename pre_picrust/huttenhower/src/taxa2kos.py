#!/usr/bin/env python

from fir import CFir
import csv
import fir
import fir_support
import logging
import sys

if len( sys.argv ) != 3:
	raise Exception( "Usage: taxa2kos.py <tree.fir> <kos.pcl> < <taxa.pcl>" )
strTree, strKOs = sys.argv[1:3]
logging.basicConfig( format = "%(asctime)s %(levelname)s:%(message)s", level = logging.INFO )

pTree = CFir( )
pTree.open( open( strTree ) )
hashTree = {}
fir_support.funcSpecies( pTree.get_roots( ), hashTree, True )

hashKOIDs = {}
hashhashKOs = {}
for astrLine in csv.reader( open( strKOs ), csv.excel_tab ):
	strID, astrData = astrLine[0], astrLine[1:]
	hashhashKOs[strID] = hashKOs = {}
	for strToken in astrData:
		strKO, strWeight = strToken.split( ":" )
		iKO = hashKOIDs.get( strKO )
		if iKO == None:
			hashKOIDs[strKO] = iKO = len( hashKOIDs )
		hashKOs[iKO] = float(strWeight)

astrKOs = [None] * len( hashKOIDs )
for strKO, iKO in hashKOIDs.items( ):
	astrKOs[iKO] = strKO
for pNode, hashKOs in hashhashKOs.items( ):
	adKOs = [0] * len( astrKOs )
	for iKO, dKO in hashKOs.items( ):
		adKOs[iKO] = dKO
	hashhashKOs[pNode] = adKOs

aadCounts = []
hashCache = {}
fFirst = True
csvw = csv.writer( sys.stdout, csv.excel_tab )
for astrLine in csv.reader( sys.stdin, csv.excel_tab ):
	strID, astrData = astrLine[0], astrLine[1:]
	if fFirst:
		fFirst = False
		csvw.writerow( astrLine )
		for i in range( len( astrKOs ) ):
			aadCounts.append( [0] * len( astrData ) )
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
			while pNode:
				if pNode.get_id( ) in hashhashKOs:
					hashNodes[pNode] = dWeight + hashNodes.get( pNode, 0 )
					break
				pNode = pNode.m_pParent
		hashCache[strCache] = hashNodes
	
	adData = [float(s) for s in astrData]
	for pNode, dNode in hashNodes.items( ):
		adKOs = hashhashKOs[pNode.get_id( )]
#		sys.stderr.write( "%s\n" % [pNode.get_id( ), dNode, adKOs[:10]] )
		for iKO in range( len( adKOs ) ):
			dKO = adKOs[iKO]
			adCounts = aadCounts[iKO]
			for i in range( len( adData ) ):
				adCounts[i] += adData[i] * dKO * dNode
for iKO in range( len( astrKOs ) ):
	adCounts = aadCounts[iKO]
	if sum( adCounts ):
		csvw.writerow( [astrKOs[iKO]] + adCounts )
