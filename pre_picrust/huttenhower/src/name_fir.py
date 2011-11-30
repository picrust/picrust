#!/usr/bin/env python

from fir import CFir
import csv
import logging
import sys

c_strName	= "name"

def funcPrefix( astr ):
	
	strRet = ""
	if astr:
		astr = sorted( astr )
		strRet = astr[0]
		for i in range( 1, len( astr ) ):
			fHit = False
			for j in range( 0, min( len( strRet ), len( astr[i] ) ) ):
				if strRet[j] != astr[i][j]:
					fHit = True
					break
			if fHit:
				strRet = strRet[:j]
				if not strRet:
					break
	return strRet

def funcName( pNode ):

	if not pNode.get( c_strName ):
		astrNames = filter( None, (p.get( c_strName ) for p in pNode.get_descendants( )) )
		strName = funcPrefix( astrNames ).strip( )
		if strName.find( " " ) > 0:
			pNode.set( c_strName, strName )
	for pChild in pNode.get_children( ):
		funcName( pChild )

def funcUnname( pNode ):
	
	strName = pNode.get( c_strName )
	for pChild in pNode.get_children( ):
		funcUnname( pChild )
		if strName:
			strChild = pChild.get( c_strName )
			if strName == strChild:
				pChild.delete( c_strName )

if len( sys.argv ) != 2:
	raise Exception( "Usage: name_fir.py <tree.fir> < map.txt" )
strTree = sys.argv[1]
logging.basicConfig( format = "%(asctime)s %(levelname)s:%(message)s", level = logging.DEBUG )

logging.info( "Opening named tree" )
pTree = CFir( )
pTree.open( open( strTree ) )

logging.info( "Naming leaves" )
for astrLine in csv.reader( sys.stdin, csv.excel_tab ):
	pNode = pTree.get( astrLine[0] )
	if not pNode:
		sys.stderr.write( "Unknown node: %s\n" % astrLine[0] )
		continue
	pNode.set( c_strName, astrLine[1] )
logging.info( "Naming ancestors" )
funcName( pTree.get_roots( ) )
logging.info( "Reducing ancestors" )
funcUnname( pTree.get_roots( ) )

pTree.save( sys.stdout )
