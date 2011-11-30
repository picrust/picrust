#!/usr/bin/env python

from fir import CFir
import csv
import logging
import re
import sys

c_strName	= "name"

def funcNames( pNode, hashhashNames, hashCache ):

	hashNames = {}
	strName = pNode.get( c_strName )
	if strName:
		hashNames[strName] = 1 + hashNames.get( strName, 0 )

	apChildren = pNode.get_children( )
	if apChildren:
		for pChild in apChildren:
			hashChild = funcNames( pChild, hashhashNames, hashCache )
			for strName, iName in hashChild.items( ):
				hashNames[strName] = iName + hashNames.get( strName, 0 )
	elif ( hashCache != None ) and strName:
		hashCache.setdefault( strName, [] ).append( pNode )

	if pNode.get_id( ):
		hashhashNames[pNode] = hashNames
	return hashNames

def funcJaccard( hashOne, hashTwo ):
	
	iInt = iUni = 0
	for pOne, iOne in hashOne.items( ):
		iTwo = hashTwo.get( pOne )
		if iTwo:
			iInt += min( iOne, iTwo )
			iUni += max( iOne, iTwo )
		else:
			iUni += iOne
	if not iInt:
		return 0
	for pTwo, iTwo in hashTwo.items( ):
		if pTwo not in hashOne:
			iUni += iTwo
	return ( float(iInt) / iUni )

def funcMatch( pTo, hashFrom, hashhashUnnamed, dMax ):

	logging.debug( "funcMatch( %s, %s, %g )" % (pTo.get_id( ), hashFrom if ( len( hashFrom ) < 10 ) else "-", dMax) )
	hashTo = hashhashUnnamed.get( pTo, {} )
	dRet = funcJaccard( hashFrom, hashTo )
	logging.debug( "Jaccard: %g" % dRet )
	pRet = None
	if ( not hashTo ) or ( dRet and ( dRet >= dMax ) ):
		pRet = pTo if hashTo else None
		for pChild in pTo.get_children( ):
			dCur, pCur = funcMatch( pChild, hashFrom, hashhashUnnamed, dRet )
			if dCur > dRet:
				dRet = dCur
				pRet = pCur
	return (dRet, pRet)

def funcMerge( pFrom, pTo, hashhashNamed, hashhashUnnamed, hashCache, aiDone ):

	if not ( aiDone[0] % 1000 ):
		logging.warning( "%d/%d" % (aiDone[0], len( hashhashNamed )) )
	aiDone[0] += 1
	logging.debug( "funcMerge( %s, %s )" % (pFrom.get_id( ), pTo.get_id( )) )
	
	pMax = None
	if not pFrom.get_children( ):
		apCache = hashCache.get( pFrom.get( c_strName ), [] )
		if len( apCache ) == 1:
			pMax = apCache[0]
			logging.debug( "Cache hit" )
	if not pMax:
		dMax, pMax = funcMatch( pTo, hashhashNamed.get( pFrom, {} ), hashhashUnnamed, 0 )
	if pMax:
		logging.info( "Matched (%s, %s) to (%s, %s)" % (pFrom.get_id( ), pFrom.get( c_strName ),
			pMax.get_id( ), pMax.get( c_strName )) )
		for strKey in pFrom.get_keys( ):
			setTo = pMax.get( strKey )
			if not setTo:
				setTo = set()
			setTo.add( pFrom.get( strKey ) )
			pMax.set( strKey, setTo )
	else:
		logging.debug( "Could not match: (%s, %s)" % (pFrom.get_id( ), pFrom.get( c_strName )) )
		pMax = pTo
	for pChild in pFrom.get_children( ):
		funcMerge( pChild, pMax, hashhashNamed, hashhashUnnamed, hashCache, aiDone )

def funcEnset( pNode ):
	
	for pKey in pNode.get_keys( ):
		pNode.set( pKey, set([pNode.get( pKey )]) )
	for pChild in pNode.get_children( ):
		funcEnset( pChild )

def funcUnset( pNode ):
	
	for pKey in pNode.get_keys( ):
		pNode.set( pKey, "\t".join( pNode.get( pKey ) ) )
	for pChild in pNode.get_children( ):
		funcUnset( pChild )

if len( sys.argv ) != 2:
	raise Exception( "Usage: merge_trees.py <named.fir> < <unnamed.fir>" )
strNamed = sys.argv[1]
logging.basicConfig( format = "%(asctime)s %(levelname)s:%(message)s", level = logging.WARNING )

logging.warning( "Opening named tree" )
pNamed = CFir( )
pNamed.open( open( strNamed ) )
logging.warning( "Opening unnamed tree" )
pUnnamed = CFir( )
pUnnamed.open( sys.stdin )

hashhashNamed = {}
hashhashUnnamed = {}
hashCache = {}
for pTree, hashhashTree, hashCur in ((pNamed, hashhashNamed, None), (pUnnamed, hashhashUnnamed, hashCache)):
	funcNames( pTree.get_roots( ), hashhashTree, hashCur )

funcEnset( pUnnamed.get_roots( ) )
funcMerge( pNamed.get_roots( ), pUnnamed.get_roots( ), hashhashNamed, hashhashUnnamed, hashCache, [0] )
funcUnset( pUnnamed.get_roots( ) )

pUnnamed.save( sys.stdout )
