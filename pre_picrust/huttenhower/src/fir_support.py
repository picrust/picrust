#!/usr/bin/env python

import logging
import re

c_strName	= "name"

def funcNames( pNode, hashTree ):
	
	strName = pNode.get( c_strName )
	if strName:
		logging.debug( "Named: %s = %s" % (strName, pNode.get_id( )) )
		for strToken in strName.split( "\t" ):
			hashTree.setdefault( strToken, [] ).append( pNode )
	for pChild in pNode.get_children( ):
		funcNames( pChild, hashTree )

def funcSpecies( pNode, hashTree, fGenera = False ):
	
	strName = pNode.get( c_strName )
	if strName:
		logging.debug( "Named: %s = %s" % (pNode.get_id( ), strName) )
		astrName = strName.split( "\t" )
		setstrNames = set(astrName)
		for strToken in astrName:
			astrToken = strToken.split( " " )
			setstrNames.add( " ".join( astrToken[:2] ) if ( len( astrToken ) > 1 ) else astrToken[0] )
			if fGenera:
				setstrNames.add( astrToken[0] )
		for strToken in setstrNames: 
			logging.debug( "Tagged: %s" % strToken )
			hashTree.setdefault( strToken, [] ).append( pNode )
	for pChild in pNode.get_children( ):
		funcSpecies( pChild, hashTree, fGenera )

def funcFind( strID, hashTree ):

	apRet = hashTree.get( strID )
	strCur = strID
	logging.debug( "Searching: \"%s\"" % strCur )
	while len( strCur ) and ( not apRet ):
		strCur = re.sub( r'\s*\S+\s*$', "", strCur )
		logging.debug( "Searching: \"%s\"" % strCur )
		apRet = hashTree.get( strCur )
	return apRet
