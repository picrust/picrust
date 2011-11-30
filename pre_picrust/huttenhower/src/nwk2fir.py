#!/usr/bin/env python

import sys
from ete2 import Tree

class CNamer():

	def __init__( pSelf ):

		pSelf.m_hashNames = {}

	def enname( pSelf, pNode ):

		strRet = pNode.name
		if strRet == "NoName":
			strRet = pSelf.m_hashNames.setdefault( pNode, "INT" +
				str(len( pSelf.m_hashNames ) + 1) )

		return strRet

if len( sys.argv ) != 2:
	raise Exception( "Usage: nwk2fir.py <tree.nwk>" )
strTree = sys.argv[1]

pNamer = CNamer( )
pTree = Tree( strTree )
for pNode in pTree.traverse( ):
	apChildren = pNode.get_children( )
	if len( apChildren ) == 0:
		continue
	sys.stdout.write( "%s" % pNamer.enname( pNode ) )
	for pChild in apChildren:
		sys.stdout.write( "	%s	%g" % (pNamer.enname( pChild ),
			pChild.dist) )
	sys.stdout.write( "\n" )
