#!/usr/bin/env python

from fir import CFir
import re
import sys

pTree = CFir( )
aastrTaxa = []
apTaxon = []
for strLine in sys.stdin:
	mtch = re.search( r'^(\s*)([^\t]+)\t(\d+)', strLine.rstrip( ) )
	if not mtch:
		continue
	strSpaces, strTaxon, strID = mtch.groups( )
	while len( apTaxon ) > len( strSpaces ):
		apTaxon.pop( )
	pNode = apTaxon[-1].add_child( strID ) if apTaxon else pTree.add_root( strID )
	pNode.set( "name", strTaxon )
	apTaxon.append( pNode )
pTree.save( sys.stdout )
