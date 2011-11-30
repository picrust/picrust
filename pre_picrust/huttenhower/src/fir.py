#!/usr/bin/env python

import csv
import StringIO
import sys
import unittest

class CFir():

	class CNode():

		def __init__( pSelf, strID ):

			pSelf.m_strID = strID
			pSelf.m_pParent = None
			pSelf.m_dLength = 0
			pSelf.m_setChildren = set()
			pSelf.m_setDescendants = None
			pSelf.m_setLeaves = None
			pSelf.m_hashProps = {}
			
		def __cmp__( pSelf, pOther ):
			
			return cmp( pSelf.get_id( ), pOther.get_id( ) )
		
		def __hash__( pSelf ):
			
			return hash(pSelf.get_id( ))
			
		def set( pSelf, pKey, pValue ):
			
			pSelf.m_hashProps[pKey] = pValue
			
		def get( pSelf, pKey ):
			
			return pSelf.m_hashProps.get( pKey )
		
		def delete( pSelf, pKey ):
			
			del pSelf.m_hashProps[pKey]
		
		def get_keys( pSelf ):
			
			return pSelf.m_hashProps.keys( )
			
		def get_parent( pSelf ):
			
			return pSelf.m_pParent
		
		def get_id( pSelf ):
			
			return pSelf.m_strID

		def add_child( pSelf, pChild, dLength = 1 ):

			if type( pChild ) == str:
				pChild = CFir.CNode( pChild )
			pChild.m_dLength = dLength
			pChild.m_pParent = pSelf
			pSelf.m_setChildren.add( pChild )
			return pChild

		def save( pSelf, ofil ):

			apChildren = sorted( pSelf.get_children( ) )
			if not( apChildren or pSelf.m_hashProps ):
				return

			if pSelf.get_id( ):
				ofil.write( "%s" % pSelf.get_id( ) )
				for pChild in apChildren:
					ofil.write( "	%s	%g" % (pChild.get_id( ), pChild.get_length( )) )
				ofil.write( "\n" )
				for pKey, pValue in pSelf.m_hashProps.items( ):
					ofil.write( "	%s	%s\n" % (str(pKey), str(pValue)) )
			for pChild in apChildren:
				pChild.save( ofil )
				
		def get_children( pSelf ):
			
			return list(pSelf.m_setChildren)
		
		def get_length( pSelf ):
			
			return pSelf.m_dLength

		def get_descendants( pSelf ):

			if not pSelf.m_setDescendants:
				pSelf.m_setDescendants = set(pSelf.m_setChildren)
				for pChild in pSelf.m_setChildren:
					pSelf.m_setDescendants |= pChild.get_descendants( )

			return set(pSelf.m_setDescendants)

		def get_leaves( pSelf ):

			if not pSelf.m_setLeaves:
				pSelf.m_setLeaves = set()
				for pDesc in pSelf.get_descendants( ):
					if len( pDesc.m_setChildren ) == 0:
						pSelf.m_setLeaves.add( pDesc )

			return set(pSelf.m_setLeaves)

		def collapse( pSelf ):

			while len( pSelf.get_children( ) ) == 1:
				for pChild in pSelf.get_children( ):
					pSelf.m_dLength += pChild.get_length( )
					pSelf.m_setChildren = pChild.get_children( )
					pSelf.m_setDescendants = pChild.get_descendants( )
					pSelf.m_setLeaves = pChild.get_leaves( )
			for pChild in pSelf.get_children( ):
				pChild.collapse( )

		def save_newick( pSelf, ofil ):

			if pSelf.get_children( ):
				ofil.write( "(" )
			fFirst = True
			for pChild in pSelf.get_children( ):
				if not fFirst:
					ofil.write( "," )
				fFirst = False
				pChild.save_newick( ofil )
			if pSelf.get_children( ):
				ofil.write( ")" )
			ofil.write( "%s:%g" % (pSelf.get_id( ), pSelf.get_length( )) )
			
		def remove( pSelf, pNode ):
			
			if pNode in pSelf.get_descendants( ):
				pSelf.m_setDescendants.remove( pNode )
				pSelf.m_setLeaves.remove( pNode )
				if pNode in pSelf.get_children( ):
					pSelf.m_setChildren.remove( pNode )
				else:
					for pChild in pSelf.get_children( ):
						pChild.remove( pNode )

	def __init__( pSelf ):

		pSelf.__clear( )
		
	def __clear( pSelf ):
		
		pSelf.m_pRoots = CFir.CNode( "" )
		pSelf.m_hashNodes = {}

	def open( pSelf, ifil ):

		pSelf.__clear( )
		for astrLine in csv.reader( ifil, csv.excel_tab ):
			strID, astrChildren = astrLine[0], astrLine[1:]
			if not strID:
				pNode.set( astrChildren[0], astrChildren[1] )
				continue
			pNode = pSelf.m_hashNodes.setdefault( strID, CFir.CNode( strID ) )
			for i in range( 0, len( astrChildren ), 2 ):
				pNode.add_child( pSelf.m_hashNodes.setdefault(
					astrChildren[i], CFir.CNode( astrChildren[i] ) ),
					float(astrChildren[i + 1]) )
		for pNode in pSelf.m_hashNodes.values( ):
			if not pNode.m_pParent:
				pSelf.m_pRoots.add_child( pNode )
				
	def add_root( pSelf, strID ):
		
		pRet = CFir.CNode( strID )
		pSelf.m_pRoots.add_child( pRet )
		return pRet

	def save( pSelf, ofil ):

		pSelf.m_pRoots.save( ofil )

	def get( pSelf, strID ):

		return pSelf.m_hashNodes.get( strID )
	
	def get_roots( pSelf ):
		
		return pSelf.m_pRoots

	def remove( pSelf, pNode ):

		if pNode.m_strID in pSelf.m_hashNodes:
			del( pSelf.m_hashNodes[pNode.m_strID] )
		pSelf.m_pRoots.remove( pNode )

	def collapse( pSelf ):

		pSelf.m_pRoots.collapse( )

	def save_newick( pSelf, ofil ):

		pSelf.m_pRoots.save_newick( ofil )
		ofil.write( ";" )

c_strNamed	= """A1	B1	1	B2	1
	name	a1
B1	C1	1	C2	1
	name	b1
C1
	name	c1
C2
	name	c2
B2	C3	1	C4	1
	name	b2
C3
	name	c3
C4
	name	c4
"""

c_strUnnamed	= """A1	B1	1	B2	1
B1	C1	1	C2	1
B2	C3	1	C4	1
C1
	name	c1
C2
	name	c2
C3
	name	c3
C4
	name	c4
"""
	
class CUTOpenSave(unittest.TestCase):
	
	def setUp( pSelf ):
		
		pSelf.m_pTree = CFir( )
	
	def test_open( pSelf ):
		
		pSelf.m_pTree.open( StringIO.StringIO( c_strNamed ) )
		stioSave = StringIO.StringIO( )
		pSelf.m_pTree.save( stioSave )
		pSelf.assertEqual( stioSave.getvalue( ), c_strNamed )
		
if __name__ == "__main__":
	unittest.main( )
