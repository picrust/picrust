import subprocess
import sys

c_strDirData	= "data/"
c_strDirSrc		= "src/"

def cmd( pE, strProg, strTo, aaArgs = [] ):
	return pipe( pE, None, strProg, strTo, aaArgs )

def d( *astrArgs ):
	return "/".join( astrArgs )

def download( pE, strURL, strT = None, fSSL = False, fGlob = True ):
	if not strT:
		strT = re.sub( '^.*\/', "", strURL )

	def funcDownload( target, source, env, strURL = strURL ):
		strT, astrSs = ts( target, source )
		iRet = ex( " ".join( ("curl", "--ftp-ssl -k" if fSSL else "",
			"" if fGlob else "-g", "-f", "-z", strT, "'" + strURL + "'") ), strT )
# 19 is curl's document-not-found code
		return ( iRet if ( iRet != 19 ) else 0 )
	return pE.Command( strT, None, funcDownload )

def ex( strCmd, strOut = None ):
	sys.stdout.write( "%s" % strCmd )
	sys.stdout.write( ( ( " > %s" % quote( strOut ) ) if strOut else "" ) + "\n" )
	if not strOut:
		return subprocess.call( strCmd, shell = True )
	pProc = subprocess.Popen( strCmd, shell = True, stdout = subprocess.PIPE )
	if not pProc:
		return 1
	strLine = pProc.stdout.readline( )
	if not strLine:
		pProc.communicate( )
		return pProc.wait( )
	with open( strOut, "w" ) as fileOut:
		fileOut.write( strLine )
		for strLine in pProc.stdout:
			fileOut.write( strLine )
	return pProc.wait( )

def pipe( pE, strFrom, strProg, strTo, aaArgs = [] ):
	strFrom, strTo, astrFiles, astrArgs = _pipeargs( strFrom, strTo, aaArgs )
	def funcPipe( target, source, env, strFrom = strFrom, astrArgs = astrArgs ):
		strT, astrSs = ts( target, source )
		return ex( " ".join( [astrSs[0]] + astrArgs +
			( ["<", quote( strFrom )] if strFrom else [] ) ), strT )
	return pE.Command( strTo, [strProg] + ( [strFrom] if strFrom else [] ) +
		astrFiles, funcPipe )

def quote( p ):
	return ( "\"" + str(p) + "\"" )

def scmd( pE, strCmd, strTo, aaArgs = [] ):
	return spipe( pE, None, strCmd, strTo, aaArgs )

def spipe( pE, strFrom, strCmd, strTo, aaArgs = [] ):
	strFrom, strTo, astrFiles, astrArgs = _pipeargs( strFrom, strTo, aaArgs )
	def funcPipe( target, source, env, strCmd = strCmd, strFrom = strFrom, astrArgs = astrArgs ):
		strT, astrSs = ts( target, source )
		return ex( " ".join( [strCmd] + astrArgs + ( ["<", strFrom] if strFrom else [] ) ),
			strT )
	return pE.Command( strTo, ( [strFrom] if strFrom else [] ) + astrFiles, funcPipe )

def ts( afileTargets, afileSources ):
	return (str(afileTargets[0]), [fileCur.get_abspath( ) for fileCur in afileSources])

def _pipeargs( strFrom, strTo, aaArgs ):
	astrFiles = []
	astrArgs = []
	for aArg in aaArgs:
		fFile, strArg = aArg[0], aArg[1]
		if fFile:
			strArg = _pipefile( strArg )
			astrFiles.append( strArg )
		astrArgs.append( quote( strArg ) )
	return ( [_pipefile( s ) for s in (strFrom, strTo)] + [astrFiles, astrArgs] )

def _pipefile( pFile ):
	return ( ( pFile.get_abspath( ) if ( "get_abspath" in dir( pFile ) ) else str(pFile) )
		if pFile else None )
