#!/usr/bin/env python

import re
import sys

for strLine in sys.stdin:
	mtch = re.search( r'prokMSA_id:(\d+)\s+\S+\s+(.+)', strLine )
	if mtch:
		print( "\t".join( mtch.groups( ) ) )
