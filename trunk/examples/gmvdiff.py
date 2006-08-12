#!/usr/bin/env python

####
#
#  Calculate the difference between two GMV files and output them to
#  a new GMV file
#
####

# import modules
import sys
import feval
import feval.fecodes.gmv.GMVFile as gmv

# get the file names
if len(sys.argv) != 4:
    print 'usage infile1 infile2 outfile'
    sys.exit(1)

infile1, infile2, outfile = sys.argv[1:]

# initalize a FE-Model (data container with "intelligence")
m1 = feval.ModelData.ModelData()
m2 = feval.ModelData.ModelData()

# read the input files
mf1 = gmv.GMVFile(m1)
mf1.readFile(infile1)
mf2 = gmv.GMVFile(m2)
mf2.readFile(infile2)

# some simple checks whether we compare similar files
if ((len(m1.Coord) != len(m2.Coord) ) or
    (len(m1.Conn)  != len(m2.Conn)  ) or
    (m1.getNodVarInfo() != m2.getNodVarInfo())):
    print 'files different'
    sys.exit(1)

for node in m1.getNodeNames():
    m1.setNodVar(node,
                 m1.getNodVar(node) - m2.getNodVar(node))
    
mf1.setWrite('gmvinput')
mf1.setWrite('nodes')
mf1.setWrite('cells')
mf1.setWrite('variable')
mf1.setWrite('endgmv')

mf1.writeFile(outfile)

