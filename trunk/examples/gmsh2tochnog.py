#!/usr/bin/env python

####
#
#  Python script for the generation of a Tochnog input file from a GMSH file
#
####

## import modules
import os
from feval.ModelData import *
from feval.fecodes.gmsh.GMSHFile import *
from feval.fecodes.tochnog.TochnogFile import *

## file names
infilename  = os.path.join(feval.__path__[0],'data','gmsh','simple.gmsh')
outfilename = os.path.join(feval.__path__[0],'data','gmsh','simple_tochnog.dat')

## initalize a FE-Model (data container with "intelligence")
m = ModelData()

## read the GMSH input file
mf = GMSHFile(m)
mf.readFile(infilename)

## initailise the Tochnog File
ff = TochnogFile(m)

## set the blocks which are written to the file 
ff.setWrite('node')
ff.setWrite('element')

## write the file
ff.writeFile(outfilename)

