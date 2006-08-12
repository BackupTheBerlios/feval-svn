#!/usr/bin/env python

####
#
#  Python script for the generation of a Tochnog input file from a Marc file
#
####

## import modules
import os
from feval.ModelData import *
from feval.fecodes.marc.MarcFile import *
from feval.fecodes.tochnog.TochnogFile import *

## file names
infilename  = os.path.join( feval.__path__[0], 'data', 'marc', 'test1.dat' )
outfilename = os.path.join( feval.__path__[0], 'data', 'marc', 'test1_tochnog.dat' )

## initalize a FE-Model (data container with "intelligence")
m = ModelData()

## read the Marc input file
mf = MarcFile(m)
mf.readFile(infilename)

## initailise the Femtool File
tn = TochnogFile(m)

## set the blocks which are written to the file 
#tn.setWrite('number_of_space_dimensions')
#tn.setWrite('echo')
tn.setWrite('element')
tn.setWrite('node')

## write the file
tn.writeFile(outfilename)

