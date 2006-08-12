# -*- coding: iso-8859-1 -*-

#============================================================
#
#           This file is part of FEval, a module for the
#           evaluation of Finite Element results
#
# Licencse: FEval is provided under the GNU General Public License (GPL)
#
# Authors:  Martin Lüthi, tinu@tnoo.net
#
# Homepage: http://feval.sourceforge.net
#
# History:  long, long, most of it in 2000
#           2001.09.21 (ml):  Code cleaned up for intial release
#
#============================================================

import string, re, os
import numpy as N
#from   UserDict import *
import feval

class TextFile:
    """A general FE input file definition.
    Input files are text files with 'magic' words,
    defining model properties""" 

    type = 'genericFETextFile'

    magicID = 'MAGIC:'
    magicIDlen = len(magicID)
    shapeFunctionDict = {}

    def __init__(self):
        self.FileName = ''
        self.MagicWords = []     # list of Magic words
        self.Comment = None
        self.verbatimData = []

    def getLine(self):
        """Read a line and exclude empty lines (and comments)
        """
        line = self.File.readline()
        if line:
            if string.strip(line) == '':
                return self.getLine()
##              # find comments if allowed in input file
##                      if self.Comment:
##                              line = self.skipComment(line)
##                      line = self.normalizeLine(line)
        return line

    def readFile(self, filename):
        """Parse an input file, call the magic handlers
        """
        self.FileName = filename
        self.File = open(filename, 'r')
        line = self.getLine()
        linelist = [line]
        while 1:
            # check whether the line starts with a letter
            # (this is not yet magic, just a word)
            cond, line = self.isKeyword(line)
            if cond:
                self.findMagic(linelist)
                linelist = [line]
            else:
                linelist.append(line)
            line = self.getLine()
            if not line: break
        self.findMagic(linelist)
        self.File.close()

    def writeFile(self, filename):
        """Parse a input file, call the magic handlers

        """
        self.FileName = filename
        self.File = open(filename, 'w')
        for data in self.verbatimData:
            if data[:self.magicIDlen] == self.magicID:
                magicKey = data[self.magicIDlen:]
##              try:
                lines=eval('self.compose_'+string.lower(magicKey))()
#               except:
#                   print 'no function compose_'+string.lower(magicKey)+' found'

                if lines:
                    if lines[-1][-1] <> '\n':
                        lines[-1] = lines[-1] +'\n'
                    self.File.writelines(lines)
        self.File.close()

    def close(self):
        self.File.close()

    def isKeyword(self, line):
        """check wether the line could contain a magic word
        return None if the line is not interesting
        """
        return (re.match(r'^[0-9]?[a-z]+[0-9]*', string.lower(line[:5]))), line

    def findMagic(self, linelist):
        """check whether a magic word occurs"""

        words = string.split(linelist[0])
        # try to match the magic word xxxx and
        # execute the corresponding handler (extract_xxxx)
        # if this is not possible: keep the input lines in a dictonary
        magicKey = string.lower(words[0])
        if magicKey in self.MagicWords:
            try:
                fct = eval('self.extract_'+string.lower(magicKey))
            except:
                fct = None
            if fct:
                map( fct, [linelist] )

    def setWrite(self, magicKey):
        if magicKey in self.MagicWords:
            self.verbatimData.append(self.magicID+magicKey)


    def skipComment(self, line):
        """This is a very simple comment handling strategy:
        If a comment-string is found then match the respective
        end-of-comment string.
        This does not support nested comments""" 
        for comment in self.Comment:
            self.verbatimData.append(line)
            # comment starts here
            cstart = string.find(line, comment[0]) 
            if cstart >= 0:
                oline = line[:cstart]
                cend = string.find(line, comment[1])
                while cend < 0:
                    line = self.File.readline()
                    cend = string.find(line, comment[1])
                    #  comment ends here
                line = oline + line[cend+1:]
        return line

    def normalizeLine(self, line):
        """Bring the line in a normalized form:
        - replace tabs and commas with whitespace
        """
        return line
        # return re.sub(r'[\t,]', ' ', line)



class FETextFile(TextFile):

    def __init__(self, model=None):
        TextFile.__init__(self)
        self.model = model
        self.FileDescriptor = FileDescriptor(self.type)
        self.MagicWords = self.FileDescriptor.Magic

    def readFile(self, filename):
        TextFile.readFile(self,filename)
        # make a new model cache
        self.model.makeModelCache()


class FileDescriptor(TextFile):

    """Read a file Descriptor
    This kind of a bootstrapping process
    o first the FileDescriptor gets its own MagicWords
    o then the actual file is scanned for MagicWords
    """
    type = 'FEFileDescriptor'

    def __init__(self, type):
        TextFile.__init__(self)
        self.MagicWords = ['translate', 'extract', 'end']
        self.Translate = {}
        self.Extract = []
        # Pairs of characters used to mark comments
        self.Comment = [( '#', '\n' )]

        ### read the description file
        path = os.path.split(feval.__file__)[0]
        filename = os.path.join( path, 'fecodes', type, type + '.fe' )
        self.readFile(filename)

        ### construct the MagicWords for the File
        self.Magic = self.Translate.keys()
        self.Magic.extend( self.Extract )

    def isKeyword(self, line):
        """check wether the line could contain a magic word
        return None if the line is not interesting
        """
        if line.startswith('['):
            return 1, string.split(line[1:], ']')[0]
        else:
            return 0, line

    def extract_translate(self, linelist):
        linelist = linelist[1:]
        for line in linelist:
            if not line.startswith('#'):
                w = string.split(line)
                self.Translate[w[0]] = w[1]

    def extract_extract(self, linelist):
        linelist = linelist[1:]
        for line in linelist:
            if not line.startswith('#'):
                w = string.split(line)
                self.Extract.append(w[0])


if __name__ == '__main__':

    f = FileDescriptor('gmv')
    print f.Magic
    print f.MagicWords
