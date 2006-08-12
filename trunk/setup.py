# -*- coding: latin-1 -*-
#!/usr/bin/env python



from distutils.core import setup

setup(name="feval", \
      version="0.2", \
      description="finite element evaluator", \
      author="Martin P. Lüthi", \
      author_email="luthi@vaw.baug.ethz.ch", \
      url="http://www.sourceforge.net/projects/feval", \
      license="GPL", \
      packages=['feval','feval.fecodes',
                'feval.fecodes.femtool',
                'feval.fecodes.gmv',
                'feval.fecodes.gmsh',
                'feval.fecodes.marc',
                'feval.fecodes.tochnog',
                'feval.fecodes.unv',
                'feval.fecodes.xdr',
                ], 
      package_data={'feval': ['fecodes/*/*.fe']}
      )

