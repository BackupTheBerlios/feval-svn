#!/usr/bin/env python
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
# Purpose:  Calculate a flowline with a simple Euler method
#
#============================================================

import feval.FEval, feval.MarcT16File
import numpy as N
import copy

__version__ = "0.1"


class Flowline:
    """A class for a flowline integrator"""

    def __init__(self, veloFunc, startPoint, timestep):
        self.veloFunc = veloFunc    # the function for calculating the velocity
        self.Points = [startPoint]  # collecting the points along the flowline
        self.timestep = timestep

    def calcFlowline(self):
        pass


class Euler(Flowline):

    def calcFlowline(self):
        p = self.Points[0]
        velo = self.veloFunc( p )
        while velo:
            p = p + velo*self.timestep
            self.Points.append(p)
            velo = self.veloFunc( p )
           

#  def main():
#      """
#      main()
#      Module mainline (for standalone execution)
#      """
#      m = fe.FEval.FEModel()
#      mf = fe.MarcT16File.MarcT16File(m, '/home/luthi/marc/colle/ca8.t16')




if __name__ == "__main__":

#      inc1, inc2 = -2, -1
    inc1, inc2 = 1, 2

#   Offset von Tinu's Colle Modell
    point = N.array([633846., 86617.,  4339.4027+2.5]) - \
            N.array([630000., 80000.0, 0.])

    marcfile = '/cfd/marc/colle/cgx8_dc0.80/cgx8dec.t16'
    #      marcfile = '/home/luthi/marc/car8.t16'
#      marcfile = '/home/luthi/marc/colle/cgx8/cgx8_gm1_z60_b0.1_n3.0_d0.85_tr_Dc/cgx8dec.t16'
    m  = fe.FEval.FEModel()

    print 'reading Model'
    mf = fe.MarcT16File.MarcT16File(m, marcfile)

    info = mf.IncrementInfo()
    t1 = info[inc1][1]
    t2 = info[inc2][1]
    dt = t2-t1
    print 'time', t2, t1, dt

    print 'reading increment %i' % (info[inc2][0])
    mf.readInc(inc2)
    nodvar = copy.copy(m.NodVar)
    print 'reading increment %i' % (info[inc1][0])
    mf.readInc(inc1)


    print 'calculating velocity'
    for k in m.NodVar.keys():
        velo = (nodvar[k][0:3]-m.NodVar[k][0:3])/dt
        #          velo = (nodvar[k][0:2]-m.NodVar[k][0:2])/dt
        m.NodVar[k] = velo
    m.NodVarInfo = ['v_x', 'v_y', 'v_z']
    
## falls man die anderen Grössen noch braucht....
    #      for k in m.NodVar.keys():
    #          velo = (nodvar[k][0:3]-m.NodVar[k][0:3])/dt
    #          m.NodVar[k] = N.concatenate( (m.NodVar[k], velo) )
    #      m.NodVarInfo.append('v_x')
    #      m.NodVarInfo.append('v_y')
    #      m.NodVarInfo.append('v_z')
    
    #      print m.getNodVar(m.Coord[20],('v_x','v_y','v_z'))
    
    def calcVelo( point ):
        return m.getNodVar(point, ('v_x','v_y','v_z'))
#          return m.getNodVar(point, ('v_x','v_y'))

import time
t1 = time.time()
print 'calculating flowline'
fl = Euler(calcVelo, point, -2.0)

import profile
profile.run('fl.calcFlowline()')

t2 = time.time()
print t2-t1

#pp = pp = N.array(fl.Points) 

#import dislin
#  dislin.disini()
#  dislin.graf(0,500,0,100, 4200, 4600, 4200,100)

#  ## plot mesh
#  for k in m.Conn.keys():
#      e = m.ElementFactory(k)
#      co = e.nodcoord
#      co = N.concatenate( (co, co[0,N.NewAxis]) )
#      for i in range(3):
#          dislin.curve( co[i:i+2,0], co[i:i+2, 1], 2)

#  dislin.curve(pp[:,0], pp[:,1], len(pp) )

#  dislin.disfin()


