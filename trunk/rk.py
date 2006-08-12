#
# The Runge-Kutta fourth order routine with 
# an adaptive timestep...
#

# first setup constants....

from Numeric import *

r0 = 2.0                 # Start orbit at x=r0, y=0
r = array([r0,0.0])

v0 = pi/2.5              # Start velocity at v_x = 0.0, v_y = pi/2
v = array([0.0,v0])

tau = 0.01               # Start with 50 steps/year

GM = 4.0*pi**2

mass = 1.0
time = 0.0
nstep = 38


#
# now for 'plain' fourth order RK
#

def rk4(x, t, tau, derivsRK, param):
    """
    rk4 takes a single Runge-Kutta timestep..
    """
    
    half_tau = 0.5*tau
    
    F1 = derivsRK( x, t, param )
    t_half = t + half_tau
    xtemp = x + half_tau*F1
    
    F2 = derivsRK( xtemp, t_half, param)
    xtemp = x + half_tau*F2
    
    F3 = derivsRK( xtemp, t_half, param)
    t_full = t + tau
    xtemp = x + tau*F3

    F4 = derivsRK( xtemp, t_full, param)

    xout = x + (tau/6.0)*(F1 + F4 + 2.0*(F2 + F3))
    return xout

#
# our 'f(x)' function to calculate derivitives
#
def gravrk(s, time, GM):
    r = s[0:2]
    rnorm = sqrt(sum(r*r))
    v = s[2:4]
    accel = -GM*r/rnorm**3
    return array([v[0], v[1], accel[0], accel[1]])

#
# Here is the 'adaptor' that dynamically adjusts tau...
#

def rka(x, t, tau, err, derivsRK, param):
    """
    Following Garcia.. p74 rka adjusts the timestep as we go.
    """

    tsave = t
    xsave = x
    safe1 = 0.9
    safe2 = 4.0
    maxtry = 100

    for itry in range(maxtry):
	half_tau = 0.5*tau
	xtemp = rk4(xsave, tsave, half_tau, derivsRK, param)  # do a half_step
	t = tsave + half_tau
	x = rk4(xtemp, t, half_tau, derivsRK, param)  # complete the next half_step
	t = tsave + tau
	xtemp = rk4(xsave, tsave, tau, derivsRK, param)
	scale = 0.5*(abs(x) + abs(xtemp))*err
	xdiff = x - xtemp
	errmax = max(abs(xdiff)/abs(scale + 1e-15))
	tau_old = tau
	tau = safe1*tau_old*errmax**(-0.20)
	tau = max((tau, tau_old/safe2))
	if errmax < 1.0:
	    break
    tau = min((tau, safe2*tau_old))
    return x, t, tau

xplot = []               # Initialize the bookeeping lists
yplot = []
tplot = []
kinetic = []
potential = []
total = []
tauList = []
rList = []

state = array([r[0], r[1], v[0], v[1]])   # init the state 
err = 1.0e-3                              # about 0.1% accuracy
    

for istep in range(nstep):
    rnorm = sqrt(sum(r*r))      # get sqrt(x*x + y*y)
 
    xplot.append(r[0])         # record all those values...
    yplot.append(r[1])
    tplot.append(time)
    kinetic.append(0.5*mass*sum(v*v))
    potential.append(-GM*mass/rnorm)
    total.append(kinetic[-1] + potential[-1])
    tauList.append(log(tau)/log(10.0))
    rList.append(log(sqrt(r[0]**2 + r[1]**2))/log(10.0))

    state, time, tau = rka( state, time, tau, err, gravrk, GM)
    r = state[0:2]
    v = state[2:4]
    
    time = time + tau

import dislin
print "x vs. y"
dislin.plot(xplot,yplot)


#  plot(xplot, yplot, width=600, yrange = (min(xplot), max(xplot)))
print
#  print "Energy: kinetic - red, potential - green, total - blue"
#  plot(tplot, (kinetic, potential, total), width=600, lineColors = (RED, GREEN, BLUE),
#       xLabel = "time", yLabel = "enery")
#  print
#  print "Log(tau (timestep)) vs time"
#  plot(tplot, tauList, width = 600, lineTypes = 1,
#       xLabel="time", yLabel="Log(tau)")
#  print
#  print "Log(r) vs Log(tau)"
#  plot(rList, tauList, width = 600, lineTypes = 1,
#       xLabel="Log(tau)", yLabel="Log(r)")

    
