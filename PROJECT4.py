# -*- coding: utf-8 -*-
"""
Spyder Editor SUPRIYA

"""
from __future__ import division
import scipy
import numpy
from scipy import optimize
import matplotlib.pyplot as plt 
"""
To find the pH of a system for the neutralisation of a weak acid (Acetic acid) and
a strong base (NaOH) by performing the following balances
1. mass balance
2. charge balance for maintaining neutral conditions
3. water equilibrium
4. equilibrium for dissociation of the weak acid
We get,
XB + 10^-pH -XA/(1+(10^-pH/kA)) -10^(pH-14)
SOlving this equation we can get the pH of the system at different instances of time
"""
class acid:
    

    # DATA
    #Equilibrium dissociation constant of the weak acid
    kA=1.75*10**(-5)
    #Volume of the reactor
    V=3.318 #l
    # acid flow rate
    Fa= 15.0/60.0 #l/sec
    #base flow rate
    Fb= 15.0/60.0 #l/sec
    """
    at any time t=0 as step change of 15 l/min is given to the base flow rate and hence the 
    pH of the system is measured as a function of time
    """
    Ca=2.0 # M or moles/l
    Cb=1.0 # M or moles/l
    """
    concentration/ strength of the two solutions coming in
    """
    def pH(self):
        Ca=self.Ca
        Cb=self.Cb
        kA=self.kA
        V=self.V
        Fa=self.Fa
        Fb=self.Fb
    
        Xa0=0.0
        Xb0=0.0
        t=0.0
        T=20.0 # sec; this is more than the residence time of the components in reactor in order to get the pH
        n=20.0
        dt=T/(n-1)
        """Now we know the inlet flow rates and volume of the Reactor is specified before
        so we know the residence time of the components with the reactor, so we assume that 
        initially it takes about 20 second for steady state to be reached and after this
        steady state is reached we give EITHER ATEP TYPE OF INPUT OR RAMP INPUT" and again 
        we can see the variation in pH due to thid disturbance."""
        
        while t<=20:
            Xa1=Xa0+dt/V*(Ca*Fa-(Fa+Fb)*Xa0)
            
            Xb1=Xb0+dt/V*(Cb*Fb-(Fa+Fb)*Xb0)
            
            def fun(x):
                return[Xb0+10**(-x[0])-Xa0/(1+(10**(-x[0])/kA))-10**(x[0]-14)]
        
            def jac(x):
                return numpy.array([-scipy.log(10)*(10**(-x[0])-(Xa0/kA)*10**(-x[0])*(1/(1+(10**(-x[0])/kA)))**2-10**(x[0]-14))])
            
            pH0=optimize.root(fun, [5], jac=jac, method='hybr')
            
            Xa0=Xa1
            Xb0=Xb1
            t=t+dt
        #therefore now we know the pH of the initial system at steady state. this is 
        #done by solving the non - polynomial equation using powells method to arrive
            # at the pH i.e. roots of the equation.
        pH=pH0
        Xa0=Xa1
        Xb0=Xb1
        
        T1=50.0 #sec; time till which we monitor the system
        # now the step input is given
        a=15.0/60.0 #l/sec @ t=0
        FB=Fb+a
        t=0
        h=pH
        x, y= [], []
        
        while t<= (T1-T):
            x += [t]
            y += [h]
            
            Xa1=Xa0+dt/V*(Ca*Fa-(Fa+FB)*Xa0)
            
            Xb1=Xb0+dt/V*(Cb*FB-(Fa+FB)*Xb0)
            
            def fun1(x):
                return[Xb1+10**(-x[0])-Xa1/(1+(10**(-x[0])/kA))-10**(x[0]-14)]
            def jac1(x):
                return numpy.array([-scipy.log(10)*(10**(-x[0])-(Xa1/kA)*10**(-x[0])*(1/(1+(10**(-x[0])/kA)))**2-10**(x[0]-14))])
            pH=optimize.root(fun1, [5], jac=jac1, method='hybr')
            h=pH
            print("h=",h)
        
            Xa0=Xa1
            Xb0=Xb1
            t=t+dt
            """this will hence show us how pH varies with time in the reactor"""
        """   
        fig = plt.figure()  #plot between pH and time after the step input is given
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(x, y, 'ro')
        ax.title.set_text('pH Vs Time')
        ax.xaxis.label.set_text('Time(min)')
        ax.yaxis.label.set_text('pH')
        fig.canvas.draw()   
        """ 
            
        return pH
#for RAMP INPUT        
    def pH1(self):
        Ca=self.Ca
        Cb=self.Cb
        kA=self.kA
        V=self.V
        Fa=self.Fa
        Fb=self.Fb
    
        Xa0=0.0
        Xb0=0.0
        t=0.0
        T=20.0 # sec; this is more than the residence time of the components in reactor in order to get the pH
        n=20.0
        dt=T/(n-1)
        while t<=20:
            Xa1=Xa0+dt/V*(Ca*Fa-(Fa+Fb)*Xa0)
            
            Xb1=Xb0+dt/V*(Cb*Fb-(Fa+Fb)*Xb0)
            
            def fun(x):
                return[Xb0+10**(-x[0])-Xa0/(1+(10**(-x[0])/kA))-10**(x[0]-14)]
        
            def jac(x):
                return numpy.array([-scipy.log(10)*(10**(-x[0])-(Xa0/kA)*10**(-x[0])*(1/(1+(10**(-x[0])/kA)))**2-10**(x[0]-14))])
            
            pH10=optimize.root(fun, [5], jac=jac, method='hybr')
            # Finds the roots of a multivariate function using MINPACK’s hybrd and hybrj routines (modified Powell method).
            
            Xa0=Xa1
            Xb0=Xb1
            t=t+dt
            
        pH1=pH10
        Xa0=Xa1
        Xb0=Xb1
        
        T1=50.0 #sec; time till which we monitor the system
        # now the ramp input is given
        #a=0.1*t #l/sec @ t=0
        
        t=0
        h=pH1
        x, y= [], []
        
        while t<= (T1-T):
            x += [t]
            y += [h]
            
            FB=Fb+0.001*t
            
            Xa1=Xa0+dt/V*(Ca*Fa-(Fa+FB)*Xa0)
            
            Xb1=Xb0+dt/V*(Cb*FB-(Fa+FB)*Xb0)
            
            def fun1(x):
                return[Xb1+10**(-x[0])-Xa1/(1+(10**(-x[0])/kA))-10**(x[0]-14)]
            def jac1(x):
                return numpy.array([-scipy.log(10)*(10**(-x[0])-(Xa1/kA)*10**(-x[0])*(1/(1+(10**(-x[0])/kA)))**2-10**(x[0]-14))])
            pH1=optimize.root(fun1, [5], jac=jac1, method='hybr')
            h=pH1
            print("h=",h)
        
            Xa0=Xa1
            Xb0=Xb1
            t=t+dt
            
        return pH1

""" for obtaining the output, we press F5 to run the python script
in the console we write:
a=acid()
a.pH() : this will give us the pH variation for the STEP input
a.pH1() : this is pH for Ramp input
"""        