
"""
@author: SUPRIYA
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
    Fa= 30.0/60.0 #l/sec
    #base flow rate
    Fb= 30.0/60.0 #l/sec
    """
    at any time t=0 as step change of 15 l/min is given to the base flow rate and hence the 
    pH of the system is measured as a function of time
    """
    Ca=2.0 # M or moles/l
    Cb=1.0 # M or moles/l
    """
    concentration/ strength of the two solutions coming in
    """
    def pH0(self):
        Ca=self.Ca
        Cb=self.Cb
        kA=self.kA
        V=self.V
        Fa=self.Fa
        Fb=self.Fb
    
        Xa0=0.0
        Xb0=0.0
        t=0.0
        T=40.0 # sec; this is more than the residence time of the components in reactor in order to get the pH
        n=40.0
        dt=T/(n-1)
        
        while t<=40:
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
        h=pH0.x
        print ("INITIAL STEADY STATE pH = ",h)
            
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
        h=4.75
        T=100.0 # sec; this is more than the residence time of the components in reactor in order to get the pH
        n=100.0
        dt=T/(n-1)
        FB=Fb
        """Now we know the inlet flow rates and volume of the Reactor is specified before
        so we know the residence time of the components with the reactor, so we assume that 
        initially it takes about 20 second for steady state to be reached and after this
        steady state is reached we give EITHER ATEP TYPE OF INPUT OR RAMP INPUT" and again 
        we can see the variation in pH due to thid disturbance."""
        
        
        a=input("STEP CHANGE TO BE GIVEN (IN LITRES/SEC)=")
        
        x, y, z= [], [], []
        
        while t<=100:
            x += [t]
            z += [h]
            y += [FB]
            
            if (t<40):
                 FB=Fb 
            else:
                FB=Fb+a
                
            Xa1=Xa0+dt/V*(Ca*Fa-(Fa+FB)*Xa0)
            
            Xb1=Xb0+dt/V*(Cb*FB-(Fa+FB)*Xb0)
            
            def fun1(x):
                return[Xb1+10**(-x[0])-Xa1/(1+(10**(-x[0])/kA))-10**(x[0]-14)]
            def jac1(x):
                return numpy.array([-scipy.log(10)*(10**(-x[0])-(Xa1/kA)*10**(-x[0])*(1/(1+(10**(-x[0])/kA)))**2-10**(x[0]-14))])
            pH=optimize.root(fun1, [5], jac=jac1, method='hybr')
            h=pH.x
            print("pH=",h)
                
            Xa0=Xa1
            Xb0=Xb1
            t=t+dt
            """this will hence show us how pH varies with time in the reactor"""
                
        
        m=float(pH.x)
        print "pH of the system at the end of the run=",m
        
        fig = plt.figure()  #plot between pH and time after the step input is given
        ax = fig.add_subplot(1, 1, 1)
        ax.plot( x, z, 'blue')
        ax.title.set_text('pH Vs Time')
        ax.xaxis.label.set_text('Time(min)')
        ax.yaxis.label.set_text('pH')
        fig.canvas.draw() 
        
        fig = plt.figure()  #plot between FB and time , shows us the type of disturbance
        ax = fig.add_subplot(1, 1, 1)
        ax.plot( x, y, 'red')
        ax.title.set_text('FB Vs Time')
        ax.xaxis.label.set_text('Time(min)')
        ax.yaxis.label.set_text('Flow rate of Base (l/sec)')
        fig.canvas.draw() 
         
            
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
        T=100.0 # sec; this is more than the residence time of the components in reactor in order to get the pH
        n=100.0
        dt=T/(n-1)
        h=4.75
        FB=Fb
        
        """
        Ramp Input basically increase linearly with time, so say after the system is 
        attaining the steady state, then we give this iput (at t=0) and after the point
        the system will recieve a linearly increasing input of the form 
        a=b*t, where b is the constant which shows us the slope of the ramp function
        """
        b=0.01
        b=input("Enter the value of slope for the ramp function=")
        
        x, y, z= [], [], []
        
        while t<=T:
            x += [t]
            z += [h]
            y += [FB]
            
            if (t<40):
                FB=Fb 
            else:
                FB=Fb+b*t
                
            Xa1=Xa0+dt/V*(Ca*Fa-(Fa+FB)*Xa0)
            
            Xb1=Xb0+dt/V*(Cb*FB-(Fa+FB)*Xb0)
            
            def fun1(x):
                return[Xb1+10**(-x[0])-Xa1/(1+(10**(-x[0])/kA))-10**(x[0]-14)]
            def jac1(x):
                return numpy.array([-scipy.log(10)*(10**(-x[0])-(Xa1/kA)*10**(-x[0])*(1/(1+(10**(-x[0])/kA)))**2-10**(x[0]-14))])
            pH1=optimize.root(fun1, [5], jac=jac1, method='hybr')
            h=pH1.x
            print("h=",h)
        
            Xa0=Xa1
            Xb0=Xb1
            t=t+dt
        
        m=float(pH1.x)
        print "pH of the system at the end of the run=",m
        
        fig = plt.figure()  #plot between pH and time after the step input is given
        ax = fig.add_subplot(1, 1, 1)
        ax.plot( x, z, 'blue')
        ax.title.set_text('pH Vs Time')
        ax.xaxis.label.set_text('Time(min)')
        ax.yaxis.label.set_text('pH')
        fig.canvas.draw()
        
        fig = plt.figure()  #plot between FB and time , shows us the type of disturbance
        ax = fig.add_subplot(1, 1, 1)
        ax.plot( x, y, 'red')
        ax.title.set_text('FB Vs Time')
        ax.xaxis.label.set_text('Time(min)')
        ax.yaxis.label.set_text('Flow rate of Base (l/sec)')
        fig.canvas.draw() 
        
        return pH1
     
# FOR PULSE INPUT     
    def pH2(self):
        Ca=self.Ca
        Cb=self.Cb
        kA=self.kA
        V=self.V
        Fa=self.Fa
        Fb=self.Fb
        
        Xa0=0.0
        Xb0=0.0
        t=0.0
        T=100.0 # sec; this is more than the residence time of the components in reactor in order to get the pH
        n=100.0
        dt=T/(n-1)
        h=4.75
        FB=Fb
        
        """ Pulse input is a deviation/ distubance that occurs for a short amount 
        of time, after which we go back again to our initial value.
        say we give a pulse input at t=40 sec for 5 sec of magnitude "c"
        """
        c=50/60 #lit/sec
        c=input("enter magnitude of the PULSE (l/sec) = ")
        
        x, y, z= [], [], []
        
        while t<=T:
            x += [t]
            z += [h]
            y += [FB]
            
            if (t<40):
                FB=Fb
            if (t>=40 and t<45):
                FB=Fb+c
            else:
                FB=Fb
            
            Xa1=Xa0+dt/V*(Ca*Fa-(Fa+FB)*Xa0)
            
            Xb1=Xb0+dt/V*(Cb*FB-(Fa+FB)*Xb0)
            
            def fun1(x):
                return[Xb1+10**(-x[0])-Xa1/(1+(10**(-x[0])/kA))-10**(x[0]-14)]
            def jac1(x):
                return numpy.array([-scipy.log(10)*(10**(-x[0])-(Xa1/kA)*10**(-x[0])*(1/(1+(10**(-x[0])/kA)))**2-10**(x[0]-14))])
            pH2=optimize.root(fun1, [5], jac=jac1, method='hybr')
            h=pH2.x
            print("h=",h)
        
            Xa0=Xa1
            Xb0=Xb1
            t=t+dt
            
            
        m=float(pH2.x)
        print "pH of the system at the end of the run=",m
        
        fig = plt.figure()  #plot between pH and time after the step input is given
        ax = fig.add_subplot(1, 1, 1)
        ax.plot( x, z, 'blue')
        ax.title.set_text('pH Vs Time')
        ax.xaxis.label.set_text('Time(min)')
        ax.yaxis.label.set_text('pH')
        fig.canvas.draw()
        
        fig = plt.figure()  #plot between FB and time , shows us the type of disturbance
        ax = fig.add_subplot(1, 1, 1)
        ax.plot( x, y, 'red')
        ax.title.set_text('FB Vs Time')
        ax.xaxis.label.set_text('Time(min)')
        ax.yaxis.label.set_text('Flow rate of Base (l/sec)')
        fig.canvas.draw() 
        
        
"""
for obtaining the output, we press F5 to run the python script
in the console we write:
a=acid()
a.pH0() : this gives the pH of the system at steady state
a.pH() : this will give us the pH variation for the STEP input
a.pH1() : this is pH for Ramp input
a.pH2() :  for Pulse Input
"""     
        
    