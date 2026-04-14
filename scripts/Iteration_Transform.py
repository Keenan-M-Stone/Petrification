# **********************************************************************************************************************
# Title:            Iteration_Transform.py
# Author:           Keenan M Stone
# Contact:          Lemma137@gmail.com
#
# Date Created:     5/12/17
# Date Modified:    5/12/17
#
# Description:      Strongly based on 'Iteration.py', a python script for 
#                   function testing, intended for visualizing fixed point 
#                   analysis. Here we explore as well the effects of a type
#                   transform on self similar functions.
# Other:
#                   "MFLAG:" for programer errors and notable issues
# **********************************************************************************************************************
# Libraries
# **********************************************************************************************************************
from pylab import *
import sympy as sym
#from matplotlib.pyplot import quiver
#MFLAG: may not work in some versions of python. Alternatives:  (numpy, scipy, matplotlib)

# **********************************************************************************************************************
#Functions
# **********************************************************************************************************************

# Logistic map or other function, I've left some options in comment
def lmap(a,x):
    #f= a*x*(1-x)                    # Traditional Logistic map
    #return f                            
    #return a*(1-exp(-x))
    #return (2.0/a)*log(x) - (1.0/x)
    #return cosh(x)/sinh(x)                         # Caternary
    #return 0.5 * (1.0 + exp(-2.0*x))
    return ((a*x*x/((1+(x*x))**2) ) )               # Migdal-Katinov Renormalization
	
# Transformation generator 
def trans(alpha, a, x):
    f = lmap(a,x)
    return (alpha*f) + (1-alpha)*x

# Symbolic evaluation. I'm new to this, hopefully it works out well.
def getsymbols():
    a = sym.symbols('a')
    x = sym.symbols('x') 
    alpha = sym.symbols('alpha')
    #f = sym.symbols('f')
    #g = sym.symbols('g')
    
    # Functions, symbolized
    ff = lmap(a,x)
    gg = trans(alpha,a,x)
    
    return ff,gg
    
    
# Global Alpha:  
#       This function findes the maximum slope of lmap and returns
#       the value of best value for alpha. Try to keep up, the double
#       letters aa is an actual numeric value.
def globify(aa):
    
    # Initialization of variables
    df_max = 0.0
    df_min = 0.0

    # Symbolic derivatives:    
    f,g = getsymbols()          # get sybolic function, not using gg
    df = sym.diff(f,'x')          # first derivative
    df = sym.simplify(df)
    ddf = sym.diff(df,'x')        # second derivative
    ddf = sym.simplify(ddf)
    
    ddF = ddf.subs('a',aa)        # replace occurances of symbol 'a' with actual value
    ddF = sym.simplify(ddF)
    dF = df.subs('a',aa)
    dF = sym.simplify(dF)
    
    Tau_x = sym.solve(ddF,'x')    # solve for roots of second derivative
    
    Tau_f = []                  # empty array of slopes
    # Loop through roots
    #-------------------------------------------------------------------------
    for t in Tau_x:
        temp = dF.subs('x',t)
        Tau_f.append(temp)
        
        # MFLAG:    In theory, the one thing that guarutees a succesful transform
        #           is that alpha < 0 -> df_max > 1. We can check this here.
        if (df_max < temp):# and (temp != inf)):
            df_max = temp
        if (df_min > temp):# and (temp != -inf)):
            df_min = temp
            
    #-------------------------------------------------------------------------        

    if(df_max > 1.0):
        print('Taking max: ' + str(df_max) )
        alpha = 2.0/(1.0-df_max)      # Value of alpha  
    elif(df_min < -1.0):
        print('Taking min: ' + str(df_min) )
        alpha = 2.0/(1-df_min)      # Value of alpha
    else:
        alpha = -1                  # Code runs but error needs to print
        print('Possible Errors!!!!!!!')
         
    return alpha 
    
    
        


# Iterative Composition    
def fpa(f,alpha,a,xi,tol):
    # Get first map value
    if(alpha == 0.0):
        x = f(a,xi)
    else:
        x = f(alpha,a,xi)
    # Initialize arrays    
    Ix = [xi,xi,x] 
    Iy = [0,x,x]
    
    # Loop variables
    maxstep = tol
    ii = 0
    #dx = tol + 1000

    # Loop
    while((maxstep > ii)):
    #while((dx>tol) and (maxstep > ii)):  # MFLAG: Alternate between looping conditions as needed; uncomment dx = ...
        xi = x
        if(alpha == 0.0):
            x = f(a,x)
        else:
            x = f(alpha,a,x)
        
        Iy = Iy+[xi] + [x]
        Ix = Ix +[xi] + [xi]
        
        #dx = abs(x - xi)
        ii = ii+1
    print(xi)
    return Ix,Iy

# **********************************************************************************************************************        
# Driver code
# **********************************************************************************************************************

# Setup 
a = 4.0                 # Population Amplitude (Fertility)
x0 = .7                 # Initial Value
maxstep = 20            # Number of itteration steps 

alpha = globify(a)      # Transformation constant, calculated from second derivative of lmap.

# Compute data for range and set domain
x,y = fpa(lmap,0.0,a,x0,maxstep)       # Used to draw the arrows in quiver
xt,yt = fpa(trans,alpha,a,x0,maxstep)  # Same thing for transformed equation

lx = linspace(-1,5,9000)            # Setting plot domain
lxt = lx

ly = lmap(a,lx)                     # Actual Function
lyt = trans(alpha,a,lx)
#ny = lmap(a,ly)                    # I forget what this was for

x = asarray(x, 'float64')                      # Just converts the data into an array. The last value is the fixed point evaluated at tolerence.
y = asarray(y, 'float64')   

xt = asarray(xt, 'float64')
yt = asarray(yt, 'float64')      
     
# **********************************************************************************************************************
# Plot tailoring
# **********************************************************************************************************************

# Set plot title based on computation results
#ti = 'Problem 8.1a: x_0 = %s, m_i = %s' %(a,x0)
#ti = ' $Fixed Point Analysis:$  \t $ y_{0} = a(1-e^{-y_{0}})$  \n $ Parameters:$ \t $ a = %s $ , \t $f^{0}(y) = %s $, \t $y_{0} = %s $' %(a,x0, x[-1])
#ti = ' $Fixed Point Analysis:$  \t $ x = -(\{ax^2/(1+x^2)^2 \} - x) + x  $  \n $ Parameters:$ \t $ a = %s $ , \t $f^{0}(x) = %s $, \t $x_{0} = %s $' %(a,x0, x[-1])
#ti = ' $Fixed Point Analysis:$  \t $ x = ax^2/((1+x^2)^2)  $  \n $ Parameters:$ \t $ a = %s $ , \t $f^{0}(x) = %s $, \t $x_{0} = %s $' %(a,x0, x[-1])
#ti = ' $Fixed Point Analysis:$  \t $ x = ax(1-x)  $  \n $ Parameters:$ \t $ a = %s $ , \t $f^{0}(x) = %s $, \t $x_{0} = %s $' %(a,x0, x[-1])
ti          =   r' $Original$  $Function$:     $ f = -ax(1-x) + 2x $  ' + '\n' + r'$  Parameters:$   $ a = %s $ ,   $f^{(0)}(x) = %s $,   $x^{*} = %s $' %(a,x0, x[-1])
ti_trans    =   r' $Transformed$  $Function$:  $ g = \alpha f + (1-\alpha)x $ ' + '\n' + r'$  Parameters:$   $ \alpha = %s $ ,   $ g^{(0)}(x) = %s $,   $ x^{*} = %s $'     %(alpha,x0, xt[-1])


# Make Plot and aesthetics

# Original
#-------------------------------------------------------------------------------------------------
figure()
lines = plot(lx,ly,lx,lx,lx,zeros(len(lx)),'r--',zeros(len(lx)),lx,'r--')
setp(lines,linewidth=2.0)
xlabel('$f^{(i-1)}(x)$',fontsize = 16) #xlabel('Domain')
ylabel('$f^{(i)}(x)$', fontsize = 16) #ylabel('Range')

# Scaling and finalizing plot
ylim([-0.2,1.2])
xlim([-0.2,1.2])
title(ti,fontsize=18)

# Draw arrows on plot to show iteration
quiver(x[:-1], y[:-1], x[1:]-x[:-1], y[1:]-y[:-1], scale_units='xy', angles='xy', scale=1.0) 

#-------------------------------------------------------------------------------------------------
# Transformed
#-------------------------------------------------------------------------------------------------
figure()
lines = plot(lxt,lyt,lxt,lxt,lxt,zeros(len(lxt)),'r--',zeros(len(lxt)),lxt,'r--')
setp(lines,linewidth=2.0)
xlabel('$g^{(i-1)}(x)$',fontsize = 16) #xlabel('Domain')
ylabel('$g^{(i)}(x)$', fontsize = 16) #ylabel('Range')

# Scaling and finalizing plot
ylim([-0.2,1.2])
xlim([-0.2,1.2])
title(ti_trans,fontsize=18)

# Draw arrows on plot to show iteration
quiver(xt[:-1], yt[:-1], xt[1:]-xt[:-1], yt[1:]-yt[:-1], scale_units='xy', angles='xy', scale=1.0) 

#-------------------------------------------------------------------------------------------------

show()
# **********************************************************************************************************************
# End of script
