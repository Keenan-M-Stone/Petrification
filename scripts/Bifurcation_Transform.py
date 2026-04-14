# **********************************************************************************************************************
# Title:            Bifurcation_Transform.py
# Author:           Keenan M Stone
# Contact:          Lemma137@gmail.com
#
# Date Created:     5/12/17
# Date Modified:    5/12/17
#
# Description:      Strongly based on 'Iteration_Transform.py', a python script  
#                   for function testing, intended for visualizing fixed point 
#                   analysis. Here we explore as well the effects of a type
#                   transform on self similar functions, not for visualizing
#                   itteration, but instead bifurcation. This one just makes
#                   one figure, but a second version which generates a .gif will
#                   be made later (may end up having less functional use though).
# Other:
#                   "MFLAG:" for programer errors and notable issues
# **********************************************************************************************************************
# Libraries
# **********************************************************************************************************************
from pylab import *
import sympy as sym
import random as rand
#from matplotlib.pyplot import quiver
#MFLAG: may not work in some versions of python. Alternatives:  (numpy, scipy, matplotlib)

# **********************************************************************************************************************
#Functions
# **********************************************************************************************************************

# Logistic map or other function, I've left some options in comment
def lmap(a,x):
    return a*x*(1.0-x)                      # Traditional Logistic map                         
    #return ((a*x*x/((1+(x*x))**2) ) )      # Migdal-Katinov Renormalization
	
# Transformation generator 
def trans(alpha, a, x):
    f = lmap(a,x)
    return (alpha*f) + (1-alpha)*x

# Symbolic evaluation. I'm new to this, hopefully it works out well.
def getsymbols():
    a = sym.symbols('a')
    x = sym.symbols('x') 
    alpha = sym.symbols('alpha')
    
    # Functions, symbolized
    ff = lmap(a,x)
    gg = trans(alpha,a,x)
    
    return ff,gg
    
    
# Global Alpha:  
#               This function findes the maximum slope of lmap and returns
#               the value of best value for alpha. Try to keep up, the double
#               letters aa is an actual numeric value.
# MFLAG:        Periodic functions may lead to errors.
def globify(aa,xx):
    
    # Initialization of variables
    df_max = 0.0
    df_min = 0.0
    error_flag = False

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
        if ((df_max < temp) and (temp != inf)):
            df_max = temp
        if ((df_min > temp) and (temp != -inf)):
            df_min = temp
            
    #-------------------------------------------------------------------------        

    if(df_max > 1.0):
        print('Taking max: ' + str(df_max) )
        alpha = 2.0/(1.0-df_max)      # Value of alpha  
    elif(df_min < -1.0):
        print('Taking min: ' + str(df_min) )
        alpha = 2.0/(1-df_min)      # Value of alpha
        print('WARNING:    Conversion of stable fixed points may fail!')
    else:
        alpha = -1                  # Code runs but error needs to print
        print('ERROR:      Conversion failure! No unstable fixed points!')
        print('max:   ' + str(df_max) +' ,   min:   ' + str(df_min))
        error_flag = True
        
    
    if (error_flag):    
        beta,error_flag = locify(aa,xx)
        if(error_flag):
            print('Error could not be resolved. Continuing with alpha = -1')
        else:
            print('Error resolved locally. Continuing with alpha = ' + str(beta) )
            alpha = beta
        
         
    return alpha,error_flag 
    
# Get alpha based on local maximum in slope    
def locify(aa,xx):
    # Initialization of variables
    df_max = 0.0
    df_min = 0.0
    alpha = -1.0
    error_flag = False

    # Symbolic derivatives:    
    f,g = getsymbols()            # get sybolic function, not using gg
    df = sym.diff(f,'x')          # first derivative
    df = sym.simplify(df)

    dF = df.subs('a',aa)
    dF = sym.simplify(dF)
    
    # Loop over range of interest.
    #-------------------------------------------------------------------------           
    for ii in xx:
        dF_temp = dF.subs('x',ii) 
        if (df_max < dF_temp):
            df_max = dF_temp
        if (df_min < dF_temp):
            df_min = dF_temp
    #-------------------------------------------------------------------------        

    if(df_max > 1.0):
        print('Taking max: ' + str(df_max) )
        alpha = 2.0/(1.0-df_max)      # Value of alpha  
    elif(df_min < -1.0):
        print('Taking min: ' + str(df_min) )
        alpha = 2.0/(1-df_min)      # Value of alpha
        print('WARNING:    Local conversion of stable fixed points may fail!')
    else:
        alpha = -1                  # Code runs but error needs to print
        print('ERROR:      Local conversion failure! No unstable fixed points!')
        print('max:   ' + str(df_max) +' ,   min:   ' + str(df_min))
        error_flag = True         
    
    return alpha, error_flag
    
        


# Iterative Composition    
def fpa(f,alpha,a,xi,tol,xmax):
    # Range condition    
    go = True
    # Get first map value
    if(alpha == 0.0):
        x = f(a,xi)
    else:
        x = f(alpha,a,xi)
    # Initialize arrays    
    Ix = [] 
    
    # Loop variables
    maxstep = tol
    ii = 0
    #dx = tol + 1000

    # Loop
    while((maxstep > ii) and (go)):
    #while((dx>tol) and (maxstep > ii)):  # MFLAG: Alternate between looping conditions as needed; uncomment dx = ...
        xi = x
        if(alpha == 0.0):
            x = f(a,x)
        else:
            x = f(alpha,a,x)
        
        #dx = abs(x - xi)
        ii = ii+1
        
        # Fill array
        if (ii > 100):
            Ix = Ix + [x]
        
        # Need to stay within range of x
        if ((xi > xmax) or (x < (-1*xmax))):        
            go = False
    #print(xi)
    return Ix

# **********************************************************************************************************************        
# Driver code
# **********************************************************************************************************************

# Setup 
maxstep = 500                  # Number of itteration steps 
xmax = 5.0
lx = linspace(-1,xmax,9000)    # Setting plot domain
lxt = lx

x_set = rand.sample(lx,100)  # Initial values
x_set = sorted(set(x_set))
a_set = linspace(2.8,3.8,500)    # Population Amplitude (Fertility)

#alpha,e_flag = globify(a,lx)      # Transformation constant, calculated from second derivative of lmap.
alpha = -1.0

fig = figure(figsize=(15.0, 9.0))
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)

# Outer Loop
for a in a_set:
    y = []                                  # Initialize return results of lmap
    yt = []
    for x0 in x_set:
        # Compute data for range and set domain
        x  = fpa(lmap,0.0,a,x0,maxstep,xmax)       # Used to draw the arrows in quiver
        xt = fpa(trans,alpha,a,x0,maxstep,xmax)   # Same thing for transformed equation
        x = sorted(set(x))
        xt = sorted(set(xt))
        y = y + x
        yt = yt + xt
    y = asarray(y,'float64')
    yt = asarray(yt,'float64')
    y = sorted(set(y))
    yt = sorted(set(yt))
    a_domain_1 = ones(len(y)) * a
    a_domain_2 = ones(len(yt)) * a
    ax1.plot(a_domain_1, y, 'ko')
    ax2.plot(a_domain_2, yt, 'ko')
        
     
# **********************************************************************************************************************
# Plot tailoring
# **********************************************************************************************************************

# Set plot title based on computation results

# Logistic Map Title
ti  =   r' $Original$  $Function$:     $ f = ax(1-x) $  '
# Migdal - Katinoff Title
#ti          =   r' $Original$  $Function$:     $ f = ax^2/((1+x^2)^2)  $  ' + '\n' + r'$  Parameters:$   $ a = %s $ ,   $f^{(0)}(x) = %s $,   $x^{*} = %s $' %(a,x0, x[-1])
#Transformed Function Title
ti_trans  =   r' $Transformed$  $Function$:  $ g = \alpha f + (1-\alpha)x $ ' + '\n' + r'$  Parameters:$   $ \alpha = %s $ '     %(alpha)


ax1.set_xlabel(r'Fertility (a)',fontsize = 16) #xlabel('Domain')
ax1.set_ylabel(r'Population', fontsize = 16) #ylabel('Range')
ax1.set_title(ti,fontsize=18)


ax2.set_xlabel(r'Fertility (a)',fontsize = 16) #xlabel('Domain')
ax2.set_ylabel(r'Population', fontsize = 16) #ylabel('Range')
ax2.set_title(ti_trans,fontsize=18)

show()
# **********************************************************************************************************************
# End of script
