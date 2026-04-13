# **********************************************************************************************************************
# Title:            Bifurcation_logistic_gif.py
# Author:           Keenan M Stone
# Contact:          Lemma137@gmail.com
#
# Date Created:     5/14/17
# Date Modified:    5/15/17
#
# Description:      Strongly based on 'Bifurcation_Transform.py'. 
#                   this script generates a .gif file for viewing output
#                   of a transform for a range of alpha values.
#
#                   Considering f'(x) = a - 2ax , a > 1 , 0 < x < 1
#                   we can expect alpha in ( 2/(1-a), 2/(1+a)). The largest
#                   range of alpha corresponds to the smallest value of a, 
#                   but here we might let x extend beyond the upper bound 1-
#                   still the range of alpha is valid. Notably, if the value at
#                   a used doesn't situate the current value of alpha one might
#                   want to make not of that point. I'm going to color them red
#                   and hope the distinction is clear.
#                   
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

from PIL import Image
from images2gif import writeGif             # Had to edit the module to make it work. 
                                            # Use >> import images2gif >> images2gif.__file__ 
                                            # Follow to the .py file and if you're me, search "MFLAG"
                                            # if you're not me, you'll see the error. 


#import matplotlib.animation as animation
#from joblib import Parallel, delayed


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


'''
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
    
'''


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
    
    
# -Convert Figure to Data/Image-
# Taken from " http://www.icare.univ-lille1.fr/node/1141 "
# CAUTION: Webutation -> 70%
# Supposedly these exist as modules somewhere, but I couldn't get that to work,
# so I just coppied them directly from the sourcecode. Not my handiwork.
# --------------------------------------------------------------------------------------------

def fig2data ( fig ):
    """
    @brief Convert a Matplotlib figure to a 4D numpy array with RGBA channels and return it
    @param fig a matplotlib figure
    @return a numpy 3D array of RGBA values
    """
    # draw the renderer
    fig.canvas.draw ( )
 
    # Get the RGBA buffer from the figure
    w,h = fig.canvas.get_width_height()
    buf = np.fromstring ( fig.canvas.tostring_argb(), dtype=np.uint8 )
    buf.shape = ( w, h,4 )
 
    # canvas.tostring_argb give pixmap in ARGB mode. Roll the ALPHA channel to have it in RGBA mode
    buf = np.roll ( buf, 3, axis = 2 )
    return buf


def fig2img ( fig ):
    """
    @brief Convert a Matplotlib figure to a PIL Image in RGBA format and return it
    @param fig a matplotlib figure
    @return a Python Imaging Library ( PIL ) image
    """
    # put the figure pixmap into a numpy array
    buf = fig2data ( fig )
    w, h, d = buf.shape
    return Image.fromstring( "RGBA", ( w ,h ), buf.tostring( ) )

# --------------------------------------------------------------------------------------------



# **********************************************************************************************************************        
# Driver code
# **********************************************************************************************************************

# Setup 
maxstep = 250                  # Number of itteration steps 
xmax = 5.0
lx = linspace(-1,xmax,9000)    # Setting plot domain
lxt = lx

x_set = rand.sample(lx,50)  # Initial values
x_set = sorted(set(x_set))
a_set = linspace(2.8,3.8,250)    # Population Amplitude (Fertility)
alpha_set = linspace(2/(1-a_set[0]), 2/(1+a_set[0]), 50) # Alpha's to consider as calculated from the derivative.

#alpha,e_flag = globify(a,lx)      # Transformation constant, calculated from second derivative of lmap.
#alpha = -1.0

fig = figure(figsize=(15.0, 9.0))
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)


def AniMe(alpha,alpha_start, a_set,x_set,maxstep,xmax,ax1,ax2):
    print('alpha = ' + str(alpha))
    ax2.clear()
    
    # Outer Loop
    for a in a_set:
        print('a = ' + str(a))
        y = []                                  # Initialize return results of lmap
        yt = []
        for x0 in x_set:
            # Compute data for range and set domain
            if (alpha == alpha_start):
                #print('this parts good')
                x  = fpa(lmap,0.0,a,x0,maxstep,xmax)       # Used to draw the arrows in quiver
                x = sorted(set(x))
                y = y + x
                y = asarray(y,'float64')
                y = sorted(set(y))
                a_domain_1 = ones(len(y)) * a
                ax1.plot(a_domain_1, y, 'ko')
            
            xt = fpa(trans,alpha,a,x0,maxstep,xmax)   # Same thing for transformed equation
            xt = sorted(set(xt))            
            yt = yt + xt
            yt = asarray(yt,'float64')            
            yt = sorted(set(yt)) 
            a_domain_2 = ones(len(yt)) * a
            
            
            # Change color if the range is getting questionable.
            if(alpha > 2.0/(1-a)) and ((alpha < 2.0/(1+a))):
                ax2.plot(a_domain_2, yt, 'ko')
            else:
                ax2.plot(a_domain_2, yt, 'ro')
    # **********************************************************************************************************************
    # Plot tailoring
    # **********************************************************************************************************************

    # Set plot title based on computation results

    if (alpha == alpha_start):
        # Logistic Map Title
        ti  =   r' $Original$  $Function$:     $ f = ax(1-x) $  '
        # Migdal - Katinoff Title
        #ti          =   r' $Original$  $Function$:     $ f = ax^2/((1+x^2)^2)  $  ' + '\n' + r'$  Parameters:$   $ a = %s $ ,   $f^{(0)}(x) = %s $,   $x^{*} = %s $' %(a,x0, x[-1])
        ax1.set_xlabel(r'Fertility (a)',fontsize = 16) #xlabel('Domain')
        ax1.set_ylabel(r'Population', fontsize = 16) #ylabel('Range')
        ax1.set_title(ti,fontsize=18)
    
    
    #Transformed Function Title
    ti_trans  =   r' $Transformed$  $Function$:  $ g = \alpha f + (1-\alpha)x $ ' + '\n' + r'$  Parameters:$   $ \alpha = %s $ '     %(alpha)
    ax2.set_xlabel(r'Fertility (a)',fontsize = 16) #xlabel('Domain')
    ax2.set_ylabel(r'Population', fontsize = 16) #ylabel('Range')
    ax2.set_title(ti_trans,fontsize=18)

    # **********************************************************************************************************************
    return 0

# Make list of images to create .gif
images = []

# 10 timesteps appears to be sufficient for generating a decent picture,
# easily modified but can take a while to run. 
for alpha in alpha_set:    
    
    # This section sets the stage on the first pass of the loop.
    # Whatever plot dimensions are initially used are reused for all
    # successive runs. Requires the loop to start at step ts = 0.
    # It keeps things from bouncind around in the animation.
    # --------------------------------------------------------------------------
    AniMe(alpha, alpha_set[0], a_set, x_set, maxstep, xmax, ax1, ax2)
    # --------------------------------------------------------------------------
    # Store image of figure (frame) to be written to the .gif    
    im = fig2img(fig)
    
    #imagefname = 'ts_%s.png' %(ts)
    #im.save(imagefname)    
    
    images.append(im)
        
    #print(ts)
    
    
# Write data to file
writeGif("images.gif",images,duration=1.5,dither=0)