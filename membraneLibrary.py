import numpy as np
import matplotlib.pylab as gr
from matplotlib.collections import LineCollection
import pyprocess as sp
print("Imported sympy as sy, numpy as np, pylab as gr, pyprocess as sp")
#

eCharge=1.60217733e-19 # Coulombs
kBoltzmann=1.38065812e-20 #mJ/K
zeroT=273.15 #deg Kelvin

numerics = {'timeMin':0.0, 'timeMax':10.0, 'timeStep':0.025,'rTol':1e-4, 'aTol':1e-9}

physicalConstants={'eCharge':1.60217733e-19, 'kBoltzmann':1.38065812e-20, 'zeroT':273.15,
    "Faraday":96485.3329, "NA":6.02214086e23}

ionsTemperature={'in_Na': 15.0, 'out_Na': 140.0, 'in_Ka': 140.0, 'out_Ka': 5.0,
                 'out_Ca': 2.5, 'in_Ca': 1e-4, 'in_Cl': 10.0, 'out_Cl': 150.0,
                 'Cm': 20.0,'tempCelcius':37.0, 'v_ATP':-430.0}

biophysCalcium={'k': 2.5e-4, 'r_Ca':1.0/25.0,
    'c_Min':1e-3, 'c_Init':1e-1, 'mMmillionC':5200.0,"freeBoundRatio":0.1}

# Transmembrane transport by channels and pumps
# Pumps
biophysNaKa={#'Name':'Natrium-Kalium ATPase',
    'eta':1.0, "v_r":-60.0,
    'a':1.0e-3, "a_Single":1.0, "N":1.0, 'rect':0.5,
    "nNa":3, "nKa":2}
# Natrium channels
biophysNaCa={#'Name':'Natrium-Calcium pump',
    'eta':-1.0, "v_r":40.0,
    'a':0.0, 'rect':0.5,
    "nNa":3, "nCa":1}
biophysNaT={#'Name':'Natrium transient channel',
    'eta':-1.0,
    'a':1.0, 'a_Single':1.2, "N":0.0, 'rect':0.5, 'v_r': -65.0,
    'v_Half_Act': -25.0, 'gain_Act': 4.0, 'rate_Act':10.0, 'exp_Act':1.0, 'bias_Act': 0.3,
    'v_Half_Inact': -67.5, 'gain_Inact':-4.3, 'rate_Inact':0.1, 'exp_Inact':1.0, 'bias_Inact': 0.2,
    }
biophysNaP={#'Name':'Natrium persistent channel',
    'eta':-1.0,
    'a':0.0, 'a_Single':1.2, "N":2000.0,  'rect':0.5,
    'v_Half_Act': -40.0, 'gain_Act': 4.0, 'rate_Act':10.0, 'exp_Act':1.0, 'bias_Act': 0.3}
# Kalium channels
biophysKaD={#'Name':'Kalium delayed rectifier channel',
    'eta':1.0,
    'a':1.0, 'a_Single':1.0, "N":0.0,  'rect':0.5, "v_r":-90.0,
    'v_Half_Act': -1.0, 'gain_Act': 3.0, 'rate_Act':2.5, 'exp_Act':1.0, 'bias_Act': 0.3,
    'v_Half_Inact': 100.0, 'gain_Inact':-4.0, 'rate_Inact':0.000001, 'exp_Inact':1.0, 'bias_Inact': 0.5}

biophysKaSK={#'Name':'Kalium SK channel',
    'eta':1.0,
    'a':1.0e-3, 'a_Single':5.0, #up to 10 pA at high Ca
    "N":0.0,  'rect':0.5, "v_r":-90.0,
    'c_Half_Act': 0.3,'gain_Act':3.0}

biophysKaIK={#'Name':'Kalium intermediate Ca-dependent K-current',
    'eta':1.0,
    'a':1.0e-3, 'a_Single':5.0, #up to 10 pA
    "N":0.0,  'rect':0.5, "v_r":-90.0,
    'c_Half_Act': 0.1,
    'v_Half_Inact': 100.0, 'gain_Inact':-4.0, 'rate_Inact':0.000001, 'exp_Inact':1.0, 'bias_Inact': 0.5}

biophysKaBK={#'Name':'Kalium BK channel',
    'eta':1.0,
    'a':1.0, 'a_Single':10.0, #up to 20 pA
    "N":0.0,  'rect':0.5,
    'cHalf_Act': 5.6e-4, 'gain_Act': 1.5e4, 'rate_Act':2.0, 'exp_Act':1.0, 'bias_Act': 0.3,
    'v_Half_Inact': 100.0, 'gain_Inact':-4.0, 'rate_Inact':0.000001, 'exp_Inact':1.0, 'bias_Inact': 0.5}

# Calcium channels
biophysCa12={#'Name':'Cav12 L-type channel',
    'eta':-2.0,
    'a':1.0e-3, 'a_Single':0.33, "N":0.0,  'rect':0.5,
    'v_Half_Act': -3.0, 'gain_Act': 4.0, 'rate_Act':10.0, 'exp_Act':1.0, 'bias_Act': 0.3,
    'v_Half_Inact': 100.0, 'gain_Inact':-4.0, 'rate_Inact':0.000001, 'exp_Inact':1.0, 'bias_Inact': 0.5}

biophysCa13={#'Name':'Cav13 L-type channel',
    'eta':-2.0,
    'a':0.0e-3, 'a_Single':0.33, "N":0.0,  'rect':0.5,
    'v_Half_Act': -20.0, 'gain_Act': 4.0, 'rate_Act':10.0, 'exp_Act':1.0, 'bias_Act': 0.3,
    'v_Half_Inact': 100.0, 'gain_Inact':-4.0, 'rate_Inact':0.000001, 'exp_Inact':1.0, 'bias_Inact': 0.5}

biophysCa31={#'Name':'T-type Cav31 channel',
    'eta':-2.0,
    'a':0.0e-3, 'a_Single':.25, "N":0.0,  'rect':0.5,
    'v_Half_Act': -50.0, 'gain_Act': 5.8, 'rate_Act':2.0, 'exp_Act':0.0, 'bias_Act': 0.3,
    'v_Half_Inact': -74.0, 'gain_Inact':-4.85, 'rate_Inact':1.0, 'exp_Inact':0.0, 'bias_Inact': 0.3}
biophysCa32={#'Name':'T-type Cav32 channel',
    'eta':-2.0,
    'a':1.0e-3, 'a_Single':.25, "N":0.0,  'rect':0.5,
    'v_Half_Act': -43.23, 'gain_Act': 5.31, 'rate_Act': 1/17.25, 'exp_Act': 0.0, 'bias_Act': 0.865,
    'v_Half_Inact': -75.02, 'gain_Inact':-4.45, 'rate_Inact': 1/550.0, 'exp_Inact': 0.0, 'bias_Inact': 0.35}
biophysCa33={#'Name':'T-type Cav33 channel',
    'eta':-2.0,
    'a':1.0e-3, 'a_Single':.25, "N":0.0,  'rect':0.5,
    'v_Half_Act': -34.81, 'gain_Act': 4.45, 'rate_Act': 1/15.0, 'exp_Act': 0.0, 'bias_Act': 0.5,
    'v_Half_Inact': -68.88, 'gain_Inact':-4.53, 'rate_Inact':1/4850.0, 'exp_Inact': 0.0, 'bias_Inact': 0.31}
biophysLFP={#'Name':'Local Field Potential forcing',
    'a_Mean':0.0,'a_SD':0.33, 'tau':1/2.0}
biophysAMPA={#'Name':'AMPA Glutamate receptor',
    'eta':-1.0, 'a_Single':1.2, "N":0.0,  'a':1.0,
    'a_Mean':0.005, 'a_SD':0.002, 'rect':0.5,'v_r':0.0, 'tau':5.0, 'nSynsPerContact':10.0}
biophysGabaA={#'Name':'GABA receptor type A',
    'eta':1.0, 'a_Single':1.2,"N":0.0, 'a':1.0,
    'a_Mean':0.006,'a_SD':0.001, 'rect':0.5,'v_r':-72.0, 'tau':10.0, 'nSynsPerContact':10.0}

# Short-term synaptic plasticity
biophysSTP={'ssPRelease':0.3, 'rate_ActRelease':0.1,
    'rate_Quanta':0.1, 'ssQuanta':0.9, "jumpPRelase":0.1, 'exp_Quanta':1.0,
    'v_Half_Act':-20.0, 'gain_Act': 5.0, 'rate_Act':1.0, 'exp_Act':1.0, 'bias_Act': 0.5}
# Synaptic input
stpGlu={'ssPRelease':0.1, 'rate_ActRelease':0.1,
    'rate_Inactuanta':0.1, 'ssQuanta':0.9, "jumpPRelase":0.1}
stpGaba={'ssPRelease':0.3, 'rate_ActRelease':0.1,
    'rate_Inactuanta':0.1, 'ssQuanta':0.9, "jumpPRelase":0.1}

# ---------------------------------------
# ---------------------------------------
# Auxiliary functions:
# Calculate Nernst potentials and transmembrane fluxes and currents
# ---------------------------------------
# ---------------------------------------
def sigmoid(a,b):
    return a/(a+b)

def MM(a,b,n):
    aa= a**n
    return aa/(aa + b**n)

def expSum_(y, s):
    ee = np.exp( y*s )
    return ee * (1 +  np.exp(-y))

def expSigmoid_(y):
    ee = np.exp(y)
    return ee/(1+ee)

def vBoltzmann_(tempCelcius=22.0):
    """Calculates the Boltzmann or thermal potential at the specified temperature"""
    return kBoltzmann *(zeroT+tempCelcius)/eCharge

def vNernst_(cOut=5.0, cIn=140.0, val=1.0, tempCelcius=22.0):
    """Calculates the Nernst potential
    vx=vNernst(cIn=140.0, cOut=5.0, val=1.0, tempCelcius=22.0)
    """
    vT= kBoltzmann*(zeroT+tempCelcius)/eCharge
    vN= vT * np.log(cOut/cIn) / val
    return vN

def vTNernst_(cIn=140.0, cOut=5.0, val=1.0):
    """Calculates the Nernst potential
    vx=vTNernst(cIn=140.0, cOut=5.0, val=1.0)
    """
    yN= np.log(cOut/cIn) / val
    return yN

def calcReversalPotentials(parDict):
    parDict['v_T'] = vBoltzmann_(parDict['tempCelcius'])
    #parDict['v_ATP']= -450.0;
    parDict['v_Cl']=vNernst_(cIn=parDict['in_Cl'], cOut=parDict['out_Cl'], val=-1.0, tempCelcius=parDict['tempCelcius'])
    parDict['v_Ca']=vNernst_(cIn=parDict['in_Ca'], cOut=parDict['out_Ca'], val=2.0, tempCelcius=parDict['tempCelcius'])
    parDict['v_Na']=vNernst_(cIn=parDict['in_Na'], cOut=parDict['out_Na'], val=1.0, tempCelcius=parDict['tempCelcius'])
    parDict['v_Ka']=vNernst_(cIn=parDict['in_Ka'], cOut=parDict['out_Ka'], val=1.0, tempCelcius=parDict['tempCelcius'])
    parDict['v_NaKa'] = 3*parDict['v_Na'] - 2*parDict['v_Ka'] + parDict['v_ATP']
    parDict['v_NaCa'] = -3 *parDict['v_Na'] + 2*parDict['v_Ca']
    parDict['v_AMPA'] = 0.0
    parDict['v_GabaA']= parDict["v_Cl"]
    return parDict

def normalizeVolts(parDict):
    """Normalize voltages with respect to the Boltzmann potential kT/q"""
    nDict={}
    vTCm= parDict['Cm'] * parDict['v_T']
    nDict['vTCm']=vTCm
    for k,i in parDict.items():
        if ((k.find('v_')==0)):
            #print(k)
            nDict["y_"+k[2:]]=parDict[k]/parDict['v_T']

    parDict.update(nDict)
    return parDict

def normalizeAmps(parDict):
    """
    Normalizes amplitudes ("a") for transmembrane currents with respect to the membrane capacitance.
    The new amplitudes are written with capital "A".
    """
    nDict={}
    vC= parDict['Cm'] * parDict['v_T']
    nDict['vTCm']=vC

    for k,i in parDict.items():
        if k.find('a_')==0:
            #print(k,i)
            str0="A_" + k[2:]
            #if type(i)==float:
            nDict[str0]= i/vC
            #print(nDict[str0])
    parDict.update(nDict)
    return parDict

def calcAmpsPA(parDict):
    """
    Calculates amplitudes ("a") for transmembrane currents with respect to the membrane capacitance.
    The new amplitudes are written with capital "A".
    """
    nDict={}
    vC= parDict['Cm'] * parDict['v_T']
    nDict['vTCm']=vC

    for k,i in parDict.items():
        if k.find('A_')==0:
            #print(k,i)
            str0="a_" + k[2:]
            #if type(i)==float:
            nDict[str0]= i*vC
            #print(nDict[str0])
    parDict.update(nDict)
    return parDict

def addSuffix(p,suff='',removeOriginal=1):
    pa= p.copy()
    for key, val in p.items():
        pa[key+suff]=val
        if removeOriginal:
            pa.pop(key)
    return pa

def prepMembraneDictionary(membraneList,suffixList):
    p= dict()
    p.update(ionsTemperature)
    p.update(numerics)
    nSuff=len(suffixList)
    for n in range(nSuff):
        p.update(addSuffix(membraneList[n],suffixList[n]))
    p= calcReversalPotentials(p)
    p= normalizeAmps(p)
    p['degree'] = np.e
    p= normalizeVolts(p)
    return p

def printAmps(dicti):
	for k,v in dicti.items():
		if (str.find(k,'A_')==0) & (str.find(k,'A_Single_') !=0):
			print('%s = %g'%(k,v))
	return

def printAmpsPA(dicti):
	for k,v in dicti.items():
		if (str.find(k,'a_')==0) & (str.find(k,'a_Single_') !=0):
			print('%s = %g'%(k,v))
	return

# --------------------------------------------------------
# Colored lines in plots
# Usage:
# lc = colored_line(x,y,c)
# ax.add_collection(lc)
# --------------------------------------------------------
def colored_line(x, y, c, **kwargs):
    """Helper function to take a set of points and turn them into a collection of
    lines colored by another array
    Parameters
    ----------
    x : array-like
        x-axis coordinates
    y : array-like
        y-axis coordinates
    c : array-like
        values used for color-mapping
    kwargs : dict
        Other keyword arguments passed to :class:`matplotlib.collections.LineCollection`
    Returns
    -------
        The created :class:`matplotlib.collections.LineCollection` instance.
    """
    # Paste values end to end
    points = np.concatenate([x, y])

    # Exploit numpy's strides to present a view of these points without copying.
    # Dimensions are (segment, start/end, x/y). Since x and y are concatenated back to back,
    # moving between segments only moves one item; moving start to end is only an item;
    # The move between x any moves from one half of the array to the other
    num_pts = points.size // 2
    final_shape = (num_pts - 1, 2, 2)
    final_strides = (points.itemsize, points.itemsize, num_pts * points.itemsize)
    segments = np.lib.stride_tricks.as_strided(points, shape=final_shape,
                                               strides=final_strides)

    # Create a LineCollection from the segments and set it to colormap based on c
    lc = LineCollection(segments, **kwargs)
    lc.set_array(c)
    return lc


# --------------------------------------------------------
# define OU forcing function
# --------------------------------------------------------
def ouProcess(delta,nSteps,diffC,tau=1,graph=0):
    X=np.zeros(nSteps)
    N = np.random.randn(nSteps)
    for n in range(1,nSteps):
        edt = np.exp(-delta/tau) 
        X[n] = X[n-1]*edt + np.sqrt((1-edt**2)*(diffC*tau)/2) * N[n]
    if graph:
        gr.figure(figsize=(17,5))
        gr.plot(X, label=r'$OU(t)$')
        gr.legend()
    return X


def ouForcing(sampTimes,mu=0,SD=1,tau=1,graph=0):
    samplePath=np.zeros(len(sampTimes))
    if SD>0:
        diffF = 2* SD / tau
        ou = sp.OU_process(theta=1/tau, mu=0, sigma=diffF)
        samplePath = mu+ ou.sample_path(sampTimes)[0]
        if graph:
            gr.figure(figsize=(17,5))
            gr.plot(sampTimes, samplePath, label=r'$OU(t)$')
            gr.legend()
    else:
        print('Could not simulate OU process')
    return samplePath

# --------------------------------------------------------
# Now we define the current injection function used to stimulate model neurons as in experiments.
# The timing of the pulse can be modified to simulate a square pulse, or ramps of different slopes.
# --------------------------------------------------------

# define current injection function
def UpTopDn(t,upStart=200.0,upStop=400, dnStart=600.0,dnStop=800.0,rampAmp=1.0):
    slope_up = rampAmp/(upStop-upStart)
    int_up = -slope_up*upStart
    slope_dn = -rampAmp/(dnStop-dnStart)
    int_dn= rampAmp-slope_dn*dnStart
    c1=np.int16((upStart<t)&(t<=upStop))
    c2=np.int16((upStop<t)&(t<=dnStart))
    c3=np.int16((t>dnStart)&(t<=dnStop))
    y=c1*(slope_up*t +int_up)+ c2*rampAmp+ c3*(slope_dn*t +int_dn)
    return y

def spikeInds(dv,dvdtThresh=100):
    i = np.where(dv>dvdtThresh)[0]-1;
    di= np.where(i[1:]-i[:-1]>1)[0]+1;
    si = list()
    si.append(i[0])
    for n in range(len(di)):
        si.append(i[di][n])
    return si

def calcISI(spikeTimes):
    isis= np.zeros(len(spikeTimes))
    isis[1:] = (spikeTimes[1:]-spikeTimes[:-1])
    return isis

def calcIFR(spikeTimes):
    ifrs= np.zeros(len(spikeTimes))
    ifrs[1:] =1/ (spikeTimes[1:]-spikeTimes[:-1])
    return ifrs

# -----------------------
# Numerics
# -----------------------
def autonomousRK2_Step(f,U,p):
    k1 = p['stepSize'] * f( U, p) / 2.0
    return U + p['stepSize'] * f( U + k1, p)

def autonomousRK2_Solver(f, p, parNames=[],parValues=[]):
    """Second-order Runge-Kutta method to solve x' = f(x) with U(t[0]) = U0.
    NOTES:
        This version is based on the algorithm presented in "Numerical
        Analysis", 6th Edition, by Burden and Faires, Brooks-Cole, 1997.
    """
    nDimensions = len(p['ic'])
    nForc=len(parNames)
    if p['nNeurons'] == 1:
        U=np.zeros((p['nSteps'], p['nDimensions']),"float64")
    else:
        U=np.zeros((p['nSteps'], p['nDimensions'],p['nNeurons']),"float64")

    U[0]=p['ic']
    if nForc>0:
        for i in range(p['nSteps']-1):
            for nn in range(nForc):
                p[parNames[nn]]=parValues[nn][i]
            U[i+1] = autonomousRK2_Step(f,U[i],p)
    else:
        for i in range(p['nSteps']-1):
            U[i+1] = autonomousRK2_Step(f,U[i],p)
    return U.transpose()





