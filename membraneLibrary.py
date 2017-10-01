# Functions and some dictionaries with parameters for different biophysical processes and ion channels found in excitable cell membranes

eCharge=1.60217733e-19
kBoltzmann=1.38065812e-20
zeroT=273.15

numerics = {'timeMin':0.0, 'timeMax':10.0, 'timeStep':3e-4,'rTol':1e-4, 'aTol':1e-9}

physicalConstants={'eCharge':1.60217733e-19, 'kBoltzmann':1.38065812e-20, 'zeroT':273.15}

biophysMembrane={'inNa': 10.0, 'outNa': 110.0, 'inKa': 100.0, 'outKa': 5.0,
                 'outCa': 2.5, 'inCa': 1e-4, 'inCl': 10.0, 'outCl': 100.0,
                 'Cm': 20.0,'tempCelsius':37.0, 'vATP':-450.0}

biophysCalcium={'kc': 2.5e-4, 'rc':1.0/250.0,'cMin':1e-5, 'cSlope': 1.5e4,'cHalf':5.6e-4,'c0':1e-4, 'mMmillionC':5200.0}

stpGlu={'ssPRelease':0.1, 'ratePRelease':0.1,'rateQuanta':0.1, 'ssQuanta':0.9}  

stpGABA={'ssPRelease':0.3, 'ratePRelease':0.1, 'rateQuanta':0.1, 'ssQuanta':0.9}  

biophysLFP={'Name':'Local Field Potential forcing','aMean':0.0,'aSD':0.33, 'symm':0.5,'aTau':1/2.0}
    
biophysGluAMPA={'Name':'AMPA Glutamate receptor','aMean':5.0,'aSD':2.0, 'symm':0.5,'v':0.0, 'aTau':5.0, 'nSynsPerContact':10.0}

biophysGABAA={'Name':'GABA receptor type A','aMean':6.0,'aSD':1.0, 'symm':0.5,'v':-72.0, 'aTau':10.0, 'nSynsPerContact':10.0}

biophysNaKa={'Name':'Natrium-Kalium ATPase','a':1.0, 'symm':0.5}  
biophysNaCa={'Name':'Natrium-Calcium pump','a':0.0, 'symm':0.5}  
biophysKaD={'Name':'Kalium delayed rectifier channel','a':1.0, 'symm':0.5,
            'vHalfP': -1.0, 'gainP': 3.0, 'rateP':2.5, 'expP':1.0, 'symmP': 0.3}  
biophysNaT={'Name':'Natrium transient channel','a':1.0, 'symm':0.5,
            'vHalfP': -25.0, 'gainP': 4.0, 'rateP':10.0, 'expP':1.0, 'symmP': 0.3,
            'vHalfQ': -67.5, 'gainQ':-4.3, 'rateQ':0.1, 'expQ':1.0, 'symmQ': 0.2}
biophysNaP={'Name':'Natrium persistent channel',
            'a':0.0, 'symm':0.5,
            'vHalfP': -40.0, 'gainP': 4.0, 'rateP':10.0, 'expP':1.0, 'symmP': 0.3}
biophysKaSK={'Name':'Kalium SK channel',
             'a':1.0, 'symm':0.5,
             'cHalfP': 5.6e-4, 'gainP': 1.5e4, 'rateP':2.0, 'expP':1.0, 'symmP': 0.3}
biophysCa12={'Name':'Cav12 channel',
             'a':1.0, 'symm':0.5,
             'vHalfP': -3.0, 'gainP': 4.0, 'rateP':10.0, 'expP':1.0, 'symmP': 0.3}
biophysCa13={'Name':'Cav13 channel',
             'a':0.0, 'symm':0.5,
             'vHalfP': -20.0, 'gainP': 4.0, 'rateP':10.0, 'expP':1.0, 'symmP': 0.3}
biophysCa31={'Name':'T-type Cav31 channel',
             'aCa31':0.0, 'sCa31':0.5,
             'pHCa31': -50.0, 'pGCa31': 5.8, 'pRCa31':2.0, 'pECa31':1.0, 'pSCa31': 0.3,
             'qHCa31': -74.0, 'qGCa31':-4.85, 'qRCa31':1.0, 'qECa31':1.0, 'qSCa31': 0.3}  
biophysCa32={'Name':'T-type Cav32 channel', 'aCa32':0.0, 'sCa32':0.5,
             'pHCa32': -40.0, 'pGCa32': 5.14, 'pRCa32':2.0, 'pECa32':1.0, 'pSCa32': 0.3,
             'qHCa32': -70.0, 'qGCa32':-4.3, 'qRCa32':1.0, 'qECa32':1.0, 'qSCa32': 0.3}  
biophysCa33={'Name':'T-type Cav33 channel', 'aCa33':0.0, 'sCa33':0.5,
             'pHCa33': -50.0, 'pGCa33': 4.3, 'pRCa33':2.0, 'pECa33':1.0, 'pSCa33': 0.3,
             'qHCa33': -75.0, 'qGCa33':-4.4, 'qRCa33':1.0, 'qECa33':1.0, 'qSCa33': 0.3}  
