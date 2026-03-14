import random
import h5py
import numpy as np
from numpy import linalg as LA
from iminuit import Minuit
from measurements import *
import matplotlib.pyplot as plt

from copy import deepcopy

#In this file we define a general fit class that allows to fit all different things we want to fit.



#Functions: different fit models, real time correlator, fourier transform, static, different channel and so on... See below how the parameters are initialized and handled by the fitting routine (constructor of Fitter class).

# If you want a new type of fit, change, add here.
def realtimeaxialcor(susc,t,parameters):
    (a,b,c) = parameters
   
    #omegap2 = v2 * (m**2 + p**2)
    #Delta = 0.5 * ( D * p**2 - Gamma * (p**2 + m**2))
    #Gammap = (D + Gamma) * p**2 + Gamma * m**2
    #Omegap= np.sqrt(np.complex(omegap2 - Delta**2))
    
    return susc * np.exp(-0.5 * a * t) * (np.cos(b * t) -c * np.sin(b* t))
    #return np.real(np.exp(-0.5 * Gammap * t) * (np.cos(Omegap * t) - Delta / Omegap * np.sin(Omegap * t)))
    
def realtimevec(susc,t,parameters):
    (a,b) = parameters
   
    return susc * np.exp(-a * t-b*t**2)
    
def realtimecondensate(t,parameters):
    (a,tau) = parameters
   
    return a * np.exp(- t / tau)
    



def axialprop(susc,omega,parameters):
    (mp,gammap) = parameters
    numerator= 2.0 * susc * gammap * mp**2
    denominator= np.square(-np.square(omega)+np.square(mp))+np.square(omega*gammap)
    return numerator/denominator

def phiprop(susc,omega, parameters):
    (mp,gammap) = parameters
    numerator= 2.0 * susc * gammap * np.square(omega) 
    denominator= np.square(-np.square(omega)+np.square(mp))+np.square(omega*gammap)
    return numerator/denominator
    
def phi0prop(omega, parameters):
    (msigma,gamma) = parameters
    numerator= 2 *gamma
    denominator= np.square(omega)+np.square(np.square(msigma)*gamma)
    return numerator/denominator

def staticPhiPropK(k, parameters):
    (chi, ms) = parameters
    numerator =  chi * np.tanh(ms / 2.0) * 2.0 * np.sinh(ms)
    denominator = 4.0 * np.sinh(ms / 2.0) ** 2 + 4.0 * np.sin(k / 2.0) ** 2
    return numerator / denominator

def staticPhiProp(L,x, parameters):
    (cs, ms) = parameters
    #numerator =  1.0 / 2.0 / ms / cs / L**2 * (np.exp(- ms * x) + np.exp(- ms * (L-x)))
    #denominator = 1.0 - np.exp(-ms * L)
    numerator =  1.0  / 2.0 / ms / cs / L**2 * (np.exp(- ms * x) + np.exp(- ms * (L-x)))
    denominator = 1 - np.exp(-ms * L)
    return numerator / denominator

def staticPhiRestoredProp(L,x, parameters):
    (cs, ms) = parameters
    #numerator =  1.0 / 2.0 / ms / cs / L**2 * (np.exp(- ms * x) + np.exp(- ms * (L-x)))
    #denominator = 1.0 - np.exp(-ms * L)
    numerator =  cs * np.exp(- ms * x) 
    return numerator

def staticPhiRestoredPropPeriodic(L,x, parameters):
    (cs, ms) = parameters
    #numerator =  1.0 / 2.0 / ms / cs / L**2 * (np.exp(- ms * x) + np.exp(- ms * (L-x)))
    #denominator = 1.0 - np.exp(-ms * L)
    numerator =  cs * (np.exp(- ms * x) + np.exp(- ms * (L-x)))
    #denominator = 1 - np.exp(-ms * L)
    return numerator 

def staticPhiPropKSusc(susc, k, parameters):
    ms = parameters[0]
    numerator =  susc * np.tanh(ms / 2.0) * 2.0 * np.sinh(ms)
    denominator = 4.0 * np.sinh(ms / 2.0) ** 2 + 4.0 * np.sin(k / 2.0) ** 2
    return numerator / denominator

def staticDSigmaPropK(k, parameters):
    (amp, ms) = parameters
    res = np.zeros(np.shape(k))
    c=0
    for kk in k:
        res[c] = amp / 2.0 / ms if np.abs(kk)<1e-10 else amp / kk * np.arctan(kk / 2.0 / ms)
        res[c] = amp - 2 * kk * np.arctan(kk / 2.0 / ms)
        c += 1
    return res
    
#classes

#Define the ch2 called by Minuit. Should not need to change. If you do simultanous fit, maybe need to double check.

class Chi2:
    def __init__(self):
        self.z = np.empty(0)
        self.resToAdd = []
        return
    def add(self, statRes):
        self.resToAdd.append(statRes)
    
    def computeWeight(self, data, err, prediction):
        return (data - prediction) / err
    
    def __call__(self,  listFunc, realOrImag = np.real):
        z = np.zeros(len(self.resToAdd[0].mean))
        for i in range(len(listFunc)):
            z += self.computeWeight(realOrImag(self.resToAdd[i].mean), realOrImag(self.resToAdd[i].err), listFunc[i])
        return np.sum(np.square(z))    

    
#Actual Fitter class.
class Fitter:
    def __init__(self, data, chi0, chiperp, L = 1, k = 0):
        obs = ["OtOttp", "OtOttpFourier", "propagatorF", "propagator"]
        
        self.V = L**3
        
        self.fits = dict()
        self.ndof = dict()
        self.averagechi2 = dict()
        self.averagechi2reduce = dict()
        self.func = dict()
        self.par = dict()
        self.parErr = dict()
        self.parName = dict()
        self.parLims = dict()
        
        self.chi0 = chi0
        self.chiperp = chiperp
        
        for o in obs:
            self.fits[o] = dict()
            self.ndof[o] = dict()
            self.averagechi2[o] = dict()
            self.averagechi2reduce[o] = dict()
            self.func[o] = dict()
            self.par[o] = dict()
            self.parErr[o] = dict()
            self.parName[o] = dict()
            self.parLims[o] = dict()

            
        self.data = data
                
        # This is a mapping between the fit key on the left and the actual quantity (ies) we are fitting. Allow to define combined fit in this way.
        self.baseKeys = dict()
        self.baseKeys["A"] = ["A"]
        self.baseKeys["V"] = ["V"]
        self.baseKeys["phi"] = ["phi"]
        self.baseKeys["phiRestored"] = ["phiRestored"]
        self.baseKeys["dsigma"] = ["dsigma"]
        self.baseKeys["dphi"] = ["dphi"]
        self.baseKeys["phi0"] = ["phi0"]
        self.baseKeys["Aphi"] = ["A", "phi"]
        
        if k > 0:
            self.baseKeys["Akk{}".format(k)] = ["Akk{}".format(k)]
            self.baseKeys["Vkk{}".format(k)] = ["Vkk{}".format(k)]
        
        self.xs = dict()
        
        
        
        # Actual definition of the function to fit, per observable and per fit key. The lambda function should receive x
        # and the fit parameters. It should return a list of array, each array corresponding to the function to be fitted
        # evaluated at x, one per observable fitted (so one for a simple fit and more for combined fit, equal to the length of baseKeys).
        # Can define a fit for any observable which is also defined in the measurements; we use meas.getObs(...) to retriee the data.
        
        #ignore
        #TODO: Define the static correlators at this level. Can define

        
        #Note to Alex: if you want to change the fit function, that is where you do it. You can also introduce new keys "Abiskk#" if for instance. If you do that, need to add the keys above as well.
        
        if k == 0:
            self.func["OtOttp"]["A"] = lambda x, par : [realtimeaxialcor(self.chi0, x, par)]
            self.func["OtOttp"]["V"] = lambda x, par : [realtimevec(self.chi0, x, par)]
        else:
            self.func["OtOttp"]["Akk{}".format(k)] = lambda x, par : [realtimeaxialcor(self.chi0, x, par)]
            self.func["OtOttp"]["Vkk{}".format(k)] = lambda x, par : [realtimevec(self.chi0, x, par)]
            
        self.func["OtOttp"]["phiRestored"] = lambda x, par : [realtimecondensate(x, par)]
        self.func["OtOttp"]["phi0"] = lambda x, par : [realtimecondensate(x, par)]

        
        self.func["OtOttpFourier"]["A"] = lambda x, par : [axialprop(self.chi0, x, par)]
        self.func["OtOttpFourier"]["phi"] = lambda x, par : [phiprop(self.chiperp, x, par)]
        self.func["OtOttpFourier"]["dphi"] = lambda x, par : [phiprop(self.chiperp, x, par)]
        self.func["OtOttpFourier"]["phi0"] = lambda x, par : [phi0prop( x, par)]
        
        self.func["OtOttpFourier"]["Aphi"] = lambda x, par : [axialprop(self.chi0, x, (par[0],par[1])), phiprop(self.chiperp, x, (par[0],par[1]))]
        
        self.func["propagatorF"]["phi"] = lambda x, par : [staticPhiPropK( x, par)]
        self.func["propagatorF"]["dsigma"] = lambda x, par : [staticDSigmaPropK( x, par)]
        
        self.func["propagator"]["phi"] = lambda x, par : [staticPhiProp(L, x, par)]
        self.func["propagator"]["phiRestored"] = lambda x, par : [staticPhiRestoredPropPeriodic(L, x, par)]
        
        
        
        #Note to Alex: That's where the parameters are defined. If you change them in the function, you may have to change here, especially if you cahnge the number of params.

        # Below we define the size, name and limits of the parameters, per observables and per fitKeys.
        
        key = "V" if k == 0 else "Vkk{}".format(k)
        self.par["OtOttp"][key] = np.zeros(2)
        self.parName["OtOttp"][key] = ["a","b"]
        self.parLims["OtOttp"][key] = [(0, 100),(-100, 100)]
        
        
        key = "A" if k == 0 else "Akk{}".format(k)
        #self.par["OtOttp"][key] = np.zeros(4)
        #self.parName["OtOttp"][key] = ["mp","Gamma", "D", "v2"]
        #self.parLims["OtOttp"][key] = [(0, None),(0, None), (0, None), (0, None)]
        self.par["OtOttp"][key] = np.zeros(3)
        self.parName["OtOttp"][key] = ["a", "b", "c"]
        self.parLims["OtOttp"][key] = [(0, 100),(0, 100), (-100,100)]
        
        
        self.par["OtOttp"]["phiRestored"] = np.zeros(2)
        self.parName["OtOttp"]["phiRestored"] = ["a","tau"]
        self.parLims["OtOttp"]["phiRestored"] = [(-1000, 1000),(0, 1000)]

        self.par["OtOttp"]["phi0"] = np.zeros(2)
        self.parName["OtOttp"]["phi0"] = ["a","tau"]
        self.parLims["OtOttp"]["phi0"] = [(-1000, 1000),(0, 1000)]
        
        
        self.par["OtOttpFourier"]["A"] = np.zeros(2)
        self.parName["OtOttpFourier"]["A"] = ["mp","gammap"]
        self.parLims["OtOttpFourier"]["A"] = [(0, None),(0, None)]
        
        self.par["OtOttpFourier"]["phi"] = np.zeros(2)
        self.parName["OtOttpFourier"]["phi"] = ["mp","gammap"]
        self.parLims["OtOttpFourier"]["phi"] = [(0, None), (0, None)]
        
        self.par["OtOttpFourier"]["dphi"] = np.zeros(2)
        self.parName["OtOttpFourier"]["dphi"] = ["mp","gammap"]
        self.parLims["OtOttpFourier"]["dphi"] = [(0, None), (0, None)]
        
        
        self.par["OtOttpFourier"]["Aphi"] = np.zeros(2)
        self.parName["OtOttpFourier"]["Aphi"] = ["mp","gammap"]
        self.parLims["OtOttpFourier"]["Aphi"] = [(0, None),(0, None)]
        


        self.par["propagatorF"]["phi"] = np.zeros(2)
        self.parName["propagatorF"]["phi"] = ["chi","m"]
        self.parLims["propagatorF"]["phi"] = [(0, None),(0, None)]

        self.par["propagatorF"]["dsigma"] = np.zeros(2)
        self.parName["propagatorF"]["dsigma"] = ["amp","m"]
        self.parLims["propagatorF"]["dsigma"] = [(0, None),(0, None)]
        
        self.par["propagator"]["phi"] = np.zeros(2)
        self.parName["propagator"]["phi"] = ["chi","m"]
        self.parLims["propagator"]["phi"] = [(0, None),(0, None)]
        
        self.par["propagator"]["phiRestored"] = np.zeros(2)
        self.parName["propagator"]["phiRestored"] = ["chi","m"]
        self.parLims["propagator"]["phiRestored"] = [(0, None),(0, None)]

        

        
    def setParValues(self, obs, key, pars):
        self.par[obs][key] = pars
      
    def fit(self, obs, key, minInd = 0, maxInd = None, prune = 1):                                      
        chi2 = Chi2()
        
        for k in self.baseKeys[key]:
            (tmpx, tmpfx) = self.data.getObs(obs, k)
            x = deepcopy(tmpx)
            fx = deepcopy(tmpfx)
            x = x[minInd:maxInd:prune]
            fx.reduce(minInd, maxInd,prune)
            #print(fx.mean)
            chi2.add(fx)
            self.xs[obs] = x
                                           
        self.fits[obs][key] = Minuit(lambda x: chi2(self.func[obs][key]( self.xs[obs], x)), self.par[obs][key], name=self.parName[obs][key])
        
        
        self.fits[obs][key].errordef = Minuit.LEAST_SQUARES
        
        self.fits[obs][key].limits = self.parLims[obs][key]
                
        self.fits[obs][key].migrad()
        self.fits[obs][key].minos()
                
        self.ndof[obs][key] = len(self.baseKeys[key]) * len(self.xs[obs]) - len(self.par[obs][key])
        self.averagechi2[obs][key] = self.fits[obs][key].fval
        self.averagechi2reduce[obs][key] = self.fits[obs][key].fval / self.ndof[obs][key]
        
        c = 0
        for v in self.fits[obs][key].values:
            self.par[obs][key][c] = v
            c += 1   
            
        c = 0
        self.parErr[obs][key] = np.zeros(len(self.par[obs][key]))
        for v in self.fits[obs][key].errors:
            self.parErr[obs][key][c] = v
            c += 1   
        
        return (self.averagechi2[obs][key], self.averagechi2reduce[obs][key], self.ndof[obs][key], self.fits[obs][key].values)

    
    def plot(self, obs, key, color="b"):
        if obs in self.xs.keys():
            x = self.xs[obs]
        else:
            x = self.data.getObs(obs, key, withY = False)
        
        #print(self.par[obs][key])
        prediction = self.func[obs][key](x, self.par[obs][key])
        for p in prediction:
            plt.plot(x, p, color)
