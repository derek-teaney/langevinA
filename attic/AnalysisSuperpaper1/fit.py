import random
import h5py
import numpy as np
from numpy import linalg as LA
from iminuit import Minuit
from measurements import *
import matplotlib.pyplot as plt


#Functions
def axialprop(omega,mp,amplitude,gammap):
    numerator= amplitude
    denominator= np.square(-np.square(omega)+np.square(mp))+np.square(omega*gammap)
    return numerator/denominator

def phiprop(omega,mp,amplitude,gammap):
    numerator= amplitude*np.square(omega)
    denominator= np.square(-np.square(omega)+np.square(mp))+np.square(omega*gammap)
    return numerator/denominator
    
def phi0prop(omega,msigma,gamma):
    numerator= 2 *gamma
    denominator= np.square(omega)+np.square(np.square(msigma)*gamma)
    return numerator/denominator
#classes


#Work in proigress, store the results and call the corersponding measurements in a user friendly way.
class FitResult:
    def __init__(self,data):
        self.data = data
        
        self.par = dict()
        
        self.par["mp"] = 0
        self.par["ampCharge"]=0
        self.par["ampPhi"]=0
        self.par["gammap"]=0
        
        
        
    def chisquareaxial(self,parameter):
        (mp,amplitudecharge,gammap)=parameter
        expdata=np.real(self.data.OtOttpFourier["A"].mean)
        prediction=axialprop(self.data.OtOttpFourier_oms["A"],mp,amplitudecharge,gammap)
        error= np.real(self.data.OtOttpFourier["A"].err)
        z = (expdata - prediction) / error
        return np.sum(np.square(z))
    
    def chisquarphi(self,parameter):
        (mp,amplitudephi,gammap)=parameter
        expdata=np.real(self.data.OtOttpFourier["phi"].mean)
        prediction=phiprop(self.data.OtOttpFourier_oms["phi"],mp,amplitudephi,gammap)
        error= np.real(self.data.OtOttpFourier["phi"].err)
        z = (expdata - prediction) / error
        return np.sum(np.square(z))
        
        
    def chisquareaxialphi(self,parameter):
        (mp,amplitudephi,amplitudecharge,gammap)=parameter
        expdata=np.real(self.data.OtOttpFourier["A"].mean)
        prediction=axialprop(self.data.OtOttpFourier_oms["A"],mp,amplitudecharge,gammap)
        error= np.real(self.data.OtOttpFourier["A"].err)
        z = (expdata - prediction) / error
        expdata=np.real(self.data.OtOttpFourier["phi"].mean)
        prediction=phiprop(self.data.OtOttpFourier_oms["phi"],mp,amplitudephi,gammap)
        error= np.real(self.data.OtOttpFourier["phi"].err)
        z+=(expdata - prediction) / error
        return np.sum(np.square(z))
    
    def modelaxial(self):
        return axialprop(self.data.OtOttpFourier_oms["A"],self.par["mp"],self.par["ampCharge"],self.par["gammap"])
        
    def modelphi(self):
        return phiprop(self.data.OtOttpFourier_oms["phi"],self.par["mp"],self.par["ampPhi"],self.par["gammap"])
    
    def fitAA(self):
        self.AAfit=Minuit(lambda x: self.chisquareaxial(x), (self.par["mp"],self.par["ampCharge"],self.par["gammap"]),name=("mp","amplitudecharge","gammap"))
        self.AAfit.errordef = Minuit.LEAST_SQUARES
        self.AAfit.limits=[(0, None), (0, None),(0, None)]
        self.AAfit.migrad()
        self.AAfit.minos()
        self.par["mp"]=self.AAfit.values[0]
        self.par["ampCharge"]=self.AAfit.values[1]
        self.par["gammap"]=self.AAfit.values[2]
        self.ndof= len(self.data.OtOttpFourier_oms["A"])-3
        self.averegechi2=self.AAfit.fval
        self.averegechi2reduce=self.AAfit.fval /self.ndof
        return (self.averegechi2,self.averegechi2reduce,self.ndof,self.AAfit.values)
        
    
    def fitPP(self):
        self.PPfit=Minuit(lambda x: self.chisquarphi(x), (self.par["mp"],self.par["ampPhi"],self.par["gammap"]),name=("mp","amplitudephi","gammap"))
        self.PPfit.errordef = Minuit.LEAST_SQUARES
        self.PPfit.limits=[(0, None), (0, None),(0, None)]
        self.PPfit.migrad()
        self.PPfit.minos()
        self.par["mp"]=self.PPfit.values[0]
        self.par["ampPhi"]=self.PPfit.values[1]
        self.par["gammap"]=self.PPfit.values[2]
        self.ndof= len(self.data.OtOttpFourier_oms["phi"])-3
        self.averegechi2=self.PPfit.fval
        self.averegechi2reduce=self.PPfit.fval/self.ndof
        return (self.averegechi2,self.averegechi2reduce,self.ndof,self.PPfit.values)
        
    def fitAAPP(self):
        self.AAPPfit=Minuit(lambda x: self.chisquareaxialphi(x), (self.par["mp"],self.par["ampPhi"],self.par["ampCharge"],self.par["gammap"]),name=("mp","amplitudephi","amplitudecharge","gammap"))
        self.AAPPfit.errordef = Minuit.LEAST_SQUARES
        self.AAPPfit.limits=[(0, None), (0, None),(0, None),(0, None)]
        self.AAPPfit.migrad()
        self.AAPPfit.minos()
        self.par["mp"]=self.AAPPfit.values[0]
        self.par["ampPhi"]=self.AAPPfit.values[1]
        self.par["ampCharge"]=self.AAPPfit.values[2]
        self.par["gammap"]=self.AAPPfit.values[3]
        self.ndof= len(self.data.OtOttpFourier_oms["phi"])+len(self.data.OtOttpFourier_oms["A"])-4
        self.averegechi2=self.AAPPfit.fval
        self.averegechi2reduce=self.AAPPfit.fval/self.ndof
        return (self.averegechi2,self.averegechi2reduce,self.ndof,self.AAPPfit.values)
        
        
    def AApplot(self,mp=None,amplitudecharge=None,gammap=None):
        if mp is not None:
            self.par["mp"]=mp
        if amplitudecharge is not None:
            self.par["ampCharge"]=amplitudecharge
        if gammap is not None:
            self.par["gammap"]=  gammap
        plt.errorbar(self.data.OtOttpFourier_oms["A"], self.data.OtOttpFourier["A"].mean, self.data.OtOttpFourier["A"].err)
        plt.plot(self.data.OtOttpFourier_oms["A"],self.modelaxial())
        
        
    def PPpplot(self,mp=None,amplitudephi=None,gammap=None):
        if mp is not None:
            self.par["mp"]=mp
        if amplitudephi is not None:
            self.par["ampPhi"]=amplitudephi
        if gammap is not None:
            self.par["gammap"]=gammap
        plt.errorbar(self.data.OtOttpFourier_oms["phi"], self.data.OtOttpFourier["phi"].mean, self.data.OtOttpFourier["phi"].err)
        plt.plot(self.data.OtOttpFourier_oms["phi"],self.modelphi())
            

    

    
