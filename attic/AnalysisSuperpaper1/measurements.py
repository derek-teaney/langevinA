import random
import h5py
import numpy as np
from numpy import linalg as LA
from copy import deepcopy
import matplotlib.pyplot as plt
import matplotlib.style as style
from scipy.interpolate import InterpolatedUnivariateSpline

import multiprocessing as mp

import functools
import pickle


class StatResult:

    def __init__(self, tup):
        self.mean = tup[0]
        self.err = tup[1]

    def rescale(self, fact):
        self.mean *= fact
        self.err *= fact

    def reduce(self, minInd=0, maxInd=None, prune=1):
        self.mean = self.mean[minInd:maxInd:prune]
        self.err = self.err[minInd:maxInd:prune]

    def save_to_txt(self,
                    fn,
                    x=None,
                    fmt="python",
                    decim=1,
                    tMin=0,
                    tMax=None,
                    yfact=1):
        if (len(np.shape(x)) == 0 and x == None):
            toSave = np.column_stack([
                (yfact * self.mean[tMin:tMax:decim] + 0.0j).view(float).reshape(
                    -1, 2),
                (yfact * self.err[tMin:tMax:decim] + 0.0j).view(float).reshape(
                    -1, 2),
            ])
            if fmt == "python" or fmt == "gnuplot":
                np.savetxt(fn, toSave)
            elif fmt == "latex":
                np.savetxt(fn,
                           toSave,
                           header="rey imy errrey errimy",
                           comments='')
        else:
            toSave = np.column_stack([
                x[tMin:tMax:decim],
                (yfact * self.mean[tMin:tMax:decim] + 0.0j).view(float).reshape(
                    -1, 2),
                (yfact * self.err[tMin:tMax:decim] + 0.0j).view(float).reshape(
                    -1, 2),
            ])
            if fmt == "python" or fmt == "gnuplot":
                np.savetxt(fn, toSave)
            elif fmt == "latex":
                np.savetxt(fn,
                           toSave,
                           header="x rey imy errrey errimy",
                           comments='')


def load_StatResult_from_txt(fn, withX=False):
    if withX:
        x, fR, fI, errR, errI = np.loadtxt(fn, unpack=True)
        return x, StatResult((fR + 1j * fI, errR + 1j * errI))
    else:
        f, err = np.loadtxt(fn, unpack=True).view(complex)
        return StatResult((f, err))


# Returns (mean, error) of an array, arr with a bootstrap resampling.
#
# A function of the array can be provided. For example to compute the mean  of
# arr^2 we have:
#
# booststrap(array, func= lambda x: x*x)
def bootstrap(arr, nSamples=100, func=lambda x: x):
    arr = np.asarray(arr)
    random.seed()
    bootstrapArr = []
    for n in range(nSamples):
        newArr = np.zeros(np.shape(arr), dtype=arr.dtype)
        for l in range(len(arr)):
            newArr[l] = arr[random.randrange(0, len(arr))]
        bootstrapArr.append(np.mean(func(newArr), axis=0))

    bootstrapArr = np.asarray(bootstrapArr)

    return (np.mean(bootstrapArr, axis=0), np.std(bootstrapArr, axis=0))
    # return (np.mean(arr,axis = 0), np.std(arr,axis = 0))


# Returns (mean, error) of an array, arr with jacknife
def jackknife(arr, nSamples=10, func=lambda x: x):
    nBlock = int(len(arr) / nSamples)
    arr = np.asarray(arr)
    random.seed()
    bootstrapArr = []
    count = 0
    thetaBar = np.mean(func(arr), axis=0)
    sigma2 = np.zeros(1 if len(np.shape(arr)) == 1 else np.shape(arr[0, :]),
                      dtype=arr.dtype)
    for n in range(nSamples):
        newArr = func(np.concatenate((arr[:count], arr[count + nBlock:])))
        bootstrapArr.append(np.mean(newArr, axis=0))
        count += nBlock
        sigma2 += (bootstrapArr[-1] - thetaBar)**2

    bootstrapArr = np.asarray(bootstrapArr)

    return (np.mean(bootstrapArr,
                    axis=0), np.sqrt((nSamples - 1.0) / nSamples * sigma2))
    # return (np.mean(arr,axis = 0), np.std(arr,axis = 0))


# Returns (mean, error) of an array, arr with blocking
def blocking(arr, nBlock=10, func=lambda x: x):
    blockSize = int(len(arr) / nBlock)
    arr = np.asarray(arr)
    blocks = []
    count = 0
    for n in range(nBlock):
        blocks.append(np.mean(arr[n * blockSize:(n + 1) * blockSize], axis=0))

    return (np.mean(np.asarray(blocks),
                    axis=0), np.std(np.asarray(blocks), axis=0), np.asarray(blocks))


# Returns (mean, error) of an array, arr with blocking and bootstrap
def blockedBootstrap(arr, nBlock=10, nSamples=100, func=lambda x: x):
    blockSize = int(len(arr) / nBlock)
    arr = np.asarray(arr)
    random.seed()
    blocks = []
    count = 0
    for n in range(nBlock):
        blocks.append(arr[n * blockSize:(n + 1) * blockSize])

    bootstrapArr = []
    for n in range(nSamples):
        newArr = []
        for l in range(nBlock):
            newArr.append(blocks[random.randrange(0, nBlock)])
        bootstrapArr.append(np.mean(np.mean(func(newArr), axis=0), axis=0))

    bootstrapArr = np.asarray(bootstrapArr)

    return (np.mean(bootstrapArr, axis=0), np.std(bootstrapArr, axis=0))


# Adapted from Hauke Sandmeyer's routines. (AF)
def autocorrelationFunction(data, tmax=1000, startvalue=0):
    int_corr = 0.5
    sq_sum = 0.0

    corr_arr = []

    norm = np.sum((data - np.mean(data, axis=0))**2)

    mean = np.mean(data, axis=0)

    for t in range(0, tmax):
        corr = 0.0
        for i in range(startvalue, len(data) - t):
            corr += (data[i] - mean) * (data[i + t] - mean)
        corr /= norm
        corr_arr.append(corr)

    return corr_arr


def intAutocorrrelationTime(arr, tmax=1000, startvalue=0, nMax=-1):
    Ct = autocorrelationFunction(arr[:nMax], tmax, startvalue)
    res = 0.5
    count = 0
    while float(count) < 5 * res and count < len(Ct) / 2.0:
        res += Ct[count]
        count += 1
    return res


# We don't need necessarily need the whole statistics to evaluate the
# autocorrelation time, nMax is the number of points we should use.
def autocorrelatedErr(arr, tmax=1000, startvalue=0, nMax=-1):
    mean = np.mean(arr, axis=0)
    naiveErr = np.std(arr, axis=0)
    tau = intAutocorrrelationTime(arr, tmax, startvalue, nMax)
    return (mean, np.sqrt(2.0 * tau) * naiveErr)


# Compute the equal time correlator between arr1 and arr2 and return it as an
# array.
#
# By default the routine will not subtract the connected part.  The routine
# will do the substraction, if passed with conn=True (conn is connected),
#
def computeOtOtpStatic(arr1, arr2, conn=False):
    Npoints = np.shape(arr1)[1]

    Cttp = np.zeros(np.shape(arr1))

    av1 = np.mean(np.mean(arr1, axis=1))
    av2 = np.mean(np.mean(arr2, axis=1))

    for x in range(Npoints):  # x is stores the result
        counter = 0
        tmpArr = arr1[:, (x) % Npoints] * arr2[:, 0]
        for x0 in range(1, Npoints):  # x0 is the origin
            tmpArr += arr1[:, (x0 + x) % Npoints] * arr2[:, x0]
            if conn:
                tmpArr -= av1 * av2
            counter += 1
        Cttp[:, x] = tmpArr / Npoints

    return Cttp


# Compute the correlator in time between arr1 and arr2.
#
# The correlation function is computed over the full range if nTMax is not
# passed as an optional argument. If nTMax is passd then the maximum range of
# the correlator is between 0...nTMax - 1.
#
# The routine will also substract the connected part by setting conn (for
# connected) to True, but this is false by default
#
# The statistical method for computing the average is passed as the third
# arguement.  The default is
#
# lambda x: (np.mean(x), 0)
#
# Examples:
#
# computeOtOtp(d, np.conj(d))
def computeOtOtp(arr1, arr2, statFunc=lambda x: (np.mean(x), 0), conn=False, nTMax=-1, decim=1):

    # Take the full possible range if nTMax is not passed
    Npoints = len(arr1)
    if nTMax == -1:  # nTMax was not passed
        nTMax = Npoints

    # Allocate space for the result
    Cttp = np.zeros(nTMax, dtype=arr1.dtype)
    CttpErr = np.zeros(nTMax, dtype=arr1.dtype)

    CttpRatio = np.zeros(nTMax-1, dtype=type(arr1[0]))
    CttpErrRatio = np.zeros(nTMax-1, dtype=type(arr1[0]))

    # Compute means if necessary for connected part
    if conn:
        av1 = np.mean(arr1)
        av2 = np.mean(arr2)

    # For each time compute the product of a1 and rolled a2
    # an average of this is the correlation ofunction
    for tt in range(nTMax):
        tmp = arr1 * np.roll(arr2, -tt)[0:len(arr1)]
        if conn:
            tmp -= av1 * av2


        Cttp[tt], CttpErr[tt] = statFunc(tmp[0:Npoints - tt:decim])
        

    return Cttp, CttpErr


# Compute the O(t) O(0) correlator over blocks in time.
#
# Inputs:
#
# array1 and array2, as well as the dimension of the correlation array, nTMax.
# steps is the nuber of array entries used to form a block average
#
# Outputs:
#
# A  python array,  containing the O(t) O(0) as an array, for each block.
#
# [ [data], [data], [data], ..... ]
def computeBlockedOtOtp(arr1, arr2, nTMax, steps, decim=1, conn=False):
    res = []
    resRatio = []
    sMax = len(arr1) - steps
    for s in range(0, sMax, steps):
        tmpR, tmpE = computeOtOtp(arr1[s:s + steps],
                                  arr2[s:s + steps],
                                  lambda x: (np.mean(x), 0),
                                  nTMax=nTMax,
                                  decim=decim,
                                  conn=conn)
        res.append(tmpR)

    return np.asarray(res)

def computeParallelBlockedOtOtpHelper(s,
                                      arr1,
                                      arr2,
                                      nTMax,
                                      steps,
                                      decim=1,
                                      conn=False):
    tmpR, tmpE  = computeOtOtp(arr1[s:s + steps],
                              arr2[s:s + steps],
                              lambda x: (np.mean(x), 0),
                              nTMax=nTMax,
                              decim=decim,
                              conn=conn)
    return tmpR


# Compute the O(t) O(0) correlator over blocks in time.
def computeParallelBlockedOtOtp(arr1,
                                arr2,
                                nTMax,
                                steps,
                                decim=1,
                                conn=False,
                                ncpus=2):
    res = []
    sMax = len(arr1) - steps
    blocks = [s for s in range(0, sMax, steps)]
    pool = mp.Pool(ncpus)
    func = functools.partial(computeParallelBlockedOtOtpHelper,
                             arr1=arr1,
                             arr2=arr2,
                             nTMax=nTMax,
                             steps=steps,
                             decim=decim,
                             conn=conn)
    res = pool.map(func, blocks)
    pool.close()
    return np.asarray(res)

# "Slow" fourier transform in time. Used only for the one in time. Can apply a
# filter
def ft_cos(om, dt, arr2, filterFunc=lambda x: 1):
    tmp = arr2[0] * filterFunc(0) * dt
    for i in range(1, len(arr2)):
        t = i * dt
        tmp += 2.0 * arr2[i] * np.cos(om * t) * filterFunc(t) * dt
    return tmp


def ft_sin(om, dt, arr2, filterFunc=lambda x: 1):
    tmp = 0.0  # arr2[0] * filterFunc(0) * dt
    for i in range(1, len(arr2)):
        t = i * dt
        tmp += 2.0 * 1j * arr2[i] * np.sin(om * t) * filterFunc(t) * dt
    return tmp


# transform a whole array to fourier space, with optimal sampling and up to
# oMax
def toFourier(arr, dt, omMax, myTmax=-1, filterFunc=lambda x: 1, symFactor=1):
    tMax = len(arr) * dt if myTmax == -1 else myTmax
    omegaIR = 2.0 * np.pi / tMax
    oms = np.arange(0, omMax, omegaIR)
    res = []
    for om in oms:
        if symFactor == 1:
            res.append(ft_cos(om, dt, arr, filterFunc))
        elif symFactor == -1:
            res.append(ft_sin(om, dt, arr, filterFunc))

    return (oms, np.asarray(res))


# Returns a symmetrized version of the original array:
#
# This assume first entry (arr[0]) is at zero and we don't want to dublicate
# it, i.e. if the array has length N, the returned array has length 2N-1
def symmetrizeArray(arr, sym=1):
    res = np.zeros(2 * len(arr) - 1, arr.dtype)
    for i in range(1, len(arr)):
        res[i - 1] = sym * arr[-i]
    for i in range(len(arr) - 1, len(res)):
        res[i] = arr[i - len(arr) + 1]
    return res


def quadratic_spline_roots(spl):
    roots = []
    knots = spl.get_knots()
    for a, b in zip(knots[:-1], knots[1:]):
        u, v, w = spl(a), spl((a + b) / 2), spl(b)
        t = np.roots([u + w - 2 * v, w - u, 2 * v])
        t = t[np.isreal(t) & (np.abs(t) <= 1)]
        roots.extend(t * (b - a) / 2 + (b + a) / 2)
    return np.array(roots)


def findMin(x_axis, y_axis):  # https://stackoverflow.com/a/50373953
    f = InterpolatedUnivariateSpline(x_axis, y_axis, k=3)
    cr_pts = quadratic_spline_roots(f.derivative())
    cr_pts = np.append(
        cr_pts,
        (x_axis[0], x_axis[-1]))  # also check the endpoints of the interval
    cr_vals = f(cr_pts)
    return cr_pts[np.argmin(cr_vals)]


def findMax(x_axis, y_axis):  # https://stackoverflow.com/a/50373953
    f = InterpolatedUnivariateSpline(x_axis, y_axis, k=3)
    cr_pts = quadratic_spline_roots(f.derivative())
    cr_pts = np.append(
        cr_pts,
        (x_axis[0], x_axis[-1]))  # also check the endpoints of the interval
    cr_vals = f(cr_pts)
    return cr_pts[np.argmax(cr_vals)]


# Compute the Fourier transform over blocks of data.
def toFourierBlocked(arr,
                     dt,
                     omMax,
                     errFunc,
                     myTmax=-1,
                     filterFunc=lambda x: 1,
                     symFactor=1):
    length = len(arr)
    res = []
    resMaxOm = []
    for i in range(length):
        oms, tmp = toFourier(arr[i], dt, omMax, myTmax, filterFunc, symFactor)
        res.append(tmp)
        # resMaxOm.append(oms[np.argmax(tmp)])
        resMaxOm.append(findMax(oms, tmp))

    f, ferr = errFunc(np.asarray(res))

    return (symmetrizeArray(oms, -1),
            StatResult(
                (symmetrizeArray(f,
                                 symFactor), symmetrizeArray(ferr))), resMaxOm)


# For testing purpose only
def manualFourier(arr):
    Nx = len(arr)
    ks = [2 * np.pi * i / float(Nx) for i in range(Nx)]
    arrF = np.zeros(Nx, dtype=complex)
    for k in range(Nx):
        for x in range(Nx):
            arrF[k] += arr[x] * np.exp(-1j * x * ks[k])

    return arrF


# static correlator in x
def twoPtAv(arr):
    N = len(arr)
    res = np.zeros(N)
    for i in range(N):
        res[i] = 0
        for j in range(N):
            res[i] += arr[j] * arr[(i + j) % N]

    return np.asarray(res)

################################################################################
# Work in progress, store the results and call the corersponding measurements
# in a "user friendly" way.


class ConfResults:

    def __init__(self,
                 fn,
                 thTime,
                 dt,
                 data_format="new",
                 processedDir='./',
                 plotDir='./',
                 loadFourierBool=False):
        # filename
        self.fn = fn
        # tag
        self.tag = self.fn.split('.')[-2].split('/')[-1]
        # data directory
        self.dataDir = self.fn[:self.fn.rfind('/')] + "/"
        # thermalization time
        self.thTime = thTime
        # dt
        self.dt = dt
        # bool to say if fourier transform of wall is saved
        self.loadFourierBool = loadFourierBool

        self.processedDir = processedDir
        self.plotDir = plotDir

        # old, semi_old, new
        self.data_format = data_format

        # A dictionary storing the names of the fields, and the associated
        # column number akeys["phi0"] is column 0 and so on.
        self.akeys = dict()
        self.wkeys = dict()

        self.akeys["phi0"] = 0
        # WARNING: dsigma is a bad name, should really be dphi0.
        self.akeys["dsigma"] = 0
        self.akeys["phi1"] = 1
        self.akeys["phi2"] = 2
        self.akeys["phi3"] = 3
        self.akeys["dphi1"] = 1
        self.akeys["dphi2"] = 2
        self.akeys["dphi3"] = 3
        self.akeys["A1"] = 4
        self.akeys["A2"] = 5
        self.akeys["A3"] = 6
        self.akeys["V1"] = 7
        self.akeys["V2"] = 8
        self.akeys["V3"] = 9
        self.akeys["phiNorm"] = 10
        self.akeys["M2"] = 11
        self.akeys["M22"] = 12

        if self.data_format == "old":
            self.wkeys["phi0"] = "wallX_phi_0"
            self.wkeys["dsigma"] = "wallX_phi_0"
            self.wkeys["phi1"] = "wallX_phi_1"
            self.wkeys["phi2"] = "wallX_phi_2"
            self.wkeys["phi3"] = "wallX_phi_3"
            self.wkeys["dphi1"] = "wallX_phi_1"
            self.wkeys["dphi2"] = "wallX_phi_2"
            self.wkeys["dphi3"] = "wallX_phi_3"
            self.wkeys["A1"] = "wallX_phi_4"
            self.wkeys["A2"] = "wallX_phi_5"
            self.wkeys["A3"] = "wallX_phi_6"
            self.wkeys["V1"] = "wallX_phi_7"
            self.wkeys["V2"] = "wallX_phi_8"
            self.wkeys["V3"] = "wallX_phi_9"
        else:
            self.wkeys = self.akeys.copy()

        # The a list of directions
        if data_format == "old" or data_format == "semi_old":
            self.directions = ["X"]
        else:
            self.directions = ["X", "Y", "Z"]

        self.scalarKeys = ["phi0", "dsigma", "phiNorm"]

        # The results are stored in dictionnary, to have a generic way to code
        # between the different fields.
        self.av = dict()
        self.meanValues = dict()
        self.variances = dict()
        self.o2 = dict()

        self.wall = dict()
        self.wallF = dict()
        self.wallCorr = dict()
        self.wallCorr_raw = dict()
        self.wallCorrEffMasses_blocks= dict()
        self.wallCorrEffMasses= dict()
        self.wallFProp = dict()
        for d in self.directions:
            self.wall[d] = dict()
            self.wallF[d] = dict()
            # self.wallCorr[d] = dict()

        self.OtOttp_blocks = dict()
        self.OtOttpRatio_blocks = dict()

        self.OtOttp = dict()
        self.OtOttpRatio = dict()
        self.OtOttp_time = dict()

        # self.OtOttpFourier_blocks = dict()
        self.OtOttpFourier = dict()
        self.OtOttpFourier_oms = dict()
        self.OtOttpFourier_omspeak_blocks = dict()

        # self.momenta_3d we also define the momenta 3d below
        # self.xs same

    def loadAv(self, key):
        if not key in self.av.keys():
            r = h5py.File(self.fn, 'r')
            self.av[key] = np.asarray(r["phi"][self.thTime:, self.akeys[key]])
            if key == "dsigma":
                bar = np.mean(self.av[key], axis=0)
                self.av[key] -= bar

    def computeFromMean(self, key, errFunc, decim=1, func=lambda x: x):
        if key == "phiH0":
            self.loadAv("phi0")
            av = []
            # av = deepcopy(self.av["phi0"])
            for k in ["phi0", "phi1", "phi2", "phi3"]:
                self.loadAv(k)
                av.append(self.av[k])
            return StatResult(errFunc(np.asarray(av).flatten(), func=func))
        else:
            self.loadAv(key)
            return StatResult(errFunc(self.av[key][:None:decim], func=func))

    def computeMean(self, key, errFunc, decim=1):
        self.meanValues[key] = self.computeFromMean(key,
                                                    errFunc,
                                                    decim=decim,
                                                    func=lambda x: x)

    def computeVar(self, key, errFunc, decim=1):
        self.variances[key] = self.computeFromMean(
            key, errFunc, decim=decim, func=lambda x: x**2 - np.mean(x)**2)

    def computeO2(self, key, errFunc, decim=1):
        self.o2[key] = self.computeFromMean(key,
                                            errFunc,
                                            decim=decim,
                                            func=lambda x: x**2)

    def loadWall(self, key, direc=None):
        if direc == None:
            direc = "X"
        r = h5py.File(self.fn, 'r')
        try:
            if not (direc == "X" or (self.data_format != "old" and
                                     self.data_format != "semi_old")):
                raise
            else:
                if self.data_format == "old":
                    self.wall[direc][key] = np.asarray(
                        r[self.wkeys[key]])[self.thTime:]
                elif self.data_format == "semi_old":
                    self.wall[direc][key] = np.asarray(
                        r["corrx"])[self.thTime:, self.wkeys[key], :]
                else:
                    if not key in self.wall[direc].keys():
                        self.wall[direc][key] = np.asarray(r["wall{}".format(
                            direc.lower())])[self.thTime:, self.wkeys[key], :]
                if key == "dsigma" or key == "dphi1" or key == "dphi2" or key == "dphi3":
                    self.wall[direc][key] -= np.mean(
                        np.mean(self.wall[direc][key], axis=0))

        except:
            print(
                "You tried to load a wall in a direction that does not exist in your data format."
            )

    # Load C generated x2f file, not in use
    def loadWallFourier_derek(self, key, direc=None):
        if direc == None:
            direc = "X"
        r = h5py.File(self.fn.split('.')[-2]+"_out.h5", 'r')
        # try:
        # if (self.data_format == "old" or
        # self.data_format != "semi_old")):
        # raise
        # else:
        if not key in self.wallF[direc].keys():
            self.wallF[direc][key] = np.asarray(r["wall{}_k".format(
                direc.lower())])[self.thTime:, self.wkeys[key], :, 0] + 1j * np.asarray(r["wall{}_k".format(
                    direc.lower())])[self.thTime:, self.wkeys[key], :, 1]
            # if key == "dsigma" or key == "dphi1" or key == "dphi2" or key == "dphi3":
            # raise

        # except:
        #    print(
        #        "You tried to load a wallF in a direction that does not exist in your data format. You may also have tried to load the fourier transform of a connected key"
        #    )

    # Load python generated fourier file.
    def loadWallFourier(self, key, direc=None):
        if direc == None:
            direc = "X"
        f = h5py.File(self.dataDir + "/" + self.tag + "_fourier.h5", 'r')
        if not key in self.wallF[direc].keys():
            self.wallF[direc][key] = np.asarray(f[direc][key])
        f.close()

    # TODO: Change the way we handle the average.
    def readAv(self):
        r = h5py.File(self.fn, 'r')
        self.phi = np.asarray(r["phi"])[self.thTime:]
        self.phiNorm = self.phi[:, 4]
        self.mag = self.phi[:, -2]
        self.mag2 = self.phi[:, -1]

    def computeMag(self, direc=0):
        self.readAv()
        self.mag, self.magErr = bootstrap(self.phi[:, direc], 100)

    # Compute the fourier transform of a given wall, specified by the
    # appropriate key.

    def computeWallFourier(self, key, direc, decim=1):
        if not direc in self.wallF.keys() or not key in self.wallF[direc].keys(
        ):
            if(self.loadFourierBool == True):
                print("Loading Fourier, no decim applied.")
                self.loadWallFourier(key, direc)
            else:
                print("Computing Fourier")
                self.loadWall(key, direc)
                Nx = len(self.wall[direc][key][0])
                ts = range(0, len(self.wall[direc][key]), decim)
                self.wallF[direc][key] = np.zeros(
                    [len(ts), len(self.wall[direc][key][0, :])], dtype=complex)
                c = 0
                for t in ts:
                    self.wallF[direc][key][c] = np.fft.fft(
                        self.wall[direc][key][t]) / Nx
                    c += 1
                self.momenta_3d = np.asarray(
                    [2.0 * np.pi / float(Nx) * n for n in range(Nx)])

    # Save fourier to hdf5

    def saveWallFourier(self, key, direc):
        f = h5py.File(self.processedDir + "/" + self.tag + "_fourier.h5", 'a')
        if not direc in f.keys():
            grp = f.create_group(direc)
        else:
            grp = f[direc]

        if not key in grp.keys():
            grp.create_dataset(key, data=self.wallF[direc][key])

    # Computes the time correlator of a given key.
    #
    # TODO: For now, compute only from one direction in fourier for all modes.
    #      Need to include other direction for non-zero k for better stat.
    # TODO: Mixed correlator and same correlator can be coded as one.
    def computeOtOtpBlocked(self,
                            key,
                            tMax,
                            nBlocks,
                            errFunc,
                            decim=1,
                            momNum=0,
                            conn=False,
                            parallel=False,
                            redo=False,
                            tEffMass = 50
                            ):

        # Number of bins
        nTMax = int(np.floor(tMax / self.dt))

        # blockSize = int(np.floor(blockSizeT / self.dt))
        if momNum != 0:
            directions = self.directions
            bkey = key + "kk{}".format(momNum)
        else:  # 0 momentum case
            directions = ["X"]
            bkey = key

        if redo or not bkey in self.OtOttp_blocks.keys():
            # Mixed correlator case, does not compute in  Y Z (because I am lazy)
            if '_' in key:
                key1, key2 = key.split('_')
                keys1 = [key1 + str(i) for i in [1, 2, 3]]
                keys2 = [key2 + str(i) for i in [1, 2, 3]]
                for k in keys1:
                    self.computeWallFourier(k, "X")
                for k in keys2:
                    self.computeWallFourier(k, "X")
                # Get volume from the length of the vector and the maximum
                # time as an integer.
                blockSize = int(len(self.wallF["X"][keys1[0]][:, 0]) / nBlocks)
                V = float(len(self.wallF["X"][keys1[0]][0, :])**3)
                # Compute the correlator per block.  First for a given flavor
                # index (so we don't need to worry about the shape of the
                # array).
                k1 = keys1[0]
                k2 = keys2[0]
                if (parallel == False):
                    self.OtOttp_blocks[bkey] =  computeBlockedOtOtp(
                        self.wallF["X"][k1][:, momNum],
                        np.conj(self.wallF["X"][k2][:, momNum]),
                        nTMax,
                        blockSize,
                        decim=decim,
                        conn=conn)
                else:
                    self.OtOttp_blocks[bkey]=  computeParallelBlockedOtOtp(
                        self.wallF["X"][k1][:, momNum],
                        np.conj(self.wallF["X"][k2][:, momNum]),
                        nTMax,
                        blockSize,
                        decim=decim,
                        conn=conn,
                        ncpus=mp.cpu_count())

                # Then for the remaining two directions.
                for i in range(1, len(keys1)):
                    # av += self.wallF["X"][k][:,momNum]
                    k1 = keys1[i]
                    k2 = keys2[i]
                    if (parallel == False):
                         self.OtOttp_blocks[bkey]  += computeBlockedOtOtp(
                            self.wallF["X"][k1][:, momNum],
                            np.conj(self.wallF["X"][k2][:, momNum]),
                            nTMax,
                            blockSize,
                            decim=decim,
                            conn=conn)
                    else:
                         self.OtOttp_blocks[bkey]  += computeParallelBlockedOtOtp(
                                self.wallF["X"][k1][:, momNum],
                                np.conj(self.wallF["X"][k2][:, momNum]),
                                nTMax,
                                blockSize,
                                decim=decim,
                                conn=conn,
                                ncpus=mp.cpu_count())
                self.OtOttp_blocks[bkey] *= (V / 3.0)
            elif not key in self.scalarKeys:
                isotropic = False
                if key == "phiRestored":
                    key = "phi"
                    isotropic == True
                keys = [key + str(i) for i in [1, 2, 3]] if not isotropic else [
                    key + str(i) for i in [0, 1, 2, 3]]
                # Compute the spatial fourier transform
                for k in keys:
                    for d in directions:
                        self.computeWallFourier(k, d)

                # Get volume from the length of the vector and the maximum
                # time as an integer.
                blockSize = int(len(self.wallF["X"][keys[0]][:, 0]) / nBlocks)
                V = float(len(self.wallF["X"][keys[0]][0, :])**3)

                # Compute the correlator per block. First for a given flavor
                # index and direction (so we don't need to worry about the
                # shape of the array).
                k = keys[0]
                if (parallel == False):
                    self.OtOttp_blocks[bkey] =  computeBlockedOtOtp(
                        self.wallF["X"][k][:, momNum],
                        np.conj(self.wallF["X"][k][:, momNum]),
                        nTMax,
                        blockSize,
                        decim=decim,
                        conn=conn)
                else:
                    self.OtOttp_blocks[bkey]=  computeParallelBlockedOtOtp(
                        self.wallF["X"][k][:, momNum],
                        np.conj(self.wallF["X"][k][:, momNum]),
                        nTMax,
                        blockSize,
                        decim=decim,
                        conn=conn,
                        ncpus=mp.cpu_count())

                # Then for the remaining cases.
                for k in keys:
                    # av += self.wallF["X"][k][:,momNum]
                    for d in directions:
                        if k != keys[0] or d != "X":  # already computed
                            if (parallel == False):
                                 self.OtOttp_blocks[bkey]  +=  computeBlockedOtOtp(
                                        self.wallF[d][k][:, momNum],
                                        np.conj(self.wallF[d][k][:, momNum]),
                                        nTMax,
                                        blockSize,
                                        decim=decim,
                                        conn=conn)
                            else:
                                 self.OtOttp_blocks[bkey]  +=  computeParallelBlockedOtOtp(
                                        self.wallF[d][k][:, momNum],
                                        np.conj(self.wallF[d][k][:, momNum]),
                                        nTMax,
                                        blockSize,
                                        decim=decim,
                                        conn=conn,
                                        ncpus=mp.cpu_count())
                self.OtOttp_blocks[bkey] *= (V / float(len(keys) * len(directions)))
            else:
                # For scalar quantities we have only one
                for d in directions:
                    self.computeWallFourier(key, d)
                blockSize = int(len(self.wallF["X"][key][:, 0]) / nBlocks)
                V = float(len(self.wallF["X"][key][0, :])**3)

                if (parallel == False):
                    self.OtOttp_blocks[bkey] =  computeParallelBlockedOtOtp(
                        self.wallF["X"][key][:, momNum],
                        np.conj(self.wallF["X"][key][:, momNum]),
                        nTMax,
                        blockSize,
                        decim=decim,
                        conn=conn)
                else:
                    self.OtOttp_blocks[bkey] =  computeParallelBlockedOtOtp(
                        self.wallF["X"][key][:, momNum],
                        np.conj(self.wallF["X"][key][:, momNum]),
                        nTMax,
                        blockSize,
                        decim=decim,
                        conn=conn,
                        ncpus=mp.cpu_count())

                for i in range(1, len(directions)):
                    d = self.directions[i]
                    if (parallel == False):
                         self.OtOttp_blocks[bkey]  += V * computeParallelBlockedOtOtp(
                            self.wallF[d][key][:, momNum],
                            np.conj(self.wallF[d][key][:, momNum]),
                            nTMax,
                            blockSize,
                            decim=decim,
                            conn=conn)
                    else:
                         self.OtOttp_blocks[bkey]  +=  computeParallelBlockedOtOtp(
                                self.wallF[d][key][:, momNum],
                                np.conj(self.wallF[d][key][:, momNum]),
                                nTMax,
                                blockSize,
                                decim=decim,
                                conn=conn,
                                ncpus=mp.cpu_count())

                self.OtOttp_blocks[bkey] *= (V / float(len(directions)))

        # We save the average wall. To get the correlator, need to multiply by
        # a volume factor.
        self.OtOttp[bkey] = StatResult(errFunc(self.OtOttp_blocks[bkey]))

        self.OtOttpRatio_blocks[bkey] = []
        for otottp in self.OtOttp_blocks[bkey]:
            self.OtOttpRatio_blocks[bkey].append(tEffMass * self.dt / np.log(otottp[:-tEffMass]/np.roll(otottp,-tEffMass)[:-tEffMass]))


        self.OtOttpRatio[bkey] = StatResult(errFunc(self.OtOttpRatio_blocks[bkey]))

        self.OtOttp_time[bkey] = np.asarray(
            [self.dt * float(i) for i in range(len(self.OtOttp[bkey].mean))])

    # Compute the statistical correlator of a given key.
    def computeStatisticalCor(self,
                              key,
                              omMax,
                              errFunc,
                              myTmax=-1,
                              filterFunc=lambda x: 1):
        try:
            if not (key in self.OtOttp_blocks.keys()):
                raise
            else:
                symFactor = 1
                if key == "A_phi":
                    symFactor = -1
                self.OtOttpFourier_oms[key], self.OtOttpFourier[
                    key], self.OtOttpFourier_omspeak_blocks[
                        key] = toFourierBlocked(self.OtOttp_blocks[key],
                                                self.dt,
                                                omMax,
                                                errFunc,
                                                myTmax=myTmax,
                                                filterFunc=filterFunc,
                                                symFactor=symFactor)
        except:
            print(
                "the two point function in time of {} was not computed. Need to call computeOtOtpBlocked first."
                .format(key))

    def computeFourierPropagator(self, key, decim, errFunc):
        # if key != "phi0" and key != "dsigma":
        if not key in self.scalarKeys:
            keys = [key + str(i) for i in [1, 2, 3]]
            for k in keys:
                for d in self.directions:
                    self.computeWallFourier(k, d, decim)
            tmp = np.zeros(np.shape(self.wallF["X"][keys[0]]),
                           dtype=self.wallF["X"][keys[0]].dtype)
            for k in keys:
                for d in self.directions:
                    tmp += self.wallF[d][k] * np.conj(self.wallF[d][k])
            tmp /= (float(len(keys)) * float(len(self.directions)))
        else:
            for d in self.directions:
                self.computeWallFourier(key, d, decim)
            tmp = np.zeros(np.shape(self.wallF["X"][key]),
                           dtype=self.wallF["X"][key].dtype)
            for d in self.directions:
                tmp += self.wallF[d][key] * np.conj(self.wallF[d][key])
            tmp /= float(len(self.directions))

        self.wallFProp[key] = StatResult(errFunc(tmp))

    def computePropagator(self, key, errFunc, decim=1, alreadyLoaded=False):
        if not alreadyLoaded:
            print("hi")
            ckey = key
            # if key != "phi0" and key != "dsigma":
            if not key in self.scalarKeys:
                isotropic = False
                if key == "phiRestored":
                    key = "phi"
                    isotropic == True
                keys = [key + str(i) for i in [1, 2, 3]] if not isotropic else [
                    key + str(i) for i in [0, 1, 2, 3]]
                for k in keys:
                    for d in self.directions:
                        print(d)
                        self.loadWall(k, d)
                V = float(len(self.wall["X"][keys[0]][0, :])**3)
                tmp = []
                for k in keys:
                    for d in self.directions:
                        tmp.append(
                            computeOtOtpStatic(self.wall[d][k],
                                               self.wall[d][k]))
                # tmp /= (float(len(keys)) * float(len(self.directions)))
            else:
                for d in self.directions:
                    self.loadWall(key, d)
                V = float(len(self.wall["X"][key][0, :])**3)
                tmp = []
                for d in self.directions:
                    tmp.append(
                        computeOtOtpStatic(self.wall[d][key][0:None:decim],
                                           self.wall[d][key][0:None:decim]))
                # tmp /= float(len(self.directions))

            tmp = np.asarray(tmp)
            tmp = np.reshape(tmp, (tmp.shape[0] * tmp.shape[1], tmp.shape[2]))

            self.wallCorr_raw[ckey] = tmp

        self.wallCorr[ckey] = StatResult(
            errFunc(self.wallCorr_raw[ckey][:None:decim, :]))
        self.xs = np.arange(len(self.wallCorr[ckey].mean))

    def computeStaticEffMasses(self, key, nBlocks, tSep, errFunc, alreadyLoaded=True):
        self.wallCorrEffMasses_blocks[key] = []
        if not alreadyLoaded:
            return
        else:
            shs = np.shape(self.wallCorr_raw[key])
            blockSize = int(shs[0] / nBlocks)
            for i in range(nBlocks):
                blocked = np.mean(self.wallCorr_raw[key][i:i+blockSize], axis = 0)
                self.wallCorrEffMasses_blocks[key].append(np.log(blocked[:-tSep]/np.roll(blocked,-tSep)[:-tSep]))


        self.wallCorrEffMasses[key] = StatResult(errFunc(self.wallCorrEffMasses_blocks[key]))

    def getObs(self, obs, key, withX=True, withY=True):
        if obs == "OtOttp":
            x = self.OtOttp_time[key]
            y = self.OtOttp[key]
        if obs == "OtOttpRatio":
            x = np.arange(len(self.OtOttpRatio[key].mean))
            y = self.OtOttpRatio[key]
        elif obs == "OtOttpFourier":
            x = self.OtOttpFourier_oms[key]
            y = self.OtOttpFourier[key]
        elif obs == "propagatorF":
            x = self.momenta_3d
            y = self.wallFProp[key]
        elif obs == "propagator":
            x = self.xs
            y = self.wallCorr[key]
        if withX == True and withY == True:
            return (x, y)
        elif withY == True and withX == False:
            return y
        else:
            return x

    def getDataFn(self, obs, key, direc, tag=''):
        return direc + "/" + self.tag + "_" + obs + "_" + key + tag + ".txt"

    def save(self,
             obs,
             key,
             direc=None,
             fmt="python",
             decim=1,
             tMin=0,
             tMax=None,
             xfact=1,
             yfact=1,
             tag=''):
        if direc == None:
            if fmt == "python":
                direc = self.processedDir
            elif fmt == "gnuplot" or fmt == "latex":
                direc = self.plotDir

        fn = self.getDataFn(obs, key, direc, tag)

        if not obs == "OtOttp_blocks" and not obs == "propagator_raw" and not obs == "wallF" and not obs == "OtOttpRatio_blocks" :
            x, y = self.getObs(obs, key)
            y.save_to_txt(fn,
                          x=xfact * x,
                          fmt=fmt,
                          decim=decim,
                          tMin=tMin,
                          tMax=tMax,
                          yfact=yfact)
        else:
            if fmt == "python":
                file = open(fn, 'wb')
                if obs == "OtOttp_blocks":
                    pickle.dump(self.OtOttp_blocks[key], file)
                elif obs == "OtOttpRatio_blocks":
                    pickle.dump(self.OtOttpRatio_blocks[key], file)
                elif obs == "propagator_raw":
                    pickle.dump(self.wallCorr_raw[key], file)
                elif obs == "wallF":
                    pickle.dump(self.wallF, file)
                file.close()
                if decim > 1:
                    print(
                        "Warning, decim was ignored for the blocked correlator."
                    )
                elif xfact != 1 or yfact != 1:
                    print("Can't rescale blocks to save, ignore.")
            else:
                print("Blocks can only be saved in fmt = python.")

    def load(self, obs, key, direc=None):
        if direc == None:
            direc = self.processedDir
        fn = self.getDataFn(obs, key, direc)
        if obs == "OtOttp":
            self.OtOttp_time[key], self.OtOttp[key] = load_StatResult_from_txt(
                fn, withX=True)
        if obs == "OtOttpRatio":
            x,self.OtOttpRatio[key] = load_StatResult_from_txt(
                fn, withX=True)
        elif obs == "OtOttpFourier":
            self.OtOttpFourier_oms[key], self.OtOttpFourier[
                key] = load_StatResult_from_txt(fn, withX=True)
        elif obs == "propagator":
            self.xs, self.wallCorr[key] = load_StatResult_from_txt(fn,
                                                                   withX=True)
        elif obs == "propagatorF":
            self.momenta_3d, self.wallFProp[key] = load_StatResult_from_txt(
                fn, withX=True)
        elif obs == "OtOttp_blocks":
            file = open(fn, 'rb')
            self.OtOttp_blocks[key] = pickle.load(file)
            file.close()
        elif obs == "propagator_raw":
            file = open(fn, 'rb')
            self.wallCorr_raw[key] = pickle.load(file)
            file.close()
        elif obs == "wallF":
            file = open(fn, 'rb')
            self.wallF = pickle.load(file)
            file.close()


class Plotter:

    def __init__(self):
        style.use('seaborn-paper')
        return

    def plot(self,
             confRes,
             obs,
             key,
             xfact=1.0,
             yfact=1.0,
             imOrReal=np.real,
             band=True,
             tMin=0,
             tMax=None,
             prune=1):
        x, y = confRes.getObs(obs, key, withX=True)
        self.errorplot(xfact * x,
                       y,
                       yfact=yfact,
                       imOrReal=imOrReal,
                       band=band,
                       tMin=tMin,
                       tMax=tMax)

    def errorplot(self,
                  x,
                  statobj,
                  yfact,
                  imOrReal=np.real,
                  band=True,
                  tMin=0,
                  tMax=None,
                  prune=1):
        if not band:
            plt.errorbar(x[tMin:tMax:prune],
                         yfact * statobj.mean[tMin:tMax:prune],
                         yfact * statobj.err[tMin:tMax:prune])
        else:
            plt.plot(x[tMin:tMax:prune],
                     yfact * imOrReal(statobj.mean[tMin:tMax:prune]), '-')
            plt.fill_between(
                x[tMin:tMax:prune],
                (np.asarray(yfact * imOrReal(statobj.mean[tMin:tMax:prune])) -
                 np.asarray(yfact * imOrReal(statobj.err[tMin:tMax:prune]))),
                (np.asarray(yfact * imOrReal(statobj.mean[tMin:tMax:prune])) +
                 np.asarray(yfact * imOrReal(statobj.err[tMin:tMax:prune]))),
                linewidth=0,
                zorder=1)


''' Deprecated, keep for the idea. Would be better o use the block structure above though.
class EnsembleResults:
    def __init__(self,fn, nEnd, nStart=0):
        self.fnList = [fn +"_"+str(i)+".h5" for i in range(nStart,nEnd+1)]
        self.nConf = nEnd-nStart+1


    def computeWallXFourier0(self, nTherm, nBoot, decim):
        resTot = []
        for fn in self.fnList:
            res = ConfResults(fn, nTherm,decim)
            print(fn, end="\r")
            res.computeWallXFourier0()
            resTot.append(res.wallXphi0F)
        self.wallXFourier0, self.wallXFourier0Err  = bootstrap(resTot, nBoot)
    def computeWallXFourierSquare0(self, nTherm, nBoot, decim):
        resTot = []
        for fn in self.fnList:
            res = ConfResults(fn, nTherm,decim)
            print(fn, end="\r")
            res.computeWallXFourier0()
            resTot.append(res.wallXphi0F * np.conj(res.wallXphi0F))
        self.wallXFourierSquare0, self.wallXFourierSquare0Err  = jackknife(resTot, nBoot)

    def corFourierPhi(self, nTherm, nBoot, decim):
        corrTot = []
        for fn in self.fnList:
            res = ConfResults(fn, nTherm,decim)
            print(fn)
            res.wallXFourierNorm()
            Nx = len(res.wallXphiNormF[0])
            corrTot.append([res.computeCttp(np.conj(res.wallXphiNormF[:,k]),res.wallXphiNormF[:,k]) for k in range(Nx)])
            print(len(corrTot[-1][0]))

        self.corFourierPhiMean, self.corFourierPhiErr  = bootstrap(corrTot, nBoot)
    def save(self, fn):
            np.savetxt(fn,[self.corFourierPhiMean, self.corFourierPhiErr])
'''
