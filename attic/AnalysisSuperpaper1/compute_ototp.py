from measurements import *
import time

if __name__ == '__main__':
    hs = ["003000"]
    fn1 = "../data/zminus2_N080_m-0501265_h003684_c00500.h5"
    fn2 = "../data/zplus_N080_m-0445648_h003684_c00500.h5"
    fn3 = "../data/zcritical_N080_m-0481100_h002000_c00500.h5"
    fn4 = "../data/zcritical_N080_m-0481100_h004000_c00500.h5"
    fn5 = "../data/zcritical_N080_m-0481100_h006000_c00500.h5"
    fn6 = "../data/zcritical_N080_m-0481100_h010000_c00500.h5"
    fn7 = "../data/zpseudocritical_N080_m-0470052_h003000_c00500.h5"
    fn8 = "../data/z1p214_N080_m-0471156_h003000_c00500.h5"
    fn9 = "../data/z1p2815_N080_m-0470604_h003000_c00500.h5"
    fn10 = "../data/z1p079_N080_m-0472261_h003000_c00500.h5"
    fn11 = "../data/z0p944_N080_m-0473366_h003000_c00500.h5"
    fn12 = "../data/zm4p75_N080_m-0520000_h003000_c00500.h5"
    fn13 = "../data/z1p753_N080_m-0466737_h003000_c00500.h5"
    fn14 = "../data/z2p023_N080_m-0464527_h003000_c00500.h5"
    fn15 = "../data/zminus2_N080_m-0499128_h003000_c00500.h5"
    fn16 = "../data/zplus_N080_m-0449406_h003000_c00500.h5"
    fn17 = "../data/zcritical_N080_m-0481100_h003000_c00500.h5"
    fn18 = "../data/pionrun8_N080_m-0501265_h003000_c00500.h5"
    fn19 = "../data/pionrun9_N080_m-0501265_h003000_c00500.h5"
    fns = [fn1, fn2,fn3,fn4,fn5,fn6,fn7,fn8,fn9,fn10, fn11, fn12]
    fns = [fn1, fn15]
    fns = [fn15]
    fns = [fn18,fn19]
    #fns = [fn1,fn2, fn7, fn8, fn9, fn10, fn11, fn12, fn13, fn14]

    dt = 0.72
    #dec = 50
    dec = 50
    
    for fn in fns:
        print(fn)
    
        t0 = time.time()

        #for k in ["phi"]:
        for k in ["dsigma", "phi", "A", "V"]:
    	     data = ConfResults(fn = fn,thTime=10000,dt=dt,  data_format="semi_old")
    	     data.computeOtOtpBlocked(k, momNum = 0, tMax = 5000.0, nBlocks = 5,  decim = dec, errFunc = lambda x : (np.mean(x, axis = 0), np.std(x, axis = 0)), parallel=True)
    	     data.computePropagator(k,errFunc = lambda x: blocking(x,5))
    	     #data.computeStatisticalCor(k, omMax=0.1, errFunc=lambda x: (np.mean(x, axis = 0), np.std(np.real(x), axis = 0) + 1j*np.std(np.imag(x), axis = 0)), filterFunc=lambda x : np.exp(-x / 2000.0))
    	     #data.computeFourierPropagator(k, decim=dec, errFunc = lambda x: blocking(x,5))
    	     #print(time.time() - t0)

    	     data.save("propagator", k)
    	     data.save("propagator_raw", k)
    	     data.save("OtOttp",k)    
    	     data.save("OtOttp_blocks",k)    
    	     #data.save("OtOttpFourier", k)
    	     #data.save("propagator", k)

#    for fn in [fn16]:
#        print(fn)
    
#        t0 = time.time()

#        for k in ["A"]:
        #for k in ["dsigma", "phi", "A"]:
#    	     data = ConfResults(fn = fn,thTime=10000,dt=dt,  data_format="new")
#    	     data.computeOtOtpBlocked(k, momNum = 0, tMax = 20000.0, nBlocks = 5,  decim = dec, errFunc = lambda x : (np.mean(x, axis = 0), np.std(x, axis = 0)), parallel=True)
    
    	     #data.computeStatisticalCor(k, omMax=0.1, errFunc=lambda x: (np.mean(x, axis = 0), np.std(x, axis = 0)), filterFunc=lambda x : np.exp(-x / 2000.0))
    	     #data.computeFourierPropagator(k, decim=dec, errFunc = lambda x: blocking(x,5))
#    	     print(time.time() - t0)

    
#    	     data.save("OtOttp",k)    
#    	     data.save("OtOttp_blocks",k)    
   	     #data.save("OtOttpFourier", k)
    	     #data.save("propagator", k)
