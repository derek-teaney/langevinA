from measurements import *
import time

if __name__ == '__main__':
    
    #L = ["16", "32", "48", "64"]
    #ms = ["473366", "470052","466000", "468000", "476000", "478000"]
    #ms = ["466000"]
    L = ["096"]
    #ms = ["471000", "472800", "473366", "470052", "476000", "466000", "474500", "472000"]
    #ms = ["472800", "473366", "472000"]
    #ms = ["470052", "468000", "476000", "478000","472800", "473366", "472000", "471000"]
    ms = ["463000"]
    fns = dict()
    for m in ms:
        fns[m] = dict()
        for l in L:
            tmp = "zeroHlongDiffuse_N" + l + "_m-0" + m + "_h000000_c00500"
            fns[m][l] = "../run_H0/" + tmp +"/" + tmp + ".h5"

    dt = 0.72
    dec = 50

    #decs = dict()

    #decs["470052"] = 100
    #decs["471000"] = 180
    #decs["472000"] = 200
    #decs["472800"] = 220
    #decs["473366"] = 240
    #decs["474500"] = 260
    #decs["466000"] = 280
    #decs["476000"] = 300

    
    for m in ms:
        for l in L:
            fn = fns[m][l]
            print(fn)
    
            t0 = time.time()
     #       dec = decs[m]

            for k in ["phiRestored"]:
                 data = ConfResults(fn = fn,thTime=2500,dt=dt,  data_format="new", loadFourierBool = True)
                 data.computeOtOtpBlocked(k, momNum = 0, tMax = 1000.0, nBlocks = 5,  decim = dec, errFunc = lambda x : (np.mean(x, axis = 0), np.std(x, axis = 0)), parallel=True,tEffMass = 100)
                 #data.computePropagator(k,errFunc = lambda x: bootstrap(x,100))

                 data.save("OtOttp",k)    
                 data.save("OtOttp_blocks",k)    
                 data.save("OtOttpRatio", k)
                 data.save("OtOttpRatio_blocks", k)

                 #data.save("propagator", k)
                 #data.save("propagator_raw", k)
