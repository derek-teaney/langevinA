from measurements import *
import time

if __name__ == '__main__':
    dataDir = "/Volumes/Pony/ModelG_data_backup/"
    fn1 = "zcritical_N080_m-0481100_h006000_c00500.h5"
    fn2 = "zcritical_N080_m-0481100_h004000_c00500.h5"
    fn3 = "zcritical_N080_m-0481100_h003000_c00500.h5"
    fn4 = "zcritical_N080_m-0481100_h002000_c00500.h5"
    fn5 = "zcritical_N080_m-0481100_h010000_c00500.h5"
    fn6 = "zminus2_N080_m-0499128_h003000_c00500.h5"
    dataDir = "/Volumes/Pony/ModelG_data_backup/superpaperIV/"
    fn = "zcritical_diffusiononly_N080_m-0481100_h003000_c00500.h5"
    dt = 0.72
    dec = 50
    
    #fns = [fn1, fn2, fn3 , fn4]
    fns = [fn]

    for fn in fns:
        print(fn)

        t0 = time.time()

        data = ConfResults(fn = dataDir + fn,thTime=10000,dt=dt,  data_format="new", loadFourierBool = False)
        #for k in ["phi"]:
        for key in ["A", "V"]:
            for k in range(0,5):
            #key = "V"
                if(k>0):
                    bkey=key+"kk{}".format(k)
                else:
                    bkey=key
                print(k)
 			#if k == 0:
			#    data.load("wallF", key, direc="../processed")
            
                data.computeOtOtpBlocked(key, momNum = k, tMax = 1000.0, nBlocks = 5,  decim = dec, errFunc = lambda x : (np.mean(x, axis = 0), np.std(x, axis = 0)), parallel=True)

                data.save("OtOttp", bkey)    
                data.save("OtOttp_blocks", bkey)    
                if k==1:
               	    data.save("wallF", key)    

