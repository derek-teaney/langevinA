import random
import numpy as np
import Loader as ld
import matplotlib.pyplot as plt

def makedata(n, mean, variance, tau):
    """
    Generates a set of test data of size n, with  specified mean, 
    variance, and correlation time
    
    Args:
        n (int) : the number of time steps
        mean (real) : the mean of the array
        variance (real) : the variance of the output array
        tau (real) : the correlation time of the array

    Returns:
        array : An array with the data

    """
    # The routine  works by implementing the langevin process:
    #
    # dp/dt = - 1/tau * p + xi
    p = 0
    data =  mean*np.ones(n)
    eta = 1./tau
    for i in range(0, n):
        xi = random.gauss(0., np.sqrt(2.*eta*variance))
        p = p - eta* p +  xi
        data[i] = data[i] + p

    return(data)


def testdata() :

    # generate test data and make a graph
    n = 40000
    mean_in = 2.
    variance_in = 2.1
    tau_in = 30.
    data = makedata(n, mean_in, variance_in, tau_in)
    #plt.plot(data)
    #plt.savefig("pltfig.png")

    # plt.show() 
    # compute the correlation function and make a graph
    correlationfunc = ld.autocorrelationFunction(data, 80, 0)
    plt.plot(correlationfunc)
    plt.savefig("pltfig.png")
    #plt.show()

    # print out the mean and error in the mean of the test data
    print("Input: ", mean_in, np.sqrt(variance_in*2.0*tau_in/n), tau_in)

    # this should agree within error 
    mean, err, tau = ld.autocorrelatedErr(data, tmax = 400)
    print("Correlation: ", mean, err, tau)

    mean, err = ld.blocking(data, nBlock=20)
    print("Blocking: ", mean, err)

    mean, err = ld.jackknife(data, nSamples=20)
    print("Jacknife: ", mean, err) 

    mean, err = ld.blockedBootstrap(data, nBlock=20)
    print("BlockedBootstrap: ", mean, err) 


if __name__ == "__main__":
    testdata()