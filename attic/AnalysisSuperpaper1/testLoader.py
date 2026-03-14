import matplotlib.pyplot as plt
import pprint
import numpy as np
import timeit
import Loader as ms


def getloader(name):
    if name == "donly":
        loader = ms.Loader(
            "/Users/derekteaney/common/superfluid/cori/run/hzero/donly_N080_m-0481100_h003000_c00500/donly_N080_m-0481100_h003000_c00500.h5",
            removemean=True,
            starttime=8000,
            nblocks=8)
    elif name == "pionrun8":
        loader = ms.Loader(
            "/Users/derekteaney/common/superfluid/cori/run/pionrun/pionrun8_N080_m-0501265_h003000_c00500/pionrun8_N080_m-0501265_h003000_c00500.h5",
            removemean=True,
            starttime=10000,
            nblocks=10)
    else:
        raise KeyError("Bad selection", name)

    return loader


def test_loader1(loader):
    data = loader.load("phi0")
    print("Test of simple load:", data.shape)


def test_loader2(loader):
    data_k = loader.load("phi0_kx", k=1)
    print("Test of k type load:", data_k.shape, data_k.dtype)

# def test_loader3(loader) :
#     data, data_k = loader.load_sets("phi0", "phi0_kx",k=1)
#     print(data.shape, data.dtype, data_k.shape, data_k.dtype)


def test_loader4(loader):
    data = loader.read_wall("phi0_x")
    print("Test of full wall load:", data.shape)


def test_loader5(loader):
    data = loader.load("phi0_x")
    block_start = loader.load("phi0_x", block=0)
    block_end = loader.load("phi0_x", block=loader.nblocks-1)
    print("Test of block loads:", data.shape,
          block_start.shape, block_end.shape)


def test_blocks(loader):
    for ib in range(0, loader.nblocks):
        print("Test Blocks: block:", ib, loader.blockstart(
            ib), loader.blockstart(ib+1))
    print("Test Blocks: ntime_total:", loader.ntime_total)

# def test_computeOtOtpStatic(loader) :
#     d = loader.load("phi1_x", block=1)
#     result1 = ms.StaticCorrelator.correlate_data(d)
#     result2 = np.mean(ms.computeOtOtpStatic(d,d), axis=0)
#     print("Test computeOtOtpStatic:", result1 - result2)


def plot_ototp(ototp, ototperr, blocks=None):
    fig, ax = plt.subplots()
    x = np.arange(len(ototp))
    ax.errorbar(x, ototp, ototperr, errorevery=100, color='black')

    if blocks is not None:
        for item in blocks:
            ax.plot(x, item)
    plt.show()


def plot_oxoxp(oxoxp, oxoxperr, blocks=None):
    fig, (ax1, ax2) = plt.subplots(1, 2)
    x = np.arange(oxoxp.shape[0])
    ax1.errorbar(x, oxoxp, oxoxperr, color='black')
    ax2.errorbar(x, oxoxp, oxoxperr, color='black')

    if blocks is not None:
        for item in blocks:
            ax1.plot(x, item)
    plt.show()


def test_correlate1(loader):
    print("Testing time Correlator with correlate_keys with one key phi0")
    tc = ms.TimeCorrelator(loader, 4, nTMax=1000)
    ototp, ototperr, blocks = tc.correlate_keys("phi0")
    plot_ototp(ototp, ototperr, blocks)


def test_correlate2(loader):
    print("Testing time Correlator with correlate_keys with keys phi1, phi2, phi3")
    tc = ms.TimeCorrelator(loader, nTMax=1000)
    ototp, ototperr, blocks = tc.correlate_keys("phi1", "phi2", "phi3")
    plot_ototp(ototp, ototperr, blocks)


def test_correlate3(loader):
    print("Testing time correlate_key with keygroup phi123")
    tc = ms.TimeCorrelator(loader, nTMax=1000)
    ototp, ototperr, blocks = tc.correlate_key("phi123")
    plot_ototp(ototp, ototperr, blocks)


def test_correlate4(loader):
    print("Testing time correlator with correlate_keys pulling out list of keys by hand")
    tc = ms.TimeCorrelator(loader, nTMax=1000)
    keys = loader.keygroups["phi123"]
    ototp, ototperr, blocks = tc.correlate_keys(*keys)
    plot_ototp(ototp, ototperr, blocks)


def test_correlate5(loader):
    print("Testing time Correlator with keygroup A123_kxyz")
    tc = ms.TimeCorrelator(loader, nTMax=4000)
    ototp, ototperr, blocks = tc.correlate_key("A123_kxyz", k=2)
    plot_ototp(ototp, ototperr, blocks)


def test_static_correlate1(loader):
    print("Testing single key static correlator")
    sc = ms.StaticCorrelator(loader)
    oxoxp, oxoxperr, blocks = sc.correlate_keys("phi1_x")
    plot_oxoxp(oxoxp, oxoxperr, blocks)


def test_static_correlate2(loader):
    print("Testing keygroup static correlator", "phi123_xyz")
    sc = ms.StaticCorrelator(loader)
    oxoxp, oxoxperr, blocks = sc.correlate_key("phi123_xyz")
    plot_oxoxp(oxoxp, oxoxperr, blocks)


def test_static_correlate3(loader):
    print("Testing multiple keys static correlator:",
          "phi1_x", "phi1_y", "phi1_z")
    sc = ms.StaticCorrelator(loader)
    oxoxp, oxoxperr, blocks = sc.correlate_keys("phi1_x", "phi1_y", "phi1_z")
    plot_oxoxp(oxoxp, oxoxperr, blocks)

def test_keygroups(loader):
    pp = pprint.PrettyPrinter()
    pp.pprint(loader.keygroups)

def test_times(loader):
    times = loader.load("times")
    print("Testing the times loader")
    print(times)

def all_tests(loader):
    test_keygroups(loader)
    test_times(pionloader)
    test_loader1(loader)
    test_loader2(loader)
    test_correlate1(loader)
    test_correlate2(loader)
    test_correlate3(loader)
    test_correlate4(loader)
    test_correlate5(pionloader)
    test_blocks(loader)
    test_loader4(loader)
    test_loader5(loader)
    test_static_correlate1(loader)
    test_static_correlate2(loader)
    test_static_correlate3(loader)

if __name__ == '__main__':
    pionloader = getloader("donly")
    all_tests(pionloader)
