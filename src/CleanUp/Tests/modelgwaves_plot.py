import h5py as h5
import numpy as np
import matplotlib.pyplot as plt
import math
import os
import pprint

even_odd = True


def plot_h5_data(file_path, time, field):
    global even_odd
    # Open the h5 file
    with h5.File(file_path, 'r') as f:
        # Read the dataset
        data = f['o4fields']

        L = data.shape[0]
        x = np.arange(L)

        for i in range(0, 1):
            v = data[int(L/2), int(L/2) + i, :, field]
            if (even_odd):
                plt.plot(x, v, '--')
                even_odd = not even_odd
            else:
                plt.plot(x, v, '-')
                even_odd = not even_odd
        plt.draw()
        plt.pause(0.4)

        # for i in range(0, 8, 4):
        #     v = data[16, :, 16 + i, 0]
        #     plt.plot(x, v, label=dataset_name + f'_y{i}')

        # for i in range(0, 8, 4):
        #     v = data[:, 16, 16 + i, 0]
        #     plt.plot(x, v, label=dataset_name + f'_z{i}')


if __name__ == "__main__":
    even_odd = False
    # Read the JSON file for additional parameters
    json_file = f'./modelgwaves/modelgwaves.json'
    if os.path.exists(json_file):
        import json
        with open(json_file, 'r') as f:
            params = json.load(f)
        print(f"Parameters for run:")
        pprint.pprint(params)
    else:
        raise FileNotFoundError(f"JSON file {json_file} not found.")

    # file = f'./modelgwaves/modelgwaves_initial_save.h5'
    # if os.path.exists(file):
    #     print(f"Plotting initial data from {file}")
    #     plot_h5_data(file, 0, 0)
    # else:
    #     raise FileNotFoundError(f"HDF5 file {file} not found.")

    # Plot time evolution of the field
    field = 0
    dt = params.get('deltat')
    finaltime = params.get('finaltime', 40)
    f = math.sqrt(params.get('f2_constant', 5.0))
    arange = np.arange(0, finaltime, dt)
    trange = int(0.125 * finaltime / dt)
    plt.ylim(-1.1*f, 1.1*f)
    for t in arange[0*trange:1*trange:int(0.5/dt)]:
        file = f'./modelgwaves/modelgwaves_t_{t:.2f}_save.h5'
        print(f"Plotting data from {file} at time {t:.2f} for field {field}")

        # check if the file exists
        if not os.path.exists(file):
            print(f"File {file} does not exist.")
            continue
        plot_h5_data(file,  t, field)
    plt.show()
