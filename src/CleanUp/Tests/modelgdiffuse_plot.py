import h5py as h5
import numpy as np
import matplotlib.pyplot as plt
import os
import pprint


def plot_h5_data(file_path, dataset_name, time, field):
    # Open the h5 file
    with h5.File(file_path, 'r') as f:
        # Read the dataset
        data = f[dataset_name]

        L = data.shape[0]
        x = np.arange(L)

        for i in range(0, 1):
            v = data[int(L/2), int(L/2) + i, :, field]
            plt.plot(x, v, 'r-', label=dataset_name + f'_x{i}')

        # for i in range(0, 8, 4):
        #     v = data[16, :, 16 + i, 0]
        #     plt.plot(x, v, label=dataset_name + f'_y{i}')

        # for i in range(0, 8, 4):
        #     v = data[:, 16, 16 + i, 0]
        #     plt.plot(x, v, label=dataset_name + f'_z{i}')


def plot_solution(time, params, field):
    print(f"Plotting solution at time {time}")
    Lx = params.get('Lx', 32)
    x = np.arange(0, Lx, 1)
    sigmax = params.get('sigmax', 3.0)
    # sigmay = params.get('sigmay', 3.0 * 10e8)
    # sigmaz = params.get('sigmaz', 3.0 * 10e8)
    amplitude = params.get('amplitude', 1.0)
    # Calculate the Gaussian function
    Gamma = params.get('Gamma', 1.0)
    D0 = params.get('D', Gamma/3.)
    if field < 4:
        D = Gamma
    else:
        D = D0

    sigmaxt = np.sqrt(sigmax**2 + 2 * D * time)

    # Make the solution periodic by summing over multiple periods
    gaussian = np.zeros_like(x, dtype=float)
    for i in range(-8, 9, 1):
        L = Lx/2 + i * Lx
        gaussian += amplitude * \
            np.exp(-((x - L)**2) / (2 * sigmaxt**2)) * sigmax/sigmaxt

    plt.plot(x, gaussian, 'b--', label=f'Gaussian at t={time:.2f}')


if __name__ == "__main__":
    # Read the JSON file for additional parameters
    json_file = f'./modelgdiffuse/modelgdiffuse.json'
    if os.path.exists(json_file):
        import json
        with open(json_file, 'r') as f:
            params = json.load(f)
        print(f"Parameters for run:")
        pprint.pprint(params)

    # Plot time evolution of the field
    field = 3
    dt = params.get('deltat')
    finaltime = params.get('finaltime', 40)
    arange = np.arange(0, finaltime, dt)
    for t in arange[::24]:
        file = f'./modelgdiffuse/modelgdiffuse_t_{t:.2f}_save.h5'

        # check if the file exists
        if not os.path.exists(file):
            print(f"File {file} does not exist.")
            continue
        plot_h5_data(file, 'o4fields', t, field)

        plot_solution(t, params, field)

    plt.tight_layout()
    plt.title(f"Field {field}")
    plt.xlabel('Position')
    plt.ylabel('Field Value')
    plt.show()
    # plt.savefig(f'./modelgdiffuse_test.png', dpi=300, bbox_inches='tight')
