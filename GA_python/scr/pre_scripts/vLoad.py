import os
import math
import random
import joblib
import numpy as np
import argparse
import datetime
import matplotlib.pyplot as plt
from scipy.integrate import odeint

from multiprocessing import Pool

output_dir = 'pre_scripts/pre_data'

beta_max = 0.0000005
p_max = 100000


def simulator_(tau_max, n_delta, beta, p):

    def ODE_(var, t_, n_delta, beta, p):

        T_0 = 133000
        c = 20
        delta_max = 1.57 * n_delta
        delta_50 = 100
        delta = delta_max*(p/(p + delta_50))

        F, V = var[0], var[1]

        dFdt = -beta * V * F
        dVdt = p * beta * T_0 * (1/c) * F * V - delta * V

        return [dFdt, dVdt]

    F0 = 1
    V0 = 0.01
    var = [F0, V0]
    dt = 0.01
    t_ = np.arange(0, tau_max, dt)

    sim_results = odeint(ODE_, var, t_, args = (n_delta, beta, p))

    vl_Log10 = sum(np.log10(sim_results[:, 1] + 1)) / (1/dt)

    return vl_Log10

def collect_args(Dmax, n_delta):

    global beta_max, p_max

    args = []

    for i in np.linspace(0.0, beta_max, 100): # beta
        for j in np.linspace(0, p_max, 100): # p
            args.append([Dmax, n_delta, i, j])

    return args

def wrapper_simulator_(args):
    return simulator_(*args)


if __name__ == '__main__':

    print('file name : ', os.path.basename(__file__))

    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--n_delta', type=float)
    parser.add_argument('-n', '--n_processes', type=int)
    args = parser.parse_args( )

    # parameters
    n_delta = args.n_delta
    n_processes = args.n_processes

    print(f'n_delta : {n_delta}, n_processes : {n_processes}')

    #simulation parameters
    Dmax = 80

    print('started at ', datetime.datetime.now())

    args = collect_args(Dmax, n_delta)

    p = Pool(processes=n_processes)
    viral_load = p.map(wrapper_simulator_, args)
    p.close()

    joblib.dump(viral_load, '{}/vLoad_Log10_delta_{}.pkl'.format(output_dir, n_delta), compress = 3)


    print('finish at ', datetime.datetime.now())
