import os
import math
import random
import joblib
import numpy as np
import argparse
import datetime
from scipy.integrate import odeint

from multiprocessing import Pool


output_dir = 'pre_scripts/pre_data'

beta_max = 0.0000005
p_max = 100000

def inf_simulator_(tau_max, n_delta, beta, p):

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

    vl_ = sim_results[:, 1]

    return vl_

def calc_every_prob(tau_max, n_delta, beta, p):

    def trans_prob(log_vl):

        if log_vl < 5:
            return 0
        elif 5 <= log_vl < 7:
            return 0.12
        elif 7 <= log_vl < 10:
            return 0.15
        elif 10 <= log_vl:
            return 0.24

    def make_data_dict(inf_id, start_time, _data):

        data_dict = {}
        data_dict['inf_id'] = inf_id
        data_dict['start_time'] = start_time

        for i in list(set(_data)):

            _indeces = [j for j, x in enumerate(_data) if x == i]

            ind_t = [_indeces[0]]

            for k in range(len(_indeces)-1):
                if _indeces[k] + 1 != _indeces[k+1]:
                    ind_t.append(_indeces[k])
                    ind_t.append(_indeces[k+1])

            ind_t.append(_indeces[-1])

            data_dict[i] = ind_t

        return data_dict


    vi = inf_simulator_(tau_max, n_delta, beta, p)

    total_prob = []

    inf_id = len(vi)/100

    for u in range(len(vi)):

        ind_prob = trans_prob(math.log10(vi[u] + 1))

        total_prob.append(ind_prob)

    if sum(total_prob) == 0:
        return {'inf_id' : inf_id}

    start_time = np.nonzero(total_prob)[0][0]
    end_time = np.nonzero(total_prob)[0][-1] + 1

    total_prob = list(np.array(total_prob)[start_time:end_time])

    data_dict = make_data_dict(inf_id, start_time, total_prob)

    return data_dict

def collect_args(duration_of_infectious_period, n_delta):

    global beta_max, p_max

    args = []

    l = 0

    for i in np.linspace(0.0, beta_max, 100): # beta
        for j in np.linspace(0, p_max, 100): # p

            k = duration_of_infectious_period[l]

            args.append([k, n_delta, i, j])

            l += 1

    return args

def wrapper_calc_every_prob(args):
    return calc_every_prob(*args)

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

    # parameters

    D50 = 198.58
    Dk = 1.17
    Dmax = 80

    print('started at ', datetime.datetime.now())

    duration_of_infectious_period = joblib.load('{}/infPeriod_delta_{}.pkl'.format(output_dir, n_delta))

    args = collect_args(duration_of_infectious_period, n_delta)

    p = Pool(processes=n_processes)
    viral_load_infPeriod = p.map(wrapper_calc_every_prob, args)
    p.close()

    joblib.dump(viral_load_infPeriod, '{}/vLoad_infPeriod_Log10_delta_{}.pkl'.format(output_dir, n_delta), compress = 8)


    print('finish at ', datetime.datetime.now())
