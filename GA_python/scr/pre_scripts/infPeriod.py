import os
import random
import joblib
import numpy as np
import argparse
import datetime

from multiprocessing import Pool

output_dir = 'pre_scripts/pre_data'

def collect_args(vi_amount, D50, Dk, Dmax):

    args = []

    for vi_i in vi_amount:
        args.append([D50, Dk, vi_i, Dmax])

    return args


def calc_infPeriod(D50, Dk, vi, Dmax):

    return Dmax*(vi**Dk/(vi**Dk+D50**Dk))

def wrapper_calc_infPeriod(args):
    return calc_infPeriod(*args)

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
    D50 = 198.58
    Dk = 1.17
    Dmax = 80

    print('started at ', datetime.datetime.now())

    viral_load = joblib.load('{}/vLoad_Log10_delta_{}.pkl'.format(output_dir, n_delta))

    args = collect_args(viral_load, D50, Dk, Dmax)

    p = Pool(processes=n_processes)
    duration_of_infectious_period = p.map(wrapper_calc_infPeriod, args)
    p.close()

    joblib.dump(duration_of_infectious_period, '{}/infPeriod_delta_{}.pkl'.format(output_dir, n_delta), compress = 3)

    print('finish at ', datetime.datetime.now())
