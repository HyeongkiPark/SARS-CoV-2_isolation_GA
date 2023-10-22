import os
import sys
import math
import joblib
import random
import argparse
import datetime
import warnings
import numpy as np 
from multiprocessing import Pool

from subFunc import save_text, calc_rtp, calc_contact_history

warnings.simplefilter('ignore')

def calc_iteration(shape, scale, duration, i_iter, f_prob, T_star_prob, total_vl_delta_1_0, total_vl_delta_1_5, vaccination_rate=0.0, imcomplete_isolation_prob=0.0):

	print("f_prob : {}, T_star_prob : {}".format(f_prob, T_star_prob), flush = True)

	incubation_sample_list = random.sample([i for i in range(i_iter)], int(i_iter * f_prob))
	vaccination_sample_list = random.sample([i for i in range(i_iter)], int(i_iter * vaccination_rate))
	T_star_day = T_star_prob * duration

	total_populations = []
	total_rtp = []

	for i in range(i_iter):

		if i % 100 == 0:
			print("iter {}".format(i), flush = True)

		ind_populations = calc_contact_history(shape, scale, duration)

		if i in vaccination_sample_list:

			if i in incubation_sample_list:
				ind_rtp = calc_rtp(ind_populations, total_vl_delta_1_5, T_star_day, imcomplete_isolation_prob)
			else:
				ind_rtp = calc_rtp(ind_populations, total_vl_delta_1_5, duration, imcomplete_isolation_prob)

		else:
			if i in incubation_sample_list:
				ind_rtp = calc_rtp(ind_populations, total_vl_delta_1_0, T_star_day, imcomplete_isolation_prob)
			else:
				ind_rtp = calc_rtp(ind_populations, total_vl_delta_1_0, duration, imcomplete_isolation_prob)

		total_populations.append(ind_populations)
		total_rtp.append(ind_rtp)

	save_joblib(shape, scale, f_prob, T_star_prob, vaccination_rate, imcomplete_isolation_prob, incubation_sample_list, vaccination_sample_list, total_populations, total_rtp)
	save_text(shape, scale, f_prob, T_star_prob, vaccination_rate, imcomplete_isolation_prob, vaccination_sample_list, total_rtp, T_star_day, duration)


def save_joblib(shape, scale, f_prob, T_star_prob, vaccination_rate, imcomplete_isolation_prob, incubation_sample_list, vaccination_sample_list, total_populations, total_rtp):

	current_dir = os.getcwd()

	new_dir_path = current_dir + '/results/vac_rate_{}_ii_prob_{}'.format(vaccination_rate, imcomplete_isolation_prob) + '/shape_{}_scale_{}/f_{}_T_{}'.format(shape, scale, f_prob, T_star_prob)

	joblib.dump(incubation_sample_list, new_dir_path + '/ind_inc_sample.pkl', compress = 3)
	joblib.dump(vaccination_sample_list, new_dir_path + '/ind_vac_sample.pkl', compress = 3)
	joblib.dump(total_populations, new_dir_path + '/ind_pop.pkl', compress = 3)
	joblib.dump(total_rtp, new_dir_path + '/ind_rtp.pkl', compress = 3)



def aggregate_info(shape_range, scale_range, duration, i_iter, f_range, T_range, total_vl_delta_1_0, total_vl_delta_1_5, vaccination_rate=0.0, imcomplete_isolation_prob=0.0):

	z = []

	for ind_shape in shape_range:

		for ind_scale in scale_range:

			for f_prob in f_range:

				for T_star_prob in T_range:

					z.append([ind_shape, ind_scale, duration, i_iter, f_prob, T_star_prob, total_vl_delta_1_0, total_vl_delta_1_5, vaccination_rate, imcomplete_isolation_prob])

	return z

def make_directory_from_current_directory(args):

	# args_np = np.array(args)

	current_dir = os.getcwd()

	data_dir = '/results'

	# shape_list = list(set(args_np[:, 0]))
	# scale_list = list(set(args_np[:, 1]))

	shape_list = list(set([args[i][0] for i in range(len(args))]))
	scale_list = list(set([args[i][1] for i in range(len(args))]))

	# f_list = list(set(args_np[:, 4]))
	# T_list = list(set(args_np[:, 5]))

	f_list = list(set([args[i][4] for i in range(len(args))]))
	T_list = list(set([args[i][5] for i in range(len(args))]))

	# vaccination_rate = list(set(args_np[:, 8]))[0]
	# imcomplete_isolation_prob = list(set(args_np[:, 9]))[0]

	vaccination_rate = args[0][8]
	imcomplete_isolation_prob = args[0][9]

	new_dir_path = current_dir + data_dir + '/vac_rate_{}_ii_prob_{}'.format(vaccination_rate, imcomplete_isolation_prob)

	try:
		os.mkdir(new_dir_path)
	except:
		print('Error : Save directory already exists; /vac_rate_{}_ii_prob_{}'.format(vaccination_rate, imcomplete_isolation_prob))
		sys.exit()

	for ind_shape in shape_list:
		for ind_scale in scale_list:

			new_dir_path = current_dir + data_dir + '/vac_rate_{}_ii_prob_{}'.format(vaccination_rate, imcomplete_isolation_prob) + '/shape_{}_scale_{}'.format(ind_shape, ind_scale)

			try:
				os.mkdir(new_dir_path)
			except:
				print('Error : Save directory already exists; /shape_{}_scale_{}'.format(ind_shape, ind_scale))
				sys.exit()

			for ind_f in f_list:
				for ind_T in T_list:

					new_dir_path = current_dir + data_dir + '/vac_rate_{}_ii_prob_{}'.format(vaccination_rate, imcomplete_isolation_prob) + '/shape_{}_scale_{}/f_{}_T_{}'.format(ind_shape, ind_scale, ind_f, ind_T)

					try:
						os.mkdir(new_dir_path)
					except:
						print('Error : Save directory already exists; /shape_{}_scale_{}/f_{}_T_{}'.format(ind_shape, ind_scale, ind_f, ind_T))
						sys.exit()


def wrapper_calc_iteration(args):
    return calc_iteration(*args)



if __name__ == '__main__':

	print('file name : ', os.path.basename(__file__))

	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--iteration', type=int)
	parser.add_argument('-f', '--f_range', required=True, nargs="*", type=float) 
	parser.add_argument('-vr', '--vaccination_rate', type=float)
	parser.add_argument('-ip', '--imcomplete_isolation_prob', type=float)
	parser.add_argument('-n', '--n_processes', type=int)
	args = parser.parse_args( )

	# parameters
	i_iter = args.iteration
	f_range = args.f_range
	vaccination_rate = args.vaccination_rate
	imcomplete_isolation_prob = args.imcomplete_isolation_prob
	n_processes = args.n_processes

	print(f'vaccination_rate : {vaccination_rate}, imcomplete_isolation_prob : {imcomplete_isolation_prob}, n_processes : {n_processes}')

	Dmax = 80

	'''

	shape is 1.5
	scale is 1.5,

	f : ratio of symptomatic individuals
	1 - f : ratio of asymptomatic individuals
	f_range 0.0, 0.2, 0.4, 0.6, 0.8, 1.0

	T* : incubation period
	T*max : 80 (Dmax) days
	T_range 0.1, 0.3, 0.5, 0.7, 0.9, 1.0

	'''
	
	try:

		total_vl_delta_1_0 = joblib.load('pre_scripts/pre_data/vLoad_infPeriod_Log10_delta_1.0.pkl')
		total_vl_delta_1_5 = joblib.load('pre_scripts/pre_data/vLoad_infPeriod_Log10_delta_1.5.pkl')

	except:
		print('Error : You need to execute "vLoad.py", "infPeriod.py", "vLoad_infPeriod.py" beforehand.')
		sys.exit()

	shape_range = [1.5]
	scale_range = [1.5]

	# default T_range
	T_range = [1/Dmax, 2/Dmax, 3/Dmax, 4/Dmax, 5/Dmax, 6/Dmax, 7/Dmax, 8/Dmax, 9/Dmax, 10/Dmax]

	print('started at ', datetime.datetime.now(), flush = True)

	args = aggregate_info(shape_range, scale_range, Dmax, i_iter, f_range, T_range, total_vl_delta_1_0, total_vl_delta_1_5, vaccination_rate, imcomplete_isolation_prob)

	print('args, ', len(args), flush = True)

	make_directory_from_current_directory(args)

	p = Pool(processes=n_processes)
	p.map(wrapper_calc_iteration, args)
	p.close()

	print('finish at ', datetime.datetime.now(), flush = True)

