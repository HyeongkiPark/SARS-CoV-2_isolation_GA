import os
import re
import glob
import random
import joblib
import numpy as np
import argparse

beta_max = 0.0000005
p_max = 100000

def save_best_txt(best_list, save_path):

    best_beta = best_list[0]
    best_p = best_list[1]
    best_i_d = best_list[2]
    best_rtp = best_list[3]

    text = 'contact\nbeta\tp\ti_d\tr0\n'

    for i in range(len(best_beta)):
        text = '{}{}\t{}\t{}\t{}\n'.format(text,best_beta[i],best_p[i],best_i_d[i], best_rtp[i])

    with open(save_path + '/best_param.txt','w') as f:
        f.write(text)

def aggregate_data(total_vl_delta_1_0, total_vl_delta_1_5, shape, scale, f_prob, T_star_prob, vaccination_prob, imcomplete_isolation_prob, vaccination_sample_list, iter_rtp):

	global beta_max, p_max

	current_dir = os.getcwd()

	save_path = current_dir + '/results/vac_prob_{}_ii_prob_{}'.format(vaccination_prob, imcomplete_isolation_prob) + '/shape_{}_scale_{}/f_{}_T_{}'.format(shape, scale, f_prob, T_star_prob)


	z = []
	for i in np.linspace(0.0, beta_max, 100): # beta
	    for j in np.linspace(0, p_max, 100): # p
	        z.append([i, j])

	best_beta = []
	best_p = []
	best_id = []
	best_rtp = []

	for i in range(len(iter_rtp)):
		
	    ind_rtp = iter_rtp[i]

	    best_r = [j for j, v in enumerate(ind_rtp) if v == max(ind_rtp)]


	    for b in best_r:
	        best_beta.append(z[b][0])
	        best_p.append(z[b][1])

	        if i in vaccination_sample_list:

	        	best_id.append(total_vl_delta_1_5[b]['inf_id'])

	        else:

	        	best_id.append(total_vl_delta_1_0[b]['inf_id'])

	        
	        best_rtp.append(max(ind_rtp))

	best_list = [best_beta, best_p, best_id, best_rtp]

	return best_list, save_path

def save_text(shape, scale, f_prob, T_star_prob, vaccination_prob, imcomplete_isolation_prob, vaccination_sample_list, iter_rtp, T_star_day, duration):

    total_vl_delta_1_0 = joblib.load('pre_scripts/pre_data/vLoad_infPeriod_Log10_delta_1.0.pkl')
    total_vl_delta_1_5 = joblib.load('pre_scripts/pre_data/vLoad_infPeriod_Log10_delta_1.5.pkl')

    best_list, save_path = aggregate_data(total_vl_delta_1_0, total_vl_delta_1_5, shape, scale, f_prob, T_star_prob, vaccination_prob, imcomplete_isolation_prob, vaccination_sample_list, iter_rtp)

    save_best_txt(best_list, save_path)


def calc_rtp(ind_populations, total_vl, T_star_day, imcomplete_isolation_prob):


    total_r_values = []

    for i in range(10000):

        if len(total_vl[i]) == 1:
            total_r_values.append(0)
            continue

        r_value = []

        total_dict = total_vl[i]

        total_prob = list(total_dict.keys())[2:]
        start_time = total_dict['start_time']

        for u in range(len(total_prob)):

            ind_prob = total_prob[u]

            ind_gens = total_dict[ind_prob]

            for e in range(int(len(ind_gens)/2)):

                ind_start_time = ind_gens[e*2 : e*2+2][0]
                ind_end_time = ind_gens[e*2 : e*2+2][1]

                for k in range(ind_start_time, ind_end_time + 1):

                    current_day = (start_time + k)//100

                    ind_pop = ind_populations[current_day]/100

                    if current_day > T_star_day:
                        
                        ind_pop *= imcomplete_isolation_prob

                    secondary_infected = ind_pop * ind_prob

                    r_value.append(secondary_infected)


        total_r_values.append(sum(r_value))

    return total_r_values


def Uniform():
    return abs(random.random())

def rand_gamma(kappa, theta):
    uni_num = 0

    int_kappa  = int(kappa)
    frac_kappa = kappa - int_kappa
    x_int = 0
    for i in range(int_kappa):
        x_int += -np.log(Uniform())


    if( abs(frac_kappa) < 0.01 ):
        x_frac = 0

    else:

        b=(np.exp(1.0) + frac_kappa)/np.exp(1.0)
        while 1 :
            u=Uniform()
            p = b * u
            uu = Uniform()
            if p <= 1.0:
                    x_frac = p**(1.0/frac_kappa)
                    if(uu <= np.exp(-x_frac)):
                        break
                
            else:
                x_frac =- np.log((b - p)/frac_kappa)
                if (uu <= x_frac**(frac_kappa - 1.0)):
                        break


            uni_num += 1

    return (x_int + x_frac) * theta


def porand(x):
    k = 0
    uni_num = 0
    a = Uniform()
    b = np.exp(-x)
    a_i = 0
    while(a >= b):
        a *= Uniform()
        k += 1
        uni_num += 1

    return k


def calc(shape, scale):
    return porand(rand_gamma(shape, scale))


def calc_contact_history(shape, scale, duration):
    p_list = []
    
    for i in range(duration):
        day_contact = calc(shape, scale)
        p_list.append(day_contact)

    
    return p_list

