import os
import re
import sys
import math
import glob
import scipy
import joblib
import random
import datetime
import warnings
import numpy as np
import pandas as pd
import seaborn as sns
from multiprocessing import Pool
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.colors import LinearSegmentedColormap

warnings.simplefilter('ignore')

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
# plt.rcParams['font.family'] = 'Arial'
plt.rcParams['xtick.direction'] = 'out'
plt.rcParams['ytick.direction'] = 'out'

beta_max = 0.0000005
p_max = 100000

beta_values = np.linspace(0.0, beta_max, 100)
p_values = np.linspace(0, p_max, 100)
P_values, Beta_values = np.meshgrid(p_values, beta_values)

colors = ['#000067', 'mediumblue', 'forestgreen', 'orange', 'red']
nodes = [0.0, 0.25, 0.5, 0.75, 1.0]
cmap_self = LinearSegmentedColormap.from_list("mycamp", list(zip(nodes, colors)))

colors = ['#CC7F02', '#d18c1d', '#d79a38', '#dda955', '#e2b770', '#e8c58d', '#eed3a8', '#f3e2c5', '#fffefc']# final color!
nodes = [0.0/10.0, 0.1/10.0, 0.2/10.0, 0.3/10.0, 0.4/10.0, 0.5/10.0, 1.0/10.0, 5.0/10.0, 10.0/10.0]
cmap_self2 = LinearSegmentedColormap.from_list("mycamp", list(zip(nodes, colors)))

colors = ['#B5B5B5', '#B5B5B5']
nodes = [0.0, 1.0]
cmap_self3 = LinearSegmentedColormap.from_list("mycamp", list(zip(nodes, colors)))

def makeFigure_delta_1_0(f_num, T_range, allResults, min_max_rtp_list, gaResultsDict, peak_day_list_delta_1_0):

	global beta_max, p_max, beta_values, p_values, Beta_values, P_values

	Dmax = 80

	_color_list = ['#ff6680', '#a0d568', '#4fc1e8', '#ac92eb']

	fig = plt.figure(figsize=(16,25))

	grid_master = GridSpec(nrows = 24, ncols = 23+3, height_ratios = [1 for _ in range(24)])

	gs_1 = GridSpecFromSubplotSpec(nrows = 1, ncols = 1, subplot_spec = grid_master[0:3, 0:5])
	axes_1 = fig.add_subplot(gs_1[:, :])
	gs_2 = GridSpecFromSubplotSpec(nrows = 1, ncols = 2, subplot_spec = grid_master[0:3, 7:12])
	axes_2 = fig.add_subplot(gs_2[:, :])
	gs_3 = GridSpecFromSubplotSpec(nrows = 1, ncols = 3, subplot_spec = grid_master[0:3, 14:19])
	axes_3 = fig.add_subplot(gs_3[:, :])
	gs_4 = GridSpecFromSubplotSpec(nrows = 1, ncols = 4, subplot_spec = grid_master[0:3, 21:26])
	axes_4 = fig.add_subplot(gs_4[:, :])
	gs_5 = GridSpecFromSubplotSpec(nrows = 2, ncols = 1, subplot_spec = grid_master[4:7, 0:5])
	axes_5 = fig.add_subplot(gs_5[:, :])
	gs_6 = GridSpecFromSubplotSpec(nrows = 2, ncols = 2, subplot_spec = grid_master[4:7, 7:12])
	axes_6 = fig.add_subplot(gs_6[:, :])
	gs_7 = GridSpecFromSubplotSpec(nrows = 2, ncols = 3, subplot_spec = grid_master[4:7, 14:19])
	axes_7 = fig.add_subplot(gs_7[:, :])
	gs_8 = GridSpecFromSubplotSpec(nrows = 2, ncols = 4, subplot_spec = grid_master[4:7, 21:26])
	axes_8 = fig.add_subplot(gs_8[:, :])
	gs_9 = GridSpecFromSubplotSpec(nrows = 3, ncols = 1, subplot_spec = grid_master[8:11, 0:5])
	axes_9 = fig.add_subplot(gs_9[:, :])
	gs_10 = GridSpecFromSubplotSpec(nrows = 3, ncols = 2, subplot_spec = grid_master[8:11, 7:12])
	axes_10 = fig.add_subplot(gs_10[:, :])
	gs_11 = GridSpecFromSubplotSpec(nrows = 3, ncols = 3, subplot_spec = grid_master[8:11, 14:19])
	axes_11 = fig.add_subplot(gs_11[:, :])
	gs_12 = GridSpecFromSubplotSpec(nrows = 3, ncols = 4, subplot_spec = grid_master[8:11, 21:26])
	axes_12 = fig.add_subplot(gs_12[:, :])
	gs_13 = GridSpecFromSubplotSpec(nrows = 4, ncols = 1, subplot_spec = grid_master[12:15, 0:5])
	axes_13 = fig.add_subplot(gs_13[:, :])
	gs_14 = GridSpecFromSubplotSpec(nrows = 4, ncols = 2, subplot_spec = grid_master[12:15, 7:12])
	axes_14 = fig.add_subplot(gs_14[:, :])
	gs_15 = GridSpecFromSubplotSpec(nrows = 4, ncols = 3, subplot_spec = grid_master[12:15, 14:19])
	axes_15 = fig.add_subplot(gs_15[:, :])
	gs_16 = GridSpecFromSubplotSpec(nrows = 4, ncols = 4, subplot_spec = grid_master[12:15, 21:26])
	axes_16 = fig.add_subplot(gs_16[:, :])
	gs_17 = GridSpecFromSubplotSpec(nrows = 5, ncols = 1, subplot_spec = grid_master[17:24, 0:12])
	axes_17 = fig.add_subplot(gs_17[:, :])
	gs_18 = GridSpecFromSubplotSpec(nrows = 5, ncols = 1, subplot_spec = grid_master[17:24, 14:26])
	axes_18 = fig.add_subplot(gs_18[:, :])

	axes_list = [
	    [axes_1, axes_2, axes_3, axes_4],
	    [axes_5, axes_6, axes_7, axes_8],
	    [axes_9, axes_10, axes_11, axes_12],
	    [axes_13, axes_14, axes_15, axes_16],
	    [axes_17, axes_18]
	]

	divider1 = make_axes_locatable(axes_list[1][-1])
	cax1 = divider1.append_axes("right", size="5%", pad=0.1)
	divider2 = make_axes_locatable(axes_list[2][-1])
	cax2 = divider2.append_axes("right", size="5%", pad=0.1)

	min_symp_rtp, max_symp_rtp, min_asymp_rtp, max_asymp_rtp = min_max_rtp_list

	min_symp_rtp, max_symp_rtp, min_asymp_rtp, max_asymp_rtp = np.log10(min_symp_rtp+1), np.log10(max_symp_rtp+1), np.log10(min_asymp_rtp+1), np.log10(max_asymp_rtp+1)

	if f_num == 1.0:
	    min_asymp_rtp, max_asymp_rtp = min_symp_rtp, max_symp_rtp

	for t_num in T_range:

		for i in range(len(allResults)):

			f_focal, t_focal = allResults[i][0], allResults[i][1]
			if (f_focal == f_num) and (t_focal == t_num):
			    ind_result = allResults[i]
			    break

		_, _, kernel_df, beta_trajectory_mean, p_trajectory_mean, symp_rtp, asymp_rtp, ga_fitness_mean = ind_result

		symp_rtp, asymp_rtp = np.log10(np.array(symp_rtp)+1), np.log10(np.array(asymp_rtp)+1)

		num_index = int(T_range.index(t_num))

		sns.kdeplot(x=kernel_df['beta'], y=kernel_df['p'], kind="kde", color=_color_list[num_index], linewidths = 3, zorder = 1, ax=axes_list[0][num_index])
		axes_list[0][num_index].scatter(np.mean(kernel_df['beta']), np.mean(kernel_df['p']), s = 180, ec = 'black', fc = _color_list[num_index], linewidths = 1.5, zorder = 4)
		axes_list[0][num_index].scatter(kernel_df['beta'], kernel_df['p'], s = 50, ec = 'black', fc = 'white', linewidths = 1.0, zorder = 3)
		axes_list[0][num_index].plot(beta_trajectory_mean, p_trajectory_mean, color = 'black', linewidth = 1.0, zorder = 2)

		try:
		    b = list(np.mean(np.array(asymp_rtp), axis = 0)).index(max(np.mean(np.array(asymp_rtp), axis = 0)))
		    b_beta = beta_values[int(b/100)]
		    b_p = p_values[b % 100]

		    axes_list[1][num_index].scatter(Beta_values, P_values, c = np.mean(np.array(asymp_rtp), axis = 0).reshape(100, 100), cmap = cmap_self, vmin=min_asymp_rtp, vmax=max_asymp_rtp)#, norm=colors.LogNorm(vmin=min_asymp_rtp,vmax=max_asymp_rtp))
		    axes_list[1][num_index].scatter(b_beta, b_p, fc = 'white', ec = 'black', s = 200, linewidth = 1.8)

		    if t_num == 10/Dmax:
		        im1 = axes_list[1][num_index].scatter(Beta_values, P_values, c = np.mean(np.array(asymp_rtp), axis = 0).reshape(100, 100), cmap = cmap_self, vmin=min_asymp_rtp, vmax=max_asymp_rtp)#, norm=colors.LogNorm(vmin=min_asymp_rtp, vmax=max_asymp_rtp))
		        axes_list[1][num_index].scatter(b_beta, b_p, fc = 'white', ec = 'black', s = 200, linewidth = 1.8)
		        cbar = fig.colorbar(im1, cax=cax1)
		        cbar.set_label('$R_{TP}$', size=20, rotation=0, labelpad=20)
		        cbar.ax.tick_params(labelsize = 17)
		        intermediate_int = []
		        intermediate_str = []
		        for j in [np.log10(i+1) for i in range(10)]:
		            if (min_asymp_rtp < j) and (j < max_asymp_rtp):
		                intermediate_int.append(j)
		                intermediate_str.append(str(int(10**j-1)) + '.0')
		        vmin_vmax_int = [round(min_asymp_rtp, 2)] + intermediate_int
		        vmin_vmax_str = [str(int(10**min_asymp_rtp-1)) + '.0'] + intermediate_str
		        cbar.ax.set_yticks(vmin_vmax_int)
		        cbar.ax.set_yticklabels(vmin_vmax_str, fontsize = 14)
		except:
		    pass
	    
		try:
		    b = list(np.mean(np.array(symp_rtp), axis = 0)).index(max(np.mean(np.array(symp_rtp), axis = 0)))
		    b_beta = beta_values[int(b/100)]
		    b_p = p_values[b % 100]

		    axes_list[2][num_index].scatter(Beta_values, P_values, c = np.mean(np.array(symp_rtp), axis = 0).reshape(100, 100), cmap = cmap_self, vmin=min_asymp_rtp, vmax=max_asymp_rtp)
		    axes_list[2][num_index].scatter(b_beta, b_p, fc = 'white', ec = 'black', s = 200, linewidth = 1.8)

		    if t_num == 10/Dmax:
		        im2 = axes_list[2][num_index].scatter(Beta_values, P_values, c = np.mean(np.array(symp_rtp), axis = 0).reshape(100, 100), cmap = cmap_self, vmin=min_asymp_rtp, vmax=max_asymp_rtp)
		        axes_list[2][num_index].scatter(b_beta, b_p, fc = 'white', ec = 'black', s = 200, linewidth = 1.8)
		        cbar = fig.colorbar(im2, cax=cax2)
		        cbar.set_label('$R_{TP}$', size=20, rotation=0, labelpad=20)
		        cbar.ax.tick_params(labelsize = 17)
		        intermediate_int = []
		        intermediate_str = []
		        
		        for j in [np.log10(i+1) for i in range(10)]:
		            if (min_asymp_rtp < j) and (j < max_asymp_rtp):
		                intermediate_int.append(j)
		                intermediate_str.append(str(int(10**j-1)) + '.0')
		        vmin_vmax_int = [round(min_asymp_rtp, 2)] + intermediate_int
		        vmin_vmax_str = [str(int(10**min_asymp_rtp-1)) + '.0'] + intermediate_str
		        cbar.ax.set_yticks(vmin_vmax_int)
		        cbar.ax.set_yticklabels(vmin_vmax_str, fontsize = 14)
		except:
		    pass

		gaIterationNum = len(ga_fitness_mean)

		for i in range(gaIterationNum):
		    axes_list[3][num_index].plot(ga_fitness_mean[i], color = _color_list[num_index], linewidth = 0.5, alpha = 0.5)
		axes_list[3][num_index].plot(np.mean(ga_fitness_mean, axis = 0), color = 'black', linewidth = 1.0)
	    
	cont = axes_list[-1][1].contour(Beta_values, P_values, np.array(peak_day_list_delta_1_0).reshape(100, 100), [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0, 5, 10],colors='black', zorder = 2, linestyles = 'dashed')
	cont = axes_list[-1][1].contourf(Beta_values, P_values, np.array(peak_day_list_delta_1_0).reshape(100, 100), [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0, 5, 10],cmap=cmap_self2, zorder = 1)
	cont = axes_list[-1][1].contourf(Beta_values, P_values, np.array(peak_day_list_delta_1_0).reshape(100, 100), [99, 100],cmap=cmap_self3, zorder = 1)


	t_list = [allResults[i][1] for i in range(len(allResults))]

	for t_num in [10/Dmax, 9/Dmax, 8/Dmax, 7/Dmax, 6/Dmax, 5/Dmax, 4/Dmax, 3/Dmax, 2/Dmax, 1/Dmax]:
	    
		t_indices = [i for i, t in enumerate(t_list) if t == t_num]

		for i in t_indices:
			
			f_focal = allResults[i][0]

			beta_lastGen_mean = np.mean(allResults[i][2]['beta']) # kernel_df['beta']
			p_lastGen_mean = np.mean(allResults[i][2]['p']) # kernel_df['p']

			if f_focal != f_num:
			    continue

			if t_num in T_range:
			    num_index = int(T_range.index(t_num))
			    axes_list[-1][1].scatter(beta_lastGen_mean, p_lastGen_mean, ec = 'black', fc = _color_list[num_index], s = 180, linewidth = 2, zorder = 3)
			else:
			    axes_list[-1][1].scatter(beta_lastGen_mean, p_lastGen_mean, ec = 'black', fc = 'white', s = 180, linewidth = 2, zorder = 2)

		gaResult_list = gaResultsDict[f'{t_num}']

		for ind_ga in gaResult_list:

			f_focal, beta_focal, p_focal = ind_ga[0], ind_ga[1], ind_ga[2]

			if f_focal != str(f_num):
			    continue

			inf_vl = simulator_(calc_infPeriod(simulator_(Dmax, beta_focal, p_focal, 1.0, log_sum = True)), beta_focal, p_focal, 1.0, log_sum = False)

			if t_num in T_range:
			    num_index = int(T_range.index(t_num))
			    axes_list[-1][0].plot(np.log10(inf_vl + 1), label = 'f_{}_T_{}'.format(f_focal, Dmax * float(t_num)), zorder = 2, color = _color_list[num_index], linewidth = 3)
			else:
			    axes_list[-1][0].plot(np.log10(inf_vl + 1), label = 'f_{}_T_{}'.format(f_focal, Dmax * float(t_num)), zorder = 1, color = 'gray', alpha = 0.5)

	for j in range(4):
	    if j < 3:
	        for i in range(4):
	            axes_list[j][i].xaxis.set_major_formatter(ScalarFormatter(useMathText=True))
	            axes_list[j][i].set_xticks(list(np.linspace(0, beta_max, 6)))
	            axes_list[j][i].set_yticks([0, 25000, 50000, 75000, 100000])
	            axes_list[j][i].yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
	            axes_list[j][i].ticklabel_format(style="sci", axis="y", scilimits=(5,5))
	            axes_list[j][i].xaxis.offsetText.set_fontsize(12)
	            axes_list[j][i].yaxis.offsetText.set_fontsize(12)
	            
	    if j < 3:
	        for i in range(4):
	            axes_list[j][i].set_xlabel(r'$\beta$', fontsize = 20)
	            axes_list[j][i].set_ylabel('$\it{ p }$', fontsize = 20)
	            axes_list[j][i].tick_params(labelsize=17)
	    if j == 3:
	        for i in range(4):
	            axes_list[j][i].set_xlabel('Generation', fontsize = 20)
	            axes_list[j][i].set_ylabel('$R_{TP}$', fontsize = 20)
	            axes_list[j][i].set_xticks([0, 100, 200, 300])
	            axes_list[j][i].set_yticks([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
	            axes_list[j][i].set_yticklabels(['0.0', '1.0', '2.0', '3.0', '4.0', '5.0', '6.0'], fontsize = 20)
	            axes_list[j][i].tick_params(labelsize=17)

	for i in range(4):
	    axes_list[0][i].set_title('T* = {} days'.format(T_range[i]*80), fontsize = 20, pad=25)
	    axes_list[3][i].set_ylim(0.0,6.0)
	    for j in range(3):
	        axes_list[j][i].set_xlim(0, beta_max)
	        axes_list[j][i].set_ylim(0, p_max)
	        axes_list[j][i].set_xlim(- beta_max * 0.05, beta_max + beta_max * 0.05)
	        axes_list[j][i].set_ylim(- p_max * 0.05, p_max + p_max * 0.05)

	axes_list[0][0].set_ylabel('$\it{ p }$', fontsize = 20)
	axes_list[1][0].set_ylabel('Asymptomatic\n$\it{ p }$', fontsize = 20)
	axes_list[2][0].set_ylabel('Symptomatic\n$\it{ p }$', fontsize = 20)
	axes_list[3][0].set_ylabel('$R_{TP}$', fontsize = 20)

	axes_list[-1][0].set_xlabel('Days after infection', fontsize = 20)
	axes_list[-1][0].set_ylabel('Log10(viral load)', fontsize = 20)
	axes_list[-1][0].tick_params(labelsize=20)
	axes_list[-1][0].set_xticks([0, 500, 1000, 1500, 2000, 2500, 3000])
	axes_list[-1][0].set_xticklabels(['0', '5', '10', '15', '20', '25', '30'], fontsize = 20)
	axes_list[-1][0].set_xlim(-150, 2150)
	axes_list[-1][1].set_xlim(0, beta_max)
	axes_list[-1][1].set_ylim(0, p_max)
	axes_list[-1][1].set_xlabel(r'$\beta$', fontsize = 20)
	axes_list[-1][1].set_ylabel('$\it{ p }$', fontsize = 20)
	axes_list[-1][1].tick_params(labelsize=19)
	axes_list[-1][1].set_yticks([0, 25000, 50000, 75000, 100000])
	axes_list[-1][1].xaxis.set_major_formatter(ScalarFormatter(useMathText=True))
	axes_list[-1][1].yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
	axes_list[-1][1].ticklabel_format(style="sci", axis="y", scilimits=(5,5))
	axes_list[-1][1].xaxis.offsetText.set_fontsize(14)
	axes_list[-1][1].yaxis.offsetText.set_fontsize(14)


	fig.align_labels([axes_list[i][0] for i in range(5)])
	plt.tight_layout()

	return fig

def makeFigure_delta_1_0_1_5(f_num, T_range, allResults, min_max_rtp_list, gaResultsDict, peak_day_list_delta_1_0, peak_day_list_delta_1_5):
    
	global beta_max, p_max, beta_values, p_values, Beta_values, P_values

	Dmax = 80

	_color_list = ['#ff6680', '#a0d568', '#4fc1e8', '#ac92eb']

	fig = plt.figure(figsize=(16,25))

	grid_master = GridSpec(nrows = 24, ncols = 23+3, height_ratios = [1 for _ in range(24)])

	gs_1 = GridSpecFromSubplotSpec(nrows = 1, ncols = 1, subplot_spec = grid_master[0:3, 0:5])
	axes_1 = fig.add_subplot(gs_1[:, :])
	gs_2 = GridSpecFromSubplotSpec(nrows = 1, ncols = 2, subplot_spec = grid_master[0:3, 7:12])
	axes_2 = fig.add_subplot(gs_2[:, :])
	gs_3 = GridSpecFromSubplotSpec(nrows = 1, ncols = 3, subplot_spec = grid_master[0:3, 14:19])
	axes_3 = fig.add_subplot(gs_3[:, :])
	gs_4 = GridSpecFromSubplotSpec(nrows = 1, ncols = 4, subplot_spec = grid_master[0:3, 21:26])
	axes_4 = fig.add_subplot(gs_4[:, :])
	gs_5 = GridSpecFromSubplotSpec(nrows = 2, ncols = 1, subplot_spec = grid_master[4:7, 0:5])
	axes_5 = fig.add_subplot(gs_5[:, :])
	gs_6 = GridSpecFromSubplotSpec(nrows = 2, ncols = 2, subplot_spec = grid_master[4:7, 7:12])
	axes_6 = fig.add_subplot(gs_6[:, :])
	gs_7 = GridSpecFromSubplotSpec(nrows = 2, ncols = 3, subplot_spec = grid_master[4:7, 14:19])
	axes_7 = fig.add_subplot(gs_7[:, :])
	gs_8 = GridSpecFromSubplotSpec(nrows = 2, ncols = 4, subplot_spec = grid_master[4:7, 21:26])
	axes_8 = fig.add_subplot(gs_8[:, :])
	gs_9 = GridSpecFromSubplotSpec(nrows = 3, ncols = 1, subplot_spec = grid_master[8:11, 0:5])
	axes_9 = fig.add_subplot(gs_9[:, :])
	gs_10 = GridSpecFromSubplotSpec(nrows = 3, ncols = 2, subplot_spec = grid_master[8:11, 7:12])
	axes_10 = fig.add_subplot(gs_10[:, :])
	gs_11 = GridSpecFromSubplotSpec(nrows = 3, ncols = 3, subplot_spec = grid_master[8:11, 14:19])
	axes_11 = fig.add_subplot(gs_11[:, :])
	gs_12 = GridSpecFromSubplotSpec(nrows = 3, ncols = 4, subplot_spec = grid_master[8:11, 21:26])
	axes_12 = fig.add_subplot(gs_12[:, :])
	gs_13 = GridSpecFromSubplotSpec(nrows = 4, ncols = 1, subplot_spec = grid_master[12:15, 0:5])
	axes_13 = fig.add_subplot(gs_13[:, :])
	gs_14 = GridSpecFromSubplotSpec(nrows = 4, ncols = 2, subplot_spec = grid_master[12:15, 7:12])
	axes_14 = fig.add_subplot(gs_14[:, :])
	gs_15 = GridSpecFromSubplotSpec(nrows = 4, ncols = 3, subplot_spec = grid_master[12:15, 14:19])
	axes_15 = fig.add_subplot(gs_15[:, :])
	gs_16 = GridSpecFromSubplotSpec(nrows = 4, ncols = 4, subplot_spec = grid_master[12:15, 21:26])
	axes_16 = fig.add_subplot(gs_16[:, :])
	gs_17 = GridSpecFromSubplotSpec(nrows = 2, ncols = 1, subplot_spec = grid_master[17:20, 0:12])
	axes_17 = fig.add_subplot(gs_17[:, :])
	gs_18 = GridSpecFromSubplotSpec(nrows = 2, ncols = 1, subplot_spec = grid_master[17:20, 14:26])
	axes_18 = fig.add_subplot(gs_18[:, :])

	gs_19 = GridSpecFromSubplotSpec(nrows = 2, ncols = 1, subplot_spec = grid_master[21:24, 0:12])
	axes_19 = fig.add_subplot(gs_19[:, :])
	gs_20 = GridSpecFromSubplotSpec(nrows = 2, ncols = 1, subplot_spec = grid_master[21:24, 14:26])
	axes_20 = fig.add_subplot(gs_20[:, :])

	axes_list = [
	    [axes_1, axes_2, axes_3, axes_4],
	    [axes_5, axes_6, axes_7, axes_8],
	    [axes_9, axes_10, axes_11, axes_12],
	    [axes_13, axes_14, axes_15, axes_16],
	    [axes_17, axes_18],
	    [axes_19, axes_20]
	]

	divider1 = make_axes_locatable(axes_list[1][-1])
	cax1 = divider1.append_axes("right", size="5%", pad=0.1)
	divider2 = make_axes_locatable(axes_list[2][-1])
	cax2 = divider2.append_axes("right", size="5%", pad=0.1)

	min_symp_rtp, max_symp_rtp, min_asymp_rtp, max_asymp_rtp = min_max_rtp_list

	min_symp_rtp, max_symp_rtp, min_asymp_rtp, max_asymp_rtp = np.log10(min_symp_rtp+1), np.log10(max_symp_rtp+1), np.log10(min_asymp_rtp+1), np.log10(max_asymp_rtp+1)

	if f_num == 1.0:
	    min_asymp_rtp, max_asymp_rtp = min_symp_rtp, max_symp_rtp

	for t_num in T_range:

		for i in range(len(allResults)):

			f_focal, t_focal = allResults[i][0], allResults[i][1]
			if (f_focal == f_num) and (t_focal == t_num):
			    ind_result = allResults[i]
			    break

		_, _, kernel_df, beta_trajectory_mean, p_trajectory_mean, symp_rtp, asymp_rtp, ga_fitness_mean = ind_result

		symp_rtp, asymp_rtp = np.log10(np.array(symp_rtp)+1), np.log10(np.array(asymp_rtp)+1)

		num_index = int(T_range.index(t_num))

		sns.kdeplot(x=kernel_df['beta'], y=kernel_df['p'], kind="kde", color=_color_list[num_index], linewidths = 3, zorder = 1, ax=axes_list[0][num_index])
		axes_list[0][num_index].scatter(np.mean(kernel_df['beta']), np.mean(kernel_df['p']), s = 180, ec = 'black', fc = _color_list[num_index], linewidths = 1.5, zorder = 4)
		axes_list[0][num_index].scatter(kernel_df['beta'], kernel_df['p'], s = 50, ec = 'black', fc = 'white', linewidths = 1.0, zorder = 3)
		axes_list[0][num_index].plot(beta_trajectory_mean, p_trajectory_mean, color = 'black', linewidth = 1.0, zorder = 2)
	    
		try:
		    b = list(np.mean(np.array(asymp_rtp), axis = 0)).index(max(np.mean(np.array(asymp_rtp), axis = 0)))
		    b_beta = beta_values[int(b/100)]
		    b_p = p_values[b % 100]

		    axes_list[1][num_index].scatter(Beta_values, P_values, c = np.mean(np.array(asymp_rtp), axis = 0).reshape(100, 100), cmap = cmap_self, vmin=min_asymp_rtp, vmax=max_asymp_rtp)#, norm=colors.LogNorm(vmin=min_asymp_rtp,vmax=max_asymp_rtp))
		    axes_list[1][num_index].scatter(b_beta, b_p, fc = 'white', ec = 'black', s = 200, linewidth = 1.8)

		    if t_num == 10/Dmax:
		        im1 = axes_list[1][num_index].scatter(Beta_values, P_values, c = np.mean(np.array(asymp_rtp), axis = 0).reshape(100, 100), cmap = cmap_self, vmin=min_asymp_rtp, vmax=max_asymp_rtp)#, norm=colors.LogNorm(vmin=min_asymp_rtp, vmax=max_asymp_rtp))
		        axes_list[1][num_index].scatter(b_beta, b_p, fc = 'white', ec = 'black', s = 200, linewidth = 1.8)
		        cbar = fig.colorbar(im1, cax=cax1)
		        cbar.set_label('$R_{TP}$', size=20, rotation=0, labelpad=20)
		        cbar.ax.tick_params(labelsize = 17)
		        intermediate_int = []
		        intermediate_str = []
		        for j in [np.log10(i+1) for i in range(10)]:
		            if (min_asymp_rtp < j) and (j < max_asymp_rtp):
		                intermediate_int.append(j)
		                intermediate_str.append(str(int(10**j-1)) + '.0')
		        vmin_vmax_int = [round(min_asymp_rtp, 2)] + intermediate_int
		        vmin_vmax_str = [str(int(10**min_asymp_rtp-1)) + '.0'] + intermediate_str
		        cbar.ax.set_yticks(vmin_vmax_int)
		        cbar.ax.set_yticklabels(vmin_vmax_str, fontsize = 14)
		except:
		    pass
	    
		try:
		    b = list(np.mean(np.array(symp_rtp), axis = 0)).index(max(np.mean(np.array(symp_rtp), axis = 0)))
		    b_beta = beta_values[int(b/100)]
		    b_p = p_values[b % 100]

		    axes_list[2][num_index].scatter(Beta_values, P_values, c = np.mean(np.array(symp_rtp), axis = 0).reshape(100, 100), cmap = cmap_self, vmin=min_asymp_rtp, vmax=max_asymp_rtp)# , vmin=min_symp_rtp, vmax=max_symp_rtp)
		    axes_list[2][num_index].scatter(b_beta, b_p, fc = 'white', ec = 'black', s = 200, linewidth = 1.8)

		    if t_num == 10/Dmax:
		        im2 = axes_list[2][num_index].scatter(Beta_values, P_values, c = np.mean(np.array(symp_rtp), axis = 0).reshape(100, 100), cmap = cmap_self, vmin=min_asymp_rtp, vmax=max_asymp_rtp)# , vmin=min_symp_rtp, vmax=max_symp_rtp)
		        axes_list[2][num_index].scatter(b_beta, b_p, fc = 'white', ec = 'black', s = 200, linewidth = 1.8)
		        cbar = fig.colorbar(im2, cax=cax2)
		        cbar.set_label('$R_{TP}$', size=20, rotation=0, labelpad=20)
		        cbar.ax.tick_params(labelsize = 17)
		        intermediate_int = []
		        intermediate_str = []
		        
		        for j in [np.log10(i+1) for i in range(10)]:
		            if (min_asymp_rtp < j) and (j < max_asymp_rtp):
		                intermediate_int.append(j)
		                intermediate_str.append(str(int(10**j-1)) + '.0')
		        vmin_vmax_int = [round(min_asymp_rtp, 2)] + intermediate_int
		        vmin_vmax_str = [str(int(10**min_asymp_rtp-1)) + '.0'] + intermediate_str
		        cbar.ax.set_yticks(vmin_vmax_int)
		        cbar.ax.set_yticklabels(vmin_vmax_str, fontsize = 14)
		except:
		    pass

		gaIterationNum = len(ga_fitness_mean)

		for i in range(gaIterationNum):
		    axes_list[3][num_index].plot(ga_fitness_mean[i], color = _color_list[num_index], linewidth = 0.5, alpha = 0.5)
		axes_list[3][num_index].plot(np.mean(ga_fitness_mean, axis = 0), color = 'black', linewidth = 1.0)


	cont = axes_list[-2][1].contour(Beta_values, P_values, np.array(peak_day_list_delta_1_0).reshape(100, 100), [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0, 5, 10],colors='black', zorder = 2, linestyles = 'dashed')
	cont = axes_list[-2][1].contourf(Beta_values, P_values, np.array(peak_day_list_delta_1_0).reshape(100, 100), [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0, 5, 10],cmap=cmap_self2, zorder = 1)
	cont = axes_list[-2][1].contourf(Beta_values, P_values, np.array(peak_day_list_delta_1_0).reshape(100, 100), [99, 100],cmap=cmap_self3, zorder = 1)

	cont = axes_list[-1][1].contour(Beta_values, P_values, np.array(peak_day_list_delta_1_5).reshape(100, 100), [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0, 5, 10],colors='black', zorder = 2, linestyles = 'dashed')
	cont = axes_list[-1][1].contourf(Beta_values, P_values, np.array(peak_day_list_delta_1_5).reshape(100, 100), [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0, 5, 10],cmap=cmap_self2, zorder = 1)
	cont = axes_list[-1][1].contourf(Beta_values, P_values, np.array(peak_day_list_delta_1_5).reshape(100, 100), [99, 100],cmap=cmap_self3, zorder = 1)

	t_list = [allResults[i][1] for i in range(len(allResults))]

	for t_num in [10/Dmax, 9/Dmax, 8/Dmax, 7/Dmax, 6/Dmax, 5/Dmax, 4/Dmax, 3/Dmax, 2/Dmax, 1/Dmax]:
	    
		t_indices = [i for i, t in enumerate(t_list) if t == t_num]

		for i in t_indices:
			
			f_focal = allResults[i][0]

			beta_lastGen_mean = np.mean(allResults[i][2]['beta']) # kernel_df['beta']
			p_lastGen_mean = np.mean(allResults[i][2]['p']) # kernel_df['p']

			if f_focal != f_num:
			    continue

			if t_num in T_range:
			    num_index = int(T_range.index(t_num))
			    axes_list[-2][1].scatter(beta_lastGen_mean, p_lastGen_mean, ec = 'black', fc = _color_list[num_index], s = 180, linewidth = 2, zorder = 3)
			    axes_list[-1][1].scatter(beta_lastGen_mean, p_lastGen_mean, ec = 'black', fc = _color_list[num_index], s = 180, linewidth = 2, zorder = 3)
			else:
			    axes_list[-2][1].scatter(beta_lastGen_mean, p_lastGen_mean, ec = 'black', fc = 'white', s = 180, linewidth = 2, zorder = 2)
			    axes_list[-1][1].scatter(beta_lastGen_mean, p_lastGen_mean, ec = 'black', fc = 'white', s = 180, linewidth = 2, zorder = 2)

		gaResult_list = gaResultsDict[f'{t_num}']

		for ind_ga in gaResult_list:

			f_focal, beta_focal, p_focal = ind_ga[0], ind_ga[1], ind_ga[2]

			if f_focal != str(f_num):
			    continue

			inf_vl = simulator_(calc_infPeriod(simulator_(Dmax, beta_focal, p_focal, 1.0, log_sum = True)), beta_focal, p_focal, 1.0, log_sum = False)

			if t_num in T_range:
			    num_index = int(T_range.index(t_num))
			    axes_list[-2][0].plot(np.log10(inf_vl + 1), label = 'f_{}_T_{}'.format(f_focal, Dmax * float(t_num)), zorder = 2, color = _color_list[num_index], linewidth = 3)
			else:
			    axes_list[-2][0].plot(np.log10(inf_vl + 1), label = 'f_{}_T_{}'.format(f_focal, Dmax * float(t_num)), zorder = 1, color = 'gray', alpha = 0.5)

			inf_vl = simulator_(calc_infPeriod(simulator_(Dmax, beta_focal, p_focal, 1.5, log_sum = True)), beta_focal, p_focal, 1.5, log_sum = False)

			if t_num in T_range:
			    num_index = int(T_range.index(t_num))
			    axes_list[-1][0].plot(np.log10(inf_vl + 1), label = 'f_{}_T_{}'.format(f_focal, Dmax * float(t_num)), zorder = 2, color = _color_list[num_index], linewidth = 3)
			else:
			    axes_list[-1][0].plot(np.log10(inf_vl + 1), label = 'f_{}_T_{}'.format(f_focal, Dmax * float(t_num)), zorder = 1, color = 'gray', alpha = 0.5)

	for j in range(4):
	    if j < 3:
	        for i in range(4):
	            axes_list[j][i].xaxis.set_major_formatter(ScalarFormatter(useMathText=True))
	            axes_list[j][i].set_xticks(list(np.linspace(0, beta_max, 6)))
	            axes_list[j][i].set_yticks([0, 25000, 50000, 75000, 100000])
	            axes_list[j][i].yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
	            axes_list[j][i].ticklabel_format(style="sci", axis="y", scilimits=(5,5))
	            axes_list[j][i].xaxis.offsetText.set_fontsize(12)
	            axes_list[j][i].yaxis.offsetText.set_fontsize(12)
	            
	    if j < 3:
	        for i in range(4):
	            axes_list[j][i].set_xlabel(r'$\beta$', fontsize = 20)
	            axes_list[j][i].set_ylabel('$\it{ p }$', fontsize = 20)
	            axes_list[j][i].tick_params(labelsize=17)
	    if j == 3:
	        for i in range(4):
	            axes_list[j][i].set_xlabel('Generation', fontsize = 20)
	            axes_list[j][i].set_ylabel('$R_{TP}$', fontsize = 20)
	            axes_list[j][i].set_xticks([0, 100, 200, 300])
	            axes_list[j][i].set_yticks([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
	            axes_list[j][i].set_yticklabels(['0.0', '1.0', '2.0', '3.0', '4.0', '5.0', '6.0'], fontsize = 20)
	            axes_list[j][i].tick_params(labelsize=17)

	for i in range(4):
	    axes_list[0][i].set_title('T* = {} days'.format(T_range[i]*80), fontsize = 20, pad=25)
	    axes_list[3][i].set_ylim(0.0,6.0)
	    for j in range(3):
	        axes_list[j][i].set_xlim(0, beta_max)
	        axes_list[j][i].set_ylim(0, p_max)
	        axes_list[j][i].set_xlim(- beta_max * 0.05, beta_max + beta_max * 0.05)
	        axes_list[j][i].set_ylim(- p_max * 0.05, p_max + p_max * 0.05)

	axes_list[0][0].set_ylabel('$\it{ p }$', fontsize = 20)
	axes_list[1][0].set_ylabel('Asymptomatic\n$\it{ p }$', fontsize = 20)
	axes_list[2][0].set_ylabel('Symptomatic\n$\it{ p }$', fontsize = 20)
	axes_list[3][0].set_ylabel('$R_{TP}$', fontsize = 20)

	axes_list[-2][0].set_xlabel('Days after infection', fontsize = 20)
	axes_list[-2][0].set_ylabel('Log10(viral load)', fontsize = 20)
	axes_list[-2][0].tick_params(labelsize=20)
	axes_list[-2][0].set_xticks([0, 500, 1000, 1500, 2000, 2500, 3000])
	axes_list[-2][0].set_xticklabels(['0', '5', '10', '15', '20', '25', '30'], fontsize = 20)
	axes_list[-2][0].set_xlim(-90, 2090)
	axes_list[-2][1].set_xlim(0, beta_max)
	axes_list[-2][1].set_ylim(0, p_max)
	axes_list[-2][1].set_xlabel(r'$\beta$', fontsize = 20)
	axes_list[-2][1].set_ylabel('$\it{ p }$', fontsize = 20)
	axes_list[-2][1].tick_params(labelsize=19)
	axes_list[-2][1].set_yticks([0, 25000, 50000, 75000, 100000])
	axes_list[-2][1].xaxis.set_major_formatter(ScalarFormatter(useMathText=True))
	axes_list[-2][1].yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
	axes_list[-2][1].ticklabel_format(style="sci", axis="y", scilimits=(5,5))
	axes_list[-2][1].xaxis.offsetText.set_fontsize(14)
	axes_list[-2][1].yaxis.offsetText.set_fontsize(14)

	axes_list[-1][0].set_xlabel('Days after infection', fontsize = 20)
	axes_list[-1][0].set_ylabel('Log10(viral load)', fontsize = 20)
	axes_list[-1][0].tick_params(labelsize=20)
	axes_list[-1][0].set_xticks([0, 500, 1000, 1500, 2000, 2500, 3000])
	axes_list[-1][0].set_xticklabels(['0', '5', '10', '15', '20', '25', '30'], fontsize = 20)
	axes_list[-1][0].set_xlim(-150, 2150)
	axes_list[-1][1].set_xlim(0, beta_max)
	axes_list[-1][1].set_ylim(0, p_max)
	axes_list[-1][1].set_xlabel(r'$\beta$', fontsize = 20)
	axes_list[-1][1].set_ylabel('$\it{ p }$', fontsize = 20)
	axes_list[-1][1].tick_params(labelsize=19)
	axes_list[-1][1].set_yticks([0, 25000, 50000, 75000, 100000])
	axes_list[-1][1].xaxis.set_major_formatter(ScalarFormatter(useMathText=True))
	axes_list[-1][1].yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
	axes_list[-1][1].ticklabel_format(style="sci", axis="y", scilimits=(5,5))
	axes_list[-1][1].xaxis.offsetText.set_fontsize(14)
	axes_list[-1][1].yaxis.offsetText.set_fontsize(14)


	fig.align_labels([axes_list[i][0] for i in range(4)])
	plt.tight_layout()

	return fig


def min_max_rtp(f_num, T_range, all_mean_rtp_min_max_dict):
    
    min_symp_rtp_list = []
    max_symp_rtp_list = []
    min_asymp_rtp_list = []
    max_asymp_rtp_list = []

    for t_num in T_range:
        for f_t in list(all_mean_rtp_min_max_dict.keys()):
            if (str(f_num) in f_t) and (str(t_num*80) in f_t):
                try:
                    min_symp_rtp_list.append(all_mean_rtp_min_max_dict[f_t]['symp_rtp'][0])
                    max_symp_rtp_list.append(all_mean_rtp_min_max_dict[f_t]['symp_rtp'][1])
                except:
                    min_symp_rtp_list.append(0)
                    max_symp_rtp_list.append(0)
                try:
                    min_asymp_rtp_list.append(all_mean_rtp_min_max_dict[f_t]['asymp_rtp'][0])
                    max_asymp_rtp_list.append(all_mean_rtp_min_max_dict[f_t]['asymp_rtp'][1])
                except:
                    min_asymp_rtp_list.append(0)
                    max_asymp_rtp_list.append(0)
                break
                
    return [min(min_symp_rtp_list), max(max_symp_rtp_list), min(min_asymp_rtp_list), max(max_asymp_rtp_list)]

def make_gaResultsDict(T_range, allResults):

	gaResultsDict = {}

	for t_num in T_range:

		gaResultsDict[f'{t_num}'] = []

		t_list = [allResults[i][1] for i in range(len(allResults))]

		t_indices = [i for i, t in enumerate(t_list) if t == t_num]

		for i in t_indices:

			f_num = allResults[i][0]

			beta_lastGen_mean = np.mean(allResults[i][2]['beta']) # kernel_df['beta']
			p_lastGen_mean = np.mean(allResults[i][2]['p']) # kernel_df['p']

			gaResultsDict[f'{t_num}'].append([f'{f_num}', beta_lastGen_mean, p_lastGen_mean])

	return gaResultsDict

def min_max_rtp(f_num, T_range, rtp_results_dict):

	Dmax = 80

	min_symp_rtp_list = []
	max_symp_rtp_list = []
	min_asymp_rtp_list = []
	max_asymp_rtp_list = []

	for t_num in T_range:
	    for f_focal in list(rtp_results_dict.keys()):
	        if (str(f_num) in f_focal) and (str(t_num*Dmax) in f_focal):

	            try:
	                min_symp_rtp_list.append(rtp_results_dict[f_focal]['symp_rtp'][0])
	                max_symp_rtp_list.append(rtp_results_dict[f_focal]['symp_rtp'][1])
	            except:
	                min_symp_rtp_list.append(0)
	                max_symp_rtp_list.append(0)
	            try:
	                min_asymp_rtp_list.append(rtp_results_dict[f_focal]['asymp_rtp'][0])
	                max_asymp_rtp_list.append(rtp_results_dict[f_focal]['asymp_rtp'][1])
	            except:
	                min_asymp_rtp_list.append(0)
	                max_asymp_rtp_list.append(0)
	            break
	            
	return [min(min_symp_rtp_list), max(max_symp_rtp_list), min(min_asymp_rtp_list), max(max_asymp_rtp_list)]



def aggregate_results(f_num, t_num, pkl_dirs, ga_dir):

    ga_pop = joblib.load(ga_dir + '/f_{}_T_{}/pop.pkl'.format(f_num, t_num))

    iterations = ga_pop.shape[0]

    lastGeneration = -1

    # Panel A

    beta_lastGen_list = [ga_pop[i][lastGeneration, :, 0].mean() for i in range(iterations)]
    p_lastGen_list = [ga_pop[i][lastGeneration, :, 1].mean() for i in range(iterations)]
    
    kernel_df = pd.DataFrame(columns = ['beta', 'p'])
    kernel_df['beta'] = beta_lastGen_list
    kernel_df['p'] = p_lastGen_list
    
    beta_trajectory = [ga_pop[i][:, :, 0].mean(axis = 1).tolist() for i in range(iterations)]
    p_trajectory = [ga_pop[i][:, :, 1].mean(axis = 1).tolist() for i in range(iterations)]

    beta_trajectory_mean = np.array(beta_trajectory).mean(axis = 0)
    p_trajectory_mean = np.array(p_trajectory).mean(axis = 0)    
    
    for _dir in pkl_dirs:

        _num = re.findall(r"\d+", _dir)[-4:]
 
        f_focal, t_focal = float(_num[0] + '.' + _num[1]), float(_num[2] + '.' + _num[3])

        if (f_focal == f_num) and (t_focal == t_num):
            break

    ind_rtp = joblib.load(_dir + '/ind_rtp.pkl')
    inc_sample = joblib.load(_dir + '/ind_inc_sample.pkl')

    # Panel B

    symp_rtp = []
    asymp_rtp = []

    for i in range(len(ind_rtp)):

        if i in inc_sample:
            symp_rtp.append(ind_rtp[i])
        else:
            asymp_rtp.append(ind_rtp[i])

    symp_rtp = np.array(symp_rtp)
    asymp_rtp = np.array(asymp_rtp)

    # Panel C

    ga_fitness_mean = np.array([ga_pop[i][:, :, 2].mean(axis = 1).tolist() for i in range(iterations)])


    return [f_num, t_num, kernel_df, beta_trajectory_mean, p_trajectory_mean, symp_rtp, asymp_rtp, ga_fitness_mean]

def collectDirs(f_range, T_range, pkl_dirs, ga_dir):

	args = []

	for f_num in f_range:
	    for t_num in T_range:

	        args.append([f_num, t_num, pkl_dirs, ga_dir])

	return args

def make_rtpResultsDict(allResults):

	Dmax = 80

	rtpResultsDict = {}

	for i in range(len(allResults)):

	    indResult = allResults[i]

	    f_focal, t_focal, symp_rtp, asymp_rtp = indResult[0], indResult[1]*Dmax, indResult[5], indResult[6]

	    rtpResultsDict[f'f : {f_focal}, t : {t_focal}'] = {'symp_rtp':[], 'asymp_rtp':[]}

	    if len(symp_rtp) != 0:

	        rtpResultsDict[f'f : {f_focal}, t : {t_focal}']['symp_rtp'] = [min(symp_rtp.mean(axis = 0)), max(symp_rtp.mean(axis = 0))]

	    if len(asymp_rtp) != 0:

	        rtpResultsDict[f'f : {f_focal}, t : {t_focal}']['asymp_rtp'] = [min(asymp_rtp.mean(axis = 0)), max(asymp_rtp.mean(axis = 0))]

	return rtpResultsDict

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

def simulator_(tau_max, beta, p, n_delta, log_sum = True):

    F0 = 1
    V0 = 0.01
    var = [F0, V0]
    dt = 0.01
    t_ = np.arange(0, tau_max, dt)
    
    sim_results = odeint(ODE_, var, t_, args = (n_delta, beta, p))

    vl_ = sim_results[:, 1]
    
    if log_sum == True:
        vl_Log10 = sum(np.log10(vl_ + 1)) / (1/dt)
        return vl_Log10
    
    else:
        return vl_
    
def calc_infPeriod(vi):
    
    D50 = 198.58
    Dk = 1.17
    Dmax = 80

    return Dmax*(vi**Dk/(vi**Dk+D50**Dk))

def calc_peak_day(dt=0.01):

	global beta_max, p_max

	Dmax = 80

	peak_day_list_delta_1_0 = []

	for i in np.linspace(0.0, beta_max, 100): # beta
	    for j in np.linspace(0, p_max, 100): # p

	        inf_vl = simulator_(calc_infPeriod(simulator_(Dmax, i, j, 1.0, log_sum = True)), i, j, 1.0, log_sum = False)

	        if sum((1.5 * 1.5)/(1/dt) * np.array(inf_vl)) < 1.0:
	            peak_day_list_delta_1_0.append(100) # white-colored
	            continue

	        peak_day = list(inf_vl).index(max(inf_vl)) * dt

	        if peak_day < 10:
	            peak_day_list_delta_1_0.append(peak_day)
	        else:
	            peak_day_list_delta_1_0.append(11) # gray-colored

	peak_day_list_delta_1_5 = []

	for i in np.linspace(0.0, beta_max, 100): # beta
	    for j in np.linspace(0, p_max, 100): # p

	        inf_vl = simulator_(calc_infPeriod(simulator_(Dmax, i, j, 1.5, log_sum = True)), i, j, 1.5, log_sum = False)

	        if sum((1.5 * 1.5)/(1/dt) * np.array(inf_vl)) < 1.0:
	            peak_day_list_delta_1_5.append(100) # white-colored
	            continue

	        peak_day = list(inf_vl).index(max(inf_vl)) * dt

	        if peak_day < 10:
	            peak_day_list_delta_1_5.append(peak_day)
	        else:
	            peak_day_list_delta_1_5.append(11) # gray-colored
	            
	return peak_day_list_delta_1_0, peak_day_list_delta_1_5

def wrapper_aggregate_results(args):
    return aggregate_results(*args)

def showFigure(vaccination_rate, imcomplete_isolation_prob, f_range, chosen_T_range, n_processes):

	global beta_max, p_max

	Dmax = 80

	current_dir = "results/vac_rate_{}_ii_prob_{}".format(vaccination_rate, imcomplete_isolation_prob)

	if not os.path.isdir(current_dir):
		print('Error : Specified directory does not exists; /vac_rate_{}_ii_prob_{}'.format(vaccination_rate, imcomplete_isolation_prob))
		print('        You need to execute "main.py" and "GA.py" with the parameter set of vaccination_rate, imcomplete_isolation_prob beforehand.')
		sys.exit()

	pkl_dirs = glob.glob('{}/*'.format(current_dir + '/shape_1.5_scale_1.5'))
	ga_dir = current_dir + '/ga'

	# default T_range
	T_range = [1/Dmax, 2/Dmax, 3/Dmax, 4/Dmax, 5/Dmax, 6/Dmax, 7/Dmax, 8/Dmax, 9/Dmax, 10/Dmax]

	resultsDirsList = collectDirs(f_range, T_range, pkl_dirs, ga_dir)

	p = Pool(processes=n_processes)
	allResults = p.map(wrapper_aggregate_results, resultsDirsList)
	p.close()

	rtpResultsDict = make_rtpResultsDict(allResults)

	gaResultsDict = make_gaResultsDict(T_range, allResults)

	peak_day_list_delta_1_0, peak_day_list_delta_1_5 = calc_peak_day()

	for f_num in f_range:

		print(f'f : {f_num}', flush = True)

		min_max_rtp_list = min_max_rtp(f_num, T_range, rtpResultsDict)

		if vaccination_rate == 0.0:

			makeFigure_delta_1_0(f_num, chosen_T_range, allResults, min_max_rtp_list, gaResultsDict, peak_day_list_delta_1_0)
	    
		else:

			makeFigure_delta_1_0_1_5(f_num, chosen_T_range, allResults, min_max_rtp_list, gaResultsDict, peak_day_list_delta_1_0, peak_day_list_delta_1_5)

