import os
import re
import sys
import glob
import random
import joblib
import warnings
import numpy as np
import argparse
import datetime

from multiprocessing import Pool

warnings.simplefilter('ignore')

beta_max = 0.0000005
p_max = 100000

class Ind():

	def __init__(self, beta_, p_):

		self.beta_ = beta_
		self.p_ = p_
		self.fitness = 0
		self.rank = 0

	def _fitness(self, fitness):

		self.fitness = fitness

	def _rank(self, rank):

		self.rank = rank

class GA():

	def __init__(self, popSize, numElitism, pselection, pcrossover, pmutation, fit_dir):

		self.popSize = popSize
		self.gen = 0
		self.pop = []
		self.numElitism = numElitism
		self.pselection = pselection
		self.pcrossover = pcrossover
		self.pmutation = pmutation
		self.fit_ = joblib.load(fit_dir + '/ind_rtp.pkl')
		self.trans_pop_ = joblib.load(fit_dir + '/ind_pop.pkl')
		self.npiNum_ = joblib.load(fit_dir + '/ind_inc_sample.pkl')
		num_ = re.findall(r"\d+", fit_dir)[-4:]
		f_, t_ = float(num_[0] + '.' + num_[1]), float(num_[2] + '.' + num_[3])
		self.f_t_ = [f_, t_]
		self.maxGen = len(self.trans_pop_)

	def popInit(self):

		self.gen = 0

		gen = self.gen

		pop = []

		for i in range(self.popSize):

			ind_beta = random.uniform(0, beta_max)
			ind_p = random.uniform(0, p_max)
			newInd = Ind(ind_beta, ind_p)

			pop.append(newInd)

		self.pop = [pop]
		self._assign_fitness(pop)
		self._assign_rank(pop)
		self.gen = 1

	def popEvolve(self):

		for j in range(1, self.maxGen):

			gen = self.gen
			genPop = self.pop[-1]

			pop = []

			# elitism
			for i in range(self.numElitism):
				for _Ind in genPop:
					if _Ind.rank == i:
						pop.append(_Ind)
						break

			while len(pop) < self.popSize:
				newInd = self._generate_newInd()
				pop.append(newInd)

			self._assign_fitness(pop)
			self._assign_rank(pop)

			self.pop.append(pop)
			self.gen += 1

	def _generate_newInd(self):

		# selection
		NumPaCandidates = int(self.popSize * self.pselection)
		totalFitness, selectedPop = self._select_pop(NumPaCandidates)

		# parent selection
		paA = self._select_roulette_wheel(totalFitness, selectedPop)
		paB = self._select_roulette_wheel(totalFitness, selectedPop)

		# crossover
		newInd = self._crossover(paA, paB)

		# mutation
		newInd = self._mutation(newInd)

		return newInd

	def _select_roulette_wheel(self, totalFitness, selectedPop):

		pick = random.random()
		offset = 0.0

		if totalFitness == 0:
			popNum = len(selectedPop)
			selectedNum = random.randint(0, popNum-1)
			pa_ = selectedPop[selectedNum]

			return pa_

		for _Ind in selectedPop:
			pa_ = _Ind

			offset += _Ind.fitness/totalFitness

			if offset >= pick:
				pa_ = _Ind
				break

		return pa_

	def _select_pop(self, NumPaCandidates):

		genPop = self.pop[-1]

		selectedPop = []

		for _Ind in genPop:

			if _Ind.rank <= NumPaCandidates:
				selectedPop.append(_Ind)

		totalFitness = sum([_Ind.fitness for _Ind in selectedPop])

		return totalFitness, selectedPop

	def _crossover(self, paA, paB):

		prob = self.pcrossover

		if prob >= random.random():

			if 0.5 <= random.random():
				beta_ = paA.beta_
			else:
				beta_ = paB.beta_

			if 0.5 <= random.random():
				p_ = paA.p_
			else:
				p_ = paB.p_

			newInd = Ind(beta_, p_)
				
		else:
			if 0.5 <= random.random():
				newInd = Ind(paA.beta_, paA.p_)
			else:
				newInd = Ind(paB.beta_, paB.p_)

		return newInd

	def _mutation(self, _Ind):

		prob = self.pmutation
		beta_ = _Ind.beta_
		p_ = _Ind.p_

		beta_mut_range = beta_max * 0.05
		p_mut_range = p_max * 0.05


		if prob >= random.random():

			if 0.5 <= random.random():
				beta_ += beta_mut_range * random.random()
			else:
				beta_ -= beta_mut_range * random.random()

		if beta_ >= beta_max:
			beta_ = beta_max

		if beta_ <= 0:
			beta_ = 0

		if prob >= random.random():

			if 0.5 <= random.random():
				p_ += p_mut_range * random.random()
			else:
				p_ -= p_mut_range * random.random()

		if p_ >= p_max:
			p_ = p_max

		if p_ <= 0:
			p_ = 0

		_Ind.beta_ = beta_
		_Ind.p_ = p_

		return _Ind

	def _assign_fitness(self, genPop):

		fit_id = self.gen
		fit_landscape = np.array(self.fit_[fit_id]).reshape(100, 100)
		trans_pop = self.trans_pop_[fit_id]
		if fit_id in self.npiNum_:
			T_star_day = 80*self.f_t_[1]
		else:
			T_star_day = 80

		def _fitFromMap(fit_landscape, trans_pop, T_star_day, _beta, _p):

			beta_rn = np.linspace(0.0, beta_max, 100) # beta
			p_rn = np.linspace(0, p_max, 100) # p

			try:
				below_b = beta_rn[_beta > beta_rn][-1]
			except:
				below_b = beta_rn[0]

			try:
				on_b = beta_rn[_beta < beta_rn][0]
			except:
				on_b = beta_rn[-1]

			try:
				below_p = p_rn[_p > p_rn][-1]
			except:
				below_p = p_rn[0]

			try:
				on_p = p_rn[_p < p_rn][0]
			except:
				on_p = p_rn[-1]

			fitList = []

			for i in [below_b, on_b]:
				b_i = beta_rn.tolist().index(i)

				for j in [below_p, on_p]:
					p_i = p_rn.tolist().index(j)

					indFit = fit_landscape[b_i, p_i]
					fitList.append(indFit)

			return np.mean(fitList)

		for _Ind in genPop:

			_beta = _Ind.beta_
			_p = _Ind.p_

			_fit = _fitFromMap(fit_landscape, trans_pop, T_star_day, _beta, _p)

			_Ind._fitness(_fit)


	def _assign_rank(self, genPop):

		fitList = [_Ind.fitness for _Ind in genPop]
		sorted_fitList = sorted(fitList, reverse = True)

		for _Ind, _Fit in zip(genPop, fitList):

			rank_ = sorted_fitList.index(_Fit)
			sorted_fitList[rank_] = -1

			_Ind._rank(rank_)

	def save_pop_as_list(self):

		pop_list = []

		for i in range(self.gen):

			ind_gen_ = []

			for j in range(self.popSize):

				beta_ = self.pop[i][j].beta_
				p_ = self.pop[i][j].p_
				fit_ = self.pop[i][j].fitness
				rank_ = self.pop[i][j].rank

				ind_gen_.append([beta_, p_, fit_, rank_])

			pop_list.append(ind_gen_)

		return pop_list

def aggregate_info(default_dir, ga_dir, f_range, T_range, popSize, numElitism, pselection, pcrossover, pmutation, n_iter):

	save_dir = default_dir + ga_dir

	try:
		os.makedirs(save_dir)
	except:
		print('Error : Specified directory does not exists; {}'.format(save_dir))
		print('        You need to execute "main.py" with the parameter set of vaccination_rate, imcomplete_isolation_prob beforehand.')
		sys.exit()

	z = []

	raw_data_dirs = glob.glob('{}/shape_1.5_scale_1.5/*'.format(default_dir))

	if len(raw_data_dirs) == 0:
		print('Error : You need to execute "main.py" with the parameter set of vaccination_rate, imcomplete_isolation_prob beforehand.')
		sys.exit()

	for fit_dir in raw_data_dirs:

		num_ = re.findall(r"\d+", fit_dir)[-4:]

		f_, t_ = float(num_[0] + '.' + num_[1]), float(num_[2] + '.' + num_[3])

		save_dir_path = save_dir + '/f_{}_T_{}'.format(f_, t_)

		try:
			os.mkdir(save_dir_path)
		except:
			print('Error : Save directory already exists; /f_{}_T_{}'.format(f_, t_))
			sys.exit()

		z.append([fit_dir, save_dir_path, f_, t_, popSize, numElitism, pselection, pcrossover, pmutation, n_iter])


	return z


def calc_iteration(fit_dir, save_dir_path, f_, t_, popSize, numElitism, pselection, pcrossover, pmutation, n_iter):

	pop_list = []

	single_ga = GA(popSize, numElitism, pselection, pcrossover, pmutation, fit_dir)

	for i in range(n_iter):

		if i % 50 == 0:
			print('f:{}, T:{}, gen:{}'.format(f_, t_, i), flush = True)

		single_ga.popInit()
		single_ga.popEvolve()
		pop_ = single_ga.save_pop_as_list()
		pop_list.append(pop_)

	joblib.dump(np.array(pop_list), '{}/pop.pkl'.format(save_dir_path), compress = True)

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
	n_iter = args.iteration
	f_range = args.f_range
	vaccination_rate = args.vaccination_rate
	imcomplete_isolation_prob = args.imcomplete_isolation_prob
	n_processes = args.n_processes

	print(f'vaccination_rate : {vaccination_rate}, imcomplete_isolation_prob : {imcomplete_isolation_prob}, n_processes : {n_processes}')

	Dmax = 80

	default_dir = 'results/vac_rate_{}_ii_prob_{}'.format(vaccination_rate, imcomplete_isolation_prob)

	ga_dir = '/ga'

	# default T_range
	T_range = [1/Dmax, 2/Dmax, 3/Dmax, 4/Dmax, 5/Dmax, 6/Dmax, 7/Dmax, 8/Dmax, 9/Dmax, 10/Dmax]

	popSize = 100
	numElitism = 2
	pselection = 1.0
	pcrossover = 0.2
	pmutation = 0.7

	args = aggregate_info(default_dir, ga_dir, f_range, T_range, popSize, numElitism, pselection, pcrossover, pmutation, n_iter)

	print('start at ', datetime.datetime.now(), flush = True)

	p = Pool(processes=n_processes)
	p.map(wrapper_calc_iteration, args)
	p.close()

	print('finish at ', datetime.datetime.now(), flush = True)
