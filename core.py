#!/usr/bin/python
import matplotlib.pyplot as plt
import numpy as np
import os
import math

class global_settings(object):

	def __init__(self):
		a = open('settings.txt', 'r')
		a = a.read()
		a = a.split('\n')
		self.Mu = float(a[1])
		self.original_ms_command = a[3]
		self.psmc_command = a[5]
		line = a[5].split(' ')
		self.psmc_N = line[line.index('-N')+1]
		axes = a[7].split(',')
		self.x_min = float(axes[0])
		self.x_max = float(axes[1])
		self.y_min = float(axes[2])
		self.y_max = float(axes[3])
		self.scaling_factor = float(a[9])
		self.s = int(a[11])
		self.number_of_experiments = int(a[13])
		self.generation_time = int(a[15])
 
class converter(object):
	def format_string(self, s, l=60):
		number_of_lines = len(s)/l
		temp = ''
		for i in range(number_of_lines):
			temp += s[(i*l):((i+1)*l)] + '\n'
		if (len(s)%l) != 0:
			temp += s[number_of_lines*l:] + '\n'
		return temp

	def write_results(self, number_of_hap, number_of_simulations, ms_command, psmc_N):
		res = open('results.txt', 'w')
		res.write('# MS command:\n')
		res.write(ms_command+'\n')
		if (number_of_hap > 2):
			for A in range(number_of_hap):
				for B in range(A+1, number_of_hap):
					res.write('# Infered History:\n')
					new_hist = os.popen('./psmc2history.pl -n '+ psmc_N +' ./output'+A.__str__()+B.__str__()+'.psmc | ./history2ms.pl').read()
					res.write(new_hist)
		if (number_of_hap==2):
			for i in range(number_of_simulations):
				res.write('# Infered History:\n')
				new_hist = os.popen('./psmc2history.pl -n '+ psmc_N +' ./output'+(i+1).__str__()+'.psmc | ./history2ms.pl').read()
				res.write(new_hist)
		res.close()

	def ms2psmcfa(self, input_filename, ouput_filename, chrA=0, chrB=1, s=100):
		a = open(input_filename)
		b = a.read()
		a.close()
		b = b.split('\n')
		psmcfa_temp = ''
		num_of_sites = int(b[3])/s + 1
		actual_pos = 0
		for line in b[4:-2]:
			line = line.split('\t')
			het_pos = int(line[0])/s + 1
			if (het_pos > actual_pos) and (line[1][chrA] <> line[1][chrB]):
				number_of_hom = het_pos - actual_pos -1
				psmcfa_temp += 'T' * number_of_hom
				psmcfa_temp += 'K'
				actual_pos = het_pos
		number_of_het = num_of_sites - actual_pos
		psmcfa_temp += 'T' * number_of_het
		new_psmcfa = self.format_string(psmcfa_temp)
		b = open(ouput_filename+chrA.__str__()+chrB.__str__()+'.psmcfa', 'w')
		b.write('>1\n')
		b.write(new_psmcfa)
		b.close()

	def ms2fasta(self, ms_output):
		ms_output = ms_output.split('\n')
		sequence_length = int(ms_output[3])
		
		# write the reference
		a = open('reference.fasta', 'w')
		a.write('>chrom1\n')
		a.write('A'*sequence_length+'\n')
		a.close()
		
		# search heterozygous positions
		het_positions = []
		ms_command = ms_output[0]
		ms_command = ms_command.split(' ')
		number_of_hap = int(ms_command[1])
		for i in range(number_of_hap):
			het_positions.append([])
		for snip in ms_output[4:-2]:
			snip = snip.split('\t')
			snip_position = int(snip[0])
			values = snip[1]
			# If there is a '1', then the haplotype has a mutation
			for j in range(number_of_hap):
				if values[j]=='1': het_positions[j].append(snip_position)
		# Now, we write the file with the haplotypes
		a = open('data.fasta', 'w')
		for i in range(number_of_hap):
			a.write('>'+i.__str__()+'\n')
			temp_hap = ['A']*sequence_length
			for position in het_positions[i]:
				temp_hap[position]='C'
			a.write(''.join(temp_hap)+'\n')
		a.close()

class plotter(object):

	def scale_vector(self, v, scaling_factor):
		temp = v
		for i in range(len(v)):
			temp[i] = float(v[i])/scaling_factor
		return temp

	def MS2fun(self, commands, x_min=1, x_max=1000000, scaling_factor=1, Mu=2.5*10**(-8), demographic_event='-eN'):
		'''
		Recieves the ms commands in a string line. Returns two arrays containing 
		the values for the demographic history. The x_min and x_max are minimum
		and maximum values for x
		'''
		commands = commands.split(' ')
		if demographic_event == '-eN': j =0
		else : j=1
		#search for theta
		theta = commands[commands.index('-t')+1]
		#search L
		L = commands[commands.index('-r')+2]
		# N0 = theta/(4*Mu*L)
		theta = float(theta)
		L = int(L)
		N0 = float(theta)/(4*Mu*L)
		# search for -eN option
		x_temp = []
		y_temp = []
		# Find all the parameters 'eN' for the demographic change
		i=-1
		# getting demographic events
		while 1:
			try:
				i = commands.index(demographic_event, i+1)
				time_of_event = float(commands[i+1])
				x_temp.append(time_of_event*4*N0)
				new_lambda = commands[i+2+j]
				new_lambda = float(new_lambda)
				new_e_size = new_lambda*N0
				y_temp.append(new_e_size)
			except:
				break
		# writing the values in the right way
		x = [x_min,]
		y = [N0,]
		for i in range(len(x_temp)):
			x.append(x_temp[i])
			y.append(y[-1])
			x.append(x[-1])
			y.append(y_temp[i])
		x = self.scale_vector(x, scaling_factor)
		y = self.scale_vector(y, scaling_factor)
		x.append(x_max)
		y.append(y[-1])		
		return [x, y, N0]

	def hist2fun(self, history, x_min=1, x_max=1000000, scaling_factor=1, Mu=2.5*10**(-8), s=100):
		history = history.split('\n')
		theta = float(history[0].split(' ')[1])
		N0 = theta/(4*Mu*s)
		x = [x_min,]
		y = [N0*float(history[3].split(' ')[1].split('\t')[1]),]
		for line in history[4:-1]:
			[new_x, new_y]  = line.split(' ')[1].split('\t')
			x.append(float(new_x)*2*N0)
			y.append(y[-1])
			x.append(x[-1])
			y.append(float(new_y)*N0)
		x = self.scale_vector(x, scaling_factor)
		y = self.scale_vector(y, scaling_factor)
		x.append(max(x[-1]*2, x_max))
		y.append(y[-1])
		return [x, y]

	def plot_fun(self, x, y, x_min, x_max, y_min, y_max, style):
		fig = plt.figure()
		plt.axis([x_min, x_max, y_min, y_max])
		ax = fig.add_subplot(1, 1, 1)
		ax.plot(x, y, 'g')
		ax.set_xscale('log')
		plt.savefig('plot.png')
	  
	def plot_results(self, results, scaling_factor, Mu, s, x_min, x_max, y_min, y_max, generation_time=1):
		'''
		This function plots the values in the file "results.txt".
		The file "results.txt" contains the ms-command and the infered
		histories.
		The time is in generations by default (i.e. 1 generation = 1 year)
		'''
		results = results.split('# Infered History:\n')
		original_ms_command = results[0].split('\n')[1]
		[x_original, y_original, N0_original] = self.MS2fun(original_ms_command, x_min, x_max, scaling_factor, Mu, '-eN')
		fig = plt.figure()
		if generation_time <> 1:
			for i in range(1,len(x_original)): x_original[i]*=generation_time
		plt.axis([x_min, x_max, y_min, y_max])
		ax = fig.add_subplot(1, 1, 1)
		for infered_history in results[1:]:
			[x, y] = self.hist2fun(infered_history, x_min, x_max, scaling_factor, Mu, s)
			if generation_time <> 1:
				for i in range(1, len(x)): x[i]*=generation_time
			ax.plot(x, y, 'g')
		ax.plot(x_original, y_original, 'r--')
		ax.set_xscale('log')
		plt.savefig('plot.png')

class accuracy(object):

	def evaluate_func(self, x, y, t):
		'''
		Computes the value f(t) when the function f is picewise constant
		and given by the pair of lists [x,y]
		'''
		x_arr = np.array(x)
		return (y[sum(x_arr<=t)-1])	
	
	def compute_errors(self, t0, t1, results, base=math.e):
		'''
		This method computes the accuracy of the PSMC estimate, according to
		the formula given in the Supplementary Information. The accuracy is
		computed in the time intervall t0, t1 and the time is in number of
		generations
		The input parameter "results" is a text containing the ms-command
		and some infered histories
		The output is a vector containing the accuracy for every infered 
		history.
		'''
		p = plotter()
		results = results.split('# Infered History:\n')
		original_ms_command = results[0].split('\n')[1]
		# The original function is N0
		[x_original, y_original, N0_original] = p.MS2fun(original_ms_command)
		x_temp = x_original[0:1] + [x_original[2*i+1] for i in range(len(x_original)/2)]
		x_original = x_temp
		y_temp = y_original[0:1] + [y_original[2*i+1] for i in range(len(y_original)/2)]
		y_original = y_temp
		errors_array= []
		for i in range(1, len(results)):
			errors_array.append(self.compute_accuracy(x_original, y_original, t0, t1, results[i], base))
		return errors_array

	def compute_accuracy(self, x_original, y_original, t0, t1, h, base=math.e):
		p = plotter()
		# The estimated function is N1
		[x, y] = p.hist2fun(history=h)
		x_temp = x[0:1] + [x[2*i+1] for i in range(len(x)/2)]
		x = x_temp
		y_temp = y[0:1] + [y[2*i+1] for i in range(len(y)/2)]
		y = y_temp
		# Now we merge the two vectors x_original and x, then we add t0 and t1
		x_temp = x_original + x
		x_temp = set(x_temp) # removing duplicates
		x_temp = list(x_temp)
		x_temp.sort()
		#print x_temp
		inf = np.array(x_temp)>t0
		sup = np.array(x_temp)<t1
		inf_sup = inf*sup
		x_temp = set(inf_sup * np.array(x_temp))
		x_temp = list(x_temp)
		x_temp[0] = t0
		x_temp.append(t1)
		x_vect = x_temp  # Now we have the vector. We can compute the formula
		sum_temp = 0
		for i in range(len(x_vect)-1):
			N0_i = float(self.evaluate_func(x_original, y_original, x_vect[i]))
			N1_i = float(self.evaluate_func(x, y, x_vect[i]))
			sum_temp += (abs(N0_i-N1_i)/(N0_i+N1_i))*(math.log(x_vect[i+1], base)-math.log(x_vect[i], base))
		result = float(sum_temp)/(math.log(t1, base)-math.log(t0, base))
		return result

  
