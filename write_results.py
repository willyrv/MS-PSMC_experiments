#!/usr/bin/python
import os
import core_no_plot
from core_no_plot import *

def main():
	g_s = global_settings()
	res = open('./temp_files/results.txt', 'w')
	res.write('# MS command:\n')
	res.write(g_s.original_ms_command+'\n')
	for i in range(1, g_s.number_of_experiments+1):
		res.write('# Infered History:\n')
		new_hist = os.popen('./utils/psmc2history.pl -n '+ g_s.psmc_N +' ./temp_files/output_'+i.__str__()+'.psmc').read()
		res.write(new_hist)
	res.close()

if __name__ == "__main__":
	main()
