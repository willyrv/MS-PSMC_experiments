#!/usr/bin/python
import os
import core_w_plot
from core_w_plot import *

def main():
	g_s = global_settings()
	
	results = open('./temp_files/results.txt', 'r')
	results = results.read()
	p = plotter()
	p.plot_results(results, g_s.scaling_factor, g_s.Mu, g_s.s, g_s.x_min, g_s.x_max, g_s.y_min, g_s.y_max)

if __name__ == "__main__":
	main()
