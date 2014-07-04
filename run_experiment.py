#!/usr/bin/python
import os
import core
from core import *
import time

def main():
	os.system('rm ./jobs/*job*')
	os.system('rm ./temp_files/*')
	os.system('rm ./joblogs/*')
	g_s = global_settings()
	all_done = False
	params = open('./jobs/cluster_params', 'r')
	params = params.read()
	for i in range(1, g_s.number_of_experiments+1):
		ms_result_filename = './temp_files/simulation_'+i.__str__()+'.ms'
		psmc_in_filename = './temp_files/input_'+i.__str__()+'.psmcfa'
		psmc_out_filename = './temp_files/output_'+i.__str__()+'.psmc'
		newjob = open('./jobs/job_'+i.__str__()+'.sh', 'w')
		newjob.write(params.replace('job_0', 'job_'+i.__str__()))
		# Run the MS command
		newjob.write('./utils/'+g_s.original_ms_command + ' > '+ms_result_filename+'\n')
		# Converting the output into psmcfa
		newjob.write('./utils/ms2psmcfa.pl ' + ms_result_filename + ' > ' + psmc_in_filename+'\n')
		# Apply the psmc
		newjob.write('./utils/'+g_s.psmc_command+' -o ' + psmc_out_filename + ' ' + psmc_in_filename+'\n')
		# Notify that it's done
		newjob.write('echo "DONE" > ./jobs/finished_job_'+i.__str__()+'\n')
		newjob.close()
		print os.popen('qsub ./jobs/job_'+i.__str__()+'.sh').read()
		time.sleep(2)
	
	while not all_done:
		job_finished = int(os.popen('ls ./jobs | grep finished_job_ | wc -l').read())
		if job_finished == g_s.number_of_experiments: all_done = True
		time.sleep(60)
	
	os.system('./plot_psmc.py')

if __name__ == "__main__":
	main()


  
