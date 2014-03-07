

import argparse
from multiprocessing import Pool

import simulated_annealing

def main():
	parser = argparse.ArgumentParser(description="Runner for TFPyModeller")
	parser.add_argument('-t','--target', help='target name to use',default='T9999')
	parser.add_argument('-d','--decoys', help='number of decoys to create',type=int,default=10)
	parser.add_argument('-n','--simulations',help='number of simulations to run',type=int,default='1000')
	t_group = parser.add_mutually_exclusive_group()
	t_group.add_argument('-l','--linear',action='store_false')
	t_group.add_argument('-s','--sigmoid',action='store_true')
	parser.add_argument('-k','--sequence',help='sequence to run')
	
	args = parser.parse_args()

	target = args.target
	decoys = args.decoys
	simulations = args.simulations
	temperature = 'sigmoid'
	if not args.sigmoid:
		temperature = 'linear'
	sequence = args.sequence

	pool = Pool(processes=decoys)

	runargs = []

	for i in range(decoys):
		itarget = "%s-%d" % (target,i)
		runargs.append((itarget,sequence,simulations,temperature))

	pool.map(simulated_annealing.runargs,runargs)
	

if __name__ == '__main__':
	main()