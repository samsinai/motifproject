# MotifMainArg.py: Main Function for comparing control and test data with command-line arguments
# Grant Kinsler and Nicholas Keone Lee
# Math 243, Spring 2014
# Authored: 2014/04/03
# Last updated: 2015/02/21, GRK

import sys
import getopt
import itertools
import collections
import numpy

def usage():
	print "Running a Motif Simulation using the parameters provided in the command-line separated by spaces\n"
	print "PARAMETERS:"
	print "number of trials; growthIterations; maxStrands; maxStrandLength; numCells; numRounds; motif; p_elong; p_elong_motif; motif_add; base_number; p_divide\n"
	print "There are two CSVs produced:"
	print "1) The resulting CSV has the first row as a list of these parameters in order, followed by data for each trial in the order: prop. of motif strands, prop. of total strands, prop. of motif cells\n"
	print "The last three rows are the averages of these throughout the number trials specified\n"
	print "2) -with CellData in the title- The CSV has the first row as a list of the parameters in order, followed by data for the contents of each cell in each round in of the first trial, and then the keys of \n"
	print "For instance, after the first row of parameters, the first numCells rows depict the cells in the first round in the first trial\n"
	print "Following this, there is a row with the each strand, under which is the the frequency of the strand/cell in each round (row1:round 1's frequencies, row2:round 2's freq, etc."

def makeKeyorder(maxStrandLength,motif): # create the keyorder for the order of our dictionary
	
	keyorder = []

	motifcounter = 0

	for n in range(maxStrandLength): # create order of keys for the ordered dictionary
		for key in itertools.product(range(2),repeat = n+1):
			mod_key = str(key).strip(" ,(),','").replace(", ", "")
			if motif in mod_key:
				keyorder.insert(motifcounter,mod_key) # place all motif strands first in the order
				motifcounter = motifcounter + 1
			else:
				keyorder.append(mod_key) # place the rest in the order they are created

	return keyorder


def flatten(items, seqtypes=(list, tuple)): # used for flattening lists
    for i, x in enumerate(items):
        while isinstance(items[i], seqtypes):
            items[i:i+1] = items[i]
    return items


def main(argv):
	import MotifFunctionDivisionUpdate as source
	import csv
	
	try:
		opts, args = getopt.getopt(argv, "h", ["help"])
	except getopt.GetoptError, error:
		sys.stderr.write(str(error)+"\n")
		usage()
		sys.exit(2)
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			usage()
			sys.exit()
		else:
			sys.stderr.write("Unknown option %s\n" %opt)
			usage()
			sys.exit(2)

	trials = int(args[0]) # number of trials we want to run

	growthIterations = int(args[1]) # number of growth interations per round
	maxStrands = int(args[2]) # max number of possible strands in cell
	maxStrandLength = int(args[3]) # max length of any individual strand
	numCells = int(args[4]) # total number of cells
	numRounds = int(args[5]) # total number of rounds

	motif = args[6]

	p_elong = float(args[7]) # probability that an existing strand will be elongated
	p_elong_motif = float(args[8]) # probability existings strand is elongated in presence of motif
	motif_add = float(args[9]) # probability that a zero is added in presence of motif (when a base is added)

	base_number = int(args[10]) # total number of nucleotides in strands of cell needed for division
	p_divide = args[11]

	parameterlist = [trials, growthIterations, maxStrands, maxStrandLength, numCells, numRounds, repr(motif), p_elong, p_elong_motif, motif_add, base_number, p_divide]

	# test_cell_ag = [None] * numRounds # initialize aggregate list
	# test_strand_ag = [None] * numRounds # initialize aggregate list
	# test_strand_count_ag = [None] * numRounds # initialize aggregate list
	# test_cell_ag = [] # initialize aggregate list
	# test_strand_ag = [] # initialize aggregate list
	# test_strand_count_ag = [] # initialize aggregate list
	test_strand_normavg = []
	test_strand_SD = []
	test_strand_count_normavg = []
	test_strand_count_SD = []
	test_cell_normavg = []
	test_cell_SD = []
	total_strand_count = [0] * numRounds

	celltracker_total = [] # intialize total list for celltracker

	for i in range(trials):
		celltracker_total.append([])

	with open('MotifMainArgDivUpdate_ParameterShotgun_Test6_MotifData_motif{motif}_len{maxStrandLength}_bias{motif_add}_elong{p_elong}_{trials}trials_numRound{numRounds}_bn{base_number}_div{p_divide}.csv'.format(motif = motif, maxStrandLength = maxStrandLength, motif_add= motif_add, p_elong=p_elong, trials=trials, numRounds=numRounds, base_number= base_number, p_divide=p_divide), 'wb') as f: 
		writer = csv.writer(f)
		writer.writerow(parameterlist)

		for i in range(trials):
			test_cell, test_strand, test_strand_count, celltracker_total[i] = source.MotifReg(motif, growthIterations, maxStrands, maxStrandLength, numCells, numRounds, p_elong, p_elong_motif, motif_add, base_number, p_divide)

			test_strand_norm = [x / float(y) for x,y in itertools.izip(test_strand,test_strand_count)] # normalize list of motif strands in trial
			test_strand_count_norm = [x / float(maxStrands*numCells) for x in test_strand_count] # normalize list of total strands in trial
			test_cell_norm = [x / float(numCells) for x in test_cell] # normalize list of cells in trial
			
			writer.writerow(test_strand_norm) # write row of the normalized motif strands for trial
			writer.writerow(test_strand_count_norm) # write row of the normalized motif strands for trial
			writer.writerow(test_cell_norm) # write row of the normalized motif strands for trial

			if i == 0:
				test_strand_ag = test_strand_norm # create a list where all of the round one data is in round one
				test_strand_count_ag = test_strand_count_norm # create a list where all of the round one data is in round one
				test_cell_ag = test_cell_norm
			else:
				test_strand_ag = [list(l) for l in zip(test_strand_ag,test_strand_norm)]  # create a list where all of the round one data is in round one
				test_strand_count_ag = [list(l) for l in zip(test_strand_count_ag,test_strand_count_norm)]  # create a list where all of the round one data is in round one
				test_cell_ag = [list(l) for l in zip(test_strand_count_ag,test_cell_norm)] 
		

			total_strand_count = [a + b for a, b in zip(total_strand_count, test_strand_count)] # add recent result to aggregate list


		for rounds in range(numRounds): # flatten the lists for each round
			test_strand_ag[rounds] = list(flatten(test_strand_ag[rounds]))
			test_strand_count_ag[rounds] = list(flatten(test_strand_count_ag[rounds]))
			test_cell_ag[rounds] = list(flatten(test_cell_ag[rounds]))


		for rounddata in test_strand_ag:
			test_strand_normavg.append(numpy.mean(rounddata))
			test_strand_SD.append(numpy.std(rounddata))

		for rounddata in test_strand_count_ag:
			test_strand_count_normavg.append(numpy.mean(rounddata))
			test_strand_count_SD.append(numpy.std(rounddata))
		
		for rounddata in test_cell_ag:
			test_cell_normavg.append(numpy.mean(rounddata))
			test_cell_SD.append(numpy.std(rounddata))



		writer.writerow(test_strand_normavg) # write row of the avg normalized motif strands for set
		writer.writerow(test_strand_count_normavg) # write row of the avg normalized number of strands 
		writer.writerow(test_cell_normavg) # write row of the avg normalized cells with motif
		
		writer.writerow(test_strand_SD) # write row of SD (per round) for motif strands
		writer.writerow(test_strand_count_SD) # write row of SD (per round) for total strands
		writer.writerow(test_cell_SD) # write row of SD (per round) for cells with motif

	f.close()

	keyorder = makeKeyorder(maxStrandLength,motif)

	dictavg = [] # create a list that will contain dictionaries
	for Round in range(numRounds): 
		# dictavg.append(collections.OrderedDict(sorted(holderdict,key = lambda i:keyorder.index(i[0])))) # create a dictionary for each round
		holderdict = {}
		for n in range(maxStrandLength):
			for key in itertools.product(range(2),repeat = n+1):
				mod_key = str(key).strip(" ,(),','").replace(", ", "")
				holderdict[mod_key] = 0 # initialize each dictionary with all the possible strands
		dictavg.append(collections.OrderedDict(sorted(holderdict.items(), key = lambda i:keyorder.index(i[0]))))

	for trial in range(trials):
		for Round in range(numRounds):
			for cell in range(len(celltracker_total[trial][Round])):
				for strand in celltracker_total[trial][Round][cell]:
					dictavg[Round][strand] =  dictavg[Round][strand] + 1 # add a value to the dictionary of the Round for each strand in each cell in each trial of this round

	for Round in range(numRounds):
		for key, value in dictavg[Round].items():
			dictavg[Round][key] = value / float(trials*total_strand_count[Round]) # normalize to get frequency of this strand/cell/trial/round (essentially the frequency in one cell of this strand in specific round)

	with open('MotifMainArgDivUpdate_ParameterShotgun_Test6_FullTrial1Data_motif{motif}_len{maxStrandLength}_bias{motif_add}_elong{p_elong}_{trials}trials_numRound{numRounds}_bn{base_number}_div{p_divide}.csv'.format(motif = motif, maxStrandLength = maxStrandLength, motif_add= motif_add, p_elong=p_elong, trials=trials, numRounds=numRounds, base_number= base_number, p_divide=p_divide), 'wb') as f: 
		writer2 = csv.writer(f, quotechar="'", quoting=csv.QUOTE_ALL) # making sure we put quotes around our strings so they're not read as numbers
		writer = csv.writer(f) # we don't need the quotes for the parameter list
		dwriter = csv.DictWriter(f,dictavg[0].keys())
		dwriter2 = csv.DictWriter(f,dictavg[0].keys(),quotechar="'", quoting=csv.QUOTE_ALL)

		writer.writerow(parameterlist) # write first row of the parameters

		for Round in range(numRounds):
			for cell in range(numCells):
				writer2.writerow(celltracker_total[0][Round][cell])  # write a row for each cell in each round in FIRST TRIAL

		dwriter2.writeheader() # write header of all the keys (the strands)
		for Round in range(numRounds):
			dwriter.writerow(dictavg[Round]) # write each row of frequencies, each row corresponds one round


	f.close()

if __name__ == "__main__":
	main(sys.argv[1:])

