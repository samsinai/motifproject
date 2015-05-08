import sys
import os

__version__ = "2015/03/18"

def readInGraphs(filename):
	csvFile = []
	graphs = []
	handle = open(filename)
	for line in handle:
		words = line.strip().split(",")
		csvFile.append(words)
	handle.close()
	lineOfInterest = int(csvFile[0][0])
	for graph in csvFile[(3*lineOfInterest)+1:(3*lineOfInterest)+4]:
		# print graph
		floats = [float(item) for item in graph] # cast each element from a string into a float
		graphs.append(floats)
	return csvFile[0], graphs # list of parameters, and 2-dimensional array (each row is one set of y-values)

def identifyEquilibriumPlateau(yvalues, variance):
	import numpy
	average = -1.00
	bend = -1
	samples = 0 # indicates how many samples taken

	for step in reversed(range(0, len(yvalues)-1)): # iterate backward through yvalues
		if samples < 150:
			if numpy.var(yvalues[step:], ddof = 1) < variance:
				samples = samples + 1 # update number of samples taken
				continue
			else:
				bend = step
				average = numpy.mean(yvalues[step:])
				break
		else:
			if numpy.var(yvalues[step:step+150], ddof = 1) < variance:
				continue
			else:
				bend = step+150
				average = numpy.mean(yvalues[step-150:])
				break
	return [average, bend]

def writeOutputFile(outputFile, paramArray, resultsArray, threshold):
	import csv
	with open(outputFile, "wb") as file:
		writer = csv.writer(file)
		writer.writerow(paramArray)
		for result in resultsArray:
			writer.writerow(result)
		writer.writerow(['threshold',threshold])
	file.close()

def version():
	sys.stderr.write("""\
Finds steady-state y-value of a given graph, given it exists.
Author: Nicholas Keone Lee and Grant Kinsler
Last updated: %s

""" % __version__)

def usage():
	version()
	sys.stderr.write("""\

Usage:

python findEquilibriumPlateauMod.py [-t <threshold>] graphs.csv

where graphs.csv is a CSV file containing the y-values of a graph (against discrete
time steps) in which a steady-state is reached by the end of the graph.

	Options:
		-t <threshold>
			Variance threshold.
""")


def main(argv):
	import getopt
	threshold = 0.0000100 # default threshold
	try:
		opts, args = getopt.getopt(argv, "v:t:", ["version", "threshold"])
	except getopt.GetoptError, error:
		sys.stderr.write(str(error)+"\n")
		usage()
		sys.exit(2)
	for opt, arg in opts:
		if opt in ("-v", "--version"):
			version()
			sys.exit()
		elif opt in ("-t", "--threshold"):
			threshold = float(arg)
		else:
			sys.stderr.write("Unknown option %s\n" % opt)
			usage()
			sys.exit(2)
	filepath = args[0]
	graphs = []
	parameters = [] # 1-dimensional
	results = [] # 2-dimensional
	parameters, graphs = readInGraphs(filepath)
	for graph in graphs:
		results.append(identifyEquilibriumPlateau(graph, threshold))
	outputname = "EquilibriumCalc_" + "thresh_" + str(threshold) + "_" + os.path.basename(filepath)
	writeOutputFile(outputname, parameters, results, threshold)


if __name__ == "__main__":
	main(sys.argv[1:])