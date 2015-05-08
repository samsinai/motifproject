import sys
import os
import collections

__version__ = "2015/03/27"


def isfloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False

def readInData(filename):
	
	with open(filename) as f:
		filelist = []
		for line in f:
			filelist.append(line.strip().split(","))
		motifs = []
		for entry in filelist[0]:
			motifs.append(entry.strip("'"))
		biases = filelist[1]
		lengths = filelist[2]

		freq_holderdict = {}
		sterr_holderdict = {}
		stdev_holderdict = {}

		for length in lengths:
			for motif in motifs:
				for bias in biases:
					freq_holderdict[(length,motif.strip("'"),bias)] = ''
					sterr_holderdict[(length,motif.strip("'"),bias)] = ''
					stdev_holderdict[(length,motif.strip("'"),bias)] = ''

		

		for i in range(3,len(filelist)):
			subfile = filelist[i][0]
			with open(subfile) as sf:
				handle = []
				for line in sf:
					handle.append(line.split(","))

				onelinelist = handle[0]
				motif = onelinelist[6].strip("'")
				bias = onelinelist[9].strip("'")
				length = onelinelist[3].strip("'")
				twolinelist = handle[1]
				fivelinelist= handle[4]
				freq_holderdict[(length,motif,bias)] = float(twolinelist[0].strip("'"))
				sterr_holderdict[(length,motif,bias)] = float(fivelinelist[2].strip("'"))
				stdev_holderdict[(length,motif,bias)] = float(fivelinelist[0].strip("'"))

			sf.close()



		freq_ordered_dict = collections.OrderedDict(sorted(freq_holderdict.items(), key = lambda x:(lengths.index(list(x)[0][0]),motifs.index(list(x)[0][1]),biases.index(list(x)[0][2]))))
		sterr_ordered_dict = collections.OrderedDict(sorted(sterr_holderdict.items(), key = lambda x:(lengths.index(list(x)[0][0]),motifs.index(list(x)[0][1]),biases.index(list(x)[0][2]))))
		stdev_ordered_dict = collections.OrderedDict(sorted(stdev_holderdict.items(), key = lambda x:(lengths.index(list(x)[0][0]),motifs.index(list(x)[0][1]),biases.index(list(x)[0][2]))))

	f.close()



	return motifs, biases, lengths, freq_ordered_dict, sterr_ordered_dict, stdev_ordered_dict



def writeOutputFile(filename, outputfile, motifs, biases, lengths, freq_data_dictionary, sterr_data_dictionary, stdev_data_dictionary):
	import csv
	
	freq_printlist = []
	freq_printlist_norm = []
	sterr_printlist = []
	stdev_printlist = []


	for length in range(len(lengths)):
		freq_printlist.append([])
		freq_printlist_norm.append([])
		sterr_printlist.append([])
		stdev_printlist.append([])
		for motif in range(len(motifs)):
			freq_printlist[length].append([])
			freq_printlist_norm[length].append([])
			sterr_printlist[length].append([])
			stdev_printlist[length].append([])
			for bias in range(len(biases)):
				freq_printlist[length][motif].append(freq_data_dictionary[(lengths[length],motifs[motif],biases[bias])])
				sterr_printlist[length][motif].append(sterr_data_dictionary[(lengths[length],motifs[motif],biases[bias])])
				stdev_printlist[length][motif].append(stdev_data_dictionary[(lengths[length],motifs[motif],biases[bias])])
				if isfloat(freq_data_dictionary[(lengths[length],motifs[motif],biases[bias])]):
					freq_printlist_norm[length][motif].append(float(freq_data_dictionary[(lengths[length],motifs[motif],biases[bias])])/float(freq_data_dictionary[(lengths[length],motifs[motif],biases[biases.index('0.5')])]))
				else:
					freq_printlist_norm[length][motif].append(freq_data_dictionary[(lengths[length],motifs[motif],biases[bias])])
			freq_printlist[length][motif].insert(0,"'" + motifs[motif] + "'" + ' motif')
			freq_printlist_norm[length][motif].insert(0,"'" + motifs[motif] + "'" + ' motif')
			sterr_printlist[length][motif].insert(0,'')
			stdev_printlist[length][motif].insert(0,'')



	with open(outputfile, "wb") as file:
		writer = csv.writer(file)
		writer.writerow(['Source file:',filename])

		writer.writerow(['Frequencies with Standard Error'])

		for length in range(len(lengths)):
			writer.writerow(['len:',int(lengths[length])])
			biaslabel = [float(bias) for bias in biases]
			biaslabel.insert(0,'')
			writer.writerow(biaslabel)
			for motif in range(len(motifs)):
				writer.writerow(freq_printlist[length][motif])
				writer.writerow(sterr_printlist[length][motif])


		writer.writerow(['Frequencies with Standard Deviation'])

		for length in range(len(lengths)):
			writer.writerow(['len:',int(lengths[length])])
			biaslabel = [float(bias) for bias in biases]
			biaslabel.insert(0,'')
			writer.writerow(biaslabel)
			for motif in range(len(motifs)):
				writer.writerow(freq_printlist[length][motif])
				writer.writerow(stdev_printlist[length][motif])

		writer.writerow(['Frequencies'])


		for length in range(len(lengths)):
			writer.writerow(['len:',int(lengths[length])])
			biaslabel = [float(bias) for bias in biases]
			biaslabel.insert(0,'')
			writer.writerow(biaslabel)
			for motif in range(len(motifs)):
				writer.writerow(freq_printlist[length][motif])

		writer.writerow(['Normalized Frequencies'])

		for length in range(len(lengths)):
			writer.writerow(['len:',int(lengths[length])])
			biaslabel = [float(bias) for bias in biases]
			biaslabel.insert(0,'')
			writer.writerow(biaslabel)
			for motif in range(len(motifs)):
				writer.writerow(freq_printlist_norm[length][motif])
	file.close()

def version():
	sys.stderr.write("""\
Compiles Equilibrim Data into a useful csv
Author: Grant Kinsler and Nicholas Keone Lee
Last updated: %s

""" % __version__)

def usage():
	version()
	sys.stderr.write("""\

Usage:

python CompileEquilibriumData.py listofcsvs.txt

where listofcsvs.txt is a txt file containing list of csv files to be read

	
""")


def main(argv):
	import getopt
	try:
		opts, args = getopt.getopt(argv, "v:", ["version"])
	except getopt.GetoptError, error:
		sys.stderr.write(str(error)+"\n")
		usage()
		sys.exit(2)
	for opt, arg in opts:
		if opt in ("-v", "--version"):
			version()
			sys.exit()
		else:
			sys.stderr.write("Unknown option %s\n" % opt)
			usage()
			sys.exit(2)
	filepath = args[0]
	motifs, biases, lengths, freq_data_dictionary, sterr_data_dictionary, stdev_data_dictionary = readInData(filepath)
	# outputname = "TESTOUTPUT_CompileEquilibriumData_" + os.path.basename(filepath).replace('.txt','.csv')
	outputname = os.path.basename(filepath).replace('.txt','.csv')
	writeOutputFile(filepath, outputname, motifs, biases, lengths, freq_data_dictionary, sterr_data_dictionary, stdev_data_dictionary)


if __name__ == "__main__":
	main(sys.argv[1:])