# MotifFunctionDivisionUpdate.py: Function running simulations where motifs increase the bias of adding 0's
# Grant Kinsler and Nicholas Keone Lee
# Math 243, Spring 2014
# Authored: 2014/04/05
# Last updated: 2015/02/19, GRK

import numpy as np 
import random as rand
import math as math
from copy import deepcopy
from copy import copy

def MotifReg(motif, growthIterations, maxStrands, maxStrandLength, numCells, numRounds, p_elong, p_elong_motif, motif_add, base_number, p_divide):

    # counters
    strandprogress = [] # initalize list of number of eligible strands
    cellprogress = [] # intialize list of number of cells with eligible strands
    strandtotalprogress = [] # initalize list of number strands

    cells = []
    for i in range(numCells):
        cells.append([])

    celltracker = []

    for rounds in range(numRounds): # for each round
        ## GROWTH STAGE
        for times in range(growthIterations): # for each growth iteration
            for i in range(numCells): # for each cell
                has_motif = 0 # intiialize indicator for having the motif

                for strand in cells[i]: # for each strand in the current cell
                    if motif in strand: # check if strand has motif
                        has_motif = 1 # update indicator
                
                for strand in range(len(cells[i])):
                    if has_motif == 1: # lengthen strands in cells with motif
                        if len(cells[i][strand]) < maxStrandLength and rand.uniform(0,1) < p_elong_motif: # if strand length less than max AND if can elongate this round
                            if rand.uniform(0,1) < motif_add: # add 1 or 0, with bias accounted for
                                cells[i][strand] = cells[i][strand] + "0"
                            else:
                                cells[i][strand] = cells[i][strand] + "1"
                    else: # lenthen strands in other cells
                        if len(cells[i][strand]) < maxStrandLength and rand.uniform(0,1) < p_elong: # if strand length less than max AND if can
                            if rand.uniform(0,1) < 0.5: # add 1 or 0, with equal probab
                                cells[i][strand] = cells[i][strand] + "0"
                            else:
                                cells[i][strand] = cells[i][strand] + "1"

                # Add new strands
                if has_motif == 1:
                    for times in range(maxStrands-len(cells[i])):  # T
                        if rand.uniform(0,1) < p_elong_motif:
                            if rand.uniform(0,1) < motif_add: # add 1 or 0 with bias accounted for
                                cells[i].append("0")
                            else:
                                cells[i].append("1")
                else:
                    for times in range(maxStrands-len(cells[i])):
                        if rand.uniform(0,1) < p_elong:
                            if rand.uniform(0,1) < 0.5: # add 1 or 0 with bias accounted for
                                cells[i].append("0")
                            else:
                                cells[i].append("1")


        # DIVISION STAGE
        for i in range(len(cells)): # for each cell
            total_len = 0
            for strand in cells[i]: # for each strand
                total_len += len(strand)
                
            if total_len > base_number and rand.uniform(0,1) < p_divide:
                daughterOfAbraham = [] # daughter cell; actually Abraham had only sons
                strand_counter = 0 # placeholder for current strand wanted
                for times in range(len(cells[i])): # for each strand in Abraham
                    if rand.uniform(0,1) > 0.5: # with probability 1/2
                        daughterOfAbraham.append(cells[i].pop(strand_counter)) # this strand goes to daughter cell
                    else:
                        strand_counter += 1 # increase location of wanted strand in Abraham list
                cells.append(daughterOfAbraham)
        
        cells = rand.sample(cells, numCells) # take random sampling of numCells to keep population size constant


        # CHECK STAGE
        cell_counter = 0
        total_strands = 0
        motif_strand_counter = 0

        for i in range(numCells): # for each cell
            cell_has_motif = 0
            for strand in cells[i]: # and for each strand in the current cell
                total_strands += 1
                # print strand
                if motif in strand:
                    # print '***** motif confirmed *****'
                    cell_has_motif = 1 # update indicator
                    motif_strand_counter += 1
            if cell_has_motif == 1:
                cell_counter += 1
        cellprogress.append(copy(cell_counter)) # for counting number of cells with motif strands
        strandprogress.append(copy(motif_strand_counter)) # for counting number of strands with motifs
        strandtotalprogress.append(copy(total_strands)) # count total number of strands in population

        celltracker.append(deepcopy(cells)) # update cell tracker with this round's cell data

    return cellprogress, strandprogress, strandtotalprogress, celltracker
