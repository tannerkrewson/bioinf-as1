from readfasta import readfasta
from genetic_code import code
from random import randint
import glob, os

def main():
    print("Bioinformatics - Assignment 1 - Group 3")

    randomWasRight = 0
    weWereRight = 0
    numberOfFilesScanned = 0
    ACTUAL_READING_FRAME = 1

    os.chdir(os.getcwd() + "/genes/")
    for file in glob.glob("*.fsa"):
        gene = readfasta(file)[0][1]
        rfs = getAllReadingFrames(gene)

        if randint(0, 5) == ACTUAL_READING_FRAME:
            randomWasRight = randomWasRight + 1

        if findBestReadingFrame(rfs) == ACTUAL_READING_FRAME:
            weWereRight = weWereRight + 1

        numberOfFilesScanned = numberOfFilesScanned + 1

    print()
    print(str(numberOfFilesScanned) + " genes scanned")
    print("Our program was right " + str(weWereRight/numberOfFilesScanned * 100) + "% of the time")
    print("Random was right " + str(randomWasRight/numberOfFilesScanned * 100) + "% of the time")

def findBestReadingFrame(rfs):
    #
    # PUT OUR CODE HERE
    #

    rf_scores = [0, 0, 0, 0, 0, 0]
    possible_orf_list = []

    # rf_number: the six possible reading frames, 0 through 5
    # rf_bases: a string with all the bases of the reading frame
    for rf_number, rf_bases in enumerate(rfs):
        thisRfScore = 0

        # gets list of pairs of the start and stop indexes of the
        # possible orfs in the rf_number-th reading frame
        porfs = possible_orfs(rf_bases)
        possible_orf_list.append(porfs)

        # if the possible orfs are in an AT rich region, 
        # increase the score
        for porf in porfs:
            if at_rich_check(rf_bases, porf[0]):
                rf_scores[rf_number] = rf_scores[rf_number] + 1

    # get the index of the best reading frame score
    bestRF = rf_scores.index(max(rf_scores))

    print(rf_scores)
    print("The best ORF is ORF" + str(bestRF) + "!")

    return bestRF

def getAllReadingFrames(gene):
    rfList = []
    rfList.append(gene)
    rfList.append(gene[1:]) # not including first character of string
    rfList.append(gene[2:]) # not including first or second

    rev = gene[::-1]
    rfList.append(rev) # reverse the string
    rfList.append(rev[1:])
    rfList.append(rev[2:])

    return rfList

#outputs possible orfs of over length 50 codons and returns them as a list
##of pairs of the index of the first and last basepair
def possible_orfs(dna):
    orf_list = []
    position_of_last_start = 0
    looking_for_start = True
    for i in range(0, len(dna), 3):
        codon = dna[i:i + 3]
        
        if codon == 'ATG' and looking_for_start == True:
            looking_for_start = False
            postion_of_last_start = i
        if is_stop_codon(codon) and looking_for_start == False:
            looking_for_start = True
            if ((i+3) - postion_of_last_start) >= 50 * 3:
                orf_list.append([postion_of_last_start, i + 3])
    return orf_list

def is_stop_codon(codon):
    return codon == 'TAA' or codon == 'TAG' or codon == 'TGA'

def at_rich_check(sequence, start_index):
    if start_index < 200:
        return False

    at_rich_region = sequence[start_index - 200:start_index - 1]
    intergenic_region = sequence[start_index + 3:start_index + 202]
    
    rich_at_count = 0
    for i in range(0, 199):
        if at_rich_region[i] == 'A' or at_rich_region[i] == 'T':
            rich_at_count += 1

    intergenic_at_count = 0
    for i in range(0, 199):
        if intergenic_region[i] == 'A' or intergenic_region[i] == 'T':
            intergenic_at_count += 1

    return rich_at_count > intergenic_at_count

#3' splice site UAG or CAG
#intron 5' sequence GTATGT
#Rian's Code
def find_intron(dna): 
    rna = dna.replace('T', 'U')
    pos_last_start = 0
    looking_for_start = True
    count = 0
    for i in range(0, len(rna), 1):
        start_seq = rna[i:i + 6]
        end_seq = rna[i:i + 3]
        if start_seq == "GUAUGU" and looking_for_start == True:
            looking_for_start = False
            pos_last_start = i
            print("\nStart of intron at " + str(pos_last_start))
        if (end_seq == "UAG" or end_seq == "CAG") and looking_for_start == False:
            looking_for_start = True
            pos_end = i + 3
            count += 1
            print("\nEnd of intron at " + str(pos_end))
            
    return count

main()
