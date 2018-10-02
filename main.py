from readfasta import readfasta
from genetic_code import code
from random import randint
import glob, os

def main():
    print("Bioinformatics - Assignment 1 - Group 3")

    randomWasRight = 0
    weWereRight = 0
    numberOfFilesScanned = 0
    ACTUAL_ORF = 1

    os.chdir(os.getcwd() + "/genes/")
    for file in glob.glob("*.fsa"):
        gene = readfasta(file)[0][1]
        orfs = getAllOrfs(gene)

        if randint(0, 5) == ACTUAL_ORF:
            randomWasRight = randomWasRight + 1

        if findBestORF(orfs) == ACTUAL_ORF:
            weWereRight = weWereRight + 1

        numberOfFilesScanned = numberOfFilesScanned + 1

    print()
    print(str(numberOfFilesScanned) + " genes scanned")
    print("Our program was right " + str(weWereRight/numberOfFilesScanned * 100) + "% of the time")
    print("Random was right " + str(randomWasRight/numberOfFilesScanned * 100) + "% of the time")

def findBestORF(orfs):
    #
    # PUT OUR CODE HERE
    #

    orf_scores = [0, 0, 0, 0, 0, 0]
    possible_orf_list = []

    for idx, reading_frame in enumerate(orfs):
        thisOrfScore = 0

        porf = possible_orfs(reading_frame)
        possible_orf_list.append(porf)

        if (len(porf) > 0):
            # if it has a orf, increase the score
            orf_scores[idx] = 1

    # get the index of the best reading frame score
    bestORF = orf_scores.index(max(orf_scores))

    print(orf_scores)
    print("The best ORF is ORF" + str(bestORF) + "!")

    return bestORF
def getAllOrfs(gene):
    orfList = []
    orfList.append(gene)
    orfList.append(gene[1:]) # not including first character of string
    orfList.append(gene[2:]) # not including first or second

    rev = gene[::-1]
    orfList.append(rev) # reverse the string
    orfList.append(rev[1:])
    orfList.append(rev[2:])

    return orfList

def translateToAminoAcid(dna):
    rna = dna.replace('T', 'U')

    aa_sequence = ''

    for i in range(0, len(rna), 3):
        codon = rna[i:i + 3]
        if len(codon) == 3:
            aa = code[codon]
            aa_sequence += aa

    return aa_sequence

def countStart(dna):
    rna = dna.replace('T', 'U')
    count = 0
    for i in range(0, len(rna), 3):
        codon = rna[i:i + 3]
        if codon == 'AUG':
            count += 1

    return count

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

main()
