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
    # PUT OUR CODE HERE
    return randint(0, 5)


def getAllOrfs(gene):
    orfList = []
    orfList.append(gene)
    orfList.append(gene[1:]) # not including first character of string
    orfList.append(gene[2:]) # not including first or second

    orfList.append(gene[::-1]) # reverse the string
    orfList.append(gene[1::-1])
    orfList.append(gene[2::-1])

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

main()