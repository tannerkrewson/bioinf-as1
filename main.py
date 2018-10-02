from readfasta import readfasta
from genetic_code import code
from random import randint
import glob, os

def main():
    print("Bioinformatics - Assignment 1 - Group 3")

    os.chdir(os.getcwd() + "/genes/")
    for file in glob.glob("*.fsa"):
        gene = readfasta(file)[0][1]
        orfs = getAllOrfs(gene)

        randomORF = randint(0, 5)
        ourORF = findBestORF(orfs)
        actualORF = 1

        print(file)
        print("Random ORF: " + str(randomORF))
        print("Our ORF Choice: " + str(ourORF))
        print("Actual ORF: " + str(actualORF))
        print()

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