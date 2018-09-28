from readfasta import readfasta
from genetic_code import code

def main():
    print("Bioinformatics - Assignment 1 - Group 3")
    print("Enter the name of the fasta file")
    filename = input()

    orf1 = readfasta(filename)[0][1]
    orf2 = orf1[1:] # not including first character of string
    orf3 = orf1[2:] # not including first or second

    print("\nForward ORF1:")
    print(translate(orf1))
    print("\nForward ORF2:")
    print(translate(orf2))
    print("\nForward ORF3:")
    print(translate(orf3))

    orf1reversed = orf1[::-1] # reverse the string

    orf4 = orf1reversed
    orf5 = orf1reversed[1:]
    orf6 = orf1reversed[2:]

    print("\nReversed ORF1:")
    print(translate(orf4))
    print("\nReversed ORF2:")
    print(translate(orf5))
    print("\nReversed ORF3:")
    print(translate(orf6))

def translate(dna):
    '''
    rna = dna.replace('T', 'U')

    aa_sequence = ''

    for i in range(0, len(rna), 3):
        codon = rna[i:i + 3]
        if len(codon) == 3:
            aa = code[codon]
            aa_sequence += aa

    return aa_sequence
    '''
    return dna

main()