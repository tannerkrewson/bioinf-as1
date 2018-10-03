from readfasta import readfasta
from genetic_code import code
from random import randint
from scipy import stats
import glob, os

def main():
    print( "***\nBioinformatics - Assignment 1 - Group 3\n***\n" )

    # for each fsa tested, a 0 will be added to the report card
    # if random/our function picked the wrong reading frame, and 
    # a 1 will be added if it picks the right reading frame
    randoms_report_card = []
    our_report_card = []

    random_was_right = 0
    we_were_right = 0
    number_of_files_scanned = 0
    ACTUAL_READING_FRAME = 1

    os.chdir( os.getcwd() + "/genes/" )
    for file in glob.glob( "*.fsa" ):
        print( file )
        gene = readfasta( file )[0][1]
        rfs = get_all_reading_frames( gene )

        # if random picks the correct reading frame
        if randint( 0, 5 ) == ACTUAL_READING_FRAME:
            random_was_right = random_was_right + 1
            randoms_report_card.append(1)
        else:
            randoms_report_card.append(0)

        # if our algorithm picks the correct reading frame
        if find_best_reading_frame( rfs ) == ACTUAL_READING_FRAME:
            we_were_right = we_were_right + 1
            our_report_card.append(1)
        else:
            our_report_card.append(0)

        number_of_files_scanned = number_of_files_scanned + 1

    print()
    print( str( number_of_files_scanned ) + " genes scanned" )
    print( "Our program was right " +
          str( we_were_right/number_of_files_scanned * 100 )
           + "% of the time" )
    print( "Random was right " +
           str( random_was_right/number_of_files_scanned * 100 )
           + "% of the time" )
    print()
    print( "Our report card: ", our_report_card )
    print( "Random's report card: ", randoms_report_card )

    ourStats = stats.ttest_ind(randoms_report_card,our_report_card)
    print("The p-value is " + str(ourStats.pvalue))
    if ourStats.pvalue < 0.05:
        print("Our program did statistically significantly better!  Woot!")
    else:
        print("Our program sucks...")

def find_best_reading_frame( rfs ):
    #
    # PUT OUR CODE HERE
    #

    rf_scores = [0, 0, 0, 0, 0, 0]
    possible_orf_list = []

    # rf_number: the six possible reading frames, 0 through 5
    # rf_bases: a string with all the bases of the reading frame
    for rf_number, rf_bases in enumerate( rfs ):

        # print( "\nThere are " + str( find_intron( rf_bases ) ) + " introns" )

        this_rf_score = 0

        # gets list of pairs of the start and stop indexes of the
        # possible orfs in the rf_number-th reading frame
        porfs = possible_orfs(rf_bases)
        possible_orf_list.append(porfs)

        # if the possible orfs are in an AT rich region, 
        # increase the score
        for porf in porfs:
            if at_rich_check(rf_bases, "start", porf[0]):
                rf_scores[rf_number] = rf_scores[rf_number] + 1

    # get the index of the best reading frame score
    best_rf = rf_scores.index( max( rf_scores ) )

    print( "Reading frame scores: ", rf_scores )
    print( "^^Best ORF is ORF" + str( best_rf ) + "!" )
    print()

    return best_rf

def get_all_reading_frames( gene ):
    rf_list = []
    rf_list.append( gene )
    rf_list.append( gene[1:] ) # not including first character of string
    rf_list.append( gene[2:] ) # not including first or second

    rev = gene[::-1]
    rf_list.append( rev ) # reverse the string
    rf_list.append( rev[1:] )
    rf_list.append( rev[2:] )

    return rf_list

#outputs possible orfs of over length 50 codons and returns them as a list
##of pairs of the index of the first and last basepair
def possible_orfs( dna ):
    orf_list = []
    position_of_last_start = 0
    looking_for_start = True
    for i in range( 0, len( dna ), 3 ):
        codon = dna[i:i + 3]
        
        if codon == 'ATG' and looking_for_start == True:
            looking_for_start = False
            postion_of_last_start = i
        if is_stop_codon( codon ) and looking_for_start == False:
            looking_for_start = True
            if ( ( i+3 ) - postion_of_last_start ) >= 50 * 3:
                orf_list.append( [postion_of_last_start, i + 3] )
    return orf_list

def is_stop_codon( codon ):
    return codon == 'TAA' or codon == 'TAG' or codon == 'TGA'

#computes the percentage of As and Ts in the region
#sequence is a string of nucleotide bases, the sequence to be analyzed
#start is the beginning of the range to check
#end is the end of the range to check
#returns a float in the range 0-100
def at_rich_percentage(sequence, start, end):

    region = sequence[start:stop]

    at_count = 0
    for i in range(0, len(region) - 1):
        if region[i] == 'A' or region[i] == 'T':
            at_count += 1

    return (at_count/len(region)) * 100

#correctly calls at_rich_percentage() for check_type
#sequence is a string of nucleotide bases, the sequence to be analyzed
#check_type is a string "intron", "start", or "stop"
#  for the type of check wanted
#start_index is the index of the start codon, stop codon, 
#  or the start of the intron
#stop_index is the index of the end of the intron,
#  only necessary for intron checks
#returns a boolean value that represents if the check type was within
#  the richness bounds
#richness bounds chosen based on upper limits of at-richness observed
#  in various tests of region types
#error input returns False
def at_rich_check(sequence, check_type, start_index, stop_index = 0):

    if(check_type == "intron")
        region_composition = at_rich_percentage(sequence, \
                                                start_index, \
                                                stop_index)
    elif(check_type == "start")
        region_composition = at_rich_percentage(sequence, \
                                                start_index - 200, \
                                                start_index - 1)
    elif(check_type == "stop")
        region_composition = at_rich_percentage(sequence, \
                                                start_index + 1, \
                                                start_index + 200)
    else
        return False

    if(check_type == "intron")
        return region_composition > 68.5
    else
        return region_composition > 63

#3' splice site UAG or CAG
#intron 5' sequence GTATGT
#Rian's Code
def find_intron( dna ): 
    rna = dna.replace( 'T', 'U' )
    pos_last_start = 0
    looking_for_start = True # while true, look for first start and ignore stops,
                             # if false, look until stop is found
    count = 0
    for i in range( 0, len( rna ), 1 ):
        start_seq = rna[ i:i + 6 ]
        end_seq = rna[ i:i + 3 ]
        if ( start_seq == "GUAUGU" or start_seq == "GUACGU" or
             start_seq == "GUAUGA" ) and looking_for_start == True:
            looking_for_start = False
            pos_last_start = i
            print( "\nStart of intron at " + str( pos_last_start ) )
        if ( end_seq == "UAG" or end_seq == "CAG" ) and looking_for_start == False:
            looking_for_start = True
            pos_end = i + 3 # Add 3 to see where the last base of the stop is
            count += 1
            print( "\nEnd of intron at " + str( pos_end ) )
            
    return count



main()
