from readfasta import readfasta
from random import randint
from scipy import stats
import glob, os

'''
main: read in the fasta files, run our algorithm to guess the reading
    frame, and statistically compare it to random
'''
def main():
    print( "*****\nBioinformatics - Assignment 1 - Group 3\n*****\n" )

    # for each fsa tested, a 0 will be added to the report card
    # if random/our function picked the wrong reading frame, and 
    # a 1 will be added if it picks the right reading frame
    randoms_report_card = []
    our_report_card = []

    random_was_right = 0
    we_were_right = 0
    number_of_files_scanned = 0

    # the actual reading frame of all fsa we input is always 2+
    # which is represented by the index 1
    ACTUAL_READING_FRAME = 1

    # scan in all fsa files in the "genes" directory
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

    print( "*****\n" )
    print( number_of_files_scanned, "genes scanned" )

    percent_we_were_right = we_were_right/number_of_files_scanned * 100
    percent_rand_was_right = random_was_right/number_of_files_scanned * 100

    percent_we_were_right = round( percent_we_were_right, 3 )
    percent_rand_was_right = round( percent_rand_was_right, 3 )

    print( "Our code was right ", percent_we_were_right, "% of the time" )
    print( "Random was right ", percent_rand_was_right, "% of the time\n" )

    print( "Our report card: ", our_report_card )
    print( "Random's report card: ", randoms_report_card )

    ourStats = stats.ttest_ind(randoms_report_card,our_report_card)
    print("The p-value is " + str(ourStats.pvalue))

    if ourStats.pvalue < 0.05:
        print("Our program did statistically significantly better!  Woot!")
    else:
        print("We did not do statistically significantly better than random")


'''
find_best_reading_frame: finds the best out of the six reading frames using 
    our algorithm

Parameter: a list of the six reading frames as strings of nucleotides

Return: the index of the reading frame our algorithm selected, 0 through 5
'''
def find_best_reading_frame( rfs ):
    rf_scores = [0, 0, 0, 0, 0, 0]

    # rf_number: the six possible reading frames, 0 through 5
    # rf_bases: a string with all the bases of the reading frame
    for rf_number, rf_bases in enumerate( rfs ):
        this_rf_score = 0

        # remove introns from the sequence
        rf_bases = remove_introns(rf_bases)

        # gets list of pairs of the start and stop indexes of the
        # possible orfs in the rf_number-th reading frame
        porfs = possible_orfs(rf_bases)

        # will hold the scores of all possible orfs
        # 0 is added to be the minimum score
        porf_scores = [ 0 ]

        # find and record the score of how confident we are in 
        # the viability of each possible reading frames
        for porf in porfs:
            this_porf_score = score_orf( rf_bases, porf[0], porf[1] )
            porf_scores.append( this_porf_score )

        best_orf_score = max( porf_scores )

        # set the score of this reading frame to the score of our guess 
        # for the most viable orf in this reading frame
        rf_scores[rf_number] = best_orf_score

    # get the index of the best reading frame score
    best_rf = rf_scores.index( max( rf_scores ) )

    print( "Reading frame scores: ", rf_scores )
    print( "So, best reading frame guess is RF" + str( convert_rf_index( best_rf ) ) )
    print()

    return best_rf


'''
get_all_reading_frames: generates all six reading frames by shifting and
    reversing the gene and returns them in a list

Parameter: the non-shifted string of nucleotides

Returns: a list of the six reading frames
'''
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


'''
possible_orfs: finds the possible orfs of the given gene

Parameter: the string of nucleotides to find orfs in

Returns: a list of pairs of the index of the first and last nucleotide
'''
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
            if ( ( i + 3 ) - postion_of_last_start ) >= 50 * 3:
                orf_list.append( [postion_of_last_start, i + 3] )
    return orf_list


'''
score_orf: makes an educated guess on the likeliness that the given 
    possible orf is the actual orf

Parameters:
    bases: a string of the nucleotides of the whole reading frame, 
        not just the possible orf
    orf_start: the beginning index of the orf
    orf_stop: the ending index of the orf

Returns: a score between 0 and 6
'''
def score_orf ( bases, orf_start, orf_stop ):
    score = 0
    orf_length = orf_stop - orf_start

    # ORF length scoring:
    # 50-99:   1 point
    # 100-149: 2 points
    # 150-199: 3 points
    # 200-249: 4 points
    # 250+:    5 points


    if orf_length >= 50*3:
        score += 1
    if orf_length >= 100*3:
        score += 1
    if orf_length >= 150*3:
        score += 1
    if orf_length >= 200*3:
        score += 1
    if orf_length >= 250*3:
        score += 1

    # if the region before the potential orf is AT rich,
    # give this orf a point
    if at_rich_check( bases, "start", orf_start ):
        score += 1

    # if the region after the potential orf is AT rich,
    # give this orf a point
    if at_rich_check( bases, "stop", orf_stop ):
        score += 1

    # if a polyadenylation sequence is found after the orf,
    # give this orf a point
    if find_polya_sequence(bases, orf_stop):
        score += 1

    return score

    
'''
is_stop_codon: returns true if the given codon is a stop codon

Parameter: a three-base codon

Returns: true if the given codon is a stop codon 
'''
def is_stop_codon( codon ):
    return codon == 'TAA' or codon == 'TAG' or codon == 'TGA'


'''
at_rich_percentage: computes the percentage of As and Ts in the region
    
Parameters:
    sequence: a string of nucleotide bases, the sequence to be analyzed
    start: the beginning of the range to check
    end: the end of the range to check

Returns: a float in the range 0-100
'''
def at_rich_percentage( sequence, start, end ):
    region = sequence[start:end]

    at_count = 0
    for i in range(0, len(region) - 1):
        if region[i] == 'A' or region[i] == 'T':
            at_count += 1
    return (at_count/len(region)) * 100


'''
at_rich_check: calls at_rich_percentage() for check_type

Parameters:
    sequence: a string of nucleotide bases, the sequence to be analyzed
    check_type: a string "intron", "start", or "stop", for the type of check 
        wanted
    start_index: the index of the start codon, stop codon, or the start of 
        the intron
    stop_index: the index of the end of the intron, only necessary for 
        intron checks

Returns: a boolean value that represents if the check type was within
    the richness bounds, which are chosen based on upper limits of 
    at-richness observed in various tests of region types
'''
def at_rich_check(sequence, check_type, start_index, stop_index = 0):

    if start_index - 200 < 0:
        return False

    if check_type == "intron":
        region_composition = at_rich_percentage(sequence, \
                                                start_index, \
                                                stop_index)
    elif check_type == "start":
        region_composition = at_rich_percentage(sequence, \
                                                start_index - 200, \
                                                start_index - 1)
    elif check_type == "stop":
        region_composition = at_rich_percentage(sequence, \
                                                start_index + 1, \
                                                start_index + 200)
    else:
        return False

    if check_type == "intron":
        return region_composition > 70
    else:
        return region_composition > 62

'''
find_polya_sequence: attempts to find a polyadenylation sequence after
    a given index

Parameters:
    sequence: a string of nucleotides to look within
    stop_index: the index of the end of the orf to look after

Returns: True if a polyadenylation sequence is found, False if not
'''
def find_polya_sequence(sequence, stop_index):

    for i in range(stop_index, stop_index + 140):
        polya_sequence = sequence[i:i + 6]

        is_polya_sequence = polya_sequence == "AATAAA" or \
                            polya_sequence == "AAAAAA" or \
                            polya_sequence == "TATGTA" or \
                            polya_sequence == "TATATA" or \
                            polya_sequence == "TACATA"

        if is_polya_sequence:
            return True

    return False


'''
remove_introns: removes introns from the given string of nucleotides

Parameter: a string of nucleotides

Returns: the same string of nucleotides, except with the potential 
    introns removed. note that the reading frame of nucleotides that 
    come after any removed introns is not preserved
'''
def remove_introns( dna ): 
    # while true, look for first start and ignore stops,
    # if false, look until stop is found
    looking_for_start = True 

    intron_list = []
    pos_last_start = 0
    for i in range( 0, len( dna ), 1 ):
        start_seq = dna[ i:i + 6 ]
        end_seq = dna[ i:i + 3 ]

        # intron 5' start sequence
        is_intron_start = start_seq == "GTATGT" \
            or start_seq == "GTACGT" \
            or start_seq == "GTATGA"

        # 3' splice site
        is_intron_end = end_seq == "TAG" or end_seq == "CAG"

        if is_intron_start and looking_for_start:
            looking_for_start = False
            pos_last_start = i
        
        if is_intron_end and not looking_for_start:
            looking_for_start = True

            # Add 3 to see where the last base of the stop is
            pos_end = i + 3

            # if the possible intron is AT rich, then count it as a
            # possible intron
            if at_rich_check(dna, "intron", pos_last_start, pos_end):
                intron_list.append([pos_last_start, pos_end])

    total_introns_length = 0
    for intron_pair in intron_list:
        intron_pair[0] -= total_introns_length
        intron_pair[1] -= total_introns_length
        dna = dna[:intron_pair[0]] + dna[intron_pair[1]:]
        total_introns_length += intron_pair[1] - intron_pair[0]
            
    return dna


'''
convert_rf_index:
    converts  0, 1, 2, 3, 4, 5
          to  1+,2+,3+,1-,2-,3-

Parameter: the index of the reading frame, 0 through 5

Returns: the corresponding human-readable reading frame number
'''
def convert_rf_index( rf_index ):
    if rf_index < 3:
        return str( rf_index+1 ) + "+"
    else:
        return str( rf_index-2 ) + "-"

main()
