from __future__ import division

import itertools
import time

from deprecation import deprecated

from collections import defaultdict

import random

import math


def format_string(string):
    string = string.replace(',', "")
    string = string.replace("'", "")
    string = string.replace("[", "")
    string = string.replace("]", "")
    print string

def patternCount(text,pattern):
    count = 0
    for i in range(len(text)):
        if text[i:i+len(pattern)] == pattern:
            count = count+1
    return count

def frequent_words(text, k):
    frequentPatterns = []
    count = []
    max_count = 0
    for i in range(len(text)-k):
        pattern = text[i:i+k]
        count.append(patternCount(text, pattern))
        if patternCount(text, pattern)>max_count:
            max_count = count[i]
    for i in range(len(text)-k):
        if count[i] == max_count:
            if text[i:i+k] not in frequentPatterns:
                frequentPatterns.append(text[i:i+k])
    return frequentPatterns

def frequent_words_mod(text, k,t):
    frequentPatterns = []
    for i in range(len(text)-k):
        pattern = text[i:i+k]
        if patternCount(text, pattern) == t and pattern not in frequentPatterns:
                frequentPatterns.append(pattern)
    return frequentPatterns


def number_to_pattern(num, length):
    pattern = ""
    reValues = {0: 'A', 3: 'T', 2: 'G', 1: 'C'}
    while num > 0:
        div = num // 4
        mod = num % 4
        pattern = reValues[mod] + pattern
        num = div

    while len(pattern) < length:
        pattern = "A" + pattern

    return (pattern)


def pattern_to_number(text):
    dict = {'A': 0,
            'C': 1,
            'G': 2,
            'T': 3}
    output = 0
    for i in range(len(text)):
        output = output + ((dict[text[i]]) * (4 ** (len(text) - 1 - i)))
    return(output)


def reverseComplement(pattern):
    output = ""
    for i in range(1,len(pattern)+1):
        if pattern[len(pattern)-i] == "A":
            output = output + "T"
        if pattern[len(pattern)-i] == "T":
            output = output + "A"
        if pattern[len(pattern)-i] == "G":
            output = output + "C"
        if pattern[len(pattern)-i] == "C":
            output = output + "G"
    return output

def pattern_matching(pattern, genome):
    index = []
    for i in range(len(genome)):
        if genome[i:i+len(pattern)] == pattern:
            index.append(i)
    return index

def clump_finding(genome, k, L, t):
    kmers = []
    for i in range(len(genome)-L):
        patterns = faster_frequent_words(genome[i:i+L],k,t)
        if patterns not in kmers:
            kmers.append(patterns)
    return kmers

def computing_frequencies(text, k):
    frequency_array = {}
    for i in range(4**k):
        frequency_array[i] = 0
    for i in range(len(text)-(k-1)):
        pattern = text[i:i+k]
        j = pattern_to_number(pattern)
        frequency_array[j] = frequency_array[j] + 1
    return frequency_array.values()

def faster_frequent_words(text, k):
    frequent_patterns = []
    frequency_array = computing_frequencies(text, k)
    max_count = max(frequency_array)
    for i in range(4**k):
        if frequency_array[i] == max_count:
            pattern = number_to_pattern(i,k)
            frequent_patterns.append(pattern)
    return frequent_patterns

def frequent_words_by_sorting(text, k,t):
    frequent_patterns = []
    index = []
    count = []
    for i in range(len(text)-k):
        index.append(pattern_to_number(text[i:i+k]))
        count.append(1)
    index.sort()
    for i in range(1,len(index)):
        if index[i] == index[i-1]:
            count[i] = count[i] + 1
    max_count = t#max(count)
    for i in range(len(count)):
        if max_count == count[i]:
            frequent_patterns.append(number_to_pattern(index[i],k))
    return frequent_patterns




@deprecated("unused")
def reverse_complement(pattern):
    output = ""
    for i in range(1, len(pattern) + 1):
        if pattern[len(pattern) - i] == "A":
            output = output + "T"
        if pattern[len(pattern) - i] == "T":
            output = output + "A"
        if pattern[len(pattern) - i] == "G":
            output = output + "C"
        if pattern[len(pattern) - i] == "C":
            output = output + "G"
    return output


@deprecated("unused")
def skew(genome):
    total_skew = [0]
    diff = 0
    for i in range(len(genome)):
        if genome[i] == "G":
            diff = diff + 1
        elif genome[i] == "C":
            diff = diff - 1
        else:
            diff = diff
        total_skew.append(diff)
    return total_skew


@deprecated("unused")
def minimum_skew(genome):
    total_skew = skew(genome)
    min_skew = total_skew[1]
    minimums = []
    for i in range(1, len(total_skew)):
        if total_skew[i] <= min_skew:
            min_skew = total_skew[i]
    for i in range(1, len(total_skew)):
        if total_skew[i] == min_skew:
            minimums.append(i)
    return minimums


@deprecated("unused")
def maximum_skew(genome):
    total_skew = skew(genome)
    max_skew = total_skew[1]
    maximums = []
    for i in range(1, len(total_skew)):
        if total_skew[i] >= max_skew:
            max_skew = total_skew[i]
    for i in range(1, len(total_skew)):
        if total_skew[i] == max_skew:
            maximums.append(i)
    return maximums


#
#
#
def hamming_distance(str1, str2):
    mismatch = 0
    if len(str1) > len(str2):
        for i in range(1, len(str2) + 1):
            if str1[len(str1) - i] != str2[len(str2) - i]:
                mismatch = mismatch + 1
    else:
        for i in range(1, len(str1) + 1):
            if str1[len(str1) - i] != str2[len(str2) - i]:
                mismatch = mismatch + 1
    return mismatch


#
#
#
def approx_match(pattern, genome, max_mm):
    matches = False
    for i in range(len(genome) - len(pattern) + 1):
        if hamming_distance(pattern, genome[i:i + len(pattern)]) <= max_mm:
            matches = True
    return matches


#
#
#
def count(genome, pattern, max_mm):
    matches = []
    for i in range(len(genome) - len(pattern) + 1):
        if hamming_distance(pattern, genome[i:i + len(pattern)]) <= max_mm:
            matches.append(i)
    return len(matches)


#
#
#
def frequent_words_with_mismatches(text, k, max_mm):
    frequencies = defaultdict(lambda: 0)
    frequent = []
    for i in range(len(text) - k + 1):
        kmer_and_friends = PermuteMotifDistanceTimes(text[i:i + k], max_mm)
        for kmer in kmer_and_friends:
            frequencies[kmer] += 1
    for kmer in frequencies:
        if frequencies[kmer] == max(frequencies.values()):
            frequent.append(kmer)
    return frequent


#
#
#
def frequent_words_with_mismatches_and_revcomp(text, k, max_mm):
    frequencies = defaultdict(lambda: 0)
    frequent = []
    for i in range(len(text) - k + 1):

        pattern = text[i:i + k]
        neighbors = PermuteMotifDistanceTimes(pattern, max_mm)
        for kmer in neighbors:
            frequencies[kmer] += 1

    text = reverseComplement(text)
    for i in range(len(text) - k + 1):
        pattern = text[i:i + k]
        neighbors = PermuteMotifDistanceTimes(pattern, max_mm)
        for kmer in neighbors:
            frequencies[kmer] += 1

    for kmer in frequencies:
        if frequencies[kmer] == max(frequencies.values()):
            frequent.append(reverseComplement(kmer))
    return frequent


def PermuteMotifDistanceTimes(motif, d):
    workingSet = {motif}
    for _ in range(d):
        workingSet = set(itertools.chain.from_iterable(map(PermuteMotifOnce, workingSet)))
    return list(workingSet)


def PermuteMotifOnce(motif, alphabet={"A", "C", "G", "T"}):
    """
    Gets all strings within hamming distance 1 of motif and returns it as a
    list.
    """
    return list(set(itertools.chain.from_iterable([[
        motif[:pos] + nucleotide + motif[pos + 1:] for
        nucleotide in alphabet] for
        pos in range(len(motif))])))


def all_kmers_over_one(d):
    workingSet = "A"
    for i in range(1, d):
        workingSet += "A"
    workingSet = {workingSet}
    for _ in range(d):
        workingSet = set(itertools.chain.from_iterable(map(PermuteMotifOnce, workingSet)))
    return list(workingSet)


def kmers_one(motif, alphabet={"A", "C", "G", "T"}):
    """
    Gets all strings within hamming distance 1 of motif and returns it as a
    list.
    """
    return list(set(itertools.chain.from_iterable([[
        motif[:pos] + nucleotide + motif[pos + 1:] for
        nucleotide in alphabet] for
        pos in range(len(motif))])))


def immediate_neighbors(pattern):
    tides = ["A", "T", "C", "G"]
    neighborhood = []
    for i in range(len(pattern)):
        symbol = pattern[i]
        for nuc in tides:
            if nuc != symbol:
                neighbor = pattern[0:i] + nuc + pattern[i + 1:]
                neighborhood.append(neighbor)
    return neighborhood


def neighbors(pattern, d):
    tides = ["A", "T", "C", "G"]
    if d == 0:
        return pattern
    if len(pattern) == 1:
        return ["A", "C", "G", "T"]
    neighborhood = [pattern]
    suffix_neighbors = neighbors(suffix(pattern), d)
    for text in suffix_neighbors:
        if hamming_distance(text, pattern) == d:
            neighborhood.append(pattern[0] + text)
        if hamming_distance(text, pattern) < d:
            for nuc in tides:
                if nuc != pattern[0]:
                    neighborhood.append(nuc + text)

    return neighborhood


def suffix(pattern):
    return pattern[1:]


# takes list of dna and length and returns all common motifs within the dna and within d distance
def motif_enumeration(Dna, k, d):
    patterns = []
    for dna in Dna:
        for i in range(len(dna) - k + 1):
            km = dna[i:i + k]
            for kmer in PermuteMotifDistanceTimes(km, d):
                match = True
                for dna in Dna:
                    if not approx_match(kmer, dna, d):
                        match = False
                if match == True:
                    patterns.append(kmer)
    patterns = list(dict.fromkeys(patterns))
    return sorted(patterns)


def minimum_hamming_distance(pattern, genome):
    min = len(pattern)
    for i in range(len(genome) - len(pattern) + 1):
        if hamming_distance(pattern, genome[i:i + len(pattern)]) < min:
            min = hamming_distance(pattern, genome[i:i + len(pattern)])
    return min


def profile_kmer(text, k, profile):
    probabilities = []
    kmers = []
    for i in range(len(text) - k + 1):
        probabilities.append(profile_score(text[i:i + k], profile))
    for i in range(len(probabilities)):
        if probabilities[i] == max(probabilities):
            kmers.append(text[i:i + k])
    return kmers


def find_consensus(motifs):
    profile = create_profile(motifs)
    consensus = ""
    nucleotides = ["A", "C", "G", "T"]
    for i in range(len(profile.get('A'))):
        column = []
        ns = []
        for n in nucleotides:
            column.append(profile.get(n)[i])
            ns.append(n)
        for j in range(len(column)):
            if column[j] == max(column):
                consensus += ns[j]
    return consensus


def find_consensus_profile(profile):
    consensus = ""
    nucleotides = ["A", "C", "G", "T"]
    for i in range(len(profile.get('A'))):
        column = []
        ns = []
        for n in nucleotides:
            column.append(profile.get(n)[i])
            ns.append(n)
        for j in range(len(column)):
            if column[j] == max(column):
                consensus += ns[j]
    return consensus


def create_profile(motifs):
    profile = {
        'A': [],
        'C': [],
        'G': [],
        'T': []
    }
    total = len(motifs)
    nucleotides = ['A', 'C', 'G', 'T']
    for j in range(len(motifs[0])):
        column = []
        for i in range(len(motifs)):
            column.append(motifs[i][j])
        for n in nucleotides:
            value = column.count(n)
            profile[n].append(float(value) / total)
    return profile



def greedy_motif_search(Dna, k, t):
    first_strand = []
    best_motifs = []  # BestMotifs to motif matrix formed by first k-mers in each string in Dna
    for i in range(t):
        best_motifs.append(Dna[i][0:k])
    for i in range(len(Dna[0]) - k + 1):
        first_strand.append(Dna[0][i:i + k])
    for kmer in first_strand:
        motifs = [kmer]
        for i in range(1, t):
            pass


# def create_profile(motifs):
#     profile = {
#     'A': [],
#     'C': [],
#     'G': [],
#     'T': []
#     }
#     nucleotides = ["A", "C", "G", "T"]
#     total = len(motifs)
#     for i in range(len(motifs)):
#         column = []
#         for j in range(len(motifs[i])):
#             column.append(motifs[i][j])
#         for n in nucleotides:
#             profile[n].append(float(column.count(n))/total)
#     return profile


def score(motifs):
    profile = create_profile(motifs)
    score = 0
    for motif in motifs:
        score += hamming_distance(motif, find_consensus_profile(profile))
    return score



def profile_score(kmer, profile):
    score = 1.0
    for i in range(len(kmer)):
        index = profile.get(kmer[i])
        score = score * index[i]
    return score


def laplace_create_profile(motifs):
    profile = {
        'A': [],
        'C': [],
        'G': [],
        'T': []
    }
    nucleotides = ['A', 'C', 'G', 'T']
    for j in range(len(motifs[0])):
        total = 0
        column = []
        for i in range(len(motifs)):
            column.append(motifs[i][j])
        for n in nucleotides:
            value = column.count(n)+1
            total += value
        for n in nucleotides:
            value = column.count(n)+1
            profile[n].append(float(value) / (total))
    return profile


def greedy_motif_search(Dna, k, t):
    first_strand = []
    best_motifs = []  # BestMotifs to motif matrix formed by first k-mers in each string in Dna
    for i in range(t):
        best_motifs.append(Dna[i][0:k])
    for i in range(len(Dna[0]) - k + 1):
        first_strand.append(Dna[0][i:i + k])
    for kmer in first_strand:
        motifs = [kmer]
        for i in range(1, t):
            profile = laplace_create_profile(motifs)
            motif_i = profile_most_probable(Dna[i], k, profile)
            motifs.append(motif_i)
        if score(motifs) < score(best_motifs):
            best_motifs = motifs
    print best_motifs



#greedy_motif_search(test, 12, 25)

def distance_pattern_strings(pattern, Dna):
    k = len(pattern)
    distance = 0
    for text in Dna:
        ham_dis = float('inf')
        for i in range(len(text)-len(pattern)+1):
            if ham_dis > hamming_distance(pattern, text[i:i+k]):
                ham_dis = hamming_distance(pattern, text[i:i+k])
        distance += ham_dis
    return distance






def score2(kmer, texts):
    score = 0
    for motif in texts:
        scores = []
        for j in range(len(motif)-len(kmer)+1):
            scores.append(hamming_distance(motif[j:j+len(kmer)], kmer))
        score += min(scores)
    return score


def median_string(Dna, k):
    distance = float('inf')
    median = ""
    for kmer in all_kmers_over_one(k):
        if distance > score2(kmer, Dna):
            distance = score2(kmer, Dna)
            median = kmer
            print median
    return median


def profile_most_probable(string, k, profile):
    kmers = []
    probabilities = []
    for i in range(len(string)-k+1):
        kmers.append(string[i:i+k])
        probabilities.append(profile_score(string[i:i+k], profile))
    for i in range(len(probabilities)):
        if probabilities[i] == max(probabilities):
            return string[i:i+k]


def motifs(profile, Dna):
    k = len(profile.get("A"))
    results = []
    for dna in Dna:
        results.append(profile_most_probable(dna, k, profile))
    return results


def randomized_motif_search(Dna, k, t):
        motif_list = []
        for dna in Dna:
            n = random.randrange(0,len(dna)-k+1)
            motif_list.append(dna[n:n+k])
        best_motifs = motif_list[:]
        while t>0:
            profile = laplace_create_profile(motif_list)
            motif_list = motifs(profile, Dna)
            if score(motif_list) < score(best_motifs):
                best_motifs = motif_list[:]
            else:
                return best_motifs


Dna = [
'AGCGGCGGATAGAATACTAAGGGGGCTCACGCTCAGGCAAGAACTGTATGTTTGACACGTCGCTCAGTACGCTTTACCCAGCCGACCAAGAAATTCGGACCAGTCGCTCCTGCCGTTTTGTGCGTAAGGGCCAATGGTAGACGTGATAGCAATTGAAAAGCAGCTTCAAGCGGCGGATAGAAT',
'ACTAAGGGGGCTCACGCTCAGGCAAGAACTGTATGTTTGACACGTCGCTCAGTACGCTTTACCCAGCCGACCAAGAAATTCGGACCAGTCGCTCCTGCCAGTCACGCCTCTACCGTTTTGTGCGTAAGGGCCAATGGTAGACGTGATAGCAATTGAAAAGCAGCTTCAAGCGGCGGATAGAAT',
'AGACCGAGCTACAGCACCGCACATACAAGATCGCTTCGGCCTTTATATAATAATACGAGTACTGTCTCCGTTGCCCGCGTTCGTTTAGTCCTAAAAGTAAATAGGCTGACCCATCTGGGCGGTCAAGGTTCCTTCGCACACACGTGGAAATTAACCGGTTCGTTATAATATGTATTCCAACCC',
'TCAGCACCGGTCCTACGGCAGTACCTCGATTCGCACGGTCCTGAGCGAGTCCCGGGGTAGCGTGCTGCTCTCCTTAGCCAGGCGACGCCTAGCCTCGCGTGACGAGCAGGGAGGTGGTGGCCCCTTAGAGCATACTGTTCGTGAACACCAAGCAACGCCCTGGCTGGCCGTACTTATTGAAAA',
'TCTGCTCAGCGTCCCGCCGTGTGGGATGTTTTGTTCATTGCCTCAGGAGCTGCGTAACGTTTATTGTTTCTCCATAGTTTAATACAGCCTGTACTCTGTACCCAGCAAGACCTCTACCAAGCCTATCGGTAGGTTTTTAGCGCAAAGTATAATTACGCCAATACGACCAGCCAGGGATTGCCC',
'ACAAGTTGCCCGGGATCAATCGCCTATAATGGTCCTAACTCTAGTTTCACCTACGCGGTTCAGCCACTAGAAGATCCACGATCTTGTTGAGAGTTCGAACAAGAGCTAGAGGGGGGGTTCATTACAGCAAGCCCTCTACTGGATGCCGTTGATGCCATAATTTCTTTAATCCATGATTCATCG',
'CCCAGCGGTAATGGACACAGCGCATCAGGCGTAACAGATTCCTAAGTGCTGAGAGAACTACCACTTGTACTTCCTATTGCGTCCTTGCATGAGGTCAAATGGATTCTCGCAAGGCGTCATTTACGGAATCGCTGTCACTAACTTTTGCACCGCCTCTAAGCACAGGGGACTTATCGATCTACC',
'TTTGGGTCCCACTGCCAGTTGAGTGCAAGCACATTGGAGCTGCCTATGCCCAACGGCAGCACCGCCATAACTTTTTTTAGCCGCTGGACCGTTTCACCCAGACGTCATGCCACTGAAAATCGTTGTGAGACGGTGGTCGCTAGACCTCTGTGTTTTCTTGTTCAGTGACTCGGCTTCCAAGGT',
'TGAACAAGATAACCTGCGAGATTCGCCGTACCCGTATTCACAGGCGACCGACGGTCTTCGAAGTGTCCTCCCCTGCAGCACCGCCTGGCCCACTTCTGTGGGCCTTATAAAGCTTACTTTTCATTAGCATAGCTTTGCCGTAAAAGTCTCGGAAAGTAGTATTACGTGTATCGATGTTTTGGA',
'TGTTGGAATTGGCGTGCTGTGCATCACCGTCTATACCGCCTCTACTATTGCAGGCGGGTAAGGGCCAAGTGTACATTGCGCGTCTCTCTAGTTCGTATGCGTTGCCTCCTAGTCGTCAGCCGGCAAGGATCAGCGCAACGGAGTATCTGCGGGTACCACGATCACTTGATGTTGGCGCTTATT',
'GGGTCAACCGACGTGCCGAGCCTTGTGATGGGCGGGAGTACACACGTCAGCGTGATAGGAGCCAAGTCAGCAACAGTTGGACACCTAGTGAAACAGCTCCCGACCCTTTTACCTGCCTAGTATATCAGAGACTAGGCAGCACAAACTCTACGGCATTGGACTTTTCCCGATTACAGTCAATTC',
'ATCTAGGGGACACTCGAAACGTTGGTAGTCTATGGTGTGCACGAGCGGTCCTGCGTGTTCGCATGAGCTTTTTTTGCCTTAGTAGGTTCTTGTAGGGGCTCACGCCTTTTGGATCTGGTAATGAGGTTCGCAGCACCAAGTCTACCCTGTACCTAAGAGATCAACTCTTCGCGTACCCTCAAA',
'GGAAAATAGAGCACCGCCTCTGTGAGACCTGCGCGGCCCTAAAAAGACTCTCCCCTATGATCAGCCGCTCGAGCTATCGAAGACGTCGGGACATTCTTGTGGCGCTTTGTTGAAATTGTTAGTCTACCTGAGACTCCGATGACACATCTGATGCCGTACCAAGAATGAATACCCTCTTCCAAC',
'CTGTCCTGCGACATGTATTGCAAGATAGAACGGCATACTATACGCCGCACTGCTGACACCGCCTCTACTTTACACGGTAGGAGCGGGTTGACGCTTGCCTACCTAACTTCTTATAATACATGCGGCCCTACAGATATTTCTCCGACTGACAAAGCGAATTAGCCCATTCTGACAACTCTAAAA',
'AATGCCAAACCGTTGCTGGCGGCCGAAGAGGGTTCCCCGATTCTCACCTACCTTTTTGTGCTGAATCAACGTCTACGGAACCGAGCTCGAGTGGGGGTCGCTTAACTTAATAGACGCTTCTGTCTCAGCTCCCACGCCGAGGGCGTTTGTGCAGGTCAGCCGAGCCTCTACATGATGTATCTA',
'AGCTTGATGGACCTCCGGGTACGTAATGTCGTTTCTCAGGGACAGGAAGGGAGCGTGCATCGCTAGAGTTTCACGTGGCGCGCGCCGGGGGTTGAAAAGACTAGGGGTCGGCGGGCTCCTAATATGTTGTTCTACACTCAGCCAAGCCTCTACTTGGTGGCGGTGGGGTGCCCTATTGCTGAA',
'TCGTAAACGAACCTCGGTGATTACATAATTGGTACCTTGCTCGTATGCGGCGGGGTTCGGCACTCAAGCCAGCGGTGAAGGAGAGAACGGAGATCTGCAAGTATGGCGCACCACTTCCGCCTCTACTGAGGGTCCGGTTCAACTAATATGTGGGAAAATAATAAAATTAACCCCTTTCAGATA',
'GTTATCCGCAGTGTCAGCACCGCCTCACTGTCTACTAGACTGCAACAGAACCTACTTTACTCCGGTAGGCTGGCAATGCACTTTAATCGCTTTGAGGGGGCGATCTTAGTACCCGGCCGACTTCTACTCTGTTTGGTCCTCCGTAGTGGATTTGTACGATATTACTTGGCCTGCGGAAAAACC',
'CCGAGTACCCCCTGTTTCCACTGAAGCAGCACTAACTCTACGTATTTGCAGCAATAGATATTTCTATAGTACTCTAACCGACATACATGCAGAACAGGTCTATTGACATCTTCTAAGTATGCCATCGTTCGTCTTCCAGACATTTTTTCGCCCACCGGGTGAAGCCTCGTGACTCTTACGGCC',
'AGTTGATGGCCAAGACTCACAGTTTCGCCTCTACCACAGGTCATCTGACTGATGTGAGCTAAGGACTCGAATATCTCGTGGGGCGTTACCAGGAAGTTTTACACGCTCGGCTCACTACACGTTAAGGACTAAAAGGCTCTTCACCTCCCCCTCGATTTGGCCGCGATCGTTTGATAAATGGAA'
]


motif_list = []
scores = []
for i in range(1000):
    current = randomized_motif_search(Dna, 15, 20)
    motif_list.append(current)
    scores.append(score(current))
for i in range(len(scores)):
    if scores[i] == min(scores):
        for motif in motif_list[i]:
            print motif
        break
