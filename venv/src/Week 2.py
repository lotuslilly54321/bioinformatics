import itertools
import time
from collections import defaultdict


def format(string):
    string = string.replace(',', "")
    string = string.replace("'","")
    string = string.replace("[","")
    string = string.replace("]", "")
    print string


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


def minimum_skew(genome):
    total_skew = skew(genome)
    min_skew = total_skew[1]
    minimums = []
    for i in range(1,len(total_skew)):
        if total_skew[i]<= min_skew:
            min_skew = total_skew[i]
    for i in range(1,len(total_skew)):
        if total_skew[i] == min_skew:
            minimums.append(i)
    return minimums


def maximum_skew(genome):
    total_skew = skew(genome)
    max_skew = total_skew[1]
    maximums = []
    for i in range(1,len(total_skew)):
        if total_skew[i]>= max_skew:
            max_skew = total_skew[i]
    for i in range(1,len(total_skew)):
        if total_skew[i] == max_skew:
            maximums.append(i)
    return maximums

def hamming_distance(str1, str2):
    mismatch = 0
    if len(str1)>len(str2):
        for i in range(1,len(str2)+1):
            if str1[len(str1)-i] != str2[len(str2)-i]:
               mismatch = mismatch + 1
    else:
        for i in range(1,len(str1)+1):
            if str1[len(str1)-i] != str2[len(str2)-i]:
               mismatch = mismatch + 1
    return mismatch


def approx_match(pattern, genome, max_mm):
    matches = []
    for i in range(len(genome)-len(pattern)+1):
        if hamming_distance(pattern, genome[i:i+len(pattern)]) <= max_mm:
            matches.append(i)
    return matches


def count(genome, pattern, max_mm):
    matches = []
    for i in range(len(genome) - len(pattern) + 1):
        if hamming_distance(pattern, genome[i:i + len(pattern)]) <= max_mm:
            matches.append(i)
    return len(matches)


def frequent_words_with_mismatches(text, k, max_mm):
    frequencies = defaultdict(lambda:0)
    frequent = []
    for i in range(len(text)-k+1):
        kmer_and_friends = PermuteMotifDistanceTimes(text[i:i+k],max_mm)
        for kmer in kmer_and_friends:
            frequencies[kmer] += 1
    for kmer in frequencies:
        if frequencies[kmer] == max(frequencies.values()):
            frequent.append(kmer)
    return frequent


def frequent_words_with_mismatches_and_revcomp(text, k, max_mm):
    frequencies = defaultdict(lambda:0)
    frequent = []
    for i in range(len(text)-k+1):

        pattern = text[i:i+k]
        neighbors = PermuteMotifDistanceTimes(pattern,max_mm)
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


def immediate_neighbors(pattern):
    tides = ["A","T","C","G"]
    neighborhood = []
    for i in range(len(pattern)):
        symbol = pattern[i]
        for nuc in tides:
            if nuc != symbol:
                neighbor = pattern[0:i]+ nuc + pattern[i+1:]
                neighborhood.append(neighbor)
    return neighborhood


def neighbors(pattern, d):
    tides = ["A","T","C","G"]
    if d == 0:
        return pattern
    if len(pattern) == 1:
        return ["A", "C", "G", "T"]
    neighborhood = [pattern]
    suffix_neighbors = neighbors(suffix(pattern),d)
    for text in suffix_neighbors:
        if hamming_distance(text, pattern) == d:
            neighborhood.append(pattern[0]+text)
        if hamming_distance(text, pattern) < d:
            for nuc in tides:
                if nuc != pattern[0]:
                    neighborhood.append(nuc+text)

    return neighborhood


def suffix(pattern):
    return pattern[1:]

answer = neighbors("CATGCTGATTAA",2)

print(answer)
format("['CATGCTGATTAA', 'AATGCTGATTAA', 'TATGCTGATTAA', 'GATGCTGATTAA', 'ATTGCTGATTAA', 'TTTGCTGATTAA', 'GTTGCTGATTAA', 'ACTGCTGATTAA', 'TCTGCTGATTAA', 'GCTGCTGATTAA', 'AGTGCTGATTAA', 'TGTGCTGATTAA', 'GGTGCTGATTAA', 'CTAGCTGATTAA', 'CCAGCTGATTAA', 'CGAGCTGATTAA', 'CTCGCTGATTAA', 'CCCGCTGATTAA', 'CGCGCTGATTAA', 'CTGGCTGATTAA', 'CCGGCTGATTAA', 'CGGGCTGATTAA', 'CAAACTGATTAA', 'CACACTGATTAA', 'CAGACTGATTAA', 'CAATCTGATTAA', 'CACTCTGATTAA', 'CAGTCTGATTAA', 'CAACCTGATTAA', 'CACCCTGATTAA', 'CAGCCTGATTAA', 'CATAATGATTAA', 'CATTATGATTAA', 'CATCATGATTAA', 'CATATTGATTAA', 'CATTTTGATTAA', 'CATCTTGATTAA', 'CATAGTGATTAA', 'CATTGTGATTAA', 'CATCGTGATTAA', 'CATGAAGATTAA', 'CATGTAGATTAA', 'CATGGAGATTAA', 'CATGACGATTAA', 'CATGTCGATTAA', 'CATGGCGATTAA', 'CATGAGGATTAA', 'CATGTGGATTAA', 'CATGGGGATTAA', 'CATGCAAATTAA', 'CATGCCAATTAA', 'CATGCGAATTAA', 'CATGCATATTAA', 'CATGCCTATTAA', 'CATGCGTATTAA', 'CATGCACATTAA', 'CATGCCCATTAA', 'CATGCGCATTAA', 'CATGCTATTTAA', 'CATGCTTTTTAA', 'CATGCTCTTTAA', 'CATGCTACTTAA', 'CATGCTTCTTAA', 'CATGCTCCTTAA', 'CATGCTAGTTAA', 'CATGCTTGTTAA', 'CATGCTCGTTAA', 'CATGCTGTATAA', 'CATGCTGCATAA', 'CATGCTGGATAA', 'CATGCTGTCTAA', 'CATGCTGCCTAA', 'CATGCTGGCTAA', 'CATGCTGTGTAA', 'CATGCTGCGTAA', 'CATGCTGGGTAA', 'CATGCTGAAAAA', 'CATGCTGACAAA', 'CATGCTGAGAAA', 'CATGCTGAACAA', 'CATGCTGACCAA', 'CATGCTGAGCAA', 'CATGCTGAAGAA', 'CATGCTGACGAA', 'CATGCTGAGGAA', 'CATGCTGATATA', 'CATGCTGATCTA', 'CATGCTGATGTA', 'CATGCTGATACA', 'CATGCTGATCCA', 'CATGCTGATGCA', 'CATGCTGATAGA', 'CATGCTGATCGA', 'CATGCTGATGGA', 'CATGCTGATTTC', 'CATGCTGATTCC', 'CATGCTGATTGC', 'CATGCTGATTTG', 'CATGCTGATTCG', 'CATGCTGATTGG', 'CATGCTGATTTT', 'CATGCTGATTCT', 'CATGCTGATTGT']")
print(len(""))
