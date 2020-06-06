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
    matches = False
    for i in range(len(genome)-len(pattern)+1):
        if hamming_distance(pattern, genome[i:i+len(pattern)]) <= max_mm:
            matches = True
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


def all_kmers_over_one(d):
    workingSet = "A"
    for i in range(1,d):
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


# takes list of dna and length and returns all common motifs within the dna and within d distance
def motif_enumeration(Dna, k, d):
    patterns = []
    for dna in Dna:
        for i in range(len(dna)-k+1):
            km = dna[i:i+k]
            for kmer in PermuteMotifDistanceTimes(km,d):
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
    for i in range(len(genome)-len(pattern)+1):
        if hamming_distance(pattern,genome[i:i+len(pattern)]) < min:
            min = hamming_distance(pattern,genome[i:i+len(pattern)])
    return min


def median_string(Dna, k):
    distance = float('inf')
    median = ""
    for kmer in all_kmers_over_one(k):
        for dna in Dna:
            if distance > score2(kmer, Dna):
                distance = score2(kmer, Dna)
                median = kmer
    return median




def profile_kmer(text, k, profile):
    probabilities = []
    kmers = []
    for i in range(len(text)-k+1):
        probabilities.append(profile_score(text[i:i+k], profile))
    for i in range(len(probabilities)):
        if probabilities[i] == max(probabilities):
            kmers.append(text[i:i+k])
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
            profile[n].append(float(value)/total)
            print float(value)/total
    return profile


def greedy_motif_search(Dna, k, t):
    first_strand = []
    best_motifs = []    #        BestMotifs to motif matrix formed by first k-mers in each string in Dna
    for i in range(t):
        best_motifs.append(Dna[i][0:k])
    for i in range(len(Dna[0])-k+1):
        first_strand.append(Dna[0][i:i+k])
    for kmer in first_strand:
        motifs = [kmer]
        for i in range(1,t):
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
    score = 1
    for motif in motifs:
        score += hamming_distance(motif,find_consensus_profile(profile))
    return score


def score2(kmer, motifs):
    score = 1
    for motif in motifs:
        score += hamming_distance(motif,kmer)
    return score


def profile_score(kmer, profile):
    score = 1.0
    for i in range(len(kmer)):
        index = profile.get(kmer[i])
        score = score * index[i]
    return score


test0 = ["GGCGTTCAGGCA",
"AAGAATCAGTCA",
"CAAGGAGTTCGC",
"CACGTCAATCAC",
"CAATAATATTCG"]


test1 = ["GCCCAA",
"GGCCTG",
"AACCTA",
"TTCCTT"]


test2 = ["GAGGCGCACATCATTATCGATAACGATTCGCCGCATTGCC",
"TCATCGAATCCGATAACTGACACCTGCTCTGGCACCGCTC",
"TCGGCGGTATAGCCAGAAAGCGTAGTGCCAATAATTTCCT",
"GAGTCGTGGTGAAGTGTGGGTTATGGGGAAAGGCAGACTG",
"GACGGCAACTACGGTTACAACGCAGCAACCGAAGAATATT",
"TCTGTTGTTGCTAACACCGTTAAAGGCGGCGACGGCAACT",
"AAGCGGCCAACGTAGGCGCGGCTTGGCATCTCGGTGTGTG",
"AATTGAAAGGCGCATCTTACTCTTTTCGCTTTCAAAAAAA"
]


test3 = ["GCAGGTTAATACCGCGGATCAGCTGAGAAACCGGAATGTGCGT",
         "CCTGCATGCCCGGTTTGAGGAACATCAGCGAAGAACTGTGCGT",
         "GCGCCAGTAACCCGTGCCAGTCAGGTTAATGGCAGTAACATTT",
         "AACCCGTGCCAGTCAGGTTAATGGCAGTAACATTTATGCCTTC",
         "ATGCCTTCCGCGCCAATTGTTCGTATCGTCGCCACTTCGAGTG"]


test4 = ["GACCTACGGTTACAACGCAGCAACCGAAGAATATTGGCAA",
"TCATTATCGATAACGATTCGCCGGAGGCCATTGCCGCACA",
"GGAGTCTGGTGAAGTGTGGGTTATGGGGCAGACTGGGAAA",
"GAATCCGATAACTGACACCTGCTCTGGCACCGCTCTCATC",
"AAGCGCGTAGGCGCGGCTTGGCATCTCGGTGTGTGGCCAA",
"AATTGAAAGGCGCATCTTACTCTTTTCGCTTAAAATCAAA",
"GGTATAGCCAGAAAGCGTAGTTAATTTCGGCTCCTGCCAA",
"TCTGTTGTTGCTAACACCGTTAAAGGCGGCGACGGCAACT"]


Dna = ["GGCGTTCAGGCA",
       "AAGAATCAGTCA",
       "CAAGGAGTTCGC",
       "CACGTCAATCAC",
       "CAATAATATTCG"]


motifs = [
    "TCGGGGGTTTTT",
    "CCGGTGACTTAC",
    "ACGGGGATTTTC",
    "TTGGGGACTTTT",
    "AAGGGGACTTCC",
    "TTGGGGACTTCC",
    "TCGGGGATTCAT",
    "TCGGGGATTCCT",
    "TAGGGGAACTAC",
    "TCGGGTATAACC",
]

profile = {'A': [0.2, 0.2, 0.0, 0.0, 0.0, 0.0, 0.9, 0.1, 0.1, 0.1, 0.3, 0.0],
           'C': [0.1, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.4, 0.1, 0.2, 0.4, 0.6],
           'T': [0.7, 0.2, 0.0, 0.0, 0.1, 0.1, 0.0, 0.5, 0.8, 0.7, 0.3, 0.4],
           'G': [0.0, 0.0, 1.0, 1.0, 0.9, 0.9, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0]}


profile = create_profile(test4)
print profile

format("")

def profile_most_probable(string, k, profile):
    scores = []
    all = []
    mins = []
    for i in range((len(string)-k)+1):
        score = profile_score(string[i:i+k],profile)
        scores.append(score)
        all.append(string[i:i+k])
    for i in range(len(scores)):
        if scores[i] == max(scores):
            mins.append(string[i:i+k])
    return scores[len(scores)-2],scores[len(scores)-1]

print profile_most_probable("AAAAAAAAGAGGC", 5, profile)
