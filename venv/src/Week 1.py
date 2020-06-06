def format(string):
    string = string.replace(',', "")
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



