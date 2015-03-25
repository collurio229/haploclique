import re
from math import sqrt

class Sequence:

    def __init__(self, identifier, bases=['A', 'G', 'C', 'T']):
        self.identifier = identifier
        self.sequence = []
        self.bases = bases

        self.base_freqs = dict()
        for b in self.bases:
            self.base_freqs[b] = 0.25

    def replicate(self):
        new_seq = Sequence(self.identifier)
        new_seq.add_sequence(''.join(self.sequence))

        return new_seq

    def add_sequence(self, seq):
        old_len = len(self.sequence)
        self.sequence.extend(list(seq))

        length = len(self.sequence)

        new_len = len(seq)

        new_freq = dict()
        for key in self.base_freqs.keys():
            new_freq[key] = 0

        for c in seq:
            new_freq[c] += 1
            
        for (base, freq) in new_freq.items():
            self.base_freqs[base] = (old_len * self.base_freqs[base] + freq) / length

    def sample_base(self):
        r = random.random()
        ct = 0

        for (base, freq) in self.base_freqs.items():
            ct += freq
            if r < ct:
                return base

    def create_snp(self, dist_fn=lambda self: self.sample_base(), pos=-1, base='N'):
        if pos < 0 or pos > len(self.sequence):
            pos = random.randint(0, len(self.sequence)-1)

        if base == 'N':
            base = dist_fn(self)

        self.sequence[pos] = base

    def create_insertion(self, mean, std_deviation, pos=-1):
        length = int(random.gauss(mean, std_deviation))

        if pos < 1 or pos > len(self.sequence):
            pos = random.randint(0, len(self.sequence)-1)

        insert = []

        for i in range(0, length):

            insert.append(self.sample_base())

        self.sequence[pos:pos] = insert

    def create_deletion(self, mean, std_deviation, pos=-1):
        length = int(random.gauss(mean, std_deviation))

        if pos < 1 or pos > len(self.sequence)-length:
            pos = random.randint(0, len(self.sequence)-length-1)

        self.sequence[pos:pos+length] = []

    def writeFASTA(self, fileobj):
        fileobj.write(bytes('>' + self.identifier + '\n', 'UTF-8'))
        ct = 0
        for base in self.sequence:
            fileobj.write(bytes(base, 'UTF-8'))
            ct += 1
            if ct == 70:
                ct = 0
                fileobj.write(b'\n')

    def edit_dist(self, compare):
        n = len(self.sequence) + 1
        m = len(compare.sequence) + 1
        
        matrix = [list(range(0, n))]
        
        for i in range(1, m):
            matrix.append(list(range(i, n+i)))

        for j in range(1, n):
            for i in range(1, m):
                if self.sequence[i-2] == compare.sequence[j-2]:
                    matrix[i][j] = matrix[i-1][j-1]
                else:
                    left = matrix[i][j-1] + 1
                    right = matrix[i-1][j] + 1
                    middle = matrix[i-1][j-1] + 1

                    matrix[i][j] = min(left, right, middle)

        for line in matrix:
            print(line)

        print(''.join(self.sequence))
        print(''.join(compare.sequence))

        return matrix[m-1][n-1]

class ParsingError(Exception):
    def __init__(self, message):
        self.message = message

    def __str__(self):
        return repr(self.message)

def compute_best_match_score(forward_set, backward_set):
    
    score = 0

    for seq in forward_set:
        result = seq.edit_dist(forward_set[0])

        for comp in backward_set[1:]:
            res = seq.edit_dist(comp)
            if res < result:
                result = res
        
        score += result*result

    return sqrt(score)

def readFASTA(filename):
    with open(filename, 'r') as f:
        header = f.readline()

        m = re.match(r'>(?P<identifier>\w*)', header)
        if not m:
            raise ParsingError('FASTA header was not valid')

        seq_results = []
        seq = Sequence(m.group('identifier'))

        for line in f:            

            if line[0] == '>':
                print('matched:', line)
                seq_results.append(seq)
                seq = Sequence(line.rstrip()[1:])
            else:

                seq.add_sequence(line.rstrip())
                print('added:', line.rstrip())

        if seq.sequence != '':
            seq_results.append(seq)

    return seq_results
