import re
import random

from math import sqrt

class Sequence:
    """Stores a sequence of characters like DNA.

    The sequence is stored as a list of characters.
    The class also includes methods to alter the sequence.
    """

    def __init__(self, identifier, bases=['A', 'G', 'C', 'T']):
        """Constructor of Sequence class."""

        self.identifier = identifier
        self.sequence = []
        self.bases = bases

    def replicate(self):
        """Returns a new copy of this sequence."""

        new_seq = Sequence(self.identifier)
        new_seq.add_sequence(''.join(self.sequence))

        return new_seq

    def add_sequence(self, seq):
        """Adds string seq to the sequence."""

        self.sequence.extend(list(seq))

    def sample_base(self):
        """Sample a random base from the sequence"""

        r = random.randint(0, len(self.sequence)-1)
        
        return self.sequence[r]

    def create_snp(self, dist_fn=lambda self: self.sample_base(), pos=-1):
        """Creates one single nucleotide polymorphism in the sequence.

        The function takes a distribution function which must return one of the characters
        of the sequence and the position where to create the snp.

        If pos=-1, then the position will be randomly sampled in the sequence. If you want
        to insert a specific base N, use lambda self: 'N'
		"""

        if pos < 0 or pos > len(self.sequence):
            pos = random.randint(0, len(self.sequence)-1)
        
        base = dist_fn(self)

        self.sequence[pos] = base

    def create_insertion(self, mean, std_dev, dist_fn=lambda self: self.sample_base(), pos=-1):
        """Create an insertion of length mean with standard deviation std_dev.

        The inserted bases are generated with the distribution function dist_fn and are
        inserted at position pos.

        If pos=-1, then the position will be randomly sampled in the sequence.
		"""

        length = int(random.gauss(mean, std_dev))

        if pos < 1 or pos > len(self.sequence):
            pos = random.randint(0, len(self.sequence)-1)

        insert = []

        for i in range(0, length):

            insert.append(dist_fn(self))

        self.sequence[pos:pos] = insert

    def create_deletion(self, mean, std_dev, pos=-1):
        """Create a deletion of length mean with standard deviation std_dev.

        The position is pos, if pos is -1, then the position will be randomly
        sampled.
		"""

        length = int(random.gauss(mean, std_dev))

        if pos < 1 or pos > len(self.sequence)-length:
            pos = random.randint(0, len(self.sequence)-length-1)

        self.sequence[pos:pos+length] = []

    def writeFASTA(self, fileobj):
        """Write the sequence on disk in FASTA format."""

        # Write header
        fileobj.write(bytes('>' + self.identifier + '\n', 'UTF-8'))

        # Write sequence
        ct = 0
        for base in self.sequence:
            fileobj.write(bytes(base, 'UTF-8'))
            ct += 1
            if ct == 70:
                ct = 0
                fileobj.write(b'\n')

    def edit_dist(self, compare):
        """Compute the edit distance between this and a given sequence."""

        n = len(self.sequence) + 1
        m = len(compare.sequence) + 1
        
        matrix = [list(range(0, n))]
        
        for i in range(1, m):
            matrix.append(list(range(i, n+i)))

        for j in range(1, n):
            for i in range(1, m):

                if self.sequence[j-1] == compare.sequence[i-1]:
                    matrix[i][j] = matrix[i-1][j-1]
                else:
                    left = matrix[i][j-1] + 1
                    right = matrix[i-1][j] + 1
                    middle = matrix[i-1][j-1] + 1

                    matrix[i][j] = min(left, right, middle)

        return matrix[m-1][n-1]

class ParsingError(Exception):
    """Use this error class, if an error while parsing a (FASTA) file occured"""

    def __init__(self, message):
        self.message = message

    def __str__(self):
        return repr(self.message)

def compute_best_match_score(forward_set, backward_set):
    """Compute the quality of match between the sequences in the two sets.

	Computes the mean root square error from the best edit distances between
    the sequences in the forward set and the backward_set.
	"""

    score = 0

    for seq in forward_set:
        result = seq.edit_dist(backward_set[0])

        for comp in backward_set[1:]:
            res = seq.edit_dist(comp)
            if res < result:
                result = res
        
        score += result*result

    return sqrt(score) / len(forward_set)

def readFASTA(filename):
    """Create all Sequences included in a FASTA file.

    Because a FASTA file can include more than one sequences, this
    function returns a list of sequences.
	"""

    with open(filename, 'r') as f:
        header = f.readline()

        m = re.match(r'>(?P<identifier>\w*)', header)
        if not m:
            raise ParsingError('FASTA header was not valid')

        seq_results = []
        seq = Sequence(m.group('identifier'))

        for line in f:            

            # Look for new sequence header
            if line[0] == '>':
                seq_results.append(seq)
                seq = Sequence(line.rstrip()[1:])
            else:
                seq.add_sequence(line.rstrip())

        # If the last generated sequence is not empty, add it to the list.
        if seq.sequence != '':
            seq_results.append(seq)

    return seq_results
