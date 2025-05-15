class seqUtils():
    '''
    A collection of sequence utility methods and attributes
    '''

    @classmethod
    def rc_class(cls, seq):
        r_string = seq[::-1]
        complement_dict = {"A": "T",
                           "T": "A",
                           "C": "G",
                           "G": "C"}
        rc_string = ''.join(complement_dict.get(char, char) for char in r_string)
        return rc_string

    # Constructor method
    def __init__(self, seq):
        # Input DNA sequence
        self.seq = seq
        self.len = len(seq)
        self.rc_called = False

    # Generate sequence reverse complement
    def rc(self):
        r_string = self.seq[::-1]
        complement_dict = {"A": "T",
                           "T": "A",
                           "C": "G",
                           "G": "C"}
        self.rc_string = ''.join(complement_dict.get(char, char) for char in r_string)
        self.rc_called = True

        return self.rc_string

    # Trim up to the first n bases
    def trim(self, n):
        # Check whether the sequence is as long as the trimming n request then trim up to n
        if self.rc_called:
            if n >= self.len:
                return self.rc_string
            else:
                return self.rc_string[:n]
        else:
            if n >= self.len:
                return self.seq
            else:
                return self.seq[:n]
            
    def rmPolyA(self):
        '''
        Remove polyA tail from the sequence
        Alternative implementation of cutadapt polyA trimming algorithm with improvements.
        '''
        n = len(self.seq)
        best_index = n
        best_score = score = 0
        total_errors = 0  # Track errors across all suffixes

        for i in reversed(range(n)):
            if self.seq[i] == "A":
                score += 1
            else:
                score -= 2
                total_errors += 1

            suffix_length = n - i
            error_rate = total_errors / suffix_length if suffix_length > 0 else 0 # Handle empty suffix

            if score > best_score and error_rate <= 0.2:
                best_index = i
                best_score = score


        return self.seq[:best_index]
    
    def is_dna(self):
        ref = {'A', 'C', 'T', 'G'}
        query = set(self.seq)
        if query.issubset(ref):
            return True
        else:
            return False

    def __str__(self):
        return f"Sequence: {self.seq}\nLength: {self.len}"
