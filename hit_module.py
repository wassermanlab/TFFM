"""
    :platform: Unix
    :synopsis: Create the class HIT which represents a TFBS hit.
"""


class HIT:
    """ Define the representation of a TFBS hit. """

    # TODO raise an error when the strand is wrong
    def __init__(self, seq_record, start, end, strand, score, tffm=None,
                 tffm_matched_state=-1):
        """
        Create an instance of the :class:`HIT`.

        :arg seq_record: Sequence containing the hit.
        :type seq_record: :class:`Bio.SeqRecord`
        :arg start: Start position of the hit.
        :type start: int
        :arg end: End position of the hit.
        :type end: int
        :arg strand: Strand of the hit on the sequence. It should be either '+'
            or '-').
        :type strand: str
        :arg score: TFFM score of the hit.
        :type score: float
        :arg tffm: TFFM used to predict the hit (default: None).
        :type tffm: :class:`TFFM`
        :arg tffm_matched_state: Matching state in the TFFM to predict the hit
            (default: -1, i.e. no TFFM).
        :type tffm_matched_state: int

        :warning: start and end are 1-based.
        :warning: The seq_record attribute is not the actual sequence of the
            hit but the whole sequence containing the hit.

        :todo: Raise an error when the strand is wrong.

        """

        self.seq_record = seq_record
        self.start = start
        self.end = end
        self.strand = strand
        self.score = score
        self.tffm_matched_state = tffm_matched_state
        self.tffm = tffm

    def __len__(self):
        """
        Give the length of the TFBS hit.

        :returns: The length of the TFBS hit.
        :rtype: int
        """

        return self.end - self.start + 1

    def __str__(self):
        """
        Give the string representation of the TFBS hit.

        :returns: The string representing the hit in the following format:
            start\tend\tstrand\tsequence\ttffm-name\ttffm-state\tscore
        :rtype: str

        """

        if self.tffm:
            name = self.tffm.name
        else:
            name = "NoName"
        string = "%s\t%d\t%d\t%s\t%s\t%s\t%d\t%s" % (self.seq_record.id,
                                                     self.start, self.end,
                                                     self.strand,
                                                     self.sequence(),
                                                     name,
                                                     self.tffm_matched_state,
                                                     repr(self.score))
        return string

    def __lt__(self, other):
        """
        Implement the **<** operator. The comparison looks at the score.
        """

        if other:
            return self.score < other.score
        else:
            return False

    def __le__(self, other):
        """
        Implement the **<=** operator. The comparison looks at the score.
        """

        if other:
            return self.score <= other.score
        else:
            return False

    def __eq__(self, other):
        """
        Implement the **==** operator. The comparison looks at the score.
        """

        if other:
            return self.score == other.score
        else:
            return False

    def __ne__(self, other):
        """
        Implement the **!=** operator. The comparison looks at the score.
        """

        if other:
            return self.score != other.score
        else:
            return True

    def __gt__(self, other):
        """
        Implement the **>** operator. The comparison looks at the score.
        """

        if other:
            return self.score > other.score
        else:
            return True

    def __ge__(self, other):
        """
        Implement the **>=** operator. The comparison looks at the score.
        """

        if other:
            return self.score >= other.score
        else:
            return True

    def sequence(self):
        """
        Give the sequence of the TFBS hit.

        :returns: The sequence of the TFBS hit.
        :rtype: str

        """

        seq = self.seq_record.seq[self.start - 1:self.end]
        if self.strand == "+" or not self.strand:
            return seq
        else:
            return seq.reverse_complement()


def get_start_end_strand(position, seq_record, tffm, negative):
    """
    Get the start and end positions of a TFBS hit given its end position on the
    positive strand.

    :arg position: End position of the TFBS hit on the positive strand of the
        sequence.
    :type position: int
    :arg seq_record: The actual sequence.
    :type seq_record: :class:`Bio.SeqRecord`
    :arg tffm: The TFFM used to predict the TFBS hit.
    :type tffm: :class:`TFFM`
    :arg negative: Boolean set to True if the TFBS hit is on the negative
        strand of the sequence and False otherwise.

    :returns: The start and end positions and the strand.
    :rtype: tuple(int, int, str)

    :note: The strand is '+' if the hit is on the positive strand and '-'
        otherwise.
    :warning: The input *position* is given 0-based as extracted from TFFM
        computations but the output start and end are 1-based since it is a
        more conventionnal way to print the TFBS hit positions.

    """

    # position is 0-based and we need 1-based coordinates
    if negative:
        start = len(seq_record) - position
        end = len(seq_record) - position + len(tffm) - 1
        strand = "-"
    else:
        end = position + 1
        start = position - len(tffm) + 2
        strand = "+"
    return start, end, strand
