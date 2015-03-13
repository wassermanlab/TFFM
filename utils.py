"""
Provide utilities.

:platform: Unix
:synopsis: Give utility function for the analysis of TFFMs.

"""
import itertools
import os

from Bio import SeqIO


def set_sequences_weight(sequences, weight):
    """
    Give the same weight *weight* to all the sequences in *sequences*.

    The sequences are constructed using the ghmm module and form a
    :class:`ghmm.SequenceSet`.

    :arg sequences: List of sequences to weight
    :type sequences: :class:`ghmm.SequenceSet`

    """

    for index in xrange(len(sequences)):
        sequences.setWeight(index, weight)


def parse_fasta(fasta_file):
    """
    Parse a fasta file and return the list of SeqRecord items.

    :arg fasta_file: Fasta file containing the sequences.
    :type fasta_file: str

    :returns: The set of Bio.SeqRecord instances found in the fasta file.
    :rtype: list of :class:`Bio.SeqRecord`

    """

    assert(os.path.isfile(fasta_file))
    # We do not use "with" for python2.4 compatibility
    with open(fasta_file) as stream:
        record_list = []
        for record in SeqIO.parse(stream, "fasta"):
            record_list.append(record)
        stream.close()
        return record_list


def roundrobin(*iterables):
    """
    Create a generator interlacing the iterables given in argument.

    For example, roundrobin('ABC', 'D', 'EF') --> A D E B F C

    :arg iterables: Pointer to the iterables.

    :returns: The generator over the intarlaced iterables.
    :rtype: :class:`Generator`

    """

    # Recipe credited to George Sakkis and found in the itertools documentation
    pending = len(iterables)
    nexts = itertools.cycle(iter(it).next for it in iterables)
    while pending:
        try:
            for the_next in nexts:
                yield the_next()
        except StopIteration:
            pending -= 1
            nexts = itertools.cycle(itertools.islice(nexts, pending))


def get_sequences_info(seq_file):
    """
    Get the number of sequences, number of residues, and starting nucleotides.

    :arg seq_file: File containing the sequences in fasta format
    :type seq_file: str

    :returns: The number of sequences, number of residues, and starting
        nucleotides occurrences
    :rtype: tuple of (int, int, dic of str->int)

    """

    # Check file and raise error if does not exist.
    # We do not use "with" for python2.4 compatibility
    if not os.path.exists(seq_file):
        # We assume MEME-ChIP was used
        return 600, 60600, {'A': 150, 'C': 150, 'G': 150, 'T': 150}
    nb_seq = 0
    nb_residues = 0
    starts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    stream = open(seq_file)
    for record in SeqIO.parse(stream, "fasta"):
        nb_seq += 1
        nb_residues += len(record.seq)
        if record.seq[0] in starts:
            starts[record.seq[0]] += 1
    stream.close()
    return nb_seq, nb_residues, starts
