"""
Module implementing the TFFMs.

:platform: Unix
:synopsis: Define the class representing the Transcription Factor Flexible
    Models and the necessary functions to manipulate them.
:todo: Allow the construction of TFFMs using a different *de novo* motif
    finding tool than MEME.

"""


# TODO Allow the construction of TFFMs using a different *de novo* motif
# finding tool.


import sys
import os
import math
import re
from Bio.Alphabet import IUPAC
from Bio import SeqIO

#
# Changed to use the newer Bio.motifs package
# DJA 2013/04/05
#
#from Bio.Motif.Parsers import MEME
from Bio import motifs

import ghmm

import drawing
import exceptions_errors
import utils
import hit_module

from constants import ALPHABET, EXTENDED_ALPHABET, TFFM_KIND, LOGO_TYPE


class TFFM(ghmm.DiscreteEmissionHMM):
    """
    Define the Transcription Factor Flexible Models.

    .. note::
        Instances of this class have to be created through the
        functions :func:`tffm_from_xml` or :func:`tffm_from_meme`.

    """

    def __init__(self, emission_domain, distribution, cmodel, kind,
                 name="TFFM"):
        """
        Construct an instance of the TFFM class.

        :arg emission_domain: The emission domain of the
            underlying :class:`ghmm.HMM`.
        :type emission_domain: ghmm.EmissionDomain
        :arg distribution: The distribution over the emission domain.
        :type distribution: :class:`ghmm.Distribution`
        :arg cmodel: The cmodel (HMM itself implemented in C) of the underlying
            :class:`ghmm.HMM`.
        :arg kind: The TFFM can be either a 0-order, a 1st-order, or a detailed
            TFFM, use `TFFM_KIND.ZERO_ORDER, or `TFFM_KIND.FIRST_ORDER`, or
            `TFFM_KIND.DETAILED` respectively.
        :type kind: Enum
        :arg name: Give the name of the TFFM. 'TFFM' is given by default.
        :type name: str
        :raises: :class:`exceptions.TFFMKindError` when the given kind is
            neither '1st-order' nor 'detailed'.

        """

        # Construct the underlying ghmm.EmissionHMM
        super(ghmm.DiscreteEmissionHMM, self).__init__(emission_domain,
                                                       distribution, cmodel)
        if(kind != TFFM_KIND.FIRST_ORDER and kind != TFFM_KIND.DETAILED and
           kind != TFFM_KIND.ZERO_ORDER):
            raise exceptions_errors.TFFMKindError(kind)
        self.kind = kind
        self.name = name

    def __del__(self):
        """
        Delete the underlying C structures.

        :note: The destruction is made using
            the :class:`ghmm.DiscreteEmissionHMM` destructor.

        """

        super(ghmm.DiscreteEmissionHMM, self).__del__()

    def __len__(self):
        """
        Give the length of the TFFM, i.e. the number of nucleotides in the
        model excluding the background.

        """

        if self.kind == TFFM_KIND.FIRST_ORDER:
            return self.N - 2
        elif self.kind == TFFM_KIND.DETAILED:
            return self.N / 4 - 1
        else:  # 0-order HMM here
            return self.N - 1

    def background_emission_proba(self):
        """
        Return the emission probabilities of the nucleotides in the background
        state.

        :returns: A dictionnary with characters 'A', 'C', 'G', and 'T' as keys
            and the corresponding probabilities as values.
        :rtype: dict

        """

        emissions = {'A': 0., 'C': 0., 'G': 0., 'T': 0.}
        if self.kind == TFFM_KIND.FIRST_ORDER:
            emissions = background_emission_proba_1storder(self)
        elif self.kind == TFFM_KIND.DETAILED:
            emissions = background_emission_proba_detailed(self)
        else:  # 0-order HMM here
            emissions = background_emission_proba_detailed(self)
            for i in xrange(4):
                emissions[ALPHABET[i]] = self.getEmission(0)[i]
        return emissions

    def print_summary_logo(self, output=sys.stdout):
        """
        Print the svg code of the corresponding summary logo (i.e. similar to a
        regular sequence logo).

        :arg output: Stream where to output the svg
            (defaut: :class:`sys.stdout`).
        :type output: file

        :note: The *output* argument is not a file name but it is an
            **already** open file stream.

        """

        drawing.draw_logo(output, self, LOGO_TYPE.SUMMARY)

    def print_dense_logo(self, output=sys.stdout):
        """
        Print the svg code of the corresponding dense logo (i.e. displaying the
        dinucleotide dependencies captured by the TFFM).

        :arg output: Stream where to output the svg
            (defaut: :class:`sys.stdout`).
        :type output: file

        :note: The *output* argument is not a file name but it is an
            **already** open file stream.

        """

        # No dependencies in 0order TFFMs
        if self.kind == TFFM_KIND.ZERO_ORDER:
            raise exceptions_errors.TFFMKindError(self.kind)
        drawing.draw_logo(output, self, LOGO_TYPE.DENSE)

    def get_positions_ic(self):
        """
        Give the information content for every positions of the motif modeled
        by the TFFM.

        :returns: A list of floats giving the information contents of the
            positions.
        :rtype: list

        :note: The output is an ordered list following the order of the
            positions within the motif.

        """

        previous_position_proba = self.background_emission_proba()
        positions_information_content = []
        # Get the position of the first matching state in the TFFM.
        start = self.get_position_start()
        # Iterate over all matching position in the TFFM.
        for pos in xrange(start, len(self) + start):
            position_proba = {'A': 0., 'C': 0., 'G': 0., 'T': 0.}
            # Compute the emission probabilities for each nucleotide at the
            # current position
            if (self.kind == TFFM_KIND.DETAILED
                    or self.kind == TFFM_KIND.FIRST_ORDER):
                for i in xrange(4):
                    emissions = self.get_emission_update_pos_proba(
                        position_proba, pos, previous_position_proba, i)
            else:  # 0-order HMM here
                emissions = self.get_emission_update_pos_proba(position_proba,
                                                               pos,
                                                               previous_position_proba,
                                                               0)
            previous_position_proba = position_proba.copy()
            if self.kind == TFFM_KIND.DETAILED:
                somme = sum(previous_position_proba.values())
                # For a detailed TFFM, one position is composed of 4 states so
                # we need to use proportions to compute proba.
                for i in xrange(4):
                    previous_position_proba[ALPHABET[i]] /= somme
            # Compute the entropy given the corresponding probabilites and then
            # the underlying information content
            values = previous_position_proba.values()
            keys = previous_position_proba.keys()
            emissions = zip(values, keys)
            entropy = compute_entropy(emissions)
            positions_information_content.append(2.0 - entropy)
        return positions_information_content

    def get_information_content(self):
        """
        Give the information content of the whole TFFM.

        :returns: A float corresponding to the information content of the TFFM.
        :rtype: float

        """

        pos_ic = self.get_positions_ic()
        return sum(pos_ic)

    def final_states(self):
        """
        Give the list of final states in the HMM (i.e. corresponding to the
        last matching position in the TFFM).

        :returns: A list of final states as int.
        :rtype: list

        """

        if self.kind == TFFM_KIND.FIRST_ORDER:
            return [len(self) + 1]
        elif self.kind == TFFM_KIND.DETAILED:
            size = len(self) * 4
            return [size, size + 1, size + 2, size + 3]
        else:  # 0-order HMM here
            return [len(self)]

    def train(self, training_file, epsilon=0.0001, max_iter=500):
        """
        Train the TFFM using the fasta sequences to learn emission and
        transition probabilities.

        :note: The training of the underlying HMM is made using the Baum-Welsh
            algorithm.

        :arg training_file: The fasta file of the sequences to train the TFFM
            on.
        :type training_file: str
        :arg epsilon: The least relative improvement cut-off in likelihood
            compared to the previous iteration of the Baum-Welsh algorithm
            (default: 0.0001).
        :type epsilon: float
        :arg max_iter: The maximum number of iteration of the Baum-Welsh
            algorithm to reestimate the probabilities (default: 500).
        :type max_iter: int

        """

        assert(os.path.isfile(training_file))
        # Only upper case is allowed in the ALPHABET, need to convert
        sequences = []
        for record in SeqIO.parse(training_file, "fasta"):
            sequence = record.seq.upper()
            # Only considering sequences with ACGTs
            if not re.search("[^AGCT]", str(sequence)):
                sequences.append(sequence)
        training_sequences = ghmm.SequenceSet(ghmm.Alphabet(ALPHABET),
                                              sequences)
        # Need to give the same weight to all the sequences since it does not
        # seem to be done by default by ghmm.
        utils.set_sequences_weight(training_sequences, 1.0)
        self.baumWelch(training_sequences, max_iter, epsilon)

    def scan_sequence(self, sequence, threshold=0.0, only_best=False):
        """
        Apply the TFFM on the fasta sequence and return the TFBS hits.

        :arg sequence: DNA sequence to apply the TFFM on.
        :type sequence: :class:`Bio.SeqRecord`
        :arg threshold: The threshold used to predict a hit (i.e. the minimal
            probability value for a position to be considered a TFBS hit)
            (default: 0.0).
        :type threshold: float
        :arg only_best: Argument to be set to :class:`True` if only the best
            TFBS hit per sequence is to be reported (default: :class:`False`)
        :type only_best: bool

        :returns: TFBS hits.
        :rtype: list of :class:`HIT`

        :note:  (**0.0<=** *threshold* **<=1.0**)

        """

        # Retrieve the hits on both the positive and the negative strands.
        hits_positive_strand = self._get_hits(sequence, threshold)
        hits_negative_strand = self._get_hits(sequence, threshold,
                                              negative=True)
        return merge_hits(hits_positive_strand, hits_negative_strand,
                          only_best)

    def scan_sequences(self, seq_file, threshold=0., only_best=False):
        """
        Apply the TFFM on the fasta sequences and return the TFBS hits.

        :arg seq_file: Fasta file giving the DNA sequences to apply the TFFM
            on.
        :type seq_file: str
        :arg threshold: The threshold used to predict a hit (i.e. the minimal
            probability value for a position to be considered a TFBS hit)
            (default: 0.0).
        :type threshold: float
        :arg only_best: Argument to be set to :class:`True` if only the best
            TFBS hit per sequence is to be reported (default: :class:`False`)
        :type only_best: bool
        :returns: TFBS hits through a generator.
        :rtype: :class:`Generator` of :class:`HIT`

        :note:  (**0.0<=** *threshold* **<=1.0**)

        """

        sequence_list = utils.parse_fasta(seq_file)
        for seq_record in sequence_list:
            hits = self.scan_sequence(seq_record, threshold, only_best)
            for hit in hits:
                yield hit

    def pocc_sequences(self, seq_file, threshold=0.):
        """
        Apply the TFFM on the fasta sequences and return the Pocc value
        (probability of occupancy) for each sequence.

        :arg seq_file: Fasta file giving the DNA sequences to apply the TFFM
            on.
        :type seq_file: str
        :arg threshold: The threshold used to predict hits that will be used to
            compute the Pocc (default: 0.0).
        :type threshold: float

        :returns: Pocc values through a generator.
        :rtype: :class:`Generator` of :class:`HIT`

        :note: (**0.0<=** *threshold* **<=1.0**)

        """

        sequence_list = utils.parse_fasta(seq_file)
        for seq_record in sequence_list:
            pocc = 1.
            for hit in self.scan_sequence(seq_record, threshold, False):
                if hit:
                    pocc *= (1. - hit.score)
            yield hit_module.HIT(seq_record, 1, len(seq_record), None, pocc,
                                 self, None)

    # Not sure it is a good habit to trim in-place
    def trim_in_place(self, threshold):
        """
        Trim the current TFFM by removing edges with low information content.

        :arg threshold: The minimal information content value for an edge TFFM
            match position to be kept.
        :type threshold: float

        :warning: Trims the TFFM in place. To preserve the TFFM, use
            the :func:`get_trimmed` method which returns a trimmed copy of the
            TFFM but does not alter this TFFM.

        :see also: :func:`get_trimmed`

        """

        first, last = self.get_significant_positions(threshold)

        # Only the underlying HMM is to be trimmed by trimming the
        # corresponding transition and emission probability matrices
        new_hmm = self._get_trimmed_hmm(first, last)
        # The trimming is made in-place so the instance attributes have to be
        # updated
        self.cmodel = new_hmm.cmodel
        self.N = new_hmm.cmodel.N
        self.M = new_hmm.cmodel.M
        self.model_type = new_hmm.cmodel.model_type
        self.maxorder = new_hmm.cmodel.maxorder

    def get_trimmed(self, threshold, new_name="TFFM"):
        """
        Trim the current TFFM by removing edges with low information content.

        :arg threshold: The minimal information content value for an edge TFFM
            match position to be kept.
        :type threshold: float
        :arg new_name: Name of the new TFFM to create (default:'TFFM').
        :type new_name: str

        :returns: A TFFM corresponding to the current TFFM trimmed.
        :rtype: :class:`TFFM`

        :see also: :func:`trim_in_place`

        """

        first, last = self.get_significant_positions(threshold)
        # Only the underlying HMM is to be trimmed by trimming the
        # corresponding transition and emission probability matrices
        new_hmm = self._get_trimmed_hmm(first, last)
        return TFFM(new_hmm.emissionDomain, new_hmm.distribution,
                    new_hmm.cmodel, new_name, self.kind)

    def _get_trimmed_hmm(self, first, last):
        """
        Return the new trimmed HMM.

        :arg first: Position of the new first matching position.
        :type first: int
        :arg last: Position of the new last matching position.
        :type last: int
        :returns: The new trimmed HMM.
        :rtype: :class:`ghmm.DiscreteEmissionHMM`

        :todo: Raise an error rather than a :func:`sys.exit` when the trimmed
            HMM becomes empty.
        """

        nb_states = last - first + 1
        if nb_states < 0:
            # TODO raise an error rather than a sys.exit
            sys.exit("\nTrimmed HMM is empty!!!\n")
        # Retrieve the matrices of the current HMM
        matrices = self.asMatrices()
        transitions = matrices[0]
        emissions = matrices[1]
        # Trim the transition matrices and update the emission matrices wrt the
        # new first and last positions by iterating through all the states
        for index in xrange(self.N):
            transitions[index] = transitions[index][0:nb_states + 2]
            if index > 1 and index < nb_states + 2:
                emissions[index] = emissions[first + 1]
                first += 1
        # Only keep the probabilities lying in the new matching positions
        transitions = transitions[0:nb_states + 2]
        transitions[nb_states + 1][nb_states + 1] = 0.0
        transitions[nb_states + 1][1] = 1.0
        emissions = emissions[0:nb_states + 2]
        initials = matrices[2][0:nb_states + 2]
        return ghmm.HMMFromMatrices(self.emissionDomain, self.distribution,
                                    transitions, emissions, initials)

    def get_significant_positions(self, threshold):
        """
        Get the first and last significant position the TFFM where the
        insignificant positions are the ones on the edges with low information
        content.

        :arg threshold: The minimal information content to consider a position
            to be significant.
        :type threshold: float

        :returns: The positions of the first and last positions that are to be
            considered significant (given in this order).
        :rtype: tuple

        """

        pos_ic = self.get_positions_ic()
        first = 1
        last = len(self)
        while first <= last and pos_ic[first - 1] < threshold:
            first += 1
        while last > 0 and pos_ic[last - 1] < threshold:
            last -= 1
        return first, last

    def _get_posterior_proba(self, sequence_split):
        """
        Get the posterior probabilities at each nucleotide position given the
        TFFM.

        :arg sequence_split: The sequence splitted in subsequences to not
            consider non ACGT nucleotides.
        :type sequence_split: list

        :returns: The posterior probabilities at each position of the sequence.
        :rtype: list of list

        :note: One example of a sequence_split is ["ACT", "N", "ATC"].

        """

        ghmm_extended_alphabet = ghmm.Alphabet(EXTENDED_ALPHABET)
        posterior_proba = []
        # null probabilities for non ACGT nucleotides.
        null_proba = [0.] * self.N
        for sequence in sequence_split:
            if re.match("[ACGT]", sequence):
                emission_sequence = ghmm.SequenceSet(ghmm_extended_alphabet,
                                                     [sequence])[0]
                posterior_proba.extend(self.posterior(emission_sequence))
            else:
                for __ in xrange(len(sequence)):
                    posterior_proba.append(null_proba)
        return posterior_proba

    def _get_hits(self, seq_record, threshold, negative=False):
        """
        Predict TFBS hits in the sequence given the TFFM.

        :arg seq_record: The sequence on which to predict TFBS hits.
        :type seq_record: :class:`Bio.SeqRecord`
        :arg threshold: The minimal probability to predict a position as a TFBS
            hit.
        :type threshold: float
        :arg negative: A boolean stating if the TFBS hits are to be predicted
            on the positive or the negative strand of the sequence. Set to True
            when on the negative strand (default: False).
        :type negative: bool

        :returns: The list of TFBS hits predicted on the sequence strand.
        :rtype: list of :class:`HIT`

        """

        sequence = seq_record.seq.upper()
        if negative:
            sequence = seq_record.reverse_complement().seq
        sequence_split = re.split("([ACGT]+)", str(sequence))
        posterior_proba = self._get_posterior_proba(sequence_split)
        hits = self._construct_hits(posterior_proba, seq_record, threshold,
                                    negative)
        return hits

    def get_emission_update_pos_proba(self, position_proba, position,
                                      previous_position_proba, index):
        """
        Get the emission probabilities of ACGT at position *position* and
        update the emission probabilities in *position_proba* given the
        emission probabilities at the previous position
        (*previous_position_proba*).

        :note: This function is used state by state and several states
            represent the same position in detailed TFFM, this is why we need
            to update the probabilities listed in position_proba.

        :arg position_proba: Probabilities of getting ACGT at the current
            position that need to be updated.
        :type position_proba: dict
        :arg position: Current position in the motif.
        :type position: int
        :arg previous_position_proba: Probabilities of getting ACGT at the
            previous position.
        :type previous_position_proba: dict
        :arg index: Represents the index of the state of the TFFM to be
            analyzed at the current position.

        :returns: The emission probabilities of ACGT by the state indexed by
            *index* at position *position* in the TFFM.
        :rtype: list

        """

        emissions_dic = {'A': 0., 'C': 0., 'G': 0., 'T': 0.}
        if self.kind == TFFM_KIND.FIRST_ORDER:
            for j in xrange(4):
                letter = ALPHABET[j]
                emissions_dic[letter] = self.getEmission(position)[index * 4 + j]
                proba = self.getEmission(position)[j * 4 + index]
                proba *= previous_position_proba[letter]
                position_proba[ALPHABET[index]] += proba
            emissions = zip(emissions_dic.values(), emissions_dic.keys())
        elif self.kind == TFFM_KIND.DETAILED:
            for j in xrange(4):
                start_state = (position - 1) * 4 + index
                end_state = position * 4 + j
                emissions_dic[ALPHABET[j]] = self.getTransition(start_state,
                                                                end_state)
                proba = previous_position_proba[ALPHABET[index]]
                proba *= self.getTransition(start_state, end_state)
                position_proba[ALPHABET[j]] += proba
            somme = sum(emissions_dic.values())
            emissions = zip(emissions_dic.values(), emissions_dic.keys())
            emissions = map(lambda (x, y): (x / somme, y), emissions)
        else:  # 0-order HMM here
            for j in xrange(4):
                letter = ALPHABET[j]
                emissions_dic[letter] = self.getEmission(position)[j]
                proba = self.getEmission(position)[j]
                position_proba[ALPHABET[j]] += proba
            emissions = zip(emissions_dic.values(), emissions_dic.keys())
        emissions.sort(reverse=True)
        return emissions

    def _construct_hits(self, posterior_proba, seq_record, threshold,
                        negative):
        """
        Compute the TFBS hits on a sequence given the posterior probabilities
        and construct the corresponding instances of
        :class:`HIT`.

        :arg posterior_proba: The posterior probabilities at each position of
            the sequence computed given the tffm.
        :type posterior_proba: list of list of float
        :arg seq_record: The sequence on which to predict TFBS hits.
        :type seq_record: :class:`Bio.SeqRecord`
        :arg threshold: The minimal probability to predict a position as a TFBS
            hit.
        :type threshold: float
        :arg negative: A boolean stating if the TFBS hits are to be predicted
            on the positive or the negative strand of the sequence. Set to True
            when on the negative strand.
        :type negative: bool

        :returns: The list of TFBS hits predicted on the sequence strand.
        :rtype: list of :class:`HIT`

        """

        hits = [None] * len(seq_record)
        # Iterate over all the positions of the sequences by starting at the
        # len(self) position since we need to read at least len(self)
        # nucleotides to predict a potential TFBS hit.
        for position in xrange(len(self) - 1, len(posterior_proba)):
            best_proba = 0.0
            # Iterate over the final states of the TFFM to predict TFBS hits by
            # looking at the posterior probabilities for these states only at
            # this specific position
            for state in self.final_states():
                proba = posterior_proba[position][state]
                # We want to only keep the final state giving the maximal
                # posterior probability.
                if proba > threshold and proba > best_proba:
                    best_proba = proba
                    start, end, strand = hit_module.get_start_end_strand(
                        position, seq_record, self, negative)
                    hits[end - 1] = hit_module.HIT(seq_record, start, end,
                                                   strand, proba, self, state)
        return hits

    def get_position_start(self):
        """
        Give the position of the first matching state.

        :returns: The position of the first matching state of the TFFM.
        :rtype: float

        :warning: The position is given 0-based.

        """

        if self.kind == TFFM_KIND.FIRST_ORDER:
            return 2
        else:  # Both detailed and 0-order
            return 1


def background_emission_proba_1storder(tffm):
    emissions = {'A': 0., 'C': 0., 'G': 0., 'T': 0.}
    last_emissions = {'A': 0., 'C': 0., 'G': 0., 'T': 0.}
    # Retrieve emission proba for the first state which is not
    # 1st-order
    for i in xrange(4):
        last_emissions[ALPHABET[i]] = tffm.getEmission(0)[i]
    # Compute emission proba for the background
    for i in xrange(4):
        for j in xrange(4):
            proba = tffm.getEmission(1)[j * 4 + i]
            proba *= last_emissions[ALPHABET[j]]
            emissions[ALPHABET[i]] += proba
    return emissions


def background_emission_proba_detailed(tffm):
    emissions = {'A': 0., 'C': 0., 'G': 0., 'T': 0.}
    # Compute emission proba for the background
    for i in xrange(4):
        for j in xrange(4):
            proba = tffm.getTransition(i, j) * 0.25
            emissions[ALPHABET[j]] += proba
    somme = sum(emissions.values())
    for i in xrange(4):
        # Need to divide by somme since four states are representing
        # the same position within a detailed TFFM
        emissions[ALPHABET[i]] /= somme
    return emissions


def best_hit(hits_positive_strand, hits_negative_strand):
    """
    Give the best hit in a sequence by considering both positive and negative
    strands.

    :arg hits_positive_strand: The list of hits on the positive strand.
    :type hits_positive_strand: list of :class:`HIT`
    :arg hits_negative_strand: The list of hits on the negative strand.
    :type hits_negative_strand: list of :class:`HIT`

    :returns: The best hit (None if no hit).
    :rtype: :class:`HIT`

    """

    max_positive_strand = max(hits_positive_strand)
    max_negative_strand = max(hits_negative_strand)
    if max_positive_strand or max_negative_strand:
        return max(max_positive_strand, max_negative_strand)
    else:
        return None


def merge_hits(hits_positive_strand, hits_negative_strand, only_best):
    """
    Merges the hits from both strands.

    :arg hits_positive_strand: The list of hits on the positive strand.
    :type hits_positive_strand: list of :class:`HIT`
    :arg hits_negative_strand: The list of hits on the negative strand.
    :type hits_negative_strand: list of :class:`HIT`
    :arg only_best: Boolean set to True only if the best TFBS hit in the
        sequence is to be kept.
    :type only_best: bool

    :returns: A list containing the TFBS hits (empty if no hit).
    :rtype: list

    :note: The two input lists are required to be ordered following the
        positions on the sequence.  The best hit per position is given. When
        no hit has been found at a position, the constant None is used.

    """

    if only_best:
        return [best_hit(hits_positive_strand, hits_negative_strand)]
    else:
        return [hit for hit in utils.roundrobin(hits_positive_strand,
                                                hits_negative_strand)]


def compute_entropy(emissions):
    """
    Compute the entropy given the emission probabilities of the ACGT
    nucleotides.

    :arg emissions: Emission probabilities of the ACGT nucleotides.
    :type emissions: list of float

    :returns: The computed entropy.
    :rtype: float

    :warning: The list gives the probabilities corresponding to A, C, G, and T
        **in this order**.

    """

    entropy = 0.
    for i in xrange(4):
        proba, __ = emissions[i]
        if proba != 0.:
            entropy += proba * math.log(proba, 2)
    return -entropy


def tffm_from_xml(xml, kind):
    """
    Construct a TFFM described in an XML file.

    :arg xml: File containing the TFFM description in XML format.
    :type xml: str
    :arg kind: Type of TFFM to construct between '1st-order' and 'detailed'.
    :type kind: str

    :returns: The TFFM described in the XML file.
    :rtype: :class:`TFFM`

    """

    hmm = ghmm.HMMOpen(xml)
    return TFFM(hmm.emissionDomain, hmm.distribution, hmm.cmodel, kind)


def create_1storder_hmm(nb_seq, nb_residues, first_letters, motif):
    """
    Create a 1st-order HMM initialized from MEME result

    :arg nb_seq: Number of sequences used by MEME
    :type nb_seq: int
    :arg nb_residues: Number of residues used by MEME
    :type nb_residues: int
    :arg first_letters: Number of occurrences of ACGT at the begining of
        sequences used by MEME
    :type first_letters: dic of str->int
    :arg motif: PFM as a Biopython motif to be used to initialize the TFFFM
    :type motif: :class:`Bio.motifs`

    :returns: The constructed HMM
    :rtype: :class:`ghmm.DiscreteEmissionHMM`

    """

    # Very first state is created with the initial nt frequencies and a '1'
    # pseudocount
    emissions = [[(first_letters['A'] + 1.) / (nb_seq + 4.),
                  (first_letters['C'] + 1.) / (nb_seq + 4.),
                  (first_letters['G'] + 1.) / (nb_seq + 4.),
                  (first_letters['T'] + 1.) / (nb_seq + 4.)]]
    # The second state is random
    emissions.append([0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25,
                      0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25])
    # Complete the emissions with the actual motif frequencies
    if motif.instances:
        # The motif.counts is computed directly when creating the motif from
        # instances
        nb_hits = len(motif.instances)
    else:
        nb_hits = nb_seq
    for position in xrange(len(motif)):
        frequencies = []
        for letter in "ACGT":
            freq = (motif.counts[letter][position] + 1.) / (nb_hits + 4.)
            frequencies.append(freq)
        emissions.append(frequencies * 4)

    # Background transitions
    transitions = [[0.0, 1.0] + [0.0] * len(motif)]
    background_to_background = 1. - float(nb_seq) / nb_residues
    background_to_foreground = 1. - background_to_background
    # Core transitions
    transitions.append(
        [0., background_to_background, background_to_foreground] + [0.] *
        (len(motif) - 1))
    for position in xrange(1, len(motif)):
        transitions.append(
            [0.] * (position + 2) + [1.] + [0.] * (len(motif) - position - 1))
    # Final transitions now
    transitions.append([0., 1.] + [0.] * len(motif))

    # Starting proba
    initials = [1.] + [0.] * (len(motif) + 1)
    return ghmm.HMMFromMatrices(ghmm.Alphabet(ALPHABET),
                                ghmm.DiscreteDistribution(
                                    ghmm.Alphabet(ALPHABET)),
                                transitions, emissions, initials)


def create_0order_hmm(nb_seq, nb_residues, first_letters, motif):
    """
    Create a 0-order HMM initialized from MEME result

    :arg nb_seq: Number of sequences used by MEME
    :type nb_seq: int
    :arg nb_residues: Number of residues used by MEME
    :type nb_residues: int
    :arg first_letters: Number of occurrences of ACGT at the begining of
        sequences used by MEME
    :type first_letters: dic of str->int
    :arg motif: PFM as a Biopython motif to be used to initialize the TFFFM
    :type motif: :class:`Bio.motifs`

    :returns: The constructed HMM
    :rtype: :class:`ghmm.DiscreteEmissionHMM`

    """

    # The first state is random
    emissions = [[0.25, 0.25, 0.25, 0.25]]
    # Complete the emissions with the actual motif frequencies
    if motif.instances:
        # The motif.counts is computed directly when creating the motif from
        # instances
        nb_hits = len(motif.instances)
    else:
        nb_hits = nb_seq
    for position in xrange(len(motif)):
        frequencies = []
        for letter in "ACGT":
            freq = (motif.counts[letter][position] + 1.) / (nb_hits + 4.)
            frequencies.append(freq)
        emissions.append(frequencies)

    # Background transitions
    transitions = []
    background_to_background = 1. - float(nb_seq) / nb_residues
    background_to_foreground = 1. - background_to_background
    transitions.append(
        [background_to_background, background_to_foreground] + [0.] *
        (len(motif) - 1))
    # Core transitions
    for position in xrange(1, len(motif)):
        transitions.append(
            [0.] * (position + 1) + [1.] + [0.] * (len(motif) - position - 1))
    # Final transitions now
    transitions.append([1.] + [0.] * len(motif))

    # Starting proba
    initials = [1.] + [0.] * len(motif)
    return ghmm.HMMFromMatrices(ghmm.Alphabet(ALPHABET),
                                ghmm.DiscreteDistribution(
                                    ghmm.Alphabet(ALPHABET)),
                                transitions, emissions, initials)


def create_detailed_hmm(nb_seq, nb_residues, first_letters, motif):
    """
    Create a detailed HMM initialized from MEME result

    :arg nb_seq: Number of sequences used by MEME
    :type nb_seq: int
    :arg nb_residues: Number of residues used by MEME
    :type nb_residues: int
    :arg first_letters: Number of occurrences of ACGT at the begining of
        sequences used by MEME
    :type first_letters: dic of str->int
    :arg motif: PFM as a Biopython motif to be used to initialize the TFFFM
    :type motif: :class:`Bio.motifs`

    :returns: The constructed HMM
    :rtype: :class:`ghmm.DiscreteEmissionHMM`

    """

    # Starting proba with the starting nucleotides of the sequences and a '1'
    # pseudocount added
    initials = [(first_letters['A'] + 1.) / (nb_seq + 4.),
                (first_letters['C'] + 1.) / (nb_seq + 4.),
                (first_letters['G'] + 1.) / (nb_seq + 4.),
                (first_letters['T'] + 1.) / (nb_seq + 4.),
                1. / nb_seq, 1. / nb_seq, 1. / nb_seq, 1. / nb_seq]
    initials += [0.] * 4 * (len(motif) - 1)

    # Emission proba
    emissions = [[1., 0., 0., 0.], [0., 1., 0., 0.], [0., 0., 1., 0.],
                 [0., 0., 0., 1.]] * (len(motif) + 1)

    # Background transitions proba
    if motif.instances:
        # The motif.counts is computed directly when creating the motif from
        # instances
        nb_hits = len(motif.instances)
    else:
        nb_hits = nb_seq
    background_to_background = 1. - float(nb_seq) / nb_residues
    background_to_foreground = 1. - background_to_background
    background_to_background /= 4.
    transi = {}
    for letter in "ACGT":
        freq = (motif.counts[letter][0] + 1.) / (nb_hits + 4.)
        transi[letter] = freq * background_to_foreground
    transitions = []
    for __ in xrange(4):
        transitions.append([background_to_background,
                            background_to_background,
                            background_to_background,
                            background_to_background,
                            transi['A'], transi['C'],
                            transi['G'], transi['T']]
                           + [0.] * 4 * (len(motif) - 1))
    pfm = [(v + 1.) / (nb_hits + 4.) for letter in 'ACGT' for v in motif.counts[letter]]
    for position in xrange(len(motif) * 4):
        transitions.append([0.] * 4 * (len(motif) + 1))
    for position in xrange(1, len(motif)):
        for line in xrange(4 * (position - 1) + 1, 4 * (position - 1) + 5):
            for column in xrange(4 * position + 1, 4 * position + 5):
                index = (column - (4 * position + 1)) * len(motif) + position
                transitions[line + 3][column + 3] = pfm[index]
    for index in xrange(4):
        state = len(motif) * 4 + index
        transitions[state][0] = 0.25
        transitions[state][1] = 0.25
        transitions[state][2] = 0.25
        transitions[state][3] = 0.25
    return ghmm.HMMFromMatrices(ghmm.Alphabet(ALPHABET),
                                ghmm.DiscreteDistribution(
                                    ghmm.Alphabet(ALPHABET)),
                                transitions, emissions, initials)


def tffm_from_meme(meme_output, kind, name="TFFM"):
    """
    Construct a TFFM from the output of MEME on ChIP-seq data.

    :arg meme_output: File containing the output of MEME.
    :type meme_output: str
    :arg kind: Type of TFFM to construct between '1st-order' and 'detailed'.
    :type kind: str
    :arg name: Name of the TFFM (default: "TFFM")
    :type name: str

    :returns: The TFFM initialized from MEME results.
    :rtype: :class:`TFFM`

    :note: As the PFM is used to initialize the TFFM, a pseudocount of 1 is
        added to all the values in the PFM

    """

    record = motifs.parse(open(meme_output), 'MEME')
    if record.alphabet != IUPAC.unambiguous_dna:
        sys.exit("### Wrong alphabet used in MEME ###")
    motif = record[0]
    nb_seq, nb_res, first_letters = utils.get_sequences_info(record.datafile)
    if kind == TFFM_KIND.FIRST_ORDER:
        hmm = create_1storder_hmm(nb_seq, nb_res, first_letters, motif)
    elif kind == TFFM_KIND.DETAILED:
        hmm = create_detailed_hmm(nb_seq, nb_res, first_letters, motif)
    else:  # 0-order HMM here
        hmm = create_0order_hmm(nb_seq, nb_res, first_letters, motif)
    return TFFM(hmm.emissionDomain, hmm.distribution, hmm.cmodel, kind, name)


def tffm_from_motif(motif, kind, name="TFFM", nb_res=None, first_letters=None):
    """
    Construct an initialized TFFM from a PFM.

    :arg motif: PFM as a Biopython motif to be used to initialize the TFFFM
    :type motif: :class:`Bio.motifs`
    :arg kind: Type of TFFM to construct between '1st-order' and 'detailed'
    :type kind: str
    :arg name: Name of the TFFM (default: "TFFM")
    :type name: str
    :arg nb_res: Number of residues to be used to compute the
        background->foreground transition probabilities
        (default: 0 meaning that we will assume a 100nt length for the ChIP-seq
        used to derive the PFM)
    :type nb_residues: int
    :arg first_letters: Number of occurrences of ACGT at the begining of
        sequences in the background (default: default values will give
        equiprobabilities for A, C, G, and T
    :type first_letters: dic of str->int

    :returns: The TFFM initialized from the PFM
    :rtype: :class:`TFFM`

    :see also: :func:`tffm_from_meme`

    :note: As the PFM is used to initialize the TFFM, a pseudocount of 1 is
        added to all the values in the PFM

    """

    # extract the counts in the first column of the PFM
    first_counts = [count[0] for count in motif.counts.values()]
    nb_seq = reduce(lambda x, y: x + y, first_counts)
    if not first_letters:
        equi = nb_seq / 4
        first_letters = {'A': equi, 'C': equi, 'G': equi, 'T': equi}
    if not nb_res:
        nb_res = nb_seq * 100
    if kind == TFFM_KIND.FIRST_ORDER:
        hmm = create_1storder_hmm(nb_seq, nb_res, first_letters, motif)
    elif kind == TFFM_KIND.DETAILED:
        hmm = create_detailed_hmm(nb_seq, nb_res, first_letters, motif)
    else:  # 0-order TFFM here
        hmm = create_0order_hmm(nb_seq, nb_res, first_letters, motif)
    return TFFM(hmm.emissionDomain, hmm.distribution, hmm.cmodel, kind, name)
