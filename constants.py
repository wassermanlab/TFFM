"""
    :platform: Unix
    :synopsis: Define different constants needed in the modules.

"""


BLACK = (0, 0, 0)
""" Defines the color black in RGB used in logos. """
BLUE = (0, 0, 128)
""" Defines the color blue in RGB used in logos. """
WHITE = (255, 255, 255)
""" Defines the color white in RGB used in logos. """
RED = (178, 34, 34)
""" Defines the color red in RGB used in logos. """
GREEN = (0, 100, 0)
""" Defines the color green in RGB used in logos. """
YELLOW = (255, 165, 0)
""" Defines the color yellow in RGB used in logos. """
ALPHABET = ['A', 'C', 'G', 'T']
""" Defines the nucleotide alphabet. """
EXTENDED_ALPHABET = ['A', 'C', 'G', 'T', 'N', 'M', 'K']
""" Defines the extended nucleotide alphabet. """


def enum(**enums):
    """ Create an enumerated type. """
    return type('Enum', (), enums)


TFFM_KIND = enum(FIRST_ORDER="1st-order", DETAILED="detailed",
                 ZERO_ORDER="0-order")
""" Defines the three kinds of TFFMs (1st-order, detailed, 0-order). """

LOGO_TYPE = enum(SUMMARY="summary", DENSE="dense")
""" Defines the two kinds of logos. """
