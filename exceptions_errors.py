"""
    :platform: Unix
    :synopsis: Define the different exceptions specific to TFFMs that might be
        raised in the process of TFFM analyses.
"""


class TFFMKindError(Exception):
    """
    Define the exception error needed to be raised when a TFFM is instanciated
    with a kind different from '1st-order' and 'detailed'.

    """

    def __init__(self, value):
        """
        Construct the TFFMKindError exception.

        :arg value: Erroneous value given as the kind of a TFFM.
        :type value: str

        """
        self.value = value

    def __str__(self):
        """ Give the string representation of the TFFMKindError exception."""

        return " ".join(["Wrong kind error:", repr(self.value)])
