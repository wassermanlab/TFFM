# TFFM-framework version 2.0 - 03/05/2013

We provide here the documentation of the TFFM-framework developed in Python.
The Transcription Factor Flexible Models (TFFMs) represent TFBSs and are based
on hidden Markov models (HMM). They are flexible and are able to model both
position interdependence within TFBSs and variable length motifs within a
single dedicated framework.

The framework also implements methods to generate a new graphical
representation of the modeled motifs that convey properties of position
interdependences.

TFFMs have been assessed on ChIP-seq data sets coming from the ENCODE project,
revealing that the new HMM-based framework performs, in most cases, better than
both PWMs and the dinucleotide weight matrix (DWM) extension in discriminating
motifs within ChIP-seq sequences from background sequences. Under the
assumption that ChIP-seq signal values are correlated with the affinity of the
TF-DNA binding, we find that TFFM scores correlate with ChIP-seq peak signals.
Moreover, using available TF-DNA affinity measurements for the Max TF, we
observe that TFFMs constructed from ChIP-seq data correlate with published
experimentally measured DNA-binding affinities. These results demonstrate the
capacity of TFFMs to accurately model DNA-protein interactions, while providing
a single unified framework suitable for the next generation of TFBS
predictions. All the details have been published in Mathelier and Wasserman,
The next generation of transcription binding site prediction.

TFFMs can be saved and opened from files using the XML format already used by
the GHMM library.

We recommend you to read the documentation to get more information at
[http://cisreg.cmmt.ubc.ca/TFFM/doc/](http://cisreg.cmmt.ubc.ca/TFFM/doc/).

# Authors

  Anthony Mathelier and Wyeth W. Wasserman

  Centre for Molecular Medicine and Therapeutics

  950 West 28th Avenue

  Vancouver, BC

  V5Z 4H4 Canada

  anthony.mathelier@ncmm.uio.no

  wyeth@cmmt.ubc.ca

  Refer to the AUTHORS file.
         

# System requirements

* The TFFM-framework 2.0 has been developed and tested under Ubuntu
    Linux operating system. It has also been tested on CentOS.
* Python should be installed (version 2.7 has been used successfully).
* Biopython (at least version 1.61) should be installed and accessible from
    your Python executable. See http://biopython.org for instructions on how to
    install it.
* The GHMM library should be installed and accessible from Python. See
    http://ghmm.org for instructions on how to install it.

# Tutorial

For a brief tutorial on how to use the TFFMs, go to
[http://cisreg.cmmt.ubc.ca/TFFM/doc/](http://cisreg.cmmt.ubc.ca/TFFM/doc/).

# Download

The TFFM-framework can be downloaded on
[github](https://github.com/wassermanlab/TFFM).

# Web application

You can generate TFFMs and scan DNA sequences with TFFMs using the dedicated
web application at
[http://cisreg.cmmt.ubc.ca/TFFM/webapp/](http://cisreg.cmmt.ubc.ca/TFFM/webapp/).

# Reference

The TFFMs and the TFFM-framework have been described in [Mathelier and
Wasserman. The next generation of transcription binding site prediction. _Plos
Computational Biology_, 2013](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003214).

Please cite this publication when using the TFFMs.


# Licence

* The TFFM-framework has been developed under the [GNU Lesser General
    Public Licence](http://www.gnu.org/copyleft/lesser.html).
* The TFFM-framework uses the GHMM library which is also licenced
    under the GNU LGPL.
