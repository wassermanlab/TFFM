.. TFFM documentation master file, created by
   sphinx-quickstart on Wed Jan  9 17:19:05 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to TFFM's documentation!
================================

We provide here the documentation of the TFFM-framework developed in `Python
<http://docs.python.org/2.7/>`_.  The **Transcription Factor Flexible Models
(TFFMs)** represent TFBSs and are based on hidden Markov models (HMM). They are
flexible and are able to model both position interdependence within TFBSs and
variable length motifs within a single dedicated framework.

The framework also implements methods to generate a new graphical
representation of the modeled motifs that convey properties of position
interdependences.

TFFMs have been assessed on ChIP-seq data sets coming from the `ENCODE project
<http://genome.ucsc.edu/ENCODE/>`_, revealing that the new HMM-based framework
performs, in most cases, better than both PWMs and the dinucleotide weight
matrix (DWM) extension in discriminating motifs within ChIP-seq sequences from
background sequences. Under the assumption that ChIP-seq signal values are
correlated with the affinity of the TF-DNA binding, we find that TFFM scores
correlate with ChIP-seq peak signals. Moreover, using available TF-DNA affinity
measurements for the Max TF, we observe that TFFMs constructed from ChIP-seq
data correlate with published experimentally measured DNA-binding affinities.
These results demonstrate the capacity of TFFMs to accurately model DNA-protein
interactions, while providing a single unified framework suitable for the next
generation of TFBS predictions. All the details have been published in
**Mathelier and Wasserman, The Next Generation of Transcription Binding Site
Prediction,** *PLOS Computational Biology* , Sept. 2013, 9(9):e1003214,
`DOI:10.1371/journal.pcbi.1003214 <http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1003214>`_.


TFFMs can be saved and opened from files using the XML format already used by
the `GHMM library <http://ghmm.org>`_.

System requirements
===================

* The TFFM-framework 2.0 has been developed and tested under `Ubuntu Linux
  <http://www.ubuntu.com/>`_ operating system. It has also been tested on
  `CentOS <http://www.centos.org>`_.
* **Python** should be installed (version 2.7 has been used
  successfully).
* **Biopython** (at least version 1.61) should be installed and accessible from
  your Python executable.  See http://biopython.org for instructions on how to
  install it.
* The **GHMM** library should be installed and accessible from Python. See
  http://ghmm.org for instructions on how to install it.

Contents
========

.. toctree::
   :maxdepth: 4

   tffm_module
   hit_module
   drawing
   utils
   constants
   exceptions_errors

Tutorial
========

Data for this tutorial can be found in the tutorial directory of the
TFFM-framework package at
https://github.com/wassermanlab/TFFM/tree/master/tutorial. To run this example,
you need to be located in the tutorial directory. The location to be added to
`sys.path` (second line of the tutorial) corresponds to the directory
containing the TFFM-framework.

>>> import sys
>>> sys.path.append(<TFFM-framework repository path>)
>>> import tffm_module
>>> from constants import TFFM_KIND
>>> 
>>> tffm_first_order = tffm_module.tffm_from_meme("meme.txt", TFFM_KIND.FIRST_ORDER)
>>> tffm_first_order.write("tffm_first_order_initial.xml")
>>> tffm_first_order.train("train.fa")
>>> tffm_first_order.write("tffm_first_order.xml")
>>> out = open("tffm_first_order_summary_logo.svg", "w")
>>> tffm_first_order.print_summary_logo(out)
>>> out.close()
>>> out = open("tffm_first_order_dense_logo.svg", "w")
>>> tffm_first_order.print_dense_logo(out)
>>> out.close()
>>> 
>>> tffm_detailed = tffm_module.tffm_from_meme("meme.txt", TFFM_KIND.DETAILED)
>>> tffm_detailed.write("tffm_detailed_initial.xml")
>>> tffm_detailed.train("train.fa")
>>> tffm_detailed.write("tffm_detailed.xml")
>>> out = open("tffm_detailed_summary_logo.svg", "w")
>>> tffm_detailed.print_summary_logo(out)
>>> out.close()
>>> out = open("tffm_detailed_dense_logo.svg", "w")
>>> tffm_detailed.print_dense_logo(out)
>>> out.close()
>>> 
>>> tffm_first_order = tffm_module.tffm_from_xml("tffm_first_order.xml",
...         TFFM_KIND.FIRST_ORDER)
>>> print "1st-order all"
1st-order all
>>> for hit in tffm_first_order.scan_sequences("test.fa"):
...     if hit:
...         print hit
... 
hg19_ct_UserTrack_3545_1911 1   15  -   CCCAATCTGGGGGAG TFFM    16  2.0412237861770574e-05
hg19_ct_UserTrack_3545_1911 2   16  -   TCCCAATCTGGGGGA TFFM    16  5.972631515438409e-10
hg19_ct_UserTrack_3545_1911 3   17  +   CCCCCAGATTGGGAG TFFM    16  7.302219996694176e-10
hg19_ct_UserTrack_3545_1911 3   17  -   CTCCCAATCTGGGGG TFFM    16  4.283052089644644e-10
hg19_ct_UserTrack_3545_1911 4   18  +   CCCCAGATTGGGAGA TFFM    16  2.3867153510744474e-11
hg19_ct_UserTrack_3545_1911 4   18  -   TCTCCCAATCTGGGG TFFM    16  1.809727472566856e-09
...
...
hg19_ct_UserTrack_3545_3739 84  98  +   AGGAAGGCAGTCTGA TFFM    16  1.4846886217808363e-13
hg19_ct_UserTrack_3545_3739 84  98  -   TCAGACTGCCTTCCT TFFM    16  1.221285402074912e-07
hg19_ct_UserTrack_3545_3739 85  99  +   GGAAGGCAGTCTGAG TFFM    16  4.6592155264443255e-12
hg19_ct_UserTrack_3545_3739 85  99  -   CTCAGACTGCCTTCC TFFM    16  9.886433151452122e-09
hg19_ct_UserTrack_3545_3739 86  100 +   GAAGGCAGTCTGAGC TFFM    16  4.299143127628496e-07
hg19_ct_UserTrack_3545_3739 87  101 +   AAGGCAGTCTGAGCC TFFM    16  3.147869508275758e-08
>>> print "1st-order best"
1st-order best
>>> for hit in tffm_first_order.scan_sequences("test.fa", only_best=True):
...     if hit:
...         print hit
... 
hg19_ct_UserTrack_3545_1911 33  47  -   TAGGCCCGGAGGGAG TFFM    16  0.7487725826821764
hg19_ct_UserTrack_3545_3739 29  43  +   CCTGCCCCCAGGGTG TFFM    16  0.9670956290505575
>>> tffm_detailed = tffm_module.tffm_from_xml("tffm_detailed.xml",
...         TFFM_KIND.DETAILED)
>>> print "detailed all"
detailed all
>>> for hit in tffm_detailed.scan_sequences("test.fa"):
...     if hit:
...         print hit
... 
hg19_ct_UserTrack_3545_1911 1   15  +   CTCCCCCAGATTGGG TFFM    62  1.2342473696772854e-08
hg19_ct_UserTrack_3545_1911 1   15  -   CCCAATCTGGGGGAG TFFM    62  6.845206246155106e-06
hg19_ct_UserTrack_3545_1911 2   16  +   TCCCCCAGATTGGGA TFFM    60  9.158988688768735e-10
hg19_ct_UserTrack_3545_1911 2   16  -   TCCCAATCTGGGGGA TFFM    60  2.6048552465107907e-11
hg19_ct_UserTrack_3545_1911 3   17  +   CCCCCAGATTGGGAG TFFM    62  1.841150116737821e-09
hg19_ct_UserTrack_3545_1911 3   17  -   CTCCCAATCTGGGGG TFFM    62  8.645334758769172e-12
...
...
hg19_ct_UserTrack_3545_3739 86  100 +   GAAGGCAGTCTGAGC TFFM    61  1.5796330106346655e-05
hg19_ct_UserTrack_3545_3739 86  100 -   GCTCAGACTGCCTTC TFFM    61  1.1487113198838786e-10
hg19_ct_UserTrack_3545_3739 87  101 +   AAGGCAGTCTGAGCC TFFM    61  4.790574256236894e-07
hg19_ct_UserTrack_3545_3739 87  101 -   GGCTCAGACTGCCTT TFFM    63  1.5113773626351342e-09
>>> print "detailed best"
detailed best
>>> for hit in tffm_detailed.scan_sequences("test.fa", only_best=True):
...     if hit:
...         print hit
... 
hg19_ct_UserTrack_3545_1911 31  45  +   CTCTCCCTCCGGGCC TFFM    61  0.8256588185762648
hg19_ct_UserTrack_3545_3739 29  43  +   CCTGCCCCCAGGGTG TFFM    62  0.9683269774058026


Download
========

The TFFM-framework can be downloaded from `GitHub <https://github.com/>`_ at 
https://github.com/wassermanlab/TFFM.

Reference
=========

The TFFMs and the TFFM-framework have been described in **Mathelier and
Wasserman, The Next Generation of Transcription Binding Site Prediction**, *PLOS
Computational Biology*, Sept. 2013, 9(9):e1003214,
`DOI:10.1371/journal.pcbi.1003214 <http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1003214>`_.
Please cite this publication when publishing results using the TFFM-framework.

Licence
=======

* The TFFM-framework has been developed under the GNU Lesser General Public
  Licence (see http://www.gnu.org/copyleft/lesser.html and LICENCE).
* The TFFM-framework uses the GHMM library which is also licenced under the GNU
  LGPL.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

