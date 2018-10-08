Project of Escherichia coli
===========================

Overview
--------

`Escherichia coli <https://en.wikipedia.org/wiki/Escherichia_coli>`__ is
a Gram-negative, rod-shaped bacterium belonging to the family
Enterobacteriaceae that was described in 1885 by a German pediatrician.
Pathogenic E.coli is versatile due to the diversity of their gene sets.
Virulence factors usually located on a virulence plasmid and can be
acquired through gene transfer. In this study, we investigate the
co-occurrance of virulence factors among all the available genome in the
`Genbank <https://www.ncbi.nlm.nih.gov/genbank/>`__

Environment
-----------

-  `Python <https://www.python.org/download/releases/2.7/>`__
-  `Anaconda <https://www.anaconda.com/>`__ If you havn’t install
   Python, I strongly recommend you to use Anaconda, which is a free and
   open source distribution of the Python.
-  Based on Macos(I will try to make it accomplishable on Windows
   System)

Genome Download
---------------

Install Wget
~~~~~~~~~~~~

On MacOs

::

   ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
   brew install wget

On Windows Download
`wget <http://gnuwin32.sourceforge.net/packages/wget.htm>`__ for Windows

Download via accession number
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

According the
`assembly_summary_genbank <ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt>`__
in NCBI ftp site, we can get the genome through the organism name

::

   import os
   os.mkdir("Ecoli_genome")
   os.chdir("Ecoli_genome")
   with open(r'/Your/File/PAth/assembly_summary_genbank.txt') as f:
       for i in f.read().split('\n')[2:-1]:
           os.system('wget '+i+'_genomic.fna.gz')

or Use `Pandas <https://pandas.pydata.org/>`__

::

   import pandas as pd
   f=pd.read_table("/Users/huhuifeng/Desktop/assembly_summary_genbank.txt",sep='\t',header=1)
   for i in f[f['organism_name']=='Escherichia coli']['ftp_path']:
       os.system('wget '+i+'_genomic.fna.gz')

Translate into amino acid sequence
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Usually I only download fna file, and use
`Prodigal <https://github.com/hyattpd/Prodigal>`__ to translate
nucleotide into amino acid sequence.

::

   for i in os.listdir('.'):
       order = 'prodigal -i '+i+' -q -a '+i'.faa'+' -d '+i+'.nucl'+' -o '+i+'.out'
       print order 
       os.system(order)

Now we can get predicted protein sequence of all E.coli genomes.
