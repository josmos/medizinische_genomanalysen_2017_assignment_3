# Medizinische Genomanalysen 2017 - Assignment 3

## dependencies:
vcf-tools:
https://vcftools.github.io/perl_module.html

    sudo python install_vcftools.py

## Overview
* Fork and clone the repository
* Complete the python program, based on the template, to calculate various properties
* Push to your repository

## Attention!
As the HGVS package does only support Python2, please run the script with Python2

## Data
* Son
  * HG002-NA24385-huAA53E0
  * ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/IonTorrent_TVC_03162015/AmpliseqExome.20141120.NA24385.vcf
* Mother
  * HG004-NA24143-hu8E87A9
  * ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/IonTorrent_TVC_03162015/AmpliseqExome.20141120.NA24143.vcf
* Father
  * HG003-NA24149-hu6E4515
  * ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/IonTorrent_TVC_03162015/AmpliseqExome.20141120.NA24149.vcf
  
## Tools
* hgvs (https://hgvs.readthedocs.io/en/master/) - only works with Python 2
* cyvcf2 (https://github.com/brentp/cyvcf2)