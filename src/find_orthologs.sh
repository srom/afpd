#!/bin/bash

set -e

conda activate afpd

cd ~/Documents/afpd

mmseqs easy-rbh \
    T_kodakarensis/T_kodakarensis.fasta \
    P_furiosus/P_furiosus.fasta \
    T_kodakarensis_P_furiosus_orthologs.tsv \
    /tmp \
    --format-mode 4
