#!/bin/bash
#
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err
#SBATCH --mail-type END
#SBATCH --mail-user pierre-guillaume.brun@sb-roscoff.fr
#
#
#SBATCH --partition long
#SBATCH --cpus-per-task 8
#SBATCH --mem 250GB

module load spades/3.15.5


spades.py -o 'w_isolate' -t 8 -m 200 --pe1 1 'fastp_Read_1.fastqsanger.gz' --pe1 2 'fastp_Read_2.fastqsanger.gz' -k 21,33,55,77,99,127 --cov-cutoff off --isolate


spades.py -o 'wo_isolate' -t 8 -m 200 --pe1 1 'fastp_Read_1.fastqsanger.gz' --pe1 2 'fastp_Read_2.fastqsanger.gz' -k 21,33,55,77,99,127 --cov-cutoff off


spades.py -o 'old_reads' -t 8 -m 200 --pe1 1 'fastp_Read_1.fastqsanger.gz' --pe1 2 'kraken_forward_read.fastq' --pe2 1 'kraken_reverse_read.fastq' --pe2 2 'fastp_Read_2.fastqsanger.gz' -k 21,33,55 --cov-cutoff off


spades.py -o 'old_reads_minion' -t 8 -m 200 --pe1 1 'fastp_Read_1.fastqsanger.gz' --pe1 2 'kraken_forward_read.fastq' --pe2 1 'kraken_reverse_read.fastq' --pe2 2 'fastp_Read_2.fastqsanger.gz' -k 21,33,55 --cov-cutoff off --nanopore Minionpassedreads.fastq


spades.py -o 'minion' -t 8 -m 200 --pe1 1 'fastp_Read_1.fastqsanger.gz' --pe1 2 'kraken_forward_read.fastq' -k 21,33,55,77,99,127 --cov-cutoff off --nanopore Minionpassedreads.fastq
