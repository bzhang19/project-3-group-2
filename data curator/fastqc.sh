#!/bin/bash -l
module load fastqc
fastqc sample/SRR1177981_1.fastq.gz sample/SRR1177981_2.fastq.gz -o fastqc_output/
fastqc sample/SRR1177982_1.fastq.gz sample/SRR1177982_2.fastq.gz -o fastqc_output/
fastqc sample/SRR1177983_1.fastq.gz sample/SRR1177983_2.fastq.gz -o fastqc_output/
fastqc sample/SRR1178008_1.fastq.gz sample/SRR1178008_2.fastq.gz -o fastqc_output/
fastqc sample/SRR1178009_1.fastq.gz sample/SRR1178009_2.fastq.gz -o fastqc_output/
fastqc sample/SRR1178010_1.fastq.gz sample/SRR1178010_2.fastq.gz -o fastqc_output/
fastqc sample/SRR1178014_1.fastq.gz sample/SRR1178014_2.fastq.gz -o fastqc_output/
fastqc sample/SRR1178021_1.fastq.gz sample/SRR1178021_2.fastq.gz -o fastqc_output/
fastqc sample/SRR1178047_1.fastq.gz sample/SRR1178047_2.fastq.gz -o fastqc_output/
