#Process sequence data from wlbu2 mutator experiment


#run trimmomatic
module load trimmomatic/trimmomatic-0.36
for i in {41..60}; do trimmomatic PE -threads 16 -trimlog trim.log -phred33 /home/data/dmux/180721/180721/CooperLabASL/072118_"$i"/*_R1_001.fastq.gz /home/data/dmux/180721/180721/CooperLabASL/072118_"$i"/*_R2_001.fastq.gz /home/cwm47/abaum/uo1/jeffrey_wlbu2/trim/"$i"_forward_paired.fq.gz /home/cwm47/abaum/uo1/jeffrey_wlbu2/trim/"$i"_forward_unpaired.fq.gz /home/cwm47/abaum/uo1/jeffrey_wlbu2/trim/"$i"_reverse_paired.fq.gz /home/cwm47/abaum/uo1/jeffrey_wlbu2/trim/"$i"_reverse_unpaired.fq.gz ILLUMINACLIP:/home/cwm47/build/Trimmomatic-0.36/adapters/NexteraPE-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:70;done

trimmomatic PE -threads 8 -trimlog trim.log -phred33 /home/data/dmux/180721/180721/CooperLabASL/072118_96/*_R1_001.fastq.gz /home/data/dmux/180721/180721/CooperLabASL/072118_96/*_R2_001.fastq.gz /home/cwm47/abaum/uo1/jeffrey_wlbu2/trim/96_forward_paired.fq.gz /home/cwm47/abaum/uo1/jeffrey_wlbu2/trim/96_forward_unpaired.fq.gz /home/cwm47/abaum/uo1/jeffrey_wlbu2/trim/96_reverse_paired.fq.gz /home/cwm47/abaum/uo1/jeffrey_wlbu2/trim/96_reverse_unpaired.fq.gz ILLUMINACLIP:/home/cwm47/build/Trimmomatic-0.36/adapters/NexteraPE-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:70

#run breseq
module load breseq/breseq-0.35.0
for i in {41..60}; do time breseq -p -r /home/cwm47/ref_genomes/Paeruginosa/PA14_8-25-20_GCF_000014625.1_ASM1462v1_genomic.gbff  /home/cwm47/abaum/uo1/jeffrey_wlbu2/trim/"$i"_forward_paired.fq.gz /home/cwm47/abaum/uo1/jeffrey_wlbu2/trim/"$i"_forward_unpaired.fq.gz /home/cwm47/abaum/uo1/jeffrey_wlbu2/trim/"$i"_reverse_paired.fq.gz -o /home/cwm47/abaum/uo1/jeffrey_wlbu2/populations_rerun/breseq_"$i" -j 16; done 

breseq -p -r /home/cwm47/ref_genomes/Paeruginosa/PA14_8-25-20_GCF_000014625.1_ASM1462v1_genomic.gbff /home/mjf123/PA14/PA14pops/PA14_b13/PA14_b13_S50_R1_001.fastq.gz /home/mjf123/PA14/PA14pops/PA14_b13/PA14_b13_S50_R2_001.fastq.gz -o /home/cwm47/abaum/uo1/jeffrey_wlbu2/populations_rerun/PA14_b13 -j 16 

breseq -p -r /home/cwm47/ref_genomes/Paeruginosa/PA14_8-25-20_GCF_000014625.1_ASM1462v1_genomic.gbff /home/mjf123/PA14/PA14pops/PA14W_a13/PAW_a13_S51_R1_001.fastq.gz /home/mjf123/PA14/PA14pops/PA14W_a13/PAW_a13_S51_R2_001.fastq.gz -o /home/cwm47/abaum/uo1/jeffrey_wlbu2/populations_rerun/PA14W_a13 -j 16 

breseq -p -r /home/cwm47/ref_genomes/Paeruginosa/PA14_8-25-20_GCF_000014625.1_ASM1462v1_genomic.gbff /home/mjf123/PA14/PA14pops/PA14W_e13/PAWe13_S52_R1_001.fastq.gz /home/mjf123/PA14/PA14pops/PA14W_e13/PAWe13_S52_R2_001.fastq.gz -o /home/cwm47/abaum/uo1/jeffrey_wlbu2/populations_rerun/PA14W_e13 -j 16 

breseq -r /home/cwm47/ref_genomes/Paeruginosa/PA14_8-25-20_GCF_000014625.1_ASM1462v1_genomic.gbff /home/mjf123/PA14/PA14clones/PAWa13s/PAWa13s_S148_R1_001.fastq.gz /home/mjf123/PA14/PA14clones/PAWa13s/PAWa13s_S148_R2_001.fastq.gz -o /home/cwm47/abaum/uo1/jeffrey_wlbu2/populations_rerun/PAWa13s_clone -j 16 

breseq -r /home/cwm47/ref_genomes/Paeruginosa/PA14_8-25-20_GCF_000014625.1_ASM1462v1_genomic.gbff /home/mjf123/PA14/PA14clones/PAWa13w/PAWa13W_S149_R1_001.fastq.gz /home/mjf123/PA14/PA14clones/PAWa13w/PAWa13W_S149_R2_001.fastq.gz -o /home/cwm47/abaum/uo1/jeffrey_wlbu2/populations_rerun/PAWa13w_clone -j 16 

breseq -r /home/cwm47/ref_genomes/Paeruginosa/PA14_8-25-20_GCF_000014625.1_ASM1462v1_genomic.gbff /home/mjf123/PA14/PA14clones/PAWe13n/PAWe13n_S146_R1_001.fastq.gz /home/mjf123/PA14/PA14clones/PAWe13n/PAWe13n_S146_R2_001.fastq.gz -o /home/cwm47/abaum/uo1/jeffrey_wlbu2/populations_rerun/PAWe13n_clone -j 16 

breseq -r /home/cwm47/ref_genomes/Paeruginosa/PA14_8-25-20_GCF_000014625.1_ASM1462v1_genomic.gbff /home/mjf123/PA14/PA14clones/PAWe13t/PAWe13t_S147_R1_001.fastq.gz /home/mjf123/PA14/PA14clones/PAWe13t/PAWe13t_S147_R2_001.fastq.gz -o /home/cwm47/abaum/uo1/jeffrey_wlbu2/populations_rerun/PAWe13t_clone -j 16 

breseq -r /home/cwm47/ref_genomes/Paeruginosa/PA14_8-25-20_GCF_000014625.1_ASM1462v1_genomic.gbff /home/cwm47/abaum/uo1/jeffrey_wlbu2/trim/55_forward_paired.fq.gz /home/cwm47/abaum/uo1/jeffrey_wlbu2/trim/55_forward_unpaired.fq.gz /home/cwm47/abaum/uo1/jeffrey_wlbu2/trim/55_reverse_paired.fq.gz -o /home/cwm47/abaum/uo1/jeffrey_wlbu2/populations_rerun/breseq_55_clone -j 16 

module load miniconda/miniconda-3
python3 ~/build/gscripts/breseq_parser.py -d /home/cwm47/abaum/uo1/jeffrey_wlbu2/populations_rerun -f xlsx -o wlbu2_breseq_rerun_pars_out

python /home/cwm47/build/gscripts/BreseqCat.py -p -d /home/cwm47/abaum/uo1/jeffrey_wlbu2/populations_rerun

#check outputs
for i in {41..60}; do echo breseq_$i && ls breseq_$i/output/; done

#Subtract the ancestral mutations (and false positives) from all evolved populations
for i in {41..60}; do gdtools SUBTRACT -o ./gdtools/"$i"_subtractref_output.gd breseq_"$i"/output/output.gd /home/cwm47/abaum/uo1/jeffrey_wlbu2/populations_rerun/breseq_55/output/output.gd; done
#annotate each gdtools file
for i in {41..60}; do gdtools ANNOTATE -o ./gdtools/"$i"_subtractref_output.txt -f TSV -r /home/cwm47/ref_genomes/Paeruginosa/PA14_8-25-20_GCF_000014625.1_ASM1462v1_genomic.gbff ./gdtools/"$i"_subtractref_output.gd; done
#count mutation statistics
gdtools COUNT -b -p -r /home/cwm47/ref_genomes/Paeruginosa/PA14_8-25-20_GCF_000014625.1_ASM1462v1_genomic.gbff -o counts_stats_all 41_subtractref_output.gd 42_subtractref_output.gd 43_subtractref_output.gd 44_subtractref_output.gd 45_subtractref_output.gd 46_subtractref_output.gd 47_subtractref_output.gd 48_subtractref_output.gd 49_subtractref_output.gd 50_subtractref_output.gd 51_subtractref_output.gd 52_subtractref_output.gd 53_subtractref_output.gd 54_subtractref_output.gd 55_subtractref_output.gd 56_subtractref_output.gd 57_subtractref_output.gd 58_subtractref_output.gd 59_subtractref_output.gd 60_subtractref_output.gd #pull files down locally
#change outputs to sample name - run locally
for i in {41..60}; do sed -i '' "s/output/sample_$i/g" "$i"_subtractref_output.txt; done
