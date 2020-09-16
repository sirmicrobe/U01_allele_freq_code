#Ciprofloxacin variant calling from Illumina sequencing reads

#Trim and filter reads
for i in {42..46}; do trimmomatic PE -phred33 /home/dmux/170215/CooperLabMRS/021417_$i/*R1_001.fastq.gz /home/dmux/170215/CooperLabMRS/021417_$i/*R2_001.fastq.gz "$i"_forward_paired.fq.gz "$i"_forward_unpaired.fq.gz "$i"_reverse_paired.fq.gz "$i"_reverse_unpaired.fq.gz ILLUMINACLIP:/home/vclocal/.linuxbrew/Cellar/trimmomatic/0.35/share/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:70;done

#Breseq on populations for variant calling
for i in {01..42}; do time breseq -p -r /home/cwm47/ref_genomes/Abaumannii17978/NZ_CP012004_GCF_001077675.1_ASM107767v1_NC009084pAB1_009083pAB2_genomic.gbff ~/abaum/uo1/alfonso_cipro_evol/trim/"$i"_forward_paired.fq.gz ~/abaum/uo1/alfonso_cipro_evol/trim/"$i"_forward_unpaired.fq.gz ~/abaum/uo1/alfonso_cipro_evol/trim/"$i"_reverse_paired.fq.gz -o ~/abaum/uo1/alfonso_cipro_evol/populations/new_ref/breseq_"$i"_v_vc17978 -j 8; done 
for i in {43..84}; do time breseq -p -r /home/cwm47/ref_genomes/Abaumannii17978/NZ_CP012004_GCF_001077675.1_ASM107767v1_NC009084pAB1_009083pAB2_genomic.gbff ~/abaum/uo1/alfonso_cipro_evol/trim/"$i"_forward_paired.fq.gz ~/abaum/uo1/alfonso_cipro_evol/trim/"$i"_forward_unpaired.fq.gz ~/abaum/uo1/alfonso_cipro_evol/trim/"$i"_reverse_paired.fq.gz -o ~/abaum/uo1/alfonso_cipro_evol/populations/new_ref/breseq_"$i"_v_vc17978 -j 8; done 

#Breseq on clones for variant calling
for i in {33..71}; do time breseq -r /home/cwm47/ref_genomes/Abaumannii17978/NZ_CP012004_GCF_001077675.1_ASM107767v1_NC009084pAB1_009083pAB2_genomic.gbff  ~/abaum/uo1/alfonso_cipro_evol/trim/clones/"$i"_forward_paired.fq.gz ~/abaum/uo1/alfonso_cipro_evol/trim/clones/"$i"_forward_unpaired.fq.gz ~/abaum/uo1/alfonso_cipro_evol/trim/clones/"$i"_reverse_paired.fq.gz -o ~/abaum/uo1/alfonso_cipro_evol/clones/new_ref/breseq_"$i"_v_vc17978 -j 8; done 

#Used GD Tools scripts in breseq for mutation type counting
#####gdtools version of variants
#subtract reference out - CIPRO
for i in {01..84}; do gdtools SUBTRACT -o /home/cwm47/abaum/uo1/alfonso_cipro_evol/populations/new_ref/gd_files/"$i"_subtractref_output.gd /home/cwm47/abaum/uo1/alfonso_cipro_evol/populations/new_ref/breseq_"$i"_v_vc17978/output/output.gd ~/abaum/uo1/alfonso_cipro_evol/breseq_atcc/NZ_CP012004_pAB123_ref/1_10p_v_vc17978/output/output.gd; done
#annotate each gdtools file
for i in {01..84}; do gdtools ANNOTATE -o "$i"_subtractref.txt -f TSV -r /home/cwm47/ref_genomes/Abaumannii17978/NZ_CP012004_GCF_001077675.1_ASM107767v1_NC009084pAB1_009083pAB2_genomic.gbff /home/cwm47/abaum/uo1/alfonso_cipro_evol/populations/new_ref/gd_files/"$i"_subtractref_output.gd; done

#change outputs to sample name - ran locally
for i in {01..84}; do sed -i '' "s/output/sample_$i/g" "$i"_subtractref.txt; done

#END - moved analysis into R
