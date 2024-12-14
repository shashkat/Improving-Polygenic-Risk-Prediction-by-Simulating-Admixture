#!/bin/bash
for i in {1..22}
do
    echo $i
    java -jar flare.jar ref=vcf_files/chromosome$i/ceu_yri_geno.vcf.gz gt=vcf_files/chromosome$i/asw_geno.vcf.gz map=map_files/plink.chr$i.GRCh37.map ref-panel=map_files/ceu_yri_sample.txt out=flare_output/asw/chromosome$i/flare.out probs=true
    java -jar flare.jar ref=vcf_files/chromosome$i/ceu_yri_geno.vcf.gz gt=vcf_files/chromosome$i/sim_geno.vcf.gz map=map_files/plink.chr$i.GRCh37.map ref-panel=map_files/ceu_yri_sample.txt out=flare_output/sim/chromosome$i/flare.out probs=true

done
