GAUDI Pipeline


1. Simulate recombination to generate simulated genotypes, and save ASW, CEU, YRI, and simulated genos to npy files. (sim_recomb_and_save_geno.ipynb)
2. Populate VCF files using SNP and GENO information. (geno_to_vcf.ipynb)
3. Predict local ancestry probabilities using generated VCF files as input. (flare.sh)
4. Parse FLARE output VCF files into npy files. (vcf_to_npy.ipynb)
5. Perform GAUDI and normal polygenic risk prediction using original geno files, and local ancestry from FLARE output. (gaudi.ipynb or gaudi_sim.ipynb)
