#!/bin/bash
#SBATCH --mail-user=  #############################
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=2
##SBATCH --constraint=haswell|skylake #sandybridge|
#SBATCH --time=72:30:00
#SBATCH --mem-per-cpu=18G
#SBATCH --job-name=PWAS_test
#SBATCH --partition=master
#SBATCH --array=1
#SBATCH --error=/methylation/script/ctp_analysis/log/CBMAP_norm_beta_OSCA_error_%A_%a.log
#SBATCH --output=/methylation/script/ctp_analysis/log/CBMAP_norm_beta_OSCA_output_%A_%a.log

taskplugin=task/affinity

#A:job id
#a:subjob id
module load osca/0.46.1

# ROSMAP
qc=/methylation/results/ctp_analysis/OSCA/format_for_osca
filen=CBMAP_norm_beta_OSCA_input
analysis=/methylation/results/ctp_analysis/OSCA/CBMAP
ref=/methylation/results/ctp_analysis/OSCA

# Output from meffil is in transposed format
osca \
--tefile ${qc}/${filen}.txt \
--methylation-beta \
--no-fid \
--make-bod \
--out ${qc}/${filen}

# Update the annotation file
osca \
--befile ${qc}/${filen} \
--update-opi ${ref}/935k.opi

#---------------------
# Remove sd02 filter
osca \
--befile ${qc}/${filen} \
--extract-probe ${ref}/935K_EPIC_aut.probe \
--exclude-probe ${ref}/935K_EPIC_mask.probe \
--make-bod \
--out ${analysis}/${filen}_aut_mask

# Output matrix
osca \
--befile ${analysis}/${filen}_aut_mask \
--make-efile \
--out ${analysis}/${filen}_aut_mask.txt 

# Output matrix
osca \
--befile ${analysis}/${filen}_aut_mask \
--make-tefile \
--out ${analysis}/${filen}_aut_mask_t.txt 
