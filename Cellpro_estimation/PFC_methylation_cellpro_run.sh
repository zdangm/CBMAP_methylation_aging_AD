#!/bin/bash
#SBATCH --mail-user=  #############################
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=2
##SBATCH --constraint=haswell|skylake #sandybridge|
#SBATCH --time=72:30:00
#SBATCH --mem-per-cpu=18G
#SBATCH --job-name=cellpro
#SBATCH --partition=master
#SBATCH --array=1
#SBATCH --error=/ctp_analysis/log/CBMAP_norm_beta_OSCA_error_%A_%a.log
#SBATCH --output=/ctp_analysis/log/CBMAP_norm_beta_OSCA_output_%A_%a.log

taskplugin=task/affinity  

#A:job id
#a:subjob id
module load osca/0.46.1


qc=/ctp_analysis/OSCA/format_for_osca
filen=CBMAP_norm_beta_OSCA_input
analysis=/ctp_analysis/OSCA/CBMAP
ref=/ctp_analysis/OSCA

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



# Apply filters
osca \
--befile ${qc}/${filen} \
--sd-min 0.02 \
--extract-probe ${ref}/935k_EPIC_aut.probe \
--exclude-probe ${ref}/935k_EPIC_mask.probe \
--make-efile \
--out ${analysis}/${filen}_sd02_aut_mask.txt

#---------------------
# Remove sd02 filter
osca \
--befile ${qc}/${filen} \
--extract-probe ${ref}/450K_aut.probe \
--exclude-probe ${ref}/450K_mask.probe \
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