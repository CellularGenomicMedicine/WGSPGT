# script to subsample bam file to specific target coverage as specified by user

# input file is a file with one column with sample name and one column with path of that sample
# target coverages 5, 10, 20, 30 were used

#!/bin/bash
input_file=$1
target_coverage=$2
seed=50
#Loop through each line of the input file
while IFS=$'\t' read -r -a tmp; do
JOBFILE="/path/to/your/job/folder/Subsampling${target_coverage}X_${tmp[0]}.job"
{
echo -e '#!/bin/bash\n'
echo "#SBATCH -J ${tmp[0]}_SubSampling"
echo "#SBATCH --ntasks=1"
echo "#SBATCH --cpus-per-task=3"
echo "#SBATCH --mem-per-cpu=16100M"
echo "#SBATCH --output=/path/to/your/output/folder/Subsampling${target_coverage}X_${tmp[0]}.out"
echo "#SBATCH --error=/path/to/your/error/folder/Subsampling${target_coverage}X_${tmp[0]}.err"
echo "#SBATCH --partition=defq"

echo "#Load Samtools"
echo "module load bioinf/samtools/1.15.1"
echo "target_coverage=$2"
echo "#Absolute value of ${target_coverage}X coverage"
echo "# Echo sample IDs and path to cram file"
echo "echo \"The sampleID is ${tmp[0]}\""
echo "echo \"The cram_file is ${tmp[1]}\""
echo "echo \"Target coverage for sampel ${tmp[0]} is \${target_coverage}\""

echo "#Calculate coverage for the sample"
echo "coverage=\$(samtools depth "${tmp[1]}" | awk '{sum+=\$3} END {print sum/NR}')"
echo "echo \"The coverage for sample ${tmp[0]} is \${coverage} X\""

echo "#Sub sample cram file to ${target_coverage} coverage"
echo "fraction=\$(bc -l <<< \${target_coverage}/\${coverage})"
echo "echo \"Fraction for sample ${tmp[0]} is \${fraction} as fraction\""
echo "samtools view -bs \${fraction} \"${tmp[1]}\" > \"/path/to/your/output/bam/folder/${target_coverage}X/${tmp[0]}.bam\""
} > "$JOBFILE"
echo ${JOBFILE}
sbatch ${JOBFILE}
done < "$input_file"