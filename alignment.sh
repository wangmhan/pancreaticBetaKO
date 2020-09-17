#!/usr/bin/env bash

#modified from Drop-seq_alignment.sh

tmpdir=
outdir=`pwd`
genomedir=
reference=
pipeline=0
echo_prefix=
dropseq_root=$(dirname $0)
star_executable=STAR
estimated_num_cells=
progname=`basename $0`

function usage () {
    cat >&2 <<EOF
USAGE: $progname [options] <unmapped-queryname-sorted.bam>
Perform Drop-seq tagging, trimming and alignment

-g <genomedir>      : Directory of STAR genome directory.  Required.
-r <referencefasta> : Reference fasta of the Drop-seq reference metadata bundle.  Required.
-n <estimated-num-cells> : estimate of number of cells in experiment.  Required.
-d <dropseq_root>   : Directory containing Drop-seq executables.  Default: directory containing this script.
-o <outputdir>      : Where to write output bam.  Default: current directory.
-t <tmpdir>         : Where to write temporary files.  Default: a new subdirectory in $TMPDIR.
-s <STAR_path>      : Full path of STAR.  Default: STAR is found via PATH environment variable.
-p                  : Reduce file I/O by pipeline commands together.  Requires more memory and processing power.
-e                  : Echo commands instead of executing them.  Cannot use with -p.
EOF
}

function error_exit() {
    echo "ERROR: $1
    " >&2
    usage
    exit 1
}

function check_set() {
    value=$1
    name=$2
    flag=$3

    if [[ -z "$value" ]]
    then error_exit "$name has not been specified.  $flag flag is required"
    fi
}

set -e
# Fail if any of the commands in a pipeline fails
set -o pipefail

while getopts ":d:t:o:pg:r:es:n:" options; do
  case $options in
    d ) dropseq_root=$OPTARG;;
    t ) tmpdir=$OPTARG;;
    o ) outdir=$OPTARG;;
    n ) estimated_num_cells=$OPTARG;;
    p ) pipeline=1;;
    g ) genomedir=$OPTARG;;
    r ) reference=$OPTARG;;
    s ) star_executable=$OPTARG;;
    e ) echo_prefix="echo";;
    h ) usage
          exit 1;;
    \? ) usage
         exit 1;;
    * ) usage
          exit 1;;
  esac
done
shift $(($OPTIND - 1))

if [[ "$pipeline" == 1 && -n "$echo_prefix" ]]
then error_exit "-p and -e cannot be used together"
fi

check_set "$dropseq_root" "Drop-seq root" "-d"
check_set "$genomedir" "Genome directory" "-g"
check_set "$reference" "Reference fasta"  "-r"
check_set "$estimated_num_cells" "Estimated num cells" "-n"

if (( $# != 1 ))
then error_exit "Incorrect number of arguments"
fi

if [[ -z "$tmpdir" ]]
then tmpdir=`mktemp -d`
     echo "Using temporary directory $tmpdir"
fi

if [[ "$star_executable" != "STAR" ]]
then if [[ ! ( -x $star_executable && -f $star_executable ) ]]
     then error_exit "STAR executable $star_executable passed via -s does not exist or is not executable"
     fi
elif which STAR > /dev/null
then echo > /dev/null
else error_exit "STAR executable must be on the path"
fi

gene_intervals=$(dirname $reference)/$(basename $reference .fasta).genes.intervals
exon_intervals=$(dirname $reference)/$(basename $reference .fasta).exon.intervals
refflat=$(dirname $reference)/$(basename $reference .fasta).refFlat
picard_jar=${dropseq_root}/3rdParty/picard/picard.jar

unmapped_bam=$1
tagged_unmapped_bam=${tmpdir}/unaligned_mc_tagged_polyA_filtered.bam
aligned_sam=${tmpdir}/star.Aligned.out.sam
aligned_sorted_bam=${tmpdir}/aligned.sorted.bam
files_to_delete="${aligned_sorted_bam} ${tagged_unmapped_bam}"
#${aligned_sam}

# Stage 1: pre-alignment tag and trim
tag_cells="${dropseq_root}/TagBamWithReadSequenceExtended SUMMARY=${outdir}/unaligned_tagged_Cellular.bam_summary.txt \
    BASE_RANGE=1-16 BASE_QUALITY=10 BARCODED_READ=2 DISCARD_READ=false TAG_NAME=XC NUM_BASES_BELOW_QUALITY=1 INPUT=${unmapped_bam}"

tag_molecules="${dropseq_root}/TagBamWithReadSequenceExtended SUMMARY=${outdir}/unaligned_tagged_Molecular.bam_summary.txt \
BASE_RANGE=17-25 BASE_QUALITY=10 BARCODED_READ=2 DISCARD_READ=true TAG_NAME=XM NUM_BASES_BELOW_QUALITY=1"

filter_bam="${dropseq_root}/FilterBAM TAG_REJECT=XQ"

trim_starting_sequence="${dropseq_root}/TrimStartingSequence OUTPUT_SUMMARY=${outdir}/adapter_trimming_report.txt \
SEQUENCE=AAGCAGTGGTATCAACGCAGAGTGAATGGG MISMATCHES=0 NUM_BASES=5"

trim_poly_a="${dropseq_root}/PolyATrimmer OUTPUT=${tagged_unmapped_bam} OUTPUT_SUMMARY=${outdir}/polyA_trimming_report.txt MISMATCHES=0 NUM_BASES=6"

# Stage 2: alignment
sam_to_fastq="java -Xmx500m -jar ${picard_jar} SamToFastq INPUT=${tmpdir}/unaligned_mc_tagged_polyA_filtered.bam"
star_align="$star_executable --genomeDir ${genomedir} --runThreadN 8 --outSAMstrandField intronMotif --sjdbGTFfile /mnt/hel/home/menwan/references/human_transgenic/hgGRCh38_transgenic.gtf --outFileNamePrefix ${tmpdir}/star."

# Stage 3: sort aligned reads (STAR does not necessarily emit reads in the same order as the input)
sort_aligned="java -Dsamjdk.buffer_size=131072 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx4000m -jar ${picard_jar} \
SortSam INPUT=${aligned_sam} OUTPUT=${aligned_sorted_bam} SORT_ORDER=queryname TMP_DIR=${tmpdir}"

# Stage 4: merge and tag aligned reads
merge_bam="java -Xmx4000m -jar ${picard_jar} MergeBamAlignment REFERENCE_SEQUENCE=${reference} UNMAPPED_BAM=${tagged_unmapped_bam} \
ALIGNED_BAM=${aligned_sorted_bam} INCLUDE_SECONDARY_ALIGNMENTS=false PAIRED_RUN=false"
tag_with_gene_exon="${dropseq_root}/TagReadWithGeneExon O=${tmpdir}/star_gene_exon_tagged.bam ANNOTATIONS_FILE=${refflat} TAG=GE"

if (( $pipeline == 1 ))
then
     # Stage 1
     $tag_cells OUTPUT=/dev/stdout COMPRESSION_LEVEL=0 | \
       $tag_molecules INPUT=/dev/stdin OUTPUT=/dev/stdout COMPRESSION_LEVEL=0 | \
       $filter_bam INPUT=/dev/stdin OUTPUT=/dev/stdout COMPRESSION_LEVEL=0 | \
       $trim_starting_sequence INPUT=/dev/stdin OUTPUT=/dev/stdout COMPRESSION_LEVEL=0 | \
       $trim_poly_a INPUT=/dev/stdin

     # Stage 2
     $sam_to_fastq FASTQ=/dev/stdout | \
       $star_align --readFilesIn /dev/stdin

     # Stage 3
     $sort_aligned

     # Stage 4
     $merge_bam OUTPUT=/dev/stdout COMPRESSION_LEVEL=0 | \
       $tag_with_gene_exon I=/dev/stdin
else
     # Stage 1
     $echo_prefix $tag_cells OUTPUT=$tmpdir/unaligned_tagged_Cell.bam
     $echo_prefix $tag_molecules INPUT=$tmpdir/unaligned_tagged_Cell.bam OUTPUT=$tmpdir/unaligned_tagged_CellMolecular.bam
     $echo_prefix $filter_bam INPUT=$tmpdir/unaligned_tagged_CellMolecular.bam OUTPUT=$tmpdir/unaligned_tagged_filtered.bam
     $echo_prefix $trim_starting_sequence INPUT=$tmpdir/unaligned_tagged_filtered.bam OUTPUT=$tmpdir/unaligned_tagged_trimmed_smart.bam
     $echo_prefix $trim_poly_a INPUT=$tmpdir/unaligned_tagged_trimmed_smart.bam
     files_to_delete="$files_to_delete $tmpdir/unaligned_tagged_Cell.bam $tmpdir/unaligned_tagged_CellMolecular.bam \
                        $tmpdir/unaligned_tagged_filtered.bam $tmpdir/unaligned_tagged_trimmed_smart.bam"


     # Stage 2
     $echo_prefix $sam_to_fastq FASTQ=$tmpdir/unaligned_mc_tagged_polyA_filtered.fastq
     $echo_prefix $star_align --readFilesIn $tmpdir/unaligned_mc_tagged_polyA_filtered.fastq
     files_to_delete="$files_to_delete $tmpdir/unaligned_mc_tagged_polyA_filtered.fastq"

     # Stage 3
     $echo_prefix $sort_aligned

     # Stage 4
     $echo_prefix $merge_bam OUTPUT=$tmpdir/merged.bam
     $echo_prefix $tag_with_gene_exon INPUT=$tmpdir/merged.bam
     files_to_delete="$files_to_delete $tmpdir/merged.bam"

fi

let num_barcodes=4*$estimated_num_cells

$echo_prefix ${dropseq_root}/DetectBeadSynthesisErrors INPUT=${tmpdir}/star_gene_exon_tagged.bam OUTPUT=${outdir}/error_detected.bam \
        CREATE_INDEX=true TMP_DIR=$tmpdir \
        SUMMARY=${outdir}/bead_synthesis_error_summary OUTPUT_STATS=${outdir}/bead_synthesis_error_detail \
        NUM_BARCODES=$num_barcodes

files_to_delete="$files_to_delete ${tmpdir}/star_gene_exon_tagged.bam"

$echo_prefix rm $files_to_delete
