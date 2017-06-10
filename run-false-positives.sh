#!/bin/bash


# Linux readlink -f alternative for Mac OS X
function readlinkUniversal() {
    targetFile=$1

    cd `dirname $targetFile`
    targetFile=`basename $targetFile`

    # iterate down a (possible) chain of symlinks
    while [ -L "$targetFile" ]
    do
        targetFile=`readlink $targetFile`
        cd `dirname $targetFile`
        targetFile=`basename $targetFile`
    done

    # compute the canonicalized name by finding the physical path 
    # for the directory we're in and appending the target file.
    phys_dir=`pwd -P`
    result=$phys_dir/$targetFile
    echo $result
}


os=`uname`
dir=""

case $os in
    Darwin)
        dir=$(dirname "$(readlinkUniversal "$0")")
    ;;
    Linux)
        dir="$(dirname "$(readlink -f "$0")")"
    ;;
    FreeBSD)
        dir=$(dirname "$(readlinkUniversal "$0")")    
    ;;
    *)
       echo "Unknown OS."
       exit 1
    ;;
esac


# Creates data for false positive estimation
function create_data() {
    numberOfClones=$1
    coverageFactor=$2
    length=$3
    sequencer=$4
    prefix="fpEstimation_clones${numberOfClones}_coverage${coverageFactor}_length${length}_seq${sequencer}_"
    # generating in silico clone set
    repseqio generateClones -a -b -c ${numberOfClones} murugan ${prefix}clones.jclns
    # normalizing abundancies
    repseqio normalizeClones ${prefix}clones.jclns ${prefix}clones_normed.jclns
    # export generated sequences to fasta
    repseqio exportCloneSequence -q ${coverageFactor} -g 'VDJTranscript+CExon1' -d NFeature[CDR3] -d AAFeature[CDR3] ${prefix}clones_normed.jclns ${prefix}clones.fasta
    # simulating RNA-Seq reads
    art_illumina --rndSeed 1234 --seqSys ${sequencer} --noALN --paired --len ${length} --mflen 200 --sdev 30 --fcov 1 --in ${prefix}clones.fasta --out ${prefix}clones_
    # aligning with MiXCR
    mixcr align -f -p rna-seq -OallowPartialAlignments=true -g -r ${prefix}align.report ${prefix}clones_1.fq ${prefix}clones_2.fq ${prefix}alignments.vdjca
    # MiXCR partial assembler
    mixcr assemblePartial -r ${prefix}assemblePartial.report ${prefix}alignments.vdjca ${prefix}alignments_rescued.vdjca
    # MiXCR alignments extender
    mixcr extendAlignments -r ${prefix}extend.report ${prefix}alignments.vdjca ${prefix}alignments_extended.vdjca
    mixcr extendAlignments -r ${prefix}extend.report ${prefix}alignments_rescued.vdjca ${prefix}alignments_rescued_extended.vdjca
    # MiXCR exporting alignments
    mixcr exportAlignments -readId -descrR1 ${prefix}alignments.vdjca ${prefix}readToDescr.txt
    mixcr exportAlignments -minFeatureQuality CDR3 -targetDescriptions -nFeature CDR3 ${prefix}alignments_rescued.vdjca ${prefix}overlaps.txt
    mixcr exportAlignments -descrR1 -defaultAnchorPoints -sequence -targetDescriptions -vHitsWithScore -jHitsWithScore -nFeature CDR3 ${prefix}alignments_extended.vdjca ${prefix}extends.txt
    # MiXCR assembling clones
    mixcr assemble -r ${prefix}rescued_extended_assemble.report ${prefix}alignments_rescued_extended.vdjca ${prefix}alignments_rescued_extended.clns
}
export -f create_data

# Simulation for various combinations of input data characteristics
parallel -j4 --line-buffer "create_data {1} {2} {3} {4}" ::: 100 1000 10000 ::: 10000 100000 ::: 50 75 100 ::: HS20 HS25
parallel -j4 --line-buffer "create_data {1} {2} {3} {4}" ::: 100 1000 10000 ::: 10000 100000 ::: 50 75 ::: NS50

# Calculating false overlap rates
ls -1 *_overlaps.txt | sed 's:_overlaps.txt::' | parallel "python $dir/getFalseOverlaps.py {}_readToDescr.txt {}_overlaps.txt {}_rescued_extended_assemble.report" > falseOverlapsResults.txt
# Calculating false extension rates
ls -1 *_overlaps.txt | sed 's:_overlaps.txt::' | parallel "python $dir/getFalseExtensions.py {}_extends.txt {}_rescued_extended_assemble.report" > falseExtensionsResults.txt

# Calculating overall statisics
python $dir/falseExtensionsStat.py falseExtensionsResults.txt
python $dir/falseOverlapsStat.py falseOverlapsResults.txt
