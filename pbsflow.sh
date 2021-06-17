#!/bin/bash

# Parse parameters
usage()
{
        printf "Usage: sicelore pipeline [options]\n\n"
        printf "Options:\n"
        printf "\t--illuminaBam\n"
        printf "\t\tillumina bam file from cellranger count\n"
        printf "\t--illuminaBC\n"
        printf "\t\tillumina barcode.tsv file from cellranger count\n"
        printf "\t--nanofq\n"
        printf "\t\tnanopore fastq file\n"
        printf "\t--juncbed\n"
        printf "\t\tjunc bed file,eg(gencode.v18.mm10.junctions.bed)\n"
        printf "\t--genome\n"
        printf "\t\tgenome fasta file\n"
        printf "\t--refFlat\n"
        printf "\t\trefFlat.txt file,can bet got form paftools.js\n"
        printf "\t--jarpath\n"
        printf "\t\tJar path for sicelore *.jar file\n"
        printf "\t-o O,--outdir OUTDIR\n"
        printf "\t\toutdir to save result\n"
        printf "\t-h, --help\n"
        printf "\t\tShow this help message and exit.\n"
}

# Default parameters
#java="/usr/bin/java"
java=`which java`
spoa=`which spoa`
minimap2=`which minimap2`
samtools=`which samtools`
racon=`which racon`
if [ -z "$java" ] || [ -z "$spoa" ] || [ -z "$samtools" ] || [ -z "$minimap2" ] || [ -z "$racon" ];
then
   echo -e "\nMissing path to required softwares:"
   echo -e "\tjava=$java"
   echo -e "\tspoa=$spoa"
   echo -e "\tsamtools=$samtools"
   echo -e "\tminimap2=$minimap2"
   echo -e "\tracon=$racon"
   echo -e "\nPlease update your \$PATH and rerun.\n\n"
   exit
fi

# Get the parameters selected by the user 
while [ "$1" != "" ]; do
        PARAM=`echo $1 | awk -F' ' '{print $1}'`
        VALUE=`echo $2 | awk -F' ' '{print $1}'`
        case $PARAM in
                -h | --help)
                        usage
                        exit
                        ;;
                -o | --outdir)
                        outdir=$VALUE
                        shift 2
                        ;;
                --illuminaBam)
                        illuminaBam=$VALUE
                        shift 2
                        ;;
                --illuminaBC)
                        illuminaBC=$VALUE
                        shift 2
                        ;;
                --nanofq)
                        nanofq=$VALUE
                        shift 2
                        ;;
                --juncbed)
                        juncbed=$VALUE
                        shift 2
                        ;;
                --genome)
                        genome=$VALUE
                        shift 2
                        ;;
                --refFlat)
                        refFlat=$VALUE
                        shift 2
                        ;;
                --jarpath)
                        jarpath=$VALUE
                        shift 2
                        ;;        
                *)
                        echo "Error: unknown parameter \"$PARAM\""
                        exit 1
                        ;;
        esac
done

tmp_dir="/tmp/sicelore"

# create output directory
if [ ! -e $outdir ];then
        mkdir -p $outdir
fi

if [ ! -e $tmp_dir ];then
        mkdir -p $tmp_dir
fi


############################   main task
# parse illumina bam file
echo "parse illumina bam file"
#$java -jar ${jarpath}/IlluminaParser-1.0.jar -i ${illuminaBam} -o $outdir/illumina.bam.obj -t ${illuminaBC} -b CB -g GN -u UB


# scan nanopore reads
echo "scan nanopore reads"
#$java -jar ${jarpath}/NanoporeReadScanner-0.5.jar -i ${nanofq}  -o ${outdir}
#filename=`basename ${nanofq}`
#fwdname="${filename%.*}FWD.fq"
# map reads to genome
echo "map reads to genome"
#$minimap2 -ax splice -uf --MD --sam-hit-only -t 4 --junc-bed ${juncbed}  ${genome}  ${outdir}/passed/${fwdname} > ${outdir}/minimap.sam
#$samtools view -Sb $outdir/minimap.sam -o $outdir/minimap.unsorted.bam
#$samtools sort $outdir/minimap.unsorted.bam -o $outdir/minimap.bam
#$samtools index $outdir/minimap.bam


# tag reads with gene name
echo "tag reads with gene name"
#$java -jar -Xmx16g ${jarpath}/Sicelore-2.0.jar AddGeneNameTag I=$outdir/minimap.bam O=$outdir/GE.bam REFFLAT=$refFlat GENETAG=GE ALLOW_MULTI_GENE_READS=true USE_STRAND_INFO=true VALIDATION_STRINGENCY=SILENT
#$samtools index $outdir/GE.bam

# tag reads with fastq sequence
echo "tag reads with fastq sequence"
#$java -jar -Xmx64g ${jarpath}/Sicelore-2.0.jar AddBamReadSequenceTag I=$outdir/GE.bam O=$outdir/GEUS.bam FASTQ=$outdir/passed/${fwdname} VALIDATION_STRINGENCY=SILENT
#$samtools index $outdir/GEUS.bam

# tag reads with cellBC/UMI barcodes
echo "tag reads with cellBC/UMI barcodes"
$java -jar -Xmx64g ${jarpath}/NanoporeBC_UMI_finder-1.0.jar -i $outdir/GEUS.bam -o $outdir/GEUS10xAttributes.bam -k $outdir/illumina.bam.obj --ncpu 8 --maxUMIfalseMatchPercent 1 --maxBCfalseMatchPercent 5 --logFile $outdir/out.log
$samtools index $outdir/GEUS10xAttributes.bam
$samtools index $outdir/GEUS10xAttributes_umifound_.bam


# compute consensus sequence
echo "compute consensus sequence"
$java -jar -Xmx64g ${jarpath}/Sicelore-2.0.jar ComputeConsensus T=10 I=$outdir/GEUS10xAttributes_umifound_.bam O=$outdir/consensus.fq TMPDIR=$tmp_dir

#tmpFas=`ls $tmp_dir/*.fa`
#for fa in $tmpFas
#do
#   name=`basename $fa`
#   nname="${name%_*}.msa"
#   $spoa $fa > $tmp_dir/$nname
#done

#$java -jar -Xmx16g ${jarpath}/Sicelore-2.0.jar ComputeConsensus T=10 I=$outdir/GEUS10xAttributes_umifound_.bam O=$outdir/consensus.fq TMPDIR=$tmp_dir
# map molecules to genome
echo "map molecules to genome"
$minimap2 -ax splice -uf --MD --sam-hit-only -t 4 --junc-bed ${juncbed} ${genome} $outdir/consensus.fq > $outdir/molecule.sam
$samtools view -Sb $outdir/molecule.sam -o $outdir/molecule.unsorted.bam
$samtools sort $outdir/molecule.unsorted.bam -o $outdir/molecule.bam
$samtools index $outdir/molecule.bam

# add cellBC/UMI tags
echo "add cellBC/UMI tags"
$java -jar -Xmx64g ${jarpath}/Sicelore-2.0.jar AddBamMoleculeTags I=$outdir/molecule.bam O=$outdir/molecule.tags.bam
$samtools index $outdir/molecule.tags.bam

# add gene name tag
echo "add gene name tag"
$java -jar -Xmx64g ${jarpath}/Sicelore-2.0.jar AddGeneNameTag I=$outdir/molecule.tags.bam O=$outdir/molecule.tags.GE.bam REFFLAT=$refFlat GENETAG=GE ALLOW_MULTI_GENE_READS=true USE_STRAND_INFO=true VALIDATION_STRINGENCY=SILENT
$samtools index $outdir/molecule.tags.GE.bam

# generate molecule isoform matrix
echo "generate molecule isoform matrix"
$java -jar -Xmx64g ${jarpath}/Sicelore-2.0.jar IsoformMatrix DELTA=2 METHOD=STRICT ISOBAM=true GENETAG=GE I=$outdir/molecule.tags.GE.bam REFFLAT=${refFlat} CSV=${illuminaBC} OUTDIR=$outdir PREFIX=sicmol VALIDATION_STRINGENCY=SILENT
$samtools index $outdir/sicmol_isobam.bam



