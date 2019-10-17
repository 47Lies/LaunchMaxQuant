AlreadyDownloaded=Channel.fromPath("/bioinfo/users/ltaing/DATA_TMP/ltaing/IsoAndSpe/Ressources/uniprot_trembl.fasta")

process GetTrEmbl{
  storeDir "Ressources"
  label "Wget"
  input:
  file "uniprot_trembl.fasta" from AlreadyDownloaded
  output:
  file "uniprot_trembl.9606.fasta" into HumanTrembl
  script:
  """
  cat -n uniprot_trembl.fasta | grep \">\" > uniprot_trembl.LinesEntryNames.txt
  FirstLineInWhole=\$(grep \"OX=9606 \" uniprot_trembl.LinesEntryNames.txt | head -1 | cut -f1 |  sed -e \"s/ //g\")
  LastLineInIndex=\$(cat -n uniprot_trembl.LinesEntryNames.txt | grep \"OX=9606 \" | tail -1 | cut -f1 |  sed -e \"s/ //g\")
  let \"FirstAfterInIndex=\$LastLineInIndex+1\"
  FirstAfterInWhole=\$(head -\$FirstAfterInIndex uniprot_trembl.LinesEntryNames.txt | tail -1 | cut -f1 |  sed -e \"s/ //g\")
  let \"LastInWhole=\$FirstAfterInWhole-1\"
  let \"NLines=\$LastInWhole-\$FirstLineInWhole+1\"
  head -\$LastInWhole uniprot_trembl.fasta | tail -\$NLines > uniprot_trembl.9606.fasta
   """
}



process GetIsoform{
  storeDir "Ressources"
  output:
  file "uniprot_sprot_varsplic.fasta" into UniprotIsoform
  
  """
 wget \"ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot_varsplic.fasta.gz\"\
 --output-document uniprot_sprot_varsplic.fasta.gz 
 gunzip uniprot_sprot_varsplic.fasta.gz
   """
}

process GetSwissProt{
  storeDir "Ressources"
  output:
  file "uniprot_sprot.fasta" into UniprotSwissProt
  """
  wget \"ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz\"\
 --output-document uniprot_sprot.fasta.gz 
 gunzip uniprot_sprot.fasta.gz
   """
}

process CreateWholeUniprot{
  beforeScript = "source /data/users/ltaing/.bashrc" 
  storeDir "Ressources"
  conda "envs/Rsystempiper.txt"
  input:
  file SP from UniprotSwissProt
  file Iso from UniprotIsoform
  file Tr from HumanTrembl
  output:
  file "UniprotComplete.9606.fasta" into Uniprot
  file "${HumanSwissProt}" into HumanSwissProt
  file "${HumanIsoProt}" into HumanIsoProt
  file "UniprotComplete.log" into UniprotCreationLog
  script:
  HumanSwissProt=SP.toString().replaceFirst(/.fasta/,".9606.fasta")
  HumanIsoProt=Iso.toString().replaceFirst(/.fasta/,".9606.fasta")

  """
    Rscript ${PWD}/Scripts/Uniprot2SpeciesUniprot.R -S ${SP} -I ${Iso} -T ${Tr} -O UniprotComplete -x 9606
  """
}

Uniprot.into{Uniprot;UniprotBLAST}

process MakeBlastDatabase{
  beforeScript = "source /data/users/ltaing/.bashrc" 
  storeDir "Ressources"
  conda "envs/BLAST.yml"
  input:
  file "UniprotComplete.9606.fasta" from UniprotBLAST
  output:
  file "UniprotComplete.9606.fasta.phr" into phr
  file "UniprotComplete.9606.fasta.pin" into pin
  file "UniprotComplete.9606.fasta.psq" into psq
  script:
  """
  makeblastdb -in UniprotComplete.9606.fasta -dbtype prot -taxid 9606
  """
}
Uniprot.into{Uniprot;UniprotRESCUE}



/**
 * name: GetRessources
 * 
 * WhatItDoes:
 * 1 download from the ebi ftp the zip gtf
 * 2 download from the ebi ftp the zip DNA genome fasta
 * 3 download from Uniprot ftp the zip Swissprot all species fasta
 * 4 download from Uniprot ftp the zip Trembl all species fasta
 * 5 download from Uniprot ftp the zip Swissprot varsplicing all species fasta
 * 4 unzip them
 * 5 create Ressources and moove the gtf into it
 * 
 * Notes: The UNIPROT files contains entries for all species we need to extract Homo Sapiens entries only
 * For the swissport and the swissport variant splicing files, the bank is order by their name entries by Alpha Numerical order
 * For the Trembl ,the bank is order by species or at least by chunck of species
 * to sort them
 * in BASH
 * cat -n uniprot_sprot_varsplic.fasta | grep ">" > uniprot_sprot_varsplic.LinesEntryNames.txt
 * FirstLineInWhole=$(grep "OX=9606 " uniprot_sprot_varsplic.LinesEntryNames.txt | head -1 | cut -f1 |  sed -e "s/ //g")
 * LastLineInIndex=$(cat -n uniprot_sprot_varsplic.LinesEntryNames.txt | grep "OX=9606 " | tail -1 | cut -f1 |  sed -e "s/ //g")
 * let "FirstAfterInIndex=$LastLineInIndex+1"
 * FirstAfterInWhole=$(head -$FirstAfterInIndex uniprot_sprot_varsplic.LinesEntryNames.txt | tail -1 | cut -f1 |  sed -e "s/ //g")
 * let "LastInWhole=$FirstAfterInWhole-1"
 * let "NLines=$LastInWhole-$FirstLineInWhole+1"
 * head -$LastInWhole uniprot_sprot_varsplic.fasta | tail -$NLines > uniprot_sprot_varsplic.InBetween.fasta
 *
 * 
 * @param storeDir
 *  the directory where the GTF file will be stored
 *  in our case "Ressources"
 * @return output GTF a file to the gencode gtf file that will be in the GTF channel
 * 
 * Risks:
 * I/O
 * error 404
 * 
 */
process GetRessources {
  storeDir "Ressources"
  label "Convert"
  output:
  file "gencode.v28.annotation.gtf" into GTF
  file "GRCh38.p12.genome.fa" into Genome  
  
  """
 wget "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gtf.gz" \
 --output-document gencode.v28.annotation.gtf.gz
 gunzip gencode.v28.annotation.gtf.gz
  wget \"ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/GRCh38.p12.genome.fa.gz\"\
 --output-document GRCh38.p12.genome.fa.gz 
 gunzip GRCh38.p12.genome.fa.gz 
  """
}

GTF.into{GTF_HST2INDEX;GTF_STRINDEX;GTF_CHAN;GTF_Isoform}
Genome.into{Genome_HST2INDEX;Genome_STRINDEX;Genome_GFFREAD}


/**
 * name: BuildHisat2Index
 * 
 * WhatItDoes:
 * 1 Extract splices sites from GTF file
 * 2 Build the Hisat 2 index
 * 
 * @param storeDir
 *  the directory where the GTF file will be stored
 *  in our case "Ressources"
 * @return output GTF a file to the gencode gtf file that will be in the GTF channel
 * 
 * Risks:
 * I/O
 * error 404
 * Memory ressources
 */
process HST2_INDEX{
  storeDir "Index/HISAT2"
  conda "envs/HISAT2.yaml"
  label "huge"
  input:
  file mygtf from GTF_HST2INDEX
  file GENOME from Genome_HST2INDEX
  
  output:
  file "${GenomeName}.SpliceSite.txt" into SSFile_MAP
  file "${GenomeName}.HISAT2.INDEX.*.ht2" into HST2INDEX
  script:
  GenomeName=GENOME[0].toString() - ~/.fa$/
  """
    hisat2_extract_splice_sites.py ${mygtf} > ${GenomeName}.SpliceSite.txt
    hisat2_extract_exons.py ${mygtf} > ${GenomeName}.Exons.txt
    hisat2-build -p ${task.cpus} --exon ${GenomeName}.Exons.txt --ss ${GenomeName}.SpliceSite.txt  --seed ${params.seed} $GENOME ${GenomeName}.HISAT2.INDEX
  """
}


/**
 * name: BuildStarIndex
 * 
 * WhatItDoes:
 * 1 Extract splices sites from GTF file
 * 2 Build the Hisat 2 index
 * 
 * @param storeDir
 *  the directory where the GTF file will be stored
 *  in our case "Ressources"
 * @return output GTF a file to the gencode gtf file that will be in the GTF channel
 * 
 * Risks:
 * I/O
 * error 404
 * memory ressources
 * 
 */
process BuildStarIndex{
  storeDir "Index/STAR"
  conda "envs/STAR.yaml"
  label "huge"
  input:
  file mygtf from GTF_STRINDEX
  file GENOME from Genome_STRINDEX
  
  output:
  file "GENOME_DIR/*" into STAR_INDEX
  script:
  GenomeName=GENOME[0].toString() - ~/.fa$/
  """
  mkdir GENOME_DIR
  STAR --runThreadN ${task.cpus} \
 --runMode genomeGenerate \
 --genomeDir GENOME_DIR \
 --genomeFastaFiles ${GENOME} \
 --sjdbGTFfile ${mygtf} \
 --sjdbOverhang 99 \
 """
}


//The string that look like the regexp of the files
FQ_REGEXP=params.RAWDIR+"*_"+params.FQ_PAIR_RADICAL+"{1,2}"+params.FQ_FILE_EXTENSION+params.FQ_COMPRESSION_EXTENSION
println("FastQ Regexp: "+FQ_REGEXP)
//Change the name of the channel to FastQ
Channel
    .fromFilePairs(FQ_REGEXP)
    .set { FastQ }
/**
 * name: CleanTrimGalore
 * 
 * WhatItDoes:
 * 0 Process the name of the file to extract something decent
 * 1 Run a FastQC
 * 2 Get the no hit, most likely the adapter
 * 3 Remove the adapter from the Fastq while preserving the pair structure
 * 
 * @param storeDir
 *  the directory where the cleaned fastq files will be stored
 *  
 * @return 
 *  cleaned fastq files to be map
 *  FASTQC rapport file to be store
 * 
 * Risks:
 * I/O
 * error 404
 * 
 */
process CleanTrimGalore{
  //a pretty name for the current process
  tag "$pair_id"
  //Where the files will be available
  storeDir "CleanedFastq"
  conda "envs/TRIMGALORE.yaml"
  label "Convert"
  input:
  set val(pair_id), file(FQ) from FastQ
  
  output:
  //Regex for the name of the output file
  set val(pair_id), file("${pair_id}_"+params.FQ_PAIR_RADICAL+"{1,2}_val_{1,2}.fq.gz") into FastQ2Map
  set val(pair_id), file("${pair_id}_"+params.FQ_PAIR_RADICAL+"{1,2}_val_{1,2}_fastqc.zip") into Nowhere
  
  
  script:
  TEST=params.FQ_FILE_EXTENSION.toString()+params.FQ_COMPRESSION_EXTENSION.toString()
  prefixR1=FQ[0].toString().replace(TEST,"")
  prefixR2=FQ[1].toString().replace(TEST,"")

  
  """
  fastqc ${FQ[0]}
  unzip ${prefixR1}_fastqc.zip
  mv ${prefixR1}_fastqc/fastqc_data.txt ${prefixR1}.fastqc.txt
  fastqc ${FQ[1]}
  unzip ${prefixR2}_fastqc.zip
  mv ${prefixR2}_fastqc/fastqc_data.txt ${prefixR2}.fastqc.txt
  ADAPT_R1=`cat ${prefixR1}.fastqc.txt | grep -e "^[ATGC]\\{50\\}" | grep -v "No Hit" | cut -f1`
  ADAPT_R2=`cat ${prefixR2}.fastqc.txt | grep -e "^[ATGC]\\{50\\}" | grep -v "No Hit" | cut -f1`
  
  if [ -z \$ADAPT_R1 ]
  then
   if [ -z \$ADAPT_R2 ]
   then
     trim_galore --paired --fastqc ${FQ}
   else
     trim_galore --paired --fastqc --adapter2 \$ADAPT_R2 ${FQ}
   fi
  else
   if [ -z \$ADAPT_R2 ]
   then
     trim_galore --paired --fastqc --adapter \$ADAPT_R1 ${FQ}
   else
     trim_galore --paired --fastqc --adapter \$ADAPT_R1 --adapter2 \$ADAPT_R2 ${FQ}
   fi
  fi
  """
}
//3 Mapping protocol 3 channels
FastQ2Map

if( params.style == "Hisat2" ){
  /**
   * name: MapWithHisat2
   * 
   * WhatItDoes:
   * 1 Map FastqFiles with Hisat2
   * 
   * @param storeDir
   *  the directory where the GTF file will be stored
   *  in our case "Ressources"
   * 
   * @input
   *  Reads in pairs 
   *  INDEX
   *  SPLICE SITE files
   * 
   * @return output Sam file into the SAM channel
   * 
   * Risks:
   * I/O
   * error 404
   * memory ressources
   * disk usage ressources Sam files are voluminous
   */
  process MapWithHisat2{
    storeDir="ALN_FILES/${params.style}"
    conda "envs/HISAT2.yaml"
    label="Mapping" 
    input:
    set val(pair_id), file(Reads) from FastQ2Map
    file INDEX from HST2INDEX.collect()
    file SS from SSFile_MAP.collect()
    
    output:
    file "${pair_id}.hisat2.sam" into RES
    
    script:
    LocalIndex=INDEX[0].toString() - ~/.1.ht2$/
    """
   hisat2 -x ${LocalIndex} -q --time --threads ${task.cpus} \
   -1 ${Reads[0]} -2 ${Reads[1]} \
   --known-splicesite-infile ${SS} \
   --novel-splicesite-outfile ${pair_id}.NSS.txt \
   --dta-cufflinks \
   --seed ${params.seed}\
   -S ${pair_id}.hisat2.sam
    
    """
  }
}else if( params.style == "OnePass" ){
  process MapWithSTR1PSS{
    tag "${PAIR_ID}"
    storeDir="ALN_FILES/${params.style}"
    conda "envs/STAR.yaml"
    label="Mapping"
    
    input:
    set val(pair_id), file(FQ) from FastQ2Map
    file RU from STAR_INDEX
    file GTF from GTF_CHAN.collect()
    
    output:
    file "${PAIR_ID}.STR.1PSS.Aligned.out.sam" into RES
    
    script:
    FILES=RU.toString().replaceAll(/[\[\]\,]/, "")
    PAIR_ID=pair_id.toString()
    
    """
    mkdir -p ./GenomeDir
    mv --target-directory=GenomeDir ${FILES}
    echo $pair_id
    ls ${FQ[0].toString()}
    ls ${FQ[1].toString()}
      
    STAR --readFilesIn ${FQ[0]} ${FQ[1]}\
    --readFilesCommand zcat\
    --runMode alignReads\
    --runThreadN ${task.cpus}\
    --runRNGseed ${params.seed}\
    --genomeLoad NoSharedMemory\
    --sjdbGTFfile ${GTF}\
    --twopassMode None\
    --outFileNamePrefix ${PAIR_ID}.STR.1PSS.
    """
  }

}else if( params.style == "TwoPass" ){
  process MapWithSTR2PSS{
    tag "${PAIR_ID}"
    storeDir="ALN_FILES/${params.style}"
    conda "envs/STAR.yaml"
    label="Mapping"
    
    input:
    set val(pair_id), file(FQ) from FastQ2Map
    file RU from STAR_INDEX.collect()
    file GTF from GTF_CHAN.collect()
    
    output:
    file "${PAIR_ID}.STR.2PSS.Aligned.out.sam" into RES
    
    script:
    FILES=RU.toString().replaceAll(/[\[\]\,]/, "")
    PAIR_ID=pair_id.toString()
    
    """
    mkdir -p ./GenomeDir
    mv --target-directory=GenomeDir ${FILES}
    echo $pair_id
    ls ${FQ[0].toString()}
    ls ${FQ[1].toString()}
      
    STAR --readFilesIn ${FQ[0]} ${FQ[1]}\
    --readFilesCommand zcat\
    --runMode alignReads\
    --runThreadN ${task.cpus}\
    --runRNGseed ${params.seed}\
    --genomeLoad NoSharedMemory\
    --sjdbGTFfile ${GTF}\
    --twopassMode Basic\
    --outFileNamePrefix ${PAIR_ID}.STR.2PSS.
    """
  }
}

/**
 * name: Sam2BAM
 * 
 * WhatItDoes:
 * 1 Convert a SAM files into a bam file
 *  
 * @input
 *  SAM .sam file flat mapping file
 * 
 * @return output .bam file into the bam_files channel
 * 
 * Risks:
 * I/O
 * error 404
 * memory ressources
 * 
 */
process Sam2BAM {
  tag "$prefix"
  storeDir="ALN_FILES/${params.style}"
  label="Convert"
  conda "envs/SAMTOOLS.yaml"
  input:
  file SAM from RES
  
  output:
  file "${prefix}.bam" into bam_files
    
  script:
  prefix=SAM.toString() - ~/.sam$/
  
  """
 samtools view -bS $SAM > ${prefix}.bam
  """
}


/**
 * name: Bam2SortedBAM
 * 
 * WhatItDoes:
 * 1 Sort a bam file on genomic coord
 *  
 * @input
 *  SAM .sam file flat mapping file
 * 
 * @return output .bam file into the bam_files channel
 * 
 * Risks:
 * I/O
 * error 404
 * memory ressources
 * 
 */
process Bam2SortedBAM {
  tag "$prefix"
  storeDir="ALN_FILES/${params.style}"
  label="Sort"
  conda "envs/SAMTOOLS.yaml"
  input:
  file BAM from bam_files
  
  output:
  file "${prefix}.sorted.bam" into Sbam_files
    
  script:
  prefix=BAM.toString() - ~/.bam$/
  
  """
samtools sort --threads ${task.cpus} -T ${prefix}.temp\
  -o ${prefix}.sorted.bam $BAM
  """
}

/**
 * name: FilterBamAndIndex
 * 
 * WhatItDoes:
 * 
 * 1 Filter the BAM
 *  based on http://broadinstitute.github.io/picard/explain-flags.html
 *   -f 3 keep read paired & mapped in proper pair
 *   -F 4 discard unmapped
 *   -F 8 discard mate ummaped
 *   -F 256 discard secondary or more alignement i.e. secondary guess for position
 *   -F 512 discard failed on platform/vendor quality control check
 *   -F 2048 discard supplementary alignement
 *   -q 60 keep reads with mapQ score higher or equal to 60
 *   -h include bam file header
 * 
 * 2 Index the bam file
 * 
 * @input
 *  SBAM .sam file flat mapping file
 * 
 * @return output .bam file into the bam_files channel
 * 
 * Risks:
 * I/O
 * error 404
 * memory ressources
 * 
 */
process FilterBamAndIndex {
  tag "$prefix"
  storeDir="ALN_FILES/${params.style}"
  label="Convert" 
  conda "envs/SAMTOOLS.yaml"
  input:
  file SBAM from Sbam_files
  
  output:
  file "${prefix}.filtered.sorted.bam" into FSbam_files
  file "${prefix}.filtered.sorted.bam.bai" into FSbamBai_files
    
  script:
  prefix=SBAM.toString() - ~/.sorted.bam$/
    
  """
 samtools view -b -h -f 3 -F 4 -F 8 -F 256 -F 512 -F 2048 \
  -q 60 ${SBAM} > ${prefix}.filtered.sorted.bam;\
  samtools index ${prefix}.filtered.sorted.bam
  """
}

/**
 * name: Isoform
 * 
 * WhatItDoes:
 * 
 * Call isoform from already known reference file and from mapping
 * It will produce a gene transfer format file GTF aka GFF2
 *  
 * @input
 *  - FSB .bam mapping file that need to be sorted by position, filter just to be 
 *  - GTF .gtf Reference file of already knowned transcripts
 *     1 for all studies
 * 
 * @return
 *  - StringtieSampleIsoforms .gtf file of the found transcript in the sample
 *      1 file per sample
 *      Names form the 
 *  - StringtieSampleIsoformsCoverage .coverage.gtf i dont know
 * 
 * Risks:
 * I/O
 * error 404
 * memory ressources
 * 
 */
process Isoform {
  tag "$prefix"
  storeDir="ALN_FILES/${params.style}"
  label="Isoform"
  conda "envs/STRINGTIE.yaml"
  input:
  file FSB from FSbam_files
  file GTF from GTF_Isoform.collect()
  
  
  output:
  file "${prefix}.Stringtie.gtf" into StringtieSampleIsoforms
  file "${prefix}.Stringtie.coverage.gtf" into StringtieSampleIsoformsCoverage
  
  script:
    prefix=FSB.toString() - ~/.filtered.sorted.bam$/
    
  """
  stringtie $FSB -G $GTF -p {task.cpus} -o ${prefix}.Stringtie.gtf -C ${prefix}.Stringtie.coverage.gtf --rf
  """
}

/**
 * name: MergeIsoforms
 * 
 * WhatItDoes:
 * 
 * Merge isoforms from all samples
 *  From https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual
 *   -c 2.5 minimum read coverage allowed for the predicted transcripts  
 *   -m 300 nucleotides minimum length
 *   -T min Trancript per million
 *   -f Minimum isoform fraction
 *   -l the name used for
 * 
 * Some parsing step to remove some comments lines
 * 
 * It will not keep the original name
 * 
 * @input
 *  - ISOs all the GTF files from all the sample files
 * 
 * @return
 *  - STRINGTIE_ISOFORMS one .gtf file
 *  one per study
 * 
 * Risks:
 * I/O
 * error 404
 * NFS delay
 */
process MergeIsoforms{
  beforeScript 'source /bioinfo/users/ltaing/.bashrc'
  storeDir "Isoforms"
  conda "envs/STRINGTIE.yaml"
  label="Isoform"
  input:
  file ISOs from StringtieSampleIsoforms.collect()
  output:
  file "Stringtie.${params.style}.Study.gtf" into STRINGTIE_ISOFORMS
  script:
  FILES=ISOs.toString().replaceAll(/[\[\]\,]/, "")
  
  """
  stringtie --merge -o merge.gtf -c 2.5 -m 300 -T 1 -f .001 -i $FILES -l strngt
 grep -v -e "^#" merge.gtf > Stringtie.${params.style}.Study.gtf
 sleep 1m
  """
}

/**
 * name: GTF_2_MRNA
 * 
 * WhatItDoes:
 * 
 * From a gtf with Transcripts and Exons
 * file extract the mrna
 *  
 * 
 * @input
 *  - STRINGTIE_ISOFORMS the GTF of the study
 *  - Genome_GFFREAD the fasta file of the genome the same that was used
 * for all process
 * 
 * @return
 *  - PutativeMRNA the fasta files of the mRNA
 * 
 * Risks:
 * I/O
 * error 404
 * memory ressources
 * 
 */
process GTF_2_MRNA{
  beforeScript 'source /bioinfo/users/ltaing/.bashrc'
  storeDir "Sequences/"
  label="Isoform"
  conda "envs/GFFREAD.txt"
  input:
  file PutativeGTF from STRINGTIE_ISOFORMS
  file GENOME from Genome_GFFREAD
  
  output:
  file "PutativeMRNA.${params.style}.nucleotides.fa" into PutativeMRNA
  script:
 """
  gffread ${PutativeGTF} -g ${GENOME} -w PutativeMRNA.${params.style}.nucleotides.fa
 """
}

Uniprot.into{Uniprot;UniprotPatch}


/**
 * name: MRNA_2_ProteinBank
 * 
 * WhatItDoes:
 * 
 * -Take some Nucleotide fasta file and translate it into a Peptide fasta
 * file using some in house R Script
 *  In the Script
 *   - 1) Translate all ORF
 *   - 2) Discard all the mRNA that have at least one perfect match in the
 * uniprot reference file
 *   - 3) Keep the longuest ORF for the database
 * 
 * @input
 *  - PutativeMRNA the Nucleotide fasta file of the study
 *  - Uniprot the fasta file of the genome the same that was used
 * for all process
 * 
 * @return
 *  - ProteinBank: the output of the Rscript
 *      the one that will be used for mass spec
 *       - the Uniprot entries that do not have any perfect match in the translated mrna
 *       - the Uniprot entries present in the translated mrna
 *       - the longuest ORF from the Mrna that does not have any match in UNIPROT
 *  - ProteinBankArchives: All the other files from the RScript files
 *     Not used for now
 *  - UniprotHuman: The Uniprot file is from the whole uniprot, it contains other species entries
 *  - UniprotMayBe: The Uniprot file with entries that are not in the mRNA
 *  - ToRescueQuery: The Longuest ORF of the mRNA that did not have perfect match in the UNIPROT
 *      - Those that might interest us because they are either unpreviously identified isoforms or SAAV or any other kind of modification
 * R
 * 
 * Risks:
 * I/O
 * error 404
 * memory ressources
 * 
 */
process MRNA_2_ProteinBank{
    beforeScript 'source /bioinfo/users/ltaing/.bashrc'
    conda "envs/Rsystempiper.txt"
    label "OneSmall"
    storeDir "Sequences/"
    input:
    file Library from PutativeMRNA
    file REFERENCE from Uniprot
    output:
    file "UniprotAndProteinBank.${params.style}.KnownAndMostLikely.fasta" into ProteinBank
    file "ProteinBank.*" into ProteinBankArchives
    file "Uniprot.9606.fasta" into UniprotHuman
    file "${params.style}.ToRescue.fasta" into ToRescueQuery
    file "SwissProt.NotInProteinBank.${params.style}MRNA.fasta" into SwissProtMayBe
    file "CreationProteinBank.${params.style}.log" into ProteinBankLog
    
    
    script:

  """
    Rscript ${PWD}/Scripts/mRNA_to_AA.R -i ${Library} -u ${REFERENCE} -o ProteinBank.${params.style}
    cp ProteinBank.${params.style}.SENSE.UNMATCH_UNIPROT.COHERENT_PROTEIN.UNIQUE_LORF.fasta ${params.style}.ToRescue.fasta
    cat SwissProt.NotInProteinBank.${params.style}MRNA.fasta ProteinBank.${params.style}.SENSE.MATCH_UNIPROT.fasta ProteinBank.${params.style}.SENSE.UNMATCH_UNIPROT.COHERENT_PROTEIN.UNIQUE_LORF.fasta > UniprotAndProteinBank.${params.style}.KnownAndMostLikely.fasta
    cp ProteinBank.${params.style}.log CreationProteinBank.${params.style}.log

  """
}



/**
 * name: BlastRescue
 * 
 * WhatItDoes:
 * 
 * -Take some Amino Acid fasta file and blast to Uniprot fasta
 * 
 * @input
 *  - ToRescueQuery the AA sequences that do not have perfect match in the Whole Uniprot
 *  - UniprotRESCUE,phr,pin,psq the fasta file and all the blast index files of all we know
 * 
 * @return
 *  - Blast output file "outfmt 6"
 *    -columns
 *      -QuerySeqId: the name of the blasted sequence
 *      -QueryLength: the length of the blasted sequence
 *      -QueryStart: the begenning of the match in the blasted sequence
 *      -QueryEnd: the end of the match in the blasted sequence
 *      -QuerySequence: the matching Amino Acid sequence of the blasted sequence
 *      -SubjectSeqId: the name the founded corresponding sequence
 *      -SubjectStart: the begenning of the match in the founded corresponding sequence
 *      -SubjectEnd: the end of the match in the founded corresponding sequence
 *      -SubjectLength: the length of the match in the founded corresponding sequence
 *      -SubjectSequence: the Amino Acid sequence of the match in the founded corresponding sequence
 *      -Length: Length of the matching part common between the blast and the found sequence
 *      -Mismatch: The number of missmatch between the 2 sequences
 *      -Gaps: The number of gap (Insertion/Deletion) between the 2 sequences
 *      -Evalue: I don't know, the lower the better
 *      -BitScore: I don't know, the higher the better
 * Risks:
 * I/O
 * error 404
 * memory ressources
 * 
 */
process BlastRescue{
  beforeScript "source /data/users/ltaing/.bashrc"
  afterScript "sleep 1m"
  storeDir "Sequences"
  conda "envs/BLAST.yml"
  label "Mapping"
  input:
  file "UniprotComplete.9606.fasta" from UniprotRESCUE
  file "UniprotComplete.9606.fasta.phr" from phr
  file "UniprotComplete.9606.fasta.pin" from pin
  file "UniprotComplete.9606.fasta.psq" from psq
  file "Query.Fasta" from ToRescueQuery
  output:
  file "${params.style}.Rescue.tsv" into Patch
  script:
  """
  echo -e "QuerySeqId\tQueryLength\tQueryStart\tQueryEnd\tQuerySequence\tSubjectSeqId\tSubjectStart\tSubjectEnd\tSubjectLength\tSubjectSequence\tLength\tMismatch\tGaps\tEvalue\tBitScore" > ${params.style}.Rescue.tsv
  blastp -db UniprotComplete.9606.fasta -query Query.Fasta -num_alignments 1 -num_threads ${task.cpus} -outfmt "6 qseqid qlen qstart qend qseq sseqid sstart send slen sseq length mismatch gaps evalue bitscore" >> ${params.style}.Rescue.tsv
  """
}

Patch.into{Patch;PatchForHTML}

/**
 * name: PatchTheProteinBank
 * 
 * WhatItDoes:
 *  take the some the
 *    protein bank,
 *    blast results and 
 *    add the sequences of the Hits in the whole UNIPROT
 * 
 * For the sequences that are identified as MRNA "Strngt" whithout
 * UNIPROT= and that have a blast hit in UNIPROT, the script will
 * annotated it with ",Blast=" and the sequence id
 * 
 * In the case where the mRNA/cDNA would be 5 prime degradated we might have
 * a Truncated form of the UNIPROT non canonical form (alternative start site)
 * 
 * I think it is also good to put the mRNA in concurrence with the UNIPROT blast hit
 * in case something went wrong with our isoform reconstruction.
 *  
 * @input
 *  - ProteinBank the AA fasta file that contains SwissProt canonical,
 *      the mRNA AA that match perfectly something in the UNIPROT, tagged by ",UNIPROT="
 *      the unmatched mRNA
 *        -those with BLAST hits
 *        -those whithout BLAST hits
 * 
 *  - UniprotPatch the fasta file of the genome the same that was used
 * for all process
 * 
 *  - Patch: the blast results for patching informations.
 * 
 * @return
 *  - PatchedProteinBank: The One protein bank
 *      Swissprot canonical
 *        if it doesn't match anything, it is there because Swissprot
 *        if it match anything in the mRNA, the mRNA is annotated with the Swissprot entry name
 *      The mRNA
 *        if it blast anything the mRNA is annotated with the UNIPROT entry name
 *        if it doesn't blast anything in UNIPROT it keeps it mRNA name
 *      The Uniprot concurrent
 *        if a mRNA blast something in UNIPROT and that is not something is swissprot, the blast hit will be added to the bank in case the 
 *  - MatchedWithBlast: The initial protein bank with the blast annotation
 *  - LookalikeProtMayBe the UNIPROT entries that had a hit from those who did not match
 *    Most of the time isoform, TR sometime with a SAAV or some minor changes
 * 
 * 
 * Risks:
 * Empty Output file
 * I/O
 * error 404
 * memory ressources
 * 
 */
process PatchTheProteinBank{
    beforeScript 'source /bioinfo/users/ltaing/.bashrc'
    conda "envs/Rsystempiper.txt"
    afterScript "sleep 1m"
    label "OneSmall"
    storeDir "Sequences/"
    input:
    file "UniprotAndProteinBank.${params.style}.KnownAndMostLikely.fasta" from ProteinBank
    file REFERENCE from UniprotPatch
    file BLAST from Patch
    output:
    file "UniprotAndProteinBank.${params.style}.KnownAndMostLikely.BlastPatched.LookalikeInUniprotNoSP.fasta" into PatchedProteinBank
    file "UniprotAndProteinBank.${params.style}.KnownAndMostLikely.BlastPatched.fasta" into MatchedWithBlast
    file "UniprotAndProteinBank.${params.style}.KnownAndMostLikely.LookalikeInUniprotNoSP.fasta" into LookalikeProtMayBe
    file "UniprotAndProteinBank.${params.style}.KnownAndMostLikely.log" into PatchLog
    
    script:

  """
    Rscript ${PWD}/Scripts/PatchProteinBankFromBlastResults.R \
    --Input UniprotAndProteinBank.${params.style}.KnownAndMostLikely.fasta\
    --Uniprot ${REFERENCE}\
    --BlastFile ${BLAST}\
    --Output UniprotAndProteinBank.${params.style}.KnownAndMostLikely
  """
}


process MakeLogs{
    afterScript "sleep 1m"
    storeDir "Sequences/"
    input:
    file UNILOG from UniprotCreationLog
    file BankLOG from ProteinBankLog
    file PATCHLOG from PatchLog
    output:
    file "${params.style}.log" into AllLog
    
    script:

  """
  cat ${UNILOG} | grep -e "^INFO" > ${params.style}.log
  cat ${BankLOG} | grep -e "^INFO" >> ${params.style}.log
  cat ${PATCHLOG} | grep -e "^INFO" >> ${params.style}.log
  
  """
}


PatchedProteinBank.into{PatchedProteinBank;PatchedProteinBankForRapport}

/**
 * name: MxQunt
 * 
 * WhatItDoes:
 * 
 * Mass Spec search withouth anything left.
 * Max Quant use to have a lot of leftovers un necessary files. Everything
 * is copied into the local scratch ${TMP_DIR}
 * 
 * Everything is process and only the results are kept
 * 
 * For now the job is functional but the way it s been written is not good with the Nextflow way
 * i.e. Only the conf need to be change according to the platform where it shoul run
 * 
 * @ Variable in use from config file
 *  -Spectro.RAW_Files_DIRECTORY thermofisher *.raw files directory all .raw files should be at the root of the directory
 *  -Spectro.SampleDescriptionFile A nice TSV sample sheet from the raw files
 *  -Spectro.TemplateMqparFile Spectro experimental design
 *  -Spectro.MonoSIMG the abs path to Mono's singularity image file to Use
 *  -Spectro.MaxQuant the abs path to MaxQuant cmd to Use
 * 
 * @input
 *  - PutativeMRNA ProteinBank of the Study
 * 
 * @return
 *  - ProteinBank: the output of the Rscript
 *      the one that will be used for mass spec
 *       - the Uniprot entries that do not have any perfect match in the translated mrna
 *       - the Uniprot entries present in the translated mrna
 *       - the longuest ORF from the Mrna that does not have any match in UNIPROT
 *  - ProteinBankArchives: All the other files from the RScript files
 *     Not used from now
 *  - UniprotHuman: The Uniprot file is from the whole uniprot, it contains other species entries
 *  - UniprotMayBe: The Uniprot file with entries that are not in the mRNA
 * 
 * Risks:
 * I/O
 * error 404
 * memory ressources
 * 
 */
process MxQunt{
  label="spectro"
  storeDir "Spectro/${params.style}"
  input:
  file DATABASE from PatchedProteinBank
  output:
  file "Used.QSUB.mqpar.xml"
  file "MQ_RES/*" into SPECTRO_RESULTS
  script:
  """
  mkdir /local/scratch/\${PBS_JOBID}/RAWFILES
  cp --dereference --target-directory=/local/scratch/\${PBS_JOBID}/RAWFILES ${params.Spectro.RAW_Files_DIRECTORY}/*.raw
  echo "raw files copy done"
  cp --dereference ${DATABASE} /local/scratch/\${PBS_JOBID}/.
  echo "database file copy done"
  python $PWD/Scripts/TemplateWithFiles.py -c\
 -f /local/scratch/\${PBS_JOBID}/${DATABASE}\
 -t /local/scratch/\${PBS_JOBID}\
 -w /local/scratch/\${PBS_JOBID}\
 -p ${task.cpus}\
 -s ${params.Spectro.SampleDescriptionFile}\
 -r /local/scratch/\${PBS_JOBID}/RAWFILES\
 -m ${params.Spectro.TemplateMqparFile}\
 -o /local/scratch/\${PBS_JOBID}/Temporary.mqpar.xml
 cp /local/scratch/\${PBS_JOBID}/Temporary.mqpar.xml Used.QSUB.mqpar.xml 
 singularity exec ${params.Spectro.MonoSIMG} mono ${params.Spectro.MaxQuant} /local/scratch/\${PBS_JOBID}/Temporary.mqpar.xml
 mkdir -p MQ_RES
 mv /local/scratch/\${PBS_JOBID}/fixedCombinedFolder/combined/txt/* MQ_RES/.
  """
}

Channel.fromPath("Spectro/${params.style}/MQ_RES/peptides.txt").set{SPECTRO_RESULTS}

process MakeRapportMarkdown{
  label="OneSmall"
  conda "envs/Rsystempiper.txt"
  storeDir "Rapport"
  input:
  file DATABASE from PatchedProteinBankForRapport
  file PEPTIDES from SPECTRO_RESULTS
  file RESCUE from PatchForHTML
  file LOG from AllLog
  output:
  file "${params.style}.html"
  file "${params.style}.md"
  file "${params.style}_files/**"
  script:
  """
  Rscript ${PWD}/Scripts/ComputeResults.R -O \${PWD}/${params.style} -M \${PWD}/${PEPTIDES} -P \${PWD}/${DATABASE} -B \${PWD}/${RESCUE} -L \${PWD}/${LOG}
  """
}
