workflow ProcessHaloplexHS {

    File SampleSheet
    # sample sheet has this structure:
    # index  sample    other stuff

    Array[Array[String]] inputData = read_tsv(SampleSheet)
    Array[String] Adapters = ["GATCGGAAGAGCACACGTCTGAACTCCAGTCAC","AGATCGGAAGAGCGTCGTGTAGGGAAA"]
    
    String IlluminaDir
    String JobGroup
    String OutputDir
    String Queue
    
    String Reference    = "/gscmnt/gc2709/info/production_reference_GRCh37-lite/CLE/Haloplex/reference/all_sequences.fa"
    String VEP          = "/gscmnt/gc2709/info/production_reference_GRCh37-lite/CLE/Haloplex/VEP_cache"
    String QcMetrics    = "/gscmnt/gc2709/info/production_reference_GRCh37-lite/CLE/Haloplex/git/cle-myeloseq/accessory_files/MyeloseqQCMetrics.txt"
    String Description  = "/gscmnt/gc2709/info/production_reference_GRCh37-lite/CLE/Haloplex/git/cle-myeloseq/accessory_files/MyeloseqDescription.1.1.txt"
    String HaplotectBed = "/gscmnt/gc2709/info/production_reference_GRCh37-lite/CLE/Haloplex/git/cle-myeloseq/accessory_files/myeloseq.haplotect_snppairs.041718.bed"
    String AmpliconBed  = "/gscmnt/gc2709/info/production_reference_GRCh37-lite/CLE/Haloplex/git/cle-myeloseq/accessory_files/04818-1516117769_Amplicon.b37.bed"
    String TargetBed    = "/gscmnt/gc2709/info/production_reference_GRCh37-lite/CLE/Haloplex/git/cle-myeloseq/accessory_files/04818-1516117769_Covered.b37.bed"
    String CoverageBed  = "/gscmnt/gc2709/info/production_reference_GRCh37-lite/CLE/Haloplex/git/cle-myeloseq/accessory_files/CoverageQC.b37.022219.bed"

    String CustomAnnotationVcf   = "/gscmnt/gc2709/info/production_reference_GRCh37-lite/CLE/Haloplex/git/cle-myeloseq/accessory_files/myeloseq_custom_annotations.annotated.011618.b37.vcf.gz"
    String CustomAnnotationIndex = "/gscmnt/gc2709/info/production_reference_GRCh37-lite/CLE/Haloplex/git/cle-myeloseq/accessory_files/myeloseq_custom_annotations.annotated.011618.b37.vcf.gz.tbi"
    
    String CustomAnnotationParameters = "MYELOSEQ,vcf,exact,0,TCGA_AC,MDS_AC,MYELOSEQBLACKLIST"
    

    call barcode_demux {
        input: Dir=IlluminaDir, #Fastqs=get_fastq_files.fastq_files,
               SampleSheet=SampleSheet,
               SampleIndexMM=0,
               queue=Queue,
               jobGroup=JobGroup
    }

    call prepare_samples {
        input: SampleSheet=SampleSheet,
               Fastq1=barcode_demux.read1_fastqs,
               Fastq2=barcode_demux.read2_fastqs,
               queue=Queue,
               jobGroup=JobGroup
    }

    scatter (samples in inputData){
        call trim_reads {
            input: Index=samples[0],
                   Name=samples[1],
                   SampleSheetFile=prepare_samples.sample_sheet,
                   Adapters=Adapters,
                   queue=Queue,
                   jobGroup=JobGroup
        }

        call align_barcode_and_sort_reads {
            input: Fastq=trim_reads.fastq_file,
                   Name=samples[1],
                   readGroup=samples[2],
                   refFasta=Reference,
                   queue=Queue,
                   jobGroup=JobGroup
        }

        call consensus_bam {
            input: Bam=align_barcode_and_sort_reads.bam_file,
                   BamIndex=align_barcode_and_sort_reads.bam_index,
                   TargetBed=TargetBed,
                   AmpliconBed=AmpliconBed,
                   Name=samples[1],
                   refFasta=Reference,
                   queue=Queue,
                   jobGroup=JobGroup
        }

        
        call run_varscan {
            input: Bam=consensus_bam.bam_file,
                   BamIndex=consensus_bam.bam_index,
                   CoverageBed=CoverageBed,
                   refFasta=Reference,
                   Name=samples[1],
                   queue=Queue,
                   jobGroup=JobGroup
        }

        call clean_variants as clean_varscan_indels {
            input: Vcf=run_varscan.varscan_indel_file,
                   Name=samples[1] + ".varscan_indel",
                   refFasta=Reference,
                   queue=Queue,
                   jobGroup=JobGroup
        }

        call run_platypus {
            input: Bam=consensus_bam.bam_file,
                   BamIndex=consensus_bam.bam_index,
                   CoverageBed=CoverageBed,
                   Name=samples[1],
                   refFasta=Reference,
                   queue=Queue,
                   jobGroup=JobGroup
        }

        call clean_variants as clean_platypus {
            input: Vcf=run_platypus.platypus_vcf_file,
                   Name=samples[1] + ".platypus",
                   refFasta=Reference,
                   queue=Queue,
                   jobGroup=JobGroup
        }

        call run_pindel_region as run_pindel_flt3itd {
            input: Bam=consensus_bam.bam_file,
                   BamIndex=consensus_bam.bam_index,
                   Reg='13:28608124-28608453',
                   refFasta=Reference,
                   Name=samples[1],
                   queue=Queue,
                   jobGroup=JobGroup
        }

        call clean_variants as clean_pindel_itd {
            input: Vcf=run_pindel_flt3itd.pindel_vcf_file,
                   Name=samples[1],
                   refFasta=Reference,
                   queue=Queue,
                   jobGroup=JobGroup
        }

        call combine_variants {
            input: VarscanSNV=run_varscan.varscan_snv_file,
                   VarscanIndel=clean_varscan_indels.cleaned_vcf_file,
                   PindelITD=clean_pindel_itd.cleaned_vcf_file,
                   Platypus=clean_platypus.cleaned_vcf_file,
                   Bam=consensus_bam.bam_file,
                   BamIndex=consensus_bam.bam_index,
                   refFasta=Reference,
                   Name=samples[1],
                   queue=Queue,
                   jobGroup=JobGroup
        }


        call run_vep {
            input: CombineVcf=combine_variants.combined_vcf_file,
                   refFasta=Reference,
                   Vepcache=VEP,
                   CustomAnnotationVcf=CustomAnnotationVcf,
                   CustomAnnotationIndex=CustomAnnotationIndex,
                   CustomAnnotationParameters=CustomAnnotationParameters,
                   Name=samples[1],
                   queue=Queue,
                   jobGroup=JobGroup
        }

        call run_haplotect {
            input: refFasta=Reference,
                   Bam=consensus_bam.bam_file,
                   Bed=HaplotectBed,
                   Name=samples[1],
                   queue=Queue,
                   jobGroup=JobGroup
        }

        
        call haloplex_qc {
            input: refFasta=Reference,
                   AlignedBam=align_barcode_and_sort_reads.bam_file,
                   ConsensusBam=consensus_bam.bam_file,
                   AnnotatedTsv=run_vep.annotated_filtered_tsv,
                   TargetBed=TargetBed,
                   AmpliconBed=AmpliconBed,
                   CoverageBed=CoverageBed,
                   QcMetrics=QcMetrics,
                   Description=Description,
                   Name=samples[1],
                   DemuxFile=prepare_samples.sample_sheet,
                   Haplotect=run_haplotect.out_file,
                   HaplotectSites=run_haplotect.sites_file,
                   queue=Queue,
                   jobGroup=JobGroup
        }

        call make_reports {
            input: Variants=run_vep.annotated_filtered_tsv,
                   QC=haloplex_qc.coverage_qc_json_file,
                   Description=Description,
                   Name=samples[1],
                   queue=Queue,
                   jobGroup=JobGroup
        }

        call igv_session {
            input: Name=samples[1],
                   queue=Queue,
                   jobGroup=JobGroup
        }

        call gather_files {
            input: OutputFiles=[trim_reads.demux_fastq1_file,
                   trim_reads.demux_fastq2_file,
                   align_barcode_and_sort_reads.bam_file,
                   align_barcode_and_sort_reads.bam_index,
                   consensus_bam.bam_file,
                   consensus_bam.bam_index,
                   haloplex_qc.coverage_qc_file,
                   haloplex_qc.coverage_qc_json_file,
                   run_varscan.varscan_snv_file,
                   clean_varscan_indels.cleaned_vcf_file,
                   clean_pindel_itd.cleaned_vcf_file,
                   clean_platypus.cleaned_vcf_file,
                   combine_variants.combined_vcf_file,
                   run_vep.annotated_vcf,
                   run_vep.annotated_filtered_vcf,
                   run_vep.annotated_filtered_tsv,
                   make_reports.variant_report,
                   igv_session.igv_xml],
                   OutputDir=OutputDir,
                   SubDir=samples[1] + "_" + samples[0],
                   queue=Queue,
                   jobGroup=JobGroup
        }
    }
}

task barcode_demux {
     String Dir
     String SampleSheet 
     Int SampleIndexMM
     String jobGroup
     String queue

     command <<<
          I1=$(ls ${Dir}/*I1*.fastq.gz) && \
          I2=$(ls ${Dir}/*I2*.fastq.gz) && \
          R1=$(ls ${Dir}/*R1*.fastq.gz) && \
          R2=$(ls ${Dir}/*R2*.fastq.gz) && \
          cut -f 1 ${SampleSheet} > index_list.txt && \
          python /usr/bin/demultiplexer_hja_bd.py \
          --mismatches ${SampleIndexMM} --index-read-file $I1 --tag-file index_list.txt \
          --read1-file $R1 --read2-file $R2 --read3-file $I2
     >>>
     runtime {
         docker_image: "registry.gsc.wustl.edu/fdu/haloplex-demux:1"
         cpu: "1"
         memory_gb: "12"
         queue: queue
         resource: "rusage[gtmp=10, mem=12000]"
         job_group: jobGroup 
}
     output {
         Array[File] read1_fastqs = glob("*.1.fq.gz") 
         Array[File] read2_fastqs = glob("*.2.fq.gz")
         File unknown_read1 = "unknown.1.fq.gz" 
         File unknown_read2 = "unknown.2.fq.gz" 
     }
}

task prepare_samples {
     File SampleSheet
     Array[File] Fastq1
     Array[File] Fastq2
     String jobGroup
     String queue

     command <<<
             /bin/cat ${write_tsv(Fastq1)} > 1.tmp.txt
             /bin/cat ${write_tsv(Fastq2)} > 2.tmp.txt
             /usr/bin/perl -e 'open(R1,"1.tmp.txt"); @r1 = <R1>; \
                 chomp @r1; close R1;\
                 open(R2,"2.tmp.txt"); @r2 = <R2>; \
                 chomp @r2; close R2; \
                 open(SS,"${SampleSheet}");
                 while(<SS>){
                     chomp;
                     my @l = split("\t",$_);
                     my $r1 = (grep /$l[0].1.fq.gz/, @r1)[0];
                     my $r2 = (grep /$l[0].2.fq.gz/, @r2)[0];
                     my $persamplereads1 = `gunzip -c $r1 | wc -l`;
                     chomp $persamplereads1;
                     my $persamplereads2 = `gunzip -c $r2 | wc -l`;
                     chomp $persamplereads2;
                     print join("\t",$l[0],$l[1],$r1,$r2,($persamplereads1 / 4) + ($persamplereads2 / 4)),"\n";
                 }
                 close SS;
                 my $r1 = (grep /unknown.1.fq.gz/, @r1)[0];
                 my $r2 = (grep /unknown.2.fq.gz/, @r2)[0];
                 my $persamplereads1 = `gunzip -c $r1 | wc -l`;
                 chomp $persamplereads1;
                 my $persamplereads2 = `gunzip -c $r2 | wc -l`;
                 chomp $persamplereads2;
                 print join("\t","unknown","unknown",$r1,$r2,($persamplereads1 / 4) + ($persamplereads2 / 4)),"\n"' > sample_sheet.txt
     >>>
     runtime {
         docker_image: "registry.gsc.wustl.edu/genome/lims-compute-xenial:1"
         cpu: "1"
         memory_gb: "4"
         queue: queue
         resource: "rusage[gtmp=10, mem=4000]"
         job_group: jobGroup
     }
     output {
         File sample_sheet = "sample_sheet.txt"
         Array[Array[String]] sample_data = read_tsv("sample_sheet.txt")
     }
}

task trim_reads {
     String Index
     File SampleSheetFile
     Array[String] Adapters
     String Name
     Int? TrimN
     String jobGroup
     String queue

     command {
         export PYTHONPATH=/opt/cutadapt/lib/python2.7/site-packages/ && \
         cp $(/bin/grep ${Index} ${SampleSheetFile} | cut -f 3) "${Name}.1.fastq.gz" && \
         cp $(/bin/grep ${Index} ${SampleSheetFile} | cut -f 4) "${Name}.2.fastq.gz" && \	
         /opt/cutadapt/bin/cutadapt --interleaved -a ${Adapters[0]} -A ${Adapters[1]} \
         $(/bin/grep ${Index} ${SampleSheetFile} | cut -f 3) $(/bin/grep ${Index} ${SampleSheetFile} | cut -f 4) | \
         /opt/cutadapt/bin/cutadapt --interleaved -u ${default=3 TrimN} -u -${default=3 TrimN} -U ${default=3 TrimN} -U -${default=3 TrimN} -m 30 -o trimmed.fq.gz - 
     }
     runtime {
         docker_image: "registry.gsc.wustl.edu/fdu/cutadapt:1"
         cpu: "1"
         memory_gb: "8"
         queue: queue
         resource: "rusage[gtmp=10, mem=8000]"
         job_group: jobGroup 
     }
     output {
         File fastq_file = "trimmed.fq.gz"
         File demux_fastq1_file = "${Name}.1.fastq.gz"
         File demux_fastq2_file = "${Name}.2.fastq.gz"
     }
}

task align_barcode_and_sort_reads {
     String Fastq
     String Name
     String refFasta
     String readGroup
     String jobGroup
     String queue
     
     command <<<     	     
         (set -eo pipefail && /usr/local/bin/bwa mem -M -t 8 -p ${refFasta} ${Fastq} -R "${readGroup}" | \
         /usr/bin/perl -ane '{ if (/^@/) { print join("\t", @F),"\n"; 
             } elsif ($F[0] =~ /:([ACGT]{8,10})$/){ 
             push @F, "X0:Z:$1"; 
             $F[0] =~ s/:[ACGT]{8,10}$//; 
             print join("\t", @F), "\n"; 
         }}' | /usr/local/bin/samtools view -b -S /dev/stdin | \
         /usr/local/bin/samtools sort -m 1G -O bam - > ${Name}.aligned_sorted.bam) && \
         /usr/local/bin/samtools index ${Name}.aligned_sorted.bam
     >>>
     runtime {
         docker_image: "registry.gsc.wustl.edu/genome/tagged-alignment:2"
         cpu: "8"
         memory_gb: "32"
         queue: queue
         resource: "rusage[gtmp=10, mem=32000]"
         job_group: jobGroup
     }
     output {
         File bam_file = "${Name}.aligned_sorted.bam"
         File bam_index = "${Name}.aligned_sorted.bam.bai"
     }
}

task haloplex_qc {
     String refFasta
     String DemuxFile
     String AlignedBam
     String ConsensusBam 
     String TargetBed
     String AmpliconBed
     String CoverageBed
     String QcMetrics
     String Description
     String AnnotatedTsv
     String Name
     String Haplotect
     String HaplotectSites
     String jobGroup
     String queue

     command {
         /usr/bin/perl /usr/local/bin/CalculateCoverageQC.pl -r ${refFasta} -d ${DemuxFile} \
         -a ${AmpliconBed} -t ${TargetBed} -b ${CoverageBed} -c ${ConsensusBam} -w ${AlignedBam} \
         -q ${QcMetrics} -i ${Description} -v ${AnnotatedTsv} -m ${Haplotect} -p ${HaplotectSites}
     }
     runtime {
         docker_image: "registry.gsc.wustl.edu/fdu/haloplex-qc:3"
         cpu: "1"
         memory_gb: "16"
         queue: queue
         resource: "rusage[gtmp=10, mem=16000]"
         job_group: jobGroup
     } 
     output {
         File coverage_qc_file = "${Name}.qc.txt"
         File coverage_qc_json_file = "${Name}.qc.json"
     }
}

task run_varscan {
     File Bam
     File BamIndex
     Int? MinCov
     Float? MinFreq
     Int? MinReads
     String CoverageBed
     String refFasta
     String Name
     String jobGroup
     String queue

     command <<<
         /usr/bin/samtools mpileup -q 1 -f ${refFasta} -l ${CoverageBed} ${Bam} > /tmp/mpileup.out && \
         java -Xmx12g -jar /opt/varscan/VarScan.jar mpileup2snp /tmp/mpileup.out --min-coverage ${default=8 MinCov} --min-reads2 ${default=5 MinReads} \
         --min-var-freq ${default="0.01" MinFreq} --output-vcf > ${Name}.snv.vcf && \
         java -Xmx12g -jar /opt/varscan/VarScan.jar mpileup2indel /tmp/mpileup.out --min-coverage ${default=8 MinCov} --min-reads2 ${default=5 MinReads} \
         --min-var-freq ${default="0.01" MinFreq} --output-vcf > ${Name}.indel.vcf
     >>>

     runtime {
         docker_image: "registry.gsc.wustl.edu/fdu/varscan-2.4.2-samtools-1.3.1:1"
         cpu: "1"
         memory_gb: "16"
         queue: queue
         resource: "rusage[gtmp=10, mem=16000]"
         job_group: jobGroup
     }
     output {
         File varscan_snv_file = "${Name}.snv.vcf"
         File varscan_indel_file = "${Name}.indel.vcf"
     }
}

task run_pindel_region {
     File Bam
     File BamIndex
     String Reg
     Int? Isize
     Int? MinReads
     String refFasta
     String Name
     String jobGroup
     String queue

    command <<<
        (set -eo pipefail && /usr/local/bin/samtools view ${Bam} ${Reg} | /opt/pindel-0.2.5b8/sam2pindel - /tmp/in.pindel ${default=250 Isize} tumor 0 Illumina-PairEnd) && \
        /usr/bin/pindel -f ${refFasta} -p /tmp/in.pindel -c ${Reg} -o /tmp/out.pindel && \
        /usr/bin/pindel2vcf -P /tmp/out.pindel -G -r ${refFasta} -e ${default=3 MinReads} -R GRCh37 -d GRCh37 -v /tmp/out.vcf && \
        sed 's/END=[0-9]*;//' /tmp/out.vcf > ${Name}.pindel.vcf
    >>>

    runtime {
        docker_image: "registry.gsc.wustl.edu/fdu/pindel2vcf-0.6.3:1"
        cpu: "1"
        memory_gb: "16"
        queue: queue
        resource: "rusage[gtmp=10, mem=16000]"
        job_group: jobGroup
    }
    output {
        File pindel_vcf_file = "${Name}.pindel.vcf"
    }
}

task consensus_bam {
     File Bam
     File BamIndex
     String TargetBed
     String AmpliconBed
     String Name
     String refFasta
     String jobGroup
     String queue

     command {
         /usr/bin/java -Xmx6g \
         -jar /opt/gatk/public/external-example/target/external-example-1.0-SNAPSHOT.jar \
         -T WalkerTRConsensus_wk5 -I ${Bam} \
         -L ${TargetBed} --ampliconBed ${AmpliconBed} -dcov 1000000 \
         -cbam ${Name}.consensus.bam -maxNM 5 -mmq 10 \
         -R ${refFasta} > condense.log.txt && \
         /usr/bin/samtools index ${Name}.consensus.bam
     }

     runtime {
         docker_image: "registry.gsc.wustl.edu/fdu/haloplex-walker:1"
         cpu: "1"
         memory_gb: "8"
         queue: queue
         resource: "rusage[gtmp=10, mem=8000]"
         job_group: jobGroup
     }
     output {
         File bam_file = "${Name}.consensus.bam"
         File bam_index = "${Name}.consensus.bam.bai"
         File log = "condense.log.txt"
     }
}

task run_platypus {
     File Bam
     File BamIndex
     String CoverageBed
     String? DocmVcf
     String Name
     String refFasta
     String jobGroup
     String queue

     command <<<
         /usr/bin/awk '{ print $1":"$2+1"-"$3; }' ${CoverageBed} > "regions.txt" && \
         /usr/bin/python /opt/platypus/Platypus_0.8.1/Platypus.py callVariants --bamFiles=${Bam} \
         --regions regions.txt --refFile=${refFasta} \
         --nCPU 1 ${"--source=" + DocmVcf} --output="${Name}.platypus.vcf" \
         --minReads 5 --filterDuplicates=0 --minFlank=10 --assemble 1 --filterReadPairsWithSmallInserts=0 \
         --trimOverlapping=0 --trimSoftClipped=0 --filterReadsWithDistantMates=0 --filterReadsWithUnmappedMates=0
     >>>
     runtime {
         docker_image: "registry.gsc.wustl.edu/fdu/platypus-0.8.1:1"
         cpu: "4"
         memory_gb: "16"
         queue: queue
         resource: "rusage[gtmp=10, mem=16000]"
         job_group: jobGroup
     }
     output {
         File platypus_vcf_file = "${Name}.platypus.vcf"
     }
}

task clean_variants {
     String Vcf
     String Name
     String refFasta
     String jobGroup
     String queue

     command {
         /usr/bin/java -Xmx16g -jar /opt/GenomeAnalysisTK.jar -T LeftAlignAndTrimVariants -R ${refFasta} --splitMultiallelics --variant ${Vcf} -o /tmp/out.vcf && \
         /usr/bin/java -Xmx16g -jar /opt/GenomeAnalysisTK.jar -T VariantsToAllelicPrimitives -R ${refFasta} -V /tmp/out.vcf -o ${Name}.cleaned.vcf
     }
     runtime {
         docker_image: "registry.gsc.wustl.edu/genome/gatk-3.6:1"
         cpu: "1"
         memory_gb: "16"
         queue: queue
         resource: "rusage[gtmp=10, mem=16000]"
         job_group: jobGroup
     }
     output {
         File cleaned_vcf_file = "${Name}.cleaned.vcf"
     }
}

task combine_variants {
     String VarscanSNV
     String VarscanIndel
     String PindelITD
     String Platypus
     String Bam
     String BamIndex
     String refFasta
     String Name
     String jobGroup
     String queue

     command {
         /usr/bin/java -Xmx8g -jar /opt/GenomeAnalysisTK.jar -T CombineVariants -R ${refFasta} --variant:varscanIndel ${VarscanIndel} \
         --variant:varscanSNV ${VarscanSNV} --variant:Platypus ${Platypus} --variant:PindelITD ${PindelITD} -o combined.vcf --genotypemergeoption UNIQUIFY && \
         /usr/bin/python /usr/bin/addAmpliconInfoAndCountReads.py -r ${refFasta} combined.vcf ${Bam} ${Name} > ${Name}.combined_and_tagged.vcf
     }
     runtime {
         docker_image: "registry.gsc.wustl.edu/fdu/gatk-3.6-biopython-pysam-scipy-amplicon-readcount:2"
         cpu: "1"
         memory_gb: "10"
         queue: queue
         resource: "rusage[gtmp=10, mem=10000]"
         job_group: jobGroup
     }
     output {
         File combined_vcf_file = "${Name}.combined_and_tagged.vcf"
     }

}

task run_vep {
     File CombineVcf
     String refFasta
     String Vepcache
     File CustomAnnotationVcf
     File CustomAnnotationIndex
     String CustomAnnotationParameters
     Float? maxAF
     String Name
     String jobGroup
     String queue

     command {
         if [ $(/bin/grep -v '^#' ${CombineVcf}|/usr/bin/wc -l) == 0 ]; then
             /bin/cp ${CombineVcf} ${Name}.annotated.vcf && \
             /bin/cp ${CombineVcf} ${Name}.annotated_filtered.vcf && \
             /opt/htslib/bin/bgzip ${Name}.annotated.vcf && /usr/bin/tabix -p vcf ${Name}.annotated.vcf.gz && \
             /opt/htslib/bin/bgzip ${Name}.annotated_filtered.vcf && /usr/bin/tabix -p vcf ${Name}.annotated_filtered.vcf.gz && \
             /usr/bin/touch ${Name}.variants_annotated.tsv
         else
             /usr/bin/perl -I /opt/lib/perl/VEP/Plugins /usr/bin/variant_effect_predictor.pl --format vcf \
             --vcf --plugin Downstream --plugin Wildtype --fasta ${refFasta} --hgvs --symbol --term SO --flag_pick \
             -i ${CombineVcf} --custom ${CustomAnnotationVcf},${CustomAnnotationParameters} --offline --cache --max_af --dir ${Vepcache} -o ${Name}.annotated.vcf && \
             /opt/htslib/bin/bgzip ${Name}.annotated.vcf && /usr/bin/tabix -p vcf ${Name}.annotated.vcf.gz && \
             /usr/bin/perl -I /opt/lib/perl/VEP/Plugins /opt/vep/ensembl-vep/filter_vep -i ${Name}.annotated.vcf.gz --format vcf \
             --filter "(MAX_AF < ${default='0.001' maxAF} or not MAX_AF) or MYELOSEQ_TCGA_AC or MYELOSEQ_MDS_AC" -o ${Name}.annotated_filtered.vcf && \
             /opt/htslib/bin/bgzip ${Name}.annotated_filtered.vcf && /usr/bin/tabix -p vcf ${Name}.annotated_filtered.vcf.gz && \
             /usr/bin/java -Xmx4g -jar /opt/GenomeAnalysisTK.jar -T VariantsToTable \
             -R ${refFasta} --showFiltered --variant ${Name}.annotated_filtered.vcf.gz -o /tmp/variants.tsv \
             -F CHROM -F POS -F ID -F FILTER -F REF -F ALT -F set -GF TAMP -GF SAMP -GF VAFTYPE -GF CVAF -GF NR -GF NV && \
             /usr/bin/python /usr/bin/add_annotations_to_table_helper.py /tmp/variants.tsv ${Name}.annotated_filtered.vcf.gz \
             Consequence,SYMBOL,EXON,INTRON,Feature_type,Feature,HGVSc,HGVSp,HGNC_ID,MAX_AF,MYELOSEQ_TCGA_AC,MYELOSEQ_MDS_AC /tmp/ && \
             mv /tmp/variants.annotated.tsv ${Name}.variants_annotated.tsv
         fi
     }
     runtime {
         docker_image: "registry.gsc.wustl.edu/fdu/vep90-gatk3.6-htslib1.3.2:1"
         cpu: "1"
         memory_gb: "10"
         queue: queue
         resource: "rusage[gtmp=10, mem=10000]"
         job_group: jobGroup
     }
     output {
         File annotated_vcf = "${Name}.annotated.vcf.gz"
         File annotated_filtered_vcf = "${Name}.annotated_filtered.vcf.gz"
         File annotated_filtered_tsv = "${Name}.variants_annotated.tsv"
     }
}

task run_haplotect {
     String Bam
     String Bed
     String Name
     String refFasta
     String jobGroup
     String queue

     Int? MinReads

     command <<<
             /usr/bin/awk -v OFS="\t" '{ $2=$2-1; print; }' ${Bed} > /tmp/pos.bed && \
             /usr/bin/java -Xmx6g \
             -jar /opt/gatk/public/external-example/target/external-example-1.0-SNAPSHOT.jar \
             -T Haplotect -R ${refFasta} \
             -mmq 20 -mbq 20 -dcov 20000 -unif 0 -gstol 0.001 \
             -mr ${default=100 MinReads} \
             -I ${Bam} \
             -htp ${Bed} \
             -L /tmp/pos.bed \
             -outPrefix ${Name} && \
             sort -u "${Name}.multihaploci.txt" > "${Name}.haplotectloci.txt" && \
             /usr/bin/perl /opt/CalculateContamination.pl "${Name}.txt" "${Name}.haplotectloci.txt" > "${Name}.haplotect.txt"
     >>>

     runtime {
             docker_image: "registry.gsc.wustl.edu/fdu/haloplex-walker:1"
             cpu: "1"
             memory_gb: "8"
             queue: queue
             resource: "rusage[gtmp=10, mem=8000]"
             job_group: jobGroup
     }
     output {
            File out_file = "${Name}.haplotect.txt"
            File sites_file = "${Name}.haplotectloci.txt"
     }
}


task make_reports {
     File Variants
     File QC
     File Description
     String Name
     String jobGroup
     String queue

     command {
         /usr/bin/perl /usr/local/bin/FormatClinicalReport.pl ${Variants} ${QC} ${Description} > ${Name}.variant_report.txt
     }

     runtime {
         docker_image: "registry.gsc.wustl.edu/fdu/haloplex-qc:2"
         queue: queue
         job_group: jobGroup
     } 

     output {
         File variant_report = "${Name}.variant_report.txt"
     }
}

task igv_session {
     String Name
     String jobGroup
     String queue

     command {
         cat <<EOF > ${Name}.igv.xml
<?xml version="1.0" encoding="UTF-8"?>
<Session genome="b37" locus="All" version="3">
    <Resources>
     <Resource name="All variants" path="${Name}.annotated.vcf.gz"/>
     <Resource name="Filtered variants" path="${Name}.annotated_filtered.vcf.gz"/>
     <Resource name="${Name}" path="${Name}.consensus.bam"/>
     <Resource name="Ensemble Genes" path="http://www.broadinstitute.org/igvdata/annotations/hg19/EnsemblGenes.ensGene"/>
    </Resources>
</Session>
EOF
     }

     runtime {
         docker_image: "registry.gsc.wustl.edu/genome/lims-compute-xenial:1"
         queue: queue
         job_group: jobGroup
     } 

     output {
         File igv_xml = "${Name}.igv.xml"
     }
}

task gather_files {
     Array[String] OutputFiles
     String OutputDir
     String? SubDir
     String jobGroup
     String queue

     command {
         if [[ ${SubDir} != "" ]] && [[ ! -e ${OutputDir}/${SubDir} ]]; then
             mkdir ${OutputDir}/${SubDir}
         fi
         /bin/mv -f -t ${OutputDir}/${SubDir} ${sep=" " OutputFiles}
     }
     runtime {
         docker_image: "registry.gsc.wustl.edu/genome/lims-compute-xenial:1"
         queue: queue
         job_group: jobGroup
     }
     output {
         String out = stdout()
     }
}

