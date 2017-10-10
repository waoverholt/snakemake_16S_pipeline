import glob
import os
import sys

#indentify configure file you'd like to use with this snakefile
configfile: "config.yaml"

#This should be the directory that contains your raw reads (expected to be .gz)
INPUT_DIR = config["read_directory"]

#Grab just the useful sample name and insert them into a list
postfix_length = len(config["raw_file_forward_extension"])
sample_names = {f[:-postfix_length]
    for f in os.listdir(config["read_directory"])
    if f.endswith(config["raw_file_forward_extension"])
    }

#FORWARD_READ = '{sample}' + config["raw_file_forward_extension"]
#REVERSE_READ = '{sample}' + config["raw_file_reverse_extension"]

rule all:
    #These are all the output files that will be created (hopefully)
    #This list is organized in reverse order
    #The {sample} wild card is defined in this rule, if you wish to
    #run only a portion of the pipeline, you may need to redefine it using
    #the expand function
    input:
        #Step7 merge the OTU maps
        "Step7_otu_table/full_qiime_otu_map.txt",
        "Step7_otu_table/otu_table.biom",
        "Step7_otu_table/otu_table.txt",

        #Step6 OTU picking with SWARM
        "Step6_swarm_otus/swarm_seeds.fa",
        "Step6_swarm_otus/swarm_otus.txt",
        "Step6_swarm_otus/swarm_qiime_otumap.txt",

        #Step5 find chimeras & Filter
        "Step5_chimera_filtering/derep_map_nodenovo.txt",
        "Step5_chimera_filtering/derep_map_nochimera.txt",
        "Step5_chimera_filtering/ref_chimera_list.txt",
        "Step5_chimera_filtering/denovo_chimera_list.txt",
        "Step5_chimera_filtering/derep_seqs_denovo_chimeras.fasta",
        "Step5_chimera_filtering/derep_seqs_nodenovochim.fasta",
        "Step5_chimera_filtering/derep_seqs_ref_chimeras.fasta",
        "Step5_chimera_filtering/derep_seqs_nochim.fasta",

        #Step4 dereplicate sequences & make qiime otu map
        "Step4_dereplicate/derep_map.txt",
        "Step4_dereplicate/derep_seqs.fasta",

        #Step3 Mothur trim primers & concatentate
        "Step3_adapttrim/concat_seqs.fasta",
        expand("Step3_adapttrim/{sample}.trim.fasta", sample = sample_names),

        #Step2 vsearch quality control & assurance
        expand("Step2_vsearch_qaqc/{sample}.assembled.qaqc.fasta", sample = sample_names),

        #Step1 merging reads with pair
        expand("Step1_pear_merged/{sample}.assembled.fastq", sample = sample_names)

#STEP1 Merge mate pairs
rule merge_pairs:
    input:
        r1 = os.path.join(INPUT_DIR, "{sample}_R1_001.fastq.gz"),
        r2 = os.path.join(INPUT_DIR, "{sample}_R2_001.fastq.gz")
    params:
        sample_base = "{sample}"
    output:
        assembled = "Step1_pear_merged/{sample}.assembled.fastq",
        discarded = "Step1_pear_merged/{sample}.discarded.fastq",
        unass_r1 = "Step1_pear_merged/{sample}.unassembled.forward.fastq",
        uass_r2 = "Step1_pear_merged/{sample}.unassembled.reverse.fastq"
    log: "logs/Step1_pear/{sample}_pear_log.txt"
    shell:
        """
        pear -f {input.r1} \
        -r {input.r2} -o Step1_pear_merged/{params.sample_base} >> {log}
        """

#Step2 Quality Control
rule vsearch_qaqc:
    input:
        seqs = "Step1_pear_merged/{sample}.assembled.fastq"
    output:
        "Step2_vsearch_qaqc/{sample}.assembled.qaqc.fasta"
    log: "logs/Step2_vsearch_qaqc/{sample}_vsearch_log.txt"
    params:
        eemax = 1,
        minlen = 290,
        maxlen = 300,
        relabel = "{sample}_"
    shell:
        """
        vsearch --fastq_filter {input.seqs} --fastq_maxee {params.eemax} \
        --fastq_minlen {params.minlen} --fastq_maxlen {params.maxlen} \
        --relabel {params.relabel} --fastaout {output} --log {log} --quiet
        """

rule trim_primers:
    input:
        seqs = "Step2_vsearch_qaqc/{sample}.assembled.qaqc.fasta"
    output:
        "Step3_adapttrim/{sample}.trim.fasta"
    log: "logs/Step3_mothur_trim/{sample}_mothur_log.txt"
    params:
        oligo = config["oligos_file"],
        sample = "{sample}"
    threads: config["mothur_threads"]
    shell:
        """
        mothur "#set.logfile(name = {log}); trim.seqs(fasta={input.seqs}, oligos={params.oligo}, processors={threads})";
        mv Step2_vsearch_qaqc/{params.sample}.assembled.qaqc.trim.fasta Step3_adapttrim/{params.sample}.trim.fasta;
        mv Step2_vsearch_qaqc/{params.sample}.assembled.qaqc.scrap.fasta Step3_adapttrim/{params.sample}.scrap.fasta;

        echo "mv Step2_vsearch_qaqc/{params.sample}.assembled.qaqc.trim.fasta Step3_adapttrim/{params.sample}.trim.fasta" >> {log};
        echo "mv Step2_vsearch_qaqc/{params.sample}.assembled.qaqc.scrap.fasta Step3_adapttrim/{params.sample}.scrap.fasta" >> {log};
        """

rule concat_seqs:
    input:
        #"Step3_adapttrim/{sample}.trim.fasta"
        expand("Step3_adapttrim/{sample}.trim.fasta", sample = sample_names)
    output:
        "Step3_adapttrim/concat_seqs.fasta"
    log: "logs/Step3_mothur_trim/concate_seqs.txt"
    shell:
        """
        cat {input} >> {output}
        echo "cat {input} >> {output}" >> {log}
        """

rule dereplicate_vsearch:
    input:
        fasta = "Step3_adapttrim/concat_seqs.fasta",
    output:
        fasta = "Step4_dereplicate/derep_seqs.fasta",
        uc = "Step4_dereplicate/derep_map.uc"
    log: "logs/Step4_dereplicate_vsearch/vsearch_log.txt"
    params:
        relabel = config["derep_otu_name"]
    shell:
        """
        vsearch --derep_fulllength {input.fasta} --sizeout \
        --fasta_width 0 --relabel {params.relabel} --relabel_keep --output {output.fasta} \
        --uc {output.uc} --log {log} --quiet
        """

rule convert_uc_to_map:
    input:
        "Step4_dereplicate/derep_map.uc"
    output:
        "Step4_dereplicate/derep_map.txt"
    params:
        derep_otu_label = config["derep_otu_name"]
    script:
        "scripts/convert_uc2map_overholt.py"

rule find_chimera_vsearch_denovo:
    input:
        "Step4_dereplicate/derep_seqs.fasta"
    output:
        chimeras = "Step5_chimera_filtering/derep_seqs_denovo_chimeras.fasta",
        nonchimeras = "Step5_chimera_filtering/derep_seqs_nodenovochim.fasta"
    log: "logs/Step5_chimera/vsearch_denovo_log.txt"
    shell:
        """
        vsearch --uchime_denovo {input} --chimeras {output.chimeras} \
        --nonchimeras {output.nonchimeras} --log {log} --quiet
        """

rule find_chimera_vsearch_reference:
    input:
        "Step5_chimera_filtering/derep_seqs_nodenovochim.fasta"
    output:
        chimeras = "Step5_chimera_filtering/derep_seqs_ref_chimeras.fasta",
        nonchimeras = "Step5_chimera_filtering/derep_seqs_nochim.fasta"
    params:
        db = config["chimera_db"]
    threads:
        config["uchime_ref_threads"]
    log: "logs/Step5_chimera/vsearch_ref_log.txt"
    shell:
        """
        vsearch --uchime_ref {input} --chimeras {output.chimeras} \
        --nonchimeras {output.nonchimeras} --db {params.db} --threads {threads} \
        --log {log} --quiet
        """

rule remove_all_chimeras:
    input:
        denovo = "Step5_chimera_filtering/derep_seqs_denovo_chimeras.fasta",
        ref = "Step5_chimera_filtering/derep_seqs_ref_chimeras.fasta"
    output:
        ref_list = "Step5_chimera_filtering/ref_chimera_list.txt",
        denovo_list = "Step5_chimera_filtering/denovo_chimera_list.txt",
        map_nodenovo = "Step5_chimera_filtering/derep_map_nodenovo.txt",
        map_nochim = "Step5_chimera_filtering/derep_map_nochimera.txt"
    params:
        otu_map = "Step4_dereplicate/derep_map.txt",
    log: "logs/Step5_chimera/removing_chimera_log.txt"
    shell:
        """
        #Grab the OTU names that are chimera & print to a list
        perl -ne 'if ($_ =~ m/^>/) {{print $_;}}' {input.denovo} > {output.denovo_list};
        perl -ne 'if ($_ =~ m/^>/) {{print $_;}}' {input.ref} > {output.ref_list};
        #Remove the sizing info from the name, to allow matching
        perl -i -pe 's/>(D.*_OTU[0-9]+);.*/$1/g' {output.denovo_list};
        perl -i -pe 's/>(D.*_OTU[0-9]+);.*/$1/g' {output.ref_list};
        #Remove chimeras from the OTU map file
        awk 'NR==FNR{{a[$1];next}};!($1 in a)' {output.denovo_list} {params.otu_map} > {output.map_nodenovo};
        awk 'NR==FNR{{a[$1];next}};!($1 in a)' {output.ref_list} {output.map_nodenovo} > {output.map_nochim};

        echo "Grabbing Chimera sequence names, removing sizing info to allow names to match, & remove chimeras from the map file\n" >> {log};
        echo "perl -ne \'if (\$_ =~ m/^>/) {{print \$_;}}\' {input.denovo} > {output.denovo_list}" >> {log};
        echo "perl -ne \'if (\$_ =~ m/^>/) {{print \$_;}}\' {input.ref} > {output.ref_list}" >> {log};
        echo "perl -i -pe \'s/>(D.*_OTU[0-9]+);.*/\$1/g\' {output.denovo_list}" >> {log};
        echo "perl -i -pe \'s/>(D.*_OTU[0-9]+);.*/\$1/g\' {output.ref_list}" >> {log};
        echo "awk \'NR==FNR{{a[\$1];next}};!(\$1 in a)\' {output.denovo_list} {params.otu_map} > {output.map_nodenovo}" >> {log};
        echo "awk \'NR==FNR{{a[\$1];next}};!(\$1 in a)\' {output.ref_list} {output.map_nodenovo} > {output.map_nochim}" >> {log};
        """


rule pick_otus_with_swarm:
    input:
        "Step5_chimera_filtering/derep_seqs_nochim.fasta"
    output:
        seeds = "Step6_swarm_otus/swarm_seeds.fa",
        seeds_rename = "Step6_swarm_otus/swarm_seeds_rename.fa",
        seeds_rename_nosize = "Step6_swarm_otus/swarm_seeds_rename_nosize.fa",
        otus = "Step6_swarm_otus/swarm_otus.txt",
        otu_map = "Step6_swarm_otus/swarm_qiime_otumap.txt"
    log: "logs/Step6_swarm/swarm_log.txt"
    threads: config["swarm_threads"]
    shell:
        """
        swarm -f -t {threads} -z -w {output.seeds} -o {output.otus} {input} -l {log};
        perl -ne 'BEGIN {{$count = 0}}; if ($_ =~ m/>/) {{print ">denovo".$count."\t".$_; $count++;}} else {{print $_}}' {output.seeds} > {output.seeds_rename};
        sed -E 's/\t.*//g' {output.seeds_rename} > {output.seeds_rename_nosize};
        perl -ne 'BEGIN {{$count = 0}}; print "denovo".$count."\t".$_; $count++;' {output.otus} > {output.otu_map};
        sed -i -E 's/;size=[0-9]+;//g' {output.otu_map};

        echo "\n\n\nReformating SWARM output\n" >> {log};
        echo "perl -ne BEGIN {{\$count = 0}}; if (\$_ =~ m/>/) {{print '>denovo'.\$count.'\t'.\$_; \$count++;}} else {{print \$_}} {output.seeds} > {output.seeds_rename}" >> {log};
        echo "sed -E \'s/\t.*//g\' {output.seeds_rename} > {output.seeds_rename_nosize}" >> {log};
        echo "perl -ne BEGIN {{\$count = 0}}; {{print 'denovo'.\$count.'\t'.\$_; \$count++;}}\' {output.otus} > {output.otu_map}" >> {log};
        echo "sed -i -E \'s/;size=[0-9]+;//g\' {output.otu_map}" >> {log};
        """

rule merge_map_files:
    input:
        derep_map = "Step4_dereplicate/derep_map.txt",
        otu_map = "Step6_swarm_otus/swarm_qiime_otumap.txt"
    output:
        "Step7_otu_table/full_qiime_otu_map.txt"
    script:
        "scripts/merge_otu_maps_overholt.py"

rule make_otu_table:
    input:
        full_map = "Step7_otu_table/full_qiime_otu_map.txt"
    output:
        biom = "Step7_otu_table/otu_table.biom",
        table = "Step7_otu_table/otu_table.txt"
    script:
        "scripts/make_biom_table_overholt.py"
