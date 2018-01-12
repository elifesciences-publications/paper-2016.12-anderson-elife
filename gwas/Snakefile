# DGRP freeze 2 variant calls, lifted to dm6
rule get_dgrp_vcf:
    output: "source_data/dgrp2_dm6.vcf.gz"
    shell:
        "wget ftp://ftp.hgsc.bcm.edu/DGRP/freeze2_Feb_2013/liftover_data_for_D.mel6.0_from_William_Gilks_Oct_2015/DGRP_liftover_Oct2015/dgrp2_dm6.vcf.gz -O {output}"

rule sort_vcf:
    input: "source_data/{dataset}.vcf.gz"
    output: "{dataset}.sorted.vcf.gz"
    shell:
        "zcat {input} | vcf-sort -c -t ./tmp -p 24 | bgzip > {output}"

# Cross entire VCF and convert to plink format
rule cross_vcf:
    message: "dgrp_cross.py {input} -> {output} with {threads} threads."
    threads: 24
    input:   "dgrp2_dm6.sorted.vcf.gz"
    output:  "dgrp2_cross_ssdef.bed"
    params:  base="dgrp2_cross_ssdef"
    shell:
        "zcat {input} | "
        " parallel --jobs {threads} --header '(\#.*\n)*' --block 100m --keep-order --pipe ./dgrp_cross.py /dev/stdin /dev/stdout | "
        " ./firstheader | "
        " plink2 --vcf /dev/stdin --allow-extra-chr --double-id --out {params.base}"

rule vcf_to_plink:
    input: "{dataset}.vcf.gz"
    output: "{dataset}.bed"
    params: base="{dataset}"
    shell:
        "zcat {input} | plink2 --vcf /dev/stdin --allow-extra-chr --double-id --out {params.base}"

rule pca:
    input: "{dataset}.bed"
    output: "{dataset}.eigenval"
    params: base="{dataset}"
    shell:
        "plink2 --bfile {params.base} --allow-extra-chr --pca --out {params.base}"

rule cluster:
    input: "{dataset}.bed"
    output: "{dataset}.mds"
    params: base="{dataset}"
    shell:
        "plink2 --bfile {params.base} --allow-extra-chr --cluster --mds-plot 3 --out {params.base}"

rule extract_pheno:
    input: "source_data/032516RevisitLines.xlsx"
    output: "rh.pheno"
    shell:
        """xlsx2csv {input} -n Sheet4 | awk -F, 'NR > 1 && $2 {{ print "line_" $2 "\tline_" $2 "\t" $19 }}' > {output}"""

rule assoc:
    input: "{dataset}.bed", ph="{pheno}.pheno"
    output: "{dataset}_{pheno}.qassoc"
    params: base="{dataset}", pheno="{pheno}"
    shell:
        "plink2 --bfile {params.base} --allow-extra-chr --maf 0.05 --pheno {input.ph} --assoc --allow-no-sex --out {params.base}_{params.pheno}"

rule linear:
    input: "{dataset}.bed", ph="{pheno}.pheno", eigenvec="{dataset}.eigenvec"
    output: "{dataset}_{pheno}.assoc.linear"
    params: base="{dataset}", pheno="{pheno}"
    shell:
        "plink2 --bfile {params.base} --allow-extra-chr --maf 0.05 --pheno {input.ph} --linear mperm=100000 hide-covar --covar {input.eigenvec} --allow-no-sex --out {params.base}_{params.pheno}"
