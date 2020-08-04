IDS, = glob_wildcards("input/{id,\d+}.bam")

rule all:
     input:
        HC = expand("output/service-data/{id}_HC",id=IDS),
        FB = expand("output/service-data/{id}_FB",id=IDS),
        SM = expand("output/service-data/{id}_SM",id=IDS),
        FB_HC = expand("output/service-data/{id}_f_h",id=IDS),
        FB_SM = expand("output/service-data/{id}_f_s",id=IDS),
        HC_SM = expand("output/service-data/{id}_h_s",id=IDS)
     output:expand("output/reports/{id}",id=IDS)
     run:
         def get_Num(x): 
             with open(x,"r") as rep:
                 vc = " ".join(rep.readlines()).split()
             for i in range(len(vc)):
                 if vc[i] == "Sites":
                     return(str(vc[i-1]))

         for j in range(len(output)):
             fb = get_Num(input.FB[j])
             hc = get_Num(input.HC[j])
             sm = get_Num(input.SM[j])
             fb_hc = get_Num(input.FB_HC[j])
             fb_sm = get_Num(input.FB_SM[j])
             hc_sm = get_Num(input.HC_SM[j])
             with open(output[j],"w") as vcf:
                 vcf.write(
                       str(fb) + "\t" + str(fb_hc) + "\t" + str(fb_sm) +"\n"+
                       str(hc) + "\tNone\t" + str(hc_sm) + "\n"+
                       str(sm) + "\tNone\tNone\n")

rule count_intersected_variants:
    input:
       vcf_f_h = expand("output/service-data/{id}_f_h.vcf.gz",id=IDS),
       vcf_f_s = expand("output/service-data/{id}_f_s.vcf.gz",id=IDS),
       vcf_h_s = expand("output/service-data/{id}_h_s.vcf.gz",id=IDS)
    output:
       r_f_h = expand("output/service-data/{id}_f_h",id=IDS),
       r_f_s = expand("output/service-data/{id}_f_s",id=IDS),
       r_h_s = expand("output/service-data/{id}_h_s",id=IDS)
    run:
       in_ = [input.vcf_f_h,input.vcf_f_s,input.vcf_h_s]
       out_ = [output.r_f_h,output.r_f_s,output.r_h_s]
       for i in range(len(in_)):
           for j in range(len(in_[i])):
               shell("vcftools --gzvcf {input_vcf} 2> {output_}".format(input_vcf=in_[i][j], output_=out_[i][j]))


rule intersect_vcf:
    input:
        vcf1=[expand("output/service-data/{id}_FB.vcf.gz",id=IDS),expand("output/service-data/{id}_FB.vcf.gz",id=IDS),expand("output/service-data/{id}_HC.vcf.gz",id=IDS)],
        vcf2=[expand("output/service-data/{id}_HC.vcf.gz",id=IDS),expand("output/service-data/{id}_SM.vcf.gz",id=IDS),expand("output/service-data/{id}_SM.vcf.gz",id=IDS)]
    output:[expand("output/service-data/{id}_f_h.vcf.gz",id=IDS),expand("output/service-data/{id}_f_s.vcf.gz",id=IDS),expand("output/service-data/{id}_h_s.vcf.gz",id=IDS)]
    run: 
        for i in range(len(input.vcf1)):
            shell("vcf-isec -f -n +2 {input_vcf1} {input_vcf2} | bgzip -c > {output_}".format(
input_vcf1=input.vcf1[i], input_vcf2=input.vcf2[i], output_=output[i]))

rule index_vcf:
     input:
        vcf=[expand("output/service-data/{id}_SM.vcf",id=IDS),expand("output/service-data/{id}_FB.vcf",id=IDS),expand("output/service-data/{id}_HC.vcf",id=IDS)]
     output:[expand("output/service-data/{id}_SM.vcf.gz",id=IDS),expand("output/service-data/{id}_FB.vcf.gz",id=IDS),expand("output/service-data/{id}_HC.vcf.gz",id=IDS)]
     run: 
        for i in range(len(input.vcf)):
            shell("bgzip {input_vcf} ; tabix -p vcf {input_vcf}.gz".format(
input_vcf=input.vcf[i]))

rule count_individual_variants:
    input:
       vcf_SM = expand("output/service-data/{id}_SM.vcf",id=IDS),
       vcf_FB = expand("output/service-data/{id}_FB.vcf",id=IDS),
       vcf_HC = expand("output/service-data/{id}_HC.vcf",id=IDS)
    output:
       r_SM = expand("output/service-data/{id}_SM",id=IDS),
       r_FB = expand("output/service-data/{id}_FB",id=IDS),
       r_HC = expand("output/service-data/{id}_HC",id=IDS)
    priority: 1
    run:
        in_ = [input.vcf_SM,input.vcf_FB,input.vcf_HC]
        out_ = [output.r_SM,output.r_FB,output.r_HC]
        for i in range(len(in_)):
            for j in range(len(in_[i])):
                shell("vcftools --vcf {input_vcf} 2> {output_}".format(input_vcf=in_[i][j], output_=out_[i][j]))

rule vc_Samtools:
    input:
        fa="data/22.fa",
        bam=expand("input/{id}.bam", id=IDS),
        r_i="data/22.fa.fai"
    output:expand("output/service-data/{id}_SM.vcf",id=IDS)
    run: 
        for i in range(len(input.bam)):
            shell("samtools mpileup -uf {fa_} {input_bam} | bcftools view -vcg - > {output_}".format(
fa_=input.fa,input_bam=input.bam[i],output_=output[i]))

rule vc_Freebayes:
    input:
        fa="data/22.fa",
        bam=expand("input/{id}.bam", id=IDS),
        r_i="data/22.fa.fai"
    output:expand("output/service-data/{id}_FB.vcf",id=IDS)
    run: 
        for i in range(len(input.bam)):
            shell("freebayes -f {fa_} {input_bam} > {output_}".format(
fa_=input.fa,input_bam=input.bam[i],output_=output[i]))

rule vc_HaplotypeCaller:
    input:
        fa="data/22.fa",
        bam=expand("input/{id}.bam", id=IDS),
        bai=expand("input/{id}.bam.bai", id=IDS),
        r_i="data/22.fa.fai",
        r_d="data/22.dict"
    output:expand("output/service-data/{id}_HC.vcf",id=IDS)
    run: 
        for i in range(len(input.bam)):
            shell("java -jar $GATK -R {fa_} -T HaplotypeCaller -I {input_bam} -o {output_}".format(
fa_=input.fa,input_bam=input.bam[i],output_=output[i]))

rule reference_index:
    input:"data/22.fa"
    output:"data/22.fa.fai"
    shell:"samtools faidx {input}"

rule reference_dict:
    input:"data/22.fa"
    output:"data/22.dict"
    shell:"java -jar $PICARD CreateSequenceDictionary R={input} O={output}"

rule index_bam:
    input:
        bam=expand("input/{id}.bam", id=IDS)
    output:expand("input/{id}.bam.bai",id=IDS)
    run: 
        for i in range(len(input.bam)):
            shell("samtools index {input_bam} {output_}".format(
input_bam=input.bam[i],output_=output[i]))
