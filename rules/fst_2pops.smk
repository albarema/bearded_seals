#conda activate marine, check if this is needed:source ~/.bash_profile
# It probably needs to be run in steps 

#module load snakemake/7.20.0
configfile: 'config.yaml'
import pandas as pd
import itertools

def get_Ninds(infile):
    df = pd.read_table(infile, header=None)
    Nhalf=len(df.index)
    return int(Nhalf/2)

SUPERPOPS=['Pacific', 'Atlantic']
OUTDIR=config['outSel']
MAPQ=config['angsd']['mapQ']

COMSUPER = ['Pacific-Atlantic']

pbsAll = ["PBS_"+ s for s in SUPERPOPS ]
wFstAll =  ["wFst_" + s[0] + "_" + s[1] for s in list(itertools.combinations(SUPERPOPS, 2))]

rule all:
    input:
        expand(os.path.join(OUTDIR,"saf-100kb/{postQC}/{pairSamples}.slidingwindow"),pairSamples=COMSUPER,postQC=['postQC']),# 'postQC-p1'

rule saf_angsd:
    input:
        bamlist=config['bamlist'],
        sites=expand(os.path.join(config['outputDir'], "final_sites/{prefix}.mapQ{mapQ}.{postQC}.loci.txt"),
            prefix='EB', mapQ=MAPQ, allow_missing=True),
        chrom = config['chrs'],
        ref=config['ref'],
    threads: 8
    params: 
        oname=os.path.join(OUTDIR,"saf/{postQC}/{superpop}"),
        gl=config['gl']['gl'],
        nInds = lambda wildcards: get_Ninds(config['bamlist'].format(superpop=wildcards.superpop))
    output:
        os.path.join(OUTDIR,"saf/{postQC}/{superpop}.saf.idx")
    shell:
        "angsd -bam {input.bamlist} "
        "-gl {params.gl} " # 1 vs 2 
        "-rf {input.chrom} " 
        "-sites {input.sites} "
        "-uniqueOnly "
        "-minMapQ 30 "
        "-minQ 30 "
        "-dosaf 1 "
        "-minInd {params.nInds} "
        "-anc {input.ref} "
        "-nThreads {threads} "
        "-out {params.oname}"


rule dsfs2:
    input:
        expand(os.path.join(OUTDIR,"saf/{postQC}/{superpop}.saf.idx"), superpop=SUPERPOPS, allow_missing=True)
    output:
        sfs2=os.path.join(OUTDIR,"saf/{postQC}/tmp/{pop1}-{pop2}.sfs")
    params: 
        pfiles=os.path.join(OUTDIR, "saf/{postQC}")
    threads: 24
    shell:
        "winsfs {params.pfiles}/{wildcards.pop1}.saf.idx " 
        "{params.pfiles}/{wildcards.pop2}.saf.idx "
        "-t {threads} "
        "> {params.pfiles}/tmp/{wildcards.pop1}-{wildcards.pop2}.sfs"

rule rm_1line:
    input:
        fils=os.path.join(OUTDIR,"saf/{postQC}/tmp/{pairSamples}.sfs")
    output:
        sfs2=os.path.join(OUTDIR,"saf/{postQC}/{pairSamples}.sfs")
    shell:
        "tail -n+2 {input} > {output} " 

rule fst:
    input:
        safs=expand(os.path.join(OUTDIR,"saf/{postQC}/{superpop}.saf.idx"), superpop=SUPERPOPS, allow_missing=True),
        ml2_call = os.path.join(OUTDIR,"saf/{postQC}/{pairSamples}.sfs"),
    output:
        sfs2=os.path.join(OUTDIR,"saf/{postQC}/{pairSamples}.fst.gz"),
        idx=os.path.join(OUTDIR,"saf/{postQC}/{pairSamples}.fst.idx")
    params:
        pfiles=os.path.join(OUTDIR, "saf/{postQC}"),
    threads: 24
    shell:
        "realSFS fst index {input.safs} "
        "-sfs {input.ml2_call} "
        "-fstout {params.pfiles}/{wildcards.pairSamples}"

rule pbs_w:
    input:
        fst=os.path.join(OUTDIR,"saf/{postQC}/{pairSamples}.fst.idx")
    output:
        stat=os.path.join(OUTDIR,"saf-100kb/{postQC}/{pairSamples}.{window}.{step}.slidingwindow")
    params:
        pfiles=os.path.join(OUTDIR, "saf-100kb/{postQC}"),
        w = lambda wildcards: wildcards.window *1e3,
        s = lambda wildcards: wildcards.step *1e3,
    shell:
        "realSFS fst stats2 {input.fst} "
        "-win {params.w} "
        "-step {params.s} "
        "> {output.stat}"

## todo
rule default_pbs_plot:
    """
    Figure and top regions for each of the pops
    """
    input:
        sel=os.path.join(OUTDIR,"saf-100kb/{postQC}/{pairSamples}.slidingwindow")
    output:
        png=os.path.join(config['outSel'],"saf-100kb/{postQC}/plots/Fst_PBS.{pairSamples}.png"),
    params:
        outdir=os.path.join(config['outSel'],"saf-100kb/{postQC}/")
    shell:
        "Rscript scripts/pbs_plot.r {input.sel} {wildcards.pairSamples} {params.outdir}"

rule topRegions_pbs:
    """
    Figure and top regions for each of the pops
    """
    input:
        sel=expand(os.path.join(OUTDIR,"saf-100kb/{postQC}/{allPairs}.slidingwindow"), allPairs=COMSUPER, allow_missing=True)
    output:
        top=os.path.join(config['outSel'],"saf-100kb/{postQC}/topRegions/{fullpopi}.{window}kb{pc}.tsv"),
    params:
        outdir=os.path.join(config['outSel'],"saf-100kb/{postQC}/"),
        popsAll = COMSUPER
    shell:
        "Rscript scripts/get_topRegions.r {input.sel} {params.popsAll} {wildcards.window} {wildcards.pc} {params.outdir}"


rule man_plot: 
    """
    Figure and top regions for each of the pops
    """
    input:
        sel=expand(os.path.join(OUTDIR,"saf-100kb/{postQC}/{allPairs}.slidingwindow"), allPairs=COMSUPER, allow_missing=True),
        top=os.path.join(config['outSel'],"saf-100kb/{postQC}/topRegions/{fullpopi}.{window}kb{pc}.tsv"),
    output:
        png=os.path.join(config['outSel'],"saf-100kb/{postQC}/plots/man.{fullpopi}.{window}kb{pc}.png"),
    params:
        outdir=os.path.join(config['outSel'],"saf-100kb/{postQC}/"),
        allp = COMSUPER
    shell:
        "Rscript scripts/pbs_manPlot.r {input.sel} {params.allp} {wildcards.fullpopi} {wildcards.window}kb {wildcards.pc} {params.outdir}"

# run topregions in my computer since the server wasn't working well. Update scripts.
rule annot_all:
    """
    Figure and top regions for each of the pops
    """
    input:
        top=os.path.join(config['outSel'],"saf-100kb/{postQC}/topRegions/{fullpopi}.{window}kb{pc}.tsv"),
        annot = "/projects/mjolnir1/people/nrb613/X204SC22021885-Z01-F004/reference.genome/autosomal.gene.annotations.txt"
    output:
        annot = os.path.join(config['outSel'],"saf-100kb/{postQC}/annRegions/{annsource}/{fullpopi}.{window}kb{pc}.genes.tsv"),
    params:
        patop=os.path.join(config['outSel'],"saf-100kb/{postQC}/topRegions/"),
        outdir = os.path.join(config['outSel'],"saf-100kb/{postQC}/annRegions/{annsource}/"),
    shell:
        "Rscript scripts/annot_dnazoo_regions.r "
        "-m {wildcards.fullpopi} "
        "-w {wildcards.window}kb "
        "-c {wildcards.pc} "
        "-p {params.patop} "
        "-a {input.annot} "
        "-o {params.outdir}"

