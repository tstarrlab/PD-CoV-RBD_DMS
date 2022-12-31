"""``snakemake`` file that runs entire analysis."""

# Imports ---------------------------------------------------------------------
import glob
import itertools
import os.path
import os
import textwrap
import urllib.request

import Bio.SeqIO

import dms_variants.codonvarianttable
import dms_variants.illuminabarcodeparser

import pandas as pd

# Configuration  --------------------------------------------------------------
configfile: 'config.yaml'

# run "quick" rules locally:
localrules: make_dag,
            make_summary

# Functions -------------------------------------------------------------------
def nb_markdown(nb):
    """Return path to Markdown results of notebook `nb`."""
    return os.path.join(config['summary_dir'],
                        os.path.basename(os.path.splitext(nb)[0]) + '.md')

# Global variables extracted from config --------------------------------------
pacbio_runs = (pd.read_csv(config['pacbio_runs'], dtype = str)
               .assign(pacbioRun=lambda x: x['library'] + '_' + x['bg'] + '_' + x['run'])
               )
assert len(pacbio_runs['pacbioRun'].unique()) == len(pacbio_runs['pacbioRun'])

# Information on samples and barcode runs -------------------------------------
barcode_runs = pd.read_csv(config['barcode_runs'])

# Rules -----------------------------------------------------------------------

# making this summary is the target rule (in place of `all`) since it
# is first rule listed.
rule make_summary:
    """Create Markdown summary of analysis."""
    input:
        dag=os.path.join(config['summary_dir'], 'dag.svg'),
        process_ccs_PDCoV=nb_markdown('process_ccs_PDCoV.ipynb'),
        barcode_variant_table_PDCoV=config['codon_variant_table_file_PDCoV'],
        variant_counts_file=config['variant_counts_file'],
        count_variants=nb_markdown('count_variants.ipynb'),
        fit_gAPN_titrations='results/summary/compute_gAPN_Kd.md',
        gAPN_Kds_file=config['gAPN_Kds_file'],
        fit_hAPN_meanF='results/summary/compute_hAPN_meanF.md',
        hAPN_meanF_file=config['hAPN_meanF_file'],
        fit_pAPN_meanF='results/summary/compute_pAPN_meanF.md',
        pAPN_meanF_file=config['pAPN_meanF_file'],        
        collapse_scores='results/summary/collapse_scores.md',
        mut_phenos_file=config['final_variant_scores_mut_file'],
        heatmap_viz=os.path.join(config['visualization_dir'], "heatmap.html")
    output:
        summary = os.path.join(config['summary_dir'], 'summary.md')
    run:
        def path(f):
            """Get path relative to `summary_dir`."""
            return os.path.relpath(f, config['summary_dir'])
        with open(output.summary, 'w') as f:
            f.write(textwrap.dedent(f"""
            # Summary

            Analysis run by [Snakefile]({path(workflow.snakefile)})
            using [this config file]({path(workflow.configfiles[0])}).
            See the [README in the top directory]({path('README.md')})
            for details.

            Here is the DAG of the computational workflow:
            ![{path(input.dag)}]({path(input.dag)})

            Here is the Markdown output of each Jupyter notebook in the
            workflow:
            
            1. Process PacBio CCSs for [PDCoV libraries]({path(input.process_ccs_PDCoV)}). Creates barcode-variant lookup table, which can be found [here]({path(input.barcode_variant_table_PDCoV)}),.
            
            2. [Count variants by barcode]({path(input.count_variants)}).
               Creates a [variant counts file]({path(input.variant_counts_file)})
               giving counts of each barcoded variant in each condition.

            3. [Fit gAPN titration curves]({path(input.fit_gAPN_titrations)}) to calculate per-barcode K<sub>D</sub>, recorded in [this file]({path(input.gAPN_Kds_file)}).
            
            4. Analyze single-concentration sort-seq experients for [hAPN]({path(input.fit_hAPN_meanF)}) and [pAPN]({path(input.fit_pAPN_meanF)}) to calculate per-barcode binding MFI, recorded in these files for [hAPN]({path(input.hAPN_meanF_file)}) and [pAPN]({path(input.pAPN_meanF_file)}).
            
            5. [Derive final genotype-level phenotypes from replicate barcoded sequences]({path(input.collapse_scores)}).
               Generates final phenotypes, recorded in [this file]({path(input.mut_phenos_file)}).
            
            7. [Analyze patterns of epistasis in the DMS data and in SARS-CoV-2 genomic data]({path(input.epistatic_shifts)}).
            
            8. Make interactive data visualizations, available [here](https://jbloomlab.github.io/SARS-CoV-2-RBD_DMS_Omicron/)

            """
            ).strip())

rule make_dag:
    # error message, but works: https://github.com/sequana/sequana/issues/115
    input:
        workflow.snakefile
    output:
        os.path.join(config['summary_dir'], 'dag.svg')
    shell:
        "snakemake --forceall --dag | dot -Tsvg > {output}"


rule interactive_heatmap_plot:
    """ Make the interactive heatmap for expression and binding.
    """
    input: 
        scores=config['final_variant_scores_mut_file']
    params:
        annotations=config['RBD_sites']
    output:
        html=os.path.join(config['visualization_dir'], "heatmap.html")
    notebook: "RBD-Heatmaps-Interactive-Visualization.ipynb"


rule epistatic_shifts:
    input:
        config['final_variant_scores_mut_file'],
    output:
        config['JSD_v_WH1_file'],
        config['JSD_v_WH1_expr_file'],
        md='results/summary/epistatic_shifts.md',
        md_files=directory('results/summary/epistatic_shifts_files')
    envmodules:
        'R/3.6.2-foss-2019b'
    params:
        nb='epistatic_shifts.Rmd',
        md='epistatic_shifts.md',
        md_files='epistatic_shifts_files'
    shell:
        """
        R -e \"rmarkdown::render(input=\'{params.nb}\')\";
        mv {params.md} {output.md};
        mv {params.md_files} {output.md_files}
        """


rule collapse_scores:
    input:
        config['gAPN_Kds_file'],
        config['hAPN_meanF_file'],
        config['pAPN_meanF_file'],
    output:
        config['final_variant_scores_mut_file'],
        md='results/summary/collapse_scores.md',
        md_files=directory('results/summary/collapse_scores_files')
    envmodules:
        'R/3.6.2-foss-2019b'
    params:
        nb='collapse_scores.Rmd',
        md='collapse_scores.md',
        md_files='collapse_scores_files'
    shell:
        """
        R -e \"rmarkdown::render(input=\'{params.nb}\')\";
        mv {params.md} {output.md};
        mv {params.md_files} {output.md_files}
        """

rule fit_pAPN_MFI:
    input:
        config['codon_variant_table_file_PDCoV'],
        config['variant_counts_file']
    output:
        config['pAPN_meanF_file'],
        md='results/summary/compute_pAPN_meanF.md',
        md_files=directory('results/summary/compute_pAPN_meanF_files')
    envmodules:
        'R/3.6.2-foss-2019b'
    params:
        nb='compute_pAPN_meanF.Rmd',
        md='compute_pAPN_meanF.md',
        md_files='compute_pAPN_meanF_files'
    shell:
        """
        R -e \"rmarkdown::render(input=\'{params.nb}\')\";
        mv {params.md} {output.md};
        mv {params.md_files} {output.md_files}
        """

rule fit_hAPN_MFI:
    input:
        config['codon_variant_table_file_PDCoV'],
        config['variant_counts_file']
    output:
        config['hAPN_meanF_file'],
        md='results/summary/compute_hAPN_meanF.md',
        md_files=directory('results/summary/compute_hAPN_meanF_files')
    envmodules:
        'R/3.6.2-foss-2019b'
    params:
        nb='compute_hAPN_meanF.Rmd',
        md='compute_hAPN_meanF.md',
        md_files='compute_hAPN_meanF_files'
    shell:
        """
        R -e \"rmarkdown::render(input=\'{params.nb}\')\";
        mv {params.md} {output.md};
        mv {params.md_files} {output.md_files}
        """

rule fit_gAPN_titrations:
    input:
        config['codon_variant_table_file_PDCoV'],
        config['variant_counts_file']
    output:
        config['gAPN_Kds_file'],
        md='results/summary/compute_gAPN_Kd.md',
        md_files=directory('results/summary/compute_gAPN_Kd_files')
    envmodules:
        'R/3.6.2-foss-2019b'
    params:
        nb='compute_gAPN_Kd.Rmd',
        md='compute_gAPN_Kd.md',
        md_files='compute_gAPN_Kd_files'
    shell:
        """
        R -e \"rmarkdown::render(input=\'{params.nb}\')\";
        mv {params.md} {output.md};
        mv {params.md_files} {output.md_files}
        """

rule count_variants:
    """Count codon variants from Illumina barcode runs."""
    input:
        config['codon_variant_table_file_PDCoV'],
        config['barcode_runs']
    output:
        config['variant_counts_file'],
        nb_markdown=nb_markdown('count_variants.ipynb')
    params:
        nb='count_variants.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"
        
rule process_ccs_PDCoV:
    """Process the PacBio CCSs for PDCoV background and build variant table."""
    input:
        expand(os.path.join(config['ccs_dir'], "{pacbioRun}_ccs.fastq.gz"),
               pacbioRun=pacbio_runs['pacbioRun']),
    output:
        config['processed_ccs_file_PDCoV'],
    	config['codon_variant_table_file_PDCoV'],
        nb_markdown=nb_markdown('process_ccs_PDCoV.ipynb')
    params:
        nb='process_ccs_PDCoV.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"


if config['seqdata_source'] == 'HutchServer':

    rule get_ccs:
        """Symbolically link CCS files."""
        input:
            ccs_fastq=lambda wildcards: (pacbio_runs
                                        .set_index('pacbioRun')
                                        .at[wildcards.pacbioRun, 'ccs']
                                        )
        output:
            ccs_fastq=os.path.join(config['ccs_dir'], "{pacbioRun}_ccs.fastq.gz")
        run:
            os.symlink(input.ccs_fastq, output.ccs_fastq)

elif config['seqdata_source'] == 'SRA':
    raise RuntimeError('getting sequence data from SRA not yet implemented')

else:
    raise ValueError(f"invalid `seqdata_source` {config['seqdata_source']}")
