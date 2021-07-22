"""``snakemake`` file that runs entire analysis."""

# Imports ---------------------------------------------------------------------
import glob
import itertools
import os.path
import textwrap
import urllib.request

import Bio.SeqIO

import dms_variants.codonvarianttable
import dms_variants.illuminabarcodeparser

import pandas as pd

# Configuration  --------------------------------------------------------------
configfile: 'config.yaml'

# run "quick" rules locally:
localrules: make_rulegraph,
            make_summary

# Functions -------------------------------------------------------------------
def nb_markdown(nb):
    """Return path to Markdown results of notebook `nb`."""
    return os.path.join(config['summary_dir'],
                        os.path.basename(os.path.splitext(nb)[0]) + '.md')

# Global variables extracted from config --------------------------------------
pacbio_runs = (pd.read_csv(config['pacbio_runs'], dtype = str)
               .assign(pacbioRun=lambda x: x['library'] + '_' + x['run'])
               )
assert len(pacbio_runs['pacbioRun'].unique()) == len(pacbio_runs['pacbioRun'])

# Information on samples and barcode runs -------------------------------------
barcode_runs = pd.read_csv(config['barcode_runs'])

# combination of the *library* and *sample* columns should be unique.
assert len(barcode_runs.groupby(['library', 'sample'])) == len(barcode_runs)

# *sample* should be the hyphen separated concatenation of
# *experiment*, *antibody*, *concentration*, and *sort_bin*.
sample_vs_expect = (
    barcode_runs
    .assign(expect=lambda x: x[['experiment', 'antibody', 'concentration',
                                'sort_bin']]
                             .apply(lambda r: '-'.join(r.values.astype(str)),
                                    axis=1),
            equal=lambda x: x['sample'] == x['expect'],
            )
    )
assert sample_vs_expect['equal'].all(), sample_vs_expect.query('equal != True')

# barcode runs with R1 files expanded by glob
barcode_runs_expandR1 = (
    barcode_runs
    .assign(R1=lambda x: x['R1'].str.split('; ').map(
                    lambda y: list(itertools.chain(*map(glob.glob, y)))),
            n_R1=lambda x: x['R1'].map(len),
            sample_lib=lambda x: x['sample'] + '_' + x['library'],
            )
    )

assert barcode_runs_expandR1['sample_lib'].nunique() == len(barcode_runs_expandR1)
if any(barcode_runs_expandR1['n_R1'] < 1):
    raise ValueError(f"no R1 for {barcode_runs_expandR1.query('n_R1 < 1')}")

# Rules -----------------------------------------------------------------------

# this is the target rule (in place of `all`) since it first rule listed
rule make_summary:
    """Create Markdown summary of analysis."""
    input:
        rulegraph=os.path.join(config['summary_dir'], 'rulegraph.svg'),
        get_mut_bind_expr=config['mut_bind_expr'],
        process_ccs=nb_markdown('process_ccs.ipynb'),
        build_variants=nb_markdown('build_variants.ipynb'),
        codon_variant_table=config['codon_variant_table'],
        aggregate_variant_counts=nb_markdown('aggregate_variant_counts.ipynb'),
        variant_counts=config['variant_counts'],
        counts_to_cells_ratio=nb_markdown('counts_to_cells_ratio.ipynb'),
        counts_to_cells_csv=config['counts_to_cells_csv'],
        # compute_meanF='results/summary/compute_expression_meanF.md',
        # expression_sortseq_file=config['expression_sortseq_file'],
        # global_epistasis_expression=nb_markdown('global_epistasis_expression.ipynb'),
        counts_to_scores=nb_markdown('counts_to_scores.ipynb'),
        scores_to_frac_escape=nb_markdown('scores_to_frac_escape.ipynb'),
        escape_fracs=config['escape_fracs'],
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

            Here is the rule graph of the computational workflow:
            ![{path(input.rulegraph)}]({path(input.rulegraph)})

            Here is the Markdown output of each notebook in the workflow:
            1. Get prior DMS mutation-level [binding and expression data]({path(input.get_mut_bind_expr)}).

            2. [Process PacBio CCSs]({path(input.process_ccs)}).

            3. [Build variants from CCSs]({path(input.build_variants)}).
               Creates a [codon variant table]({path(input.codon_variant_table)})
               linking barcodes to the mutations in the variants.

            4. Count variants and then
               [aggregate counts]({path(input.aggregate_variant_counts)}) to create
               to create [variant counts file]({path(input.variant_counts)}).
            
            5. [Analyze sequencing counts to cells ratio]({path(input.counts_to_cells_ratio)});
               this prints a list of any samples where this ratio too low. Also
               creates [a CSV]({path(input.counts_to_cells_csv)}) with the
               sequencing counts, number of sorted cells, and ratios for
               all samples.
            
            6. [Escape scores from variant counts]({path(input.counts_to_scores)}).
            
            7. [Global epistasis deconvolution of escape fractions for single mutations]
               ({path(input.scores_to_frac_escape)}); creating
               [mutation escape fraction file]({path(input.escape_fracs)}).


            """
            ).strip())


rule make_rulegraph:
    # error message, but works: https://github.com/sequana/sequana/issues/115
    input:
        workflow.snakefile
    output:
        os.path.join(config['summary_dir'], 'rulegraph.svg')
    shell:
        "snakemake --forceall --rulegraph | dot -Tsvg > {output}"

rule scores_to_frac_escape:
    """Estimate mutation-level escape scores."""
    input:
        escape_score_samples=config['escape_score_samples'],
        escape_scores=config['escape_scores'],
    output:
        nb_markdown=nb_markdown('scores_to_frac_escape.ipynb'),
        escape_fracs=config['escape_fracs'],
    params:
        nb='scores_to_frac_escape.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"

rule counts_to_scores:
    """Analyze variant counts to compute escape scores."""
    input:
        config['variant_counts'],
        config['wildtype_sequence'],
        # config['mut_bind_expr'],
        # config['variant_expr'],
        # config['variant_bind'],
    output:
        nb_markdown=nb_markdown('counts_to_scores.ipynb'),
        escape_scores=config['escape_scores'],
        escape_score_samples=config['escape_score_samples'],
    params:
        nb='counts_to_scores.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"

rule counts_to_cells_ratio:
    input:
        config['variant_counts'],
        config['barcode_runs'],
        config['wildtype_sequence'],
    output:
        nb_markdown=nb_markdown('counts_to_cells_ratio.ipynb'),
        counts_to_cells_csv=config['counts_to_cells_csv'],
    params:
        nb='counts_to_cells_ratio.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"

rule aggregate_variant_counts:
    input:
        counts=expand(os.path.join(config['counts_dir'],
                                   "{sample_lib}_counts.csv"),
                      sample_lib=barcode_runs_expandR1['sample_lib']),
        fates=expand(os.path.join(config['counts_dir'],
                                  "{sample_lib}_fates.csv"),
                     sample_lib=barcode_runs_expandR1['sample_lib']),
        variant_table=config['codon_variant_table'],
        wt_seq=config['wildtype_sequence'],
        barcode_runs=config['barcode_runs'],
    output:
        config['variant_counts'],
        nb_markdown=nb_markdown('aggregate_variant_counts.ipynb')
    params:
        nb='aggregate_variant_counts.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"

rule count_variants:
    """Count variants for a specific sample."""
    input:
        variant_table=config['codon_variant_table'],
        wt_seq=config['wildtype_sequence'],
        r1s=lambda wildcards: (barcode_runs_expandR1
                               .set_index('sample_lib')
                               .at[wildcards.sample_lib, 'R1']
                               ),
    output:
        counts=os.path.join(config['counts_dir'], "{sample_lib}_counts.csv"),
        fates=os.path.join(config['counts_dir'], "{sample_lib}_fates.csv"),
    params:
        sample_lib="{sample_lib}"
    run:
        # parse sample and library from `sample_lib` wildcard
        lib = params.sample_lib.split('_')[-1]
        sample = params.sample_lib[: -len(lib) - 1]
        assert sample == (barcode_runs_expandR1
                          .set_index('sample_lib')
                          .at[params.sample_lib, 'sample']
                          )
        assert lib == (barcode_runs_expandR1
                       .set_index('sample_lib')
                       .at[params.sample_lib, 'library']
                       )
        # initialize `CodonVariantTable` (used to get valid barcodes)
        wt_seqrecord = Bio.SeqIO.read(input.wt_seq, 'fasta')
        geneseq = str(wt_seqrecord.seq)
        primary_target = wt_seqrecord.name
        variants=dms_variants.codonvarianttable.CodonVariantTable(
                    geneseq=geneseq,
                    barcode_variant_file=input.variant_table,
                    substitutions_are_codon=True,
                    substitutions_col='codon_substitutions',
                    primary_target=primary_target)
        # initialize `IlluminaBarcodeParser`
        parser = dms_variants.illuminabarcodeparser.IlluminaBarcodeParser(
                    valid_barcodes=variants.valid_barcodes(lib),
                    **config['illumina_barcode_parser_params'])
        # parse barcodes
        counts, fates = parser.parse(input.r1s,
                                     add_cols={'library': lib,
                                               'sample': sample})
        # write files
        counts.to_csv(output.counts, index=False)
        fates.to_csv(output.fates, index=False)

rule build_variants:
    """Build variant table from processed CCSs."""
    input:
        config['processed_ccs_file']
    output:
        config['codon_variant_table'],
        nb_markdown=nb_markdown('build_variants.ipynb')
    params:
        nb='build_variants.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"

rule get_mut_bind_expr:
    """Download SARS-CoV-2 mutation ACE2-binding and expression from URL."""
    output:
        file=config['mut_bind_expr']
    run:
        urllib.request.urlretrieve(config['mut_bind_expr_url'], output.file)

rule process_ccs:
    """Process the PacBio CCSs."""
    input:
        expand(os.path.join(config['ccs_dir'], "{pacbioRun}_ccs.fastq.gz"),
               pacbioRun=pacbio_runs['pacbioRun']),
        config['amplicons'],
    output:
        config['processed_ccs_file'],
        nb_markdown=nb_markdown('process_ccs.ipynb')
    params:
        nb='process_ccs.ipynb'
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
