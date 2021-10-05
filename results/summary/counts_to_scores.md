# Analyze counts and compute escape scores
This Python Jupyter notebook analyzes the variant counts and looks at mutation coverage and jackpotting.
It then computes an "escape scores" for each variant after grouping by barcode or substitutions as specified in the configuration.

## Set up analysis

This notebook primarily makes use of the Bloom lab's [dms_variants](https://jbloomlab.github.io/dms_variants) package, and uses [plotnine](https://github.com/has2k1/plotnine) for ggplot2-like plotting syntax:


```python
import collections
import math
import os
import warnings

import Bio.SeqIO

import dms_variants.codonvarianttable
from dms_variants.constants import CBPALETTE
import dms_variants.plotnine_themes

from IPython.display import display, HTML

import matplotlib.pyplot as plt

import numpy

import pandas as pd

from plotnine import *

import seaborn

import yaml
```

Set [plotnine](https://github.com/has2k1/plotnine) theme to the gray-grid one defined in [dms_variants](https://jbloomlab.github.io/dms_variants):


```python
theme_set(dms_variants.plotnine_themes.theme_graygrid())
```

Versions of key software:


```python
print(f"Using dms_variants version {dms_variants.__version__}")
```

    Using dms_variants version 0.8.10


Ignore warnings that clutter output:


```python
warnings.simplefilter('ignore')
```

Read the configuration file:


```python
with open('config.yaml') as f:
    config = yaml.safe_load(f)
```

Create output directory:


```python
os.makedirs(config['escape_scores_dir'], exist_ok=True)
```

Read information about the samples:


```python
samples_df = pd.read_csv(config['barcode_runs'])
```

## Initialize codon-variant table
Initialize [CodonVariantTable](https://jbloomlab.github.io/dms_variants/dms_variants.codonvarianttable.html#dms_variants.codonvarianttable.CodonVariantTable) from wildtype gene sequence and the variant counts CSV file.
We will then use the plotting functions of this variant table to analyze the counts per sample:


```python
wt_seqrecord = Bio.SeqIO.read(config['wildtype_sequence'], 'fasta')
geneseq = str(wt_seqrecord.seq)
primary_target = wt_seqrecord.name
print(f"Read sequence of {len(geneseq)} nt for {primary_target} from {config['wildtype_sequence']}")
      
print(f"Initializing CodonVariantTable from gene sequence and {config['variant_counts']}")
      
variants = dms_variants.codonvarianttable.CodonVariantTable.from_variant_count_df(
                geneseq=geneseq,
                variant_count_df_file=config['variant_counts'],
                primary_target=primary_target)
      
print('Done initializing CodonVariantTable.')
```

    Read sequence of 603 nt for B1351 from data/wildtype_sequence.fasta
    Initializing CodonVariantTable from gene sequence and results/counts/variant_counts.csv.gz
    Done initializing CodonVariantTable.


## Sequencing counts per sample
Average counts per variant for each sample.
Note that these are the **sequencing** counts, in some cases they may outstrip the actual number of sorted cells:


```python
p = variants.plotAvgCountsPerVariant(libraries=variants.libraries,
                                     by_target=False,
                                     orientation='v')
p = p + theme(panel_grid_major_x=element_blank())  # no vertical grid lines
_ = p.draw()
```


    
![png](counts_to_scores_files/counts_to_scores_19_0.png)
    


And the numerical values plotted above:


```python
display(HTML(
 variants.avgCountsPerVariant(libraries=variants.libraries,
                               by_target=False)
 .pivot_table(index='sample',
              columns='library',
              values='avg_counts_per_variant')
 .round(1)
 .to_html()
 ))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>lib1</th>
      <th>lib2</th>
    </tr>
    <tr>
      <th>sample</th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>bg4e_6-9-none-0-ref</th>
      <td>495.2</td>
      <td>560.9</td>
    </tr>
    <tr>
      <th>bg4e_8-K115_old-80-abneg</th>
      <td>33.9</td>
      <td>33.8</td>
    </tr>
    <tr>
      <th>bg4e_10-K007-500-abneg</th>
      <td>77.0</td>
      <td>75.0</td>
    </tr>
    <tr>
      <th>bg4e_11-K031-500-abneg</th>
      <td>108.5</td>
      <td>97.5</td>
    </tr>
    <tr>
      <th>bg4e_12-K033-500-abneg</th>
      <td>118.8</td>
      <td>117.0</td>
    </tr>
    <tr>
      <th>bg4e_13-K040-500-abneg</th>
      <td>101.1</td>
      <td>112.5</td>
    </tr>
    <tr>
      <th>bg4e_14-K041-500-abneg</th>
      <td>100.3</td>
      <td>99.5</td>
    </tr>
    <tr>
      <th>bg4e_10-14-reference-0-ref</th>
      <td>1474.4</td>
      <td>1522.3</td>
    </tr>
    <tr>
      <th>bg4e_15-18-reference-0-ref</th>
      <td>2165.6</td>
      <td>2552.5</td>
    </tr>
    <tr>
      <th>bg4e_15-K046-200-abneg</th>
      <td>130.6</td>
      <td>116.3</td>
    </tr>
    <tr>
      <th>bg4e_16-K114-200-abneg</th>
      <td>150.7</td>
      <td>118.1</td>
    </tr>
    <tr>
      <th>bg4e_17-K115-80-abneg</th>
      <td>311.6</td>
      <td>426.6</td>
    </tr>
    <tr>
      <th>bg4e_18-K119-200-abneg</th>
      <td>129.5</td>
      <td>185.4</td>
    </tr>
  </tbody>
</table>


## Mutations per variant
Average number of mutations per gene among all variants of the primary target, separately for each date:


```python
#this plotting is very slow when lots of samples, so for now plots are commented out

for date, date_df in samples_df.groupby('date', sort=False):
   p = variants.plotNumCodonMutsByType(variant_type='all',
                                       orientation='v',
                                       libraries=variants.libraries,
                                       samples=date_df['sample'].unique().tolist(),
                                       widthscale=2)
   p = p + theme(panel_grid_major_x=element_blank())  # no vertical grid lines
   fig = p.draw()
   display(fig)
   plt.close(fig)
```


    
![png](counts_to_scores_files/counts_to_scores_23_0.png)
    



    
![png](counts_to_scores_files/counts_to_scores_23_1.png)
    



    
![png](counts_to_scores_files/counts_to_scores_23_2.png)
    


Now similar plots but showing mutation frequency across the gene:


```python
# this plotting is very slow when lots of samples, so for now code commented out

for date, date_df in samples_df.groupby('date', sort=False):
   p = variants.plotMutFreqs(variant_type='all',
                             mut_type='codon',
                             orientation='v',
                             libraries=variants.libraries,
                             samples=date_df['sample'].unique().tolist(),
                             widthscale=1.5)
   fig = p.draw()
   display(fig)
   plt.close(fig)
```


    
![png](counts_to_scores_files/counts_to_scores_25_0.png)
    



    
![png](counts_to_scores_files/counts_to_scores_25_1.png)
    



    
![png](counts_to_scores_files/counts_to_scores_25_2.png)
    


## Jackpotting and mutation coverage in pre-selection libraries
We look at the distribution of counts in the "reference" (pre-selection) libraries to see if they seem jackpotted (a few variants at very high frequency):


```python
pre_samples_df = samples_df.query('selection == "reference"')
```

Distribution of mutations along the gene for the pre-selection samples; big spikes may indicate jackpotting:


```python
# this plotting is very slow when lots of samples, so for now code commented out

p = variants.plotMutFreqs(variant_type='all',
                         mut_type='codon',
                         orientation='v',
                         libraries=variants.libraries,
                         samples=pre_samples_df['sample'].unique().tolist(),
                         widthscale=1.5)
_ = p.draw()
```


    
![png](counts_to_scores_files/counts_to_scores_29_0.png)
    


How many mutations are observed frequently in pre-selection libraries?
Note that the libraries have been pre-selected for ACE2 binding, so we expect stop variants to mostly be missing.
Make the plot both for all variants and just single-mutant variants:


```python
# this plotting is very slow when lots of samples, so for now code commented out

for variant_type in ['all', 'single']:
   p = variants.plotCumulMutCoverage(
                         variant_type=variant_type,
                         mut_type='aa',
                         orientation='v',
                         libraries=variants.libraries,
                         samples=pre_samples_df['sample'].unique().tolist(),
                         widthscale=1.8,
                         heightscale=1.2)
   _ = p.draw()
```


    
![png](counts_to_scores_files/counts_to_scores_31_0.png)
    



    
![png](counts_to_scores_files/counts_to_scores_31_1.png)
    


Now make a plot showing what number and fraction of counts are for each variant in each pre-selection sample / library.
If some variants constitute a very high fraction, then that indicates extensive bottlenecking:


```python
for ystat in ['frac_counts', 'count']:
    p = variants.plotCountsPerVariant(ystat=ystat,
                                      libraries=variants.libraries,
                                      samples=pre_samples_df['sample'].unique().tolist(),
                                      orientation='v',
                                      widthscale=1.75,
                                      )
    _ = p.draw()
```


    
![png](counts_to_scores_files/counts_to_scores_33_0.png)
    



    
![png](counts_to_scores_files/counts_to_scores_33_1.png)
    


Now make the same plot breaking down by variant class, which enables determination of which types of variants are at high and low frequencies.
For this plot (unlike one above not classified by category) we only show variants of the primary target (not the homologs), and also group synonymous with wildtype in order to reduce number of plotted categories to make more interpretable:


```python
for ystat in ['frac_counts', 'count']:
    p = variants.plotCountsPerVariant(ystat=ystat,
                                      libraries=variants.libraries,
                                      samples=pre_samples_df['sample'].unique().tolist(),
                                      orientation='v',
                                      widthscale=1.75,
                                      by_variant_class=True,
                                      classifyVariants_kwargs={'syn_as_wt': True},
                                      primary_target_only=True,
                                      )
    _ = p.draw()
```


    
![png](counts_to_scores_files/counts_to_scores_35_0.png)
    



    
![png](counts_to_scores_files/counts_to_scores_35_1.png)
    


We also directly look to see what is the variant in each reference library / sample with the highest fraction counts.
Knowing if the highest frequency variant is shared helps determine **where** in the experiment the jackpotting happened:


```python
frac_counts_per_variant = (
        variants.add_frac_counts(variants.variant_count_df)
        .query(f"sample in {pre_samples_df['sample'].unique().tolist()}")
        )

display(HTML(
    frac_counts_per_variant
    .sort_values('frac_counts', ascending=False)
    .groupby(['library', 'sample'])
    .head(n=1)
    .sort_values(['library', 'sample'])
    .set_index(['library', 'sample'])
    [['frac_counts', 'target', 'barcode', 'aa_substitutions', 'codon_substitutions']]
    .round(4)
    .to_html()
    ))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th></th>
      <th>frac_counts</th>
      <th>target</th>
      <th>barcode</th>
      <th>aa_substitutions</th>
      <th>codon_substitutions</th>
    </tr>
    <tr>
      <th>library</th>
      <th>sample</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th rowspan="3" valign="top">lib1</th>
      <th>bg4e_6-9-none-0-ref</th>
      <td>0.0047</td>
      <td>B1351</td>
      <td>TAACAGGGAAAGACGA</td>
      <td>G9S I142E</td>
      <td>GGT9TCT ATT142GAA</td>
    </tr>
    <tr>
      <th>bg4e_10-14-reference-0-ref</th>
      <td>0.0003</td>
      <td>B1351</td>
      <td>GCATATGCTAGTAATG</td>
      <td>I104M</td>
      <td>ATA104ATG</td>
    </tr>
    <tr>
      <th>bg4e_15-18-reference-0-ref</th>
      <td>0.0003</td>
      <td>B1351</td>
      <td>AAAGATACTACATGGT</td>
      <td>P7R A190S</td>
      <td>CCT7AGA GCG190TCT</td>
    </tr>
    <tr>
      <th rowspan="3" valign="top">lib2</th>
      <th>bg4e_6-9-none-0-ref</th>
      <td>0.0039</td>
      <td>B1351</td>
      <td>TTTATAGCAAACCAAA</td>
      <td>P96E</td>
      <td>CCT96GAA</td>
    </tr>
    <tr>
      <th>bg4e_10-14-reference-0-ref</th>
      <td>0.0003</td>
      <td>B1351</td>
      <td>ACGGTGTTCTATTCGA</td>
      <td>Y143H</td>
      <td>TAT143CAT</td>
    </tr>
    <tr>
      <th>bg4e_15-18-reference-0-ref</th>
      <td>0.0003</td>
      <td>B1351</td>
      <td>ACTAACAAACGAGACA</td>
      <td>Q176S</td>
      <td>CAG176TCT</td>
    </tr>
  </tbody>
</table>


To further where the jackpotting relative to the generation of the reference samples, we plot the correlation among the fraction of counts for the different reference samples.
If the fractions are highly correlated, that indicates that the jackpotting occurred in some upstream step common to the reference samples:


```python
# this code makes a full matrix of scatter plots, but is REALLY SLOW. So for now,
# it is commented out in favor of code that just makes correlation matrix.
for lib, lib_df in frac_counts_per_variant.groupby('library'):
   wide_lib_df = lib_df.pivot_table(index=['target', 'barcode'],
                                    columns='sample',
                                    values='frac_counts')
   g = seaborn.pairplot(wide_lib_df, corner=True, plot_kws={'alpha': 0.5}, diag_kind='kde')
   _ = g.fig.suptitle(lib, size=18)
   plt.show()
```


    
![png](counts_to_scores_files/counts_to_scores_39_0.png)
    



    
![png](counts_to_scores_files/counts_to_scores_39_1.png)
    


## Examine counts for wildtype variants
The type of score we use to quantify escape depends on how well represented wildtype is in the selected libraries.
If wildtype is still well represented, we can use a more conventional functional score that gives differential selection relative to wildtype.
If wildtype is not well represented, then we need an alternative score that does not involve normalizing frequencies to wildtype.

First get average fraction of counts per variant for each variant class:


```python
counts_by_class = (
    variants.variant_count_df
    .pipe(variants.add_frac_counts)
    .pipe(variants.classifyVariants,
          primary_target=variants.primary_target,
          non_primary_target_class='homolog',
          class_as_categorical=True)
    .groupby(['library', 'sample', 'variant_class'])
    .aggregate(avg_frac_counts=pd.NamedAgg('frac_counts', 'mean'))
    .reset_index()
    .merge(samples_df[['sample', 'library', 'date', 'antibody', 'concentration', 'selection']],
           on=['sample', 'library'], validate='many_to_one')
    )

display(HTML(counts_by_class.head().to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>sample</th>
      <th>variant_class</th>
      <th>avg_frac_counts</th>
      <th>date</th>
      <th>antibody</th>
      <th>concentration</th>
      <th>selection</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>lib1</td>
      <td>bg4e_6-9-none-0-ref</td>
      <td>wildtype</td>
      <td>0.000029</td>
      <td>210622</td>
      <td>none</td>
      <td>0</td>
      <td>reference</td>
    </tr>
    <tr>
      <td>lib1</td>
      <td>bg4e_6-9-none-0-ref</td>
      <td>synonymous</td>
      <td>0.000024</td>
      <td>210622</td>
      <td>none</td>
      <td>0</td>
      <td>reference</td>
    </tr>
    <tr>
      <td>lib1</td>
      <td>bg4e_6-9-none-0-ref</td>
      <td>1 nonsynonymous</td>
      <td>0.000022</td>
      <td>210622</td>
      <td>none</td>
      <td>0</td>
      <td>reference</td>
    </tr>
    <tr>
      <td>lib1</td>
      <td>bg4e_6-9-none-0-ref</td>
      <td>&gt;1 nonsynonymous</td>
      <td>0.000015</td>
      <td>210622</td>
      <td>none</td>
      <td>0</td>
      <td>reference</td>
    </tr>
    <tr>
      <td>lib1</td>
      <td>bg4e_6-9-none-0-ref</td>
      <td>stop</td>
      <td>0.000003</td>
      <td>210622</td>
      <td>none</td>
      <td>0</td>
      <td>reference</td>
    </tr>
  </tbody>
</table>


Plot average fraction of all counts per variant for each variant class.
If the values for wildtype are low for the non-reference samples (such as more similar to stop the nonsynonymous), then normalizing by wildtype in calculating scores will probably not work well as wildtype is too depleted:


```python
min_frac = 1e-7  # plot values < this as this

p = (ggplot(counts_by_class
            .assign(avg_frac_counts=lambda x: numpy.clip(x['avg_frac_counts'], min_frac, None))
            ) +
     aes('avg_frac_counts', 'sample', color='selection') +
     geom_point(size=2) +
     scale_color_manual(values=CBPALETTE[1:]) +
     facet_grid('library ~ variant_class') +
     scale_x_log10() +
     theme(axis_text_x=element_text(angle=90),
           figure_size=(2.5 * counts_by_class['variant_class'].nunique(),
                        0.2 * counts_by_class['library'].nunique() * 
                        counts_by_class['sample'].nunique())
           ) +
     geom_vline(xintercept=min_frac, linetype='dotted', color=CBPALETTE[3])
     )

_ = p.draw()
```


    
![png](counts_to_scores_files/counts_to_scores_43_0.png)
    


## Compute escape scores
We use the escape score metric, which does **not** involve normalizing to wildtype and so isn't strongly affected by low wildtype counts.
We compute the scores using the method [dms_variants.codonvarianttable.CodonVariantTable.escape_scores](https://jbloomlab.github.io/dms_variants/dms_variants.codonvarianttable.html?highlight=escape_scores#dms_variants.codonvarianttable.CodonVariantTable.escape_scores).

First, define what samples to compare for each calculation, matching each post-selection (escape) to the pre-selection (reference) sample on the same date:


```python
score_sample_df = (
    samples_df
    .query('selection == "escape"')
    .rename(columns={'sample': 'post_sample',
                     'number_cells': 'pre_cells_sorted'})
    .merge(samples_df
           .query('selection == "reference"')
           [['sample', 'library', 'date', 'number_cells']]
           .rename(columns={'sample': 'pre_sample',
                            'number_cells': 'post_cells_sorted'}),
           how='left', on=['date', 'library'], validate='many_to_one',
           )
    .assign(name=lambda x: x['antibody'] + '_' + x['concentration'].astype(str))
    # add dates to names where needed to make unique
    .assign(n_libs=lambda x: x.groupby(['name', 'date'])['pre_sample'].transform('count'))
    .sort_values(['name', 'date', 'n_libs'], ascending=False)
    .assign(i_name=lambda x: x.groupby(['library', 'name'], sort=False)['pre_sample'].cumcount(),
            name=lambda x: x.apply(lambda r: r['name'] + '_' + str(r['date']) if r['i_name'] > 0 else r['name'],
                                   axis=1),
            )
    .sort_values(['antibody', 'concentration', 'library', 'i_name'])
    # get columns of interest
    [['name', 'library', 'antibody', 'concentration', 'date',
      'pre_sample', 'post_sample', 'frac_escape', 'pre_cells_sorted', 'post_cells_sorted']]
    )

assert len(score_sample_df.groupby(['name', 'library'])) == len(score_sample_df)

print(f"Writing samples used to compute escape scores to {config['escape_score_samples']}\n")
score_sample_df.to_csv(config['escape_score_samples'], index=False)

display(HTML(score_sample_df.to_html(index=False)))
```

    Writing samples used to compute escape scores to results/escape_scores/samples.csv
    



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>name</th>
      <th>library</th>
      <th>antibody</th>
      <th>concentration</th>
      <th>date</th>
      <th>pre_sample</th>
      <th>post_sample</th>
      <th>frac_escape</th>
      <th>pre_cells_sorted</th>
      <th>post_cells_sorted</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>K007_500</td>
      <td>lib1</td>
      <td>K007</td>
      <td>500</td>
      <td>210702</td>
      <td>bg4e_10-14-reference-0-ref</td>
      <td>bg4e_10-K007-500-abneg</td>
      <td>0.038</td>
      <td>510107.0</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>K007_500</td>
      <td>lib2</td>
      <td>K007</td>
      <td>500</td>
      <td>210702</td>
      <td>bg4e_10-14-reference-0-ref</td>
      <td>bg4e_10-K007-500-abneg</td>
      <td>0.033</td>
      <td>494538.0</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>K031_500</td>
      <td>lib1</td>
      <td>K031</td>
      <td>500</td>
      <td>210702</td>
      <td>bg4e_10-14-reference-0-ref</td>
      <td>bg4e_11-K031-500-abneg</td>
      <td>0.072</td>
      <td>730562.0</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>K031_500</td>
      <td>lib2</td>
      <td>K031</td>
      <td>500</td>
      <td>210702</td>
      <td>bg4e_10-14-reference-0-ref</td>
      <td>bg4e_11-K031-500-abneg</td>
      <td>0.069</td>
      <td>712255.0</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>K033_500</td>
      <td>lib1</td>
      <td>K033</td>
      <td>500</td>
      <td>210702</td>
      <td>bg4e_10-14-reference-0-ref</td>
      <td>bg4e_12-K033-500-abneg</td>
      <td>0.077</td>
      <td>802239.0</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>K033_500</td>
      <td>lib2</td>
      <td>K033</td>
      <td>500</td>
      <td>210702</td>
      <td>bg4e_10-14-reference-0-ref</td>
      <td>bg4e_12-K033-500-abneg</td>
      <td>0.074</td>
      <td>773070.0</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>K040_500</td>
      <td>lib1</td>
      <td>K040</td>
      <td>500</td>
      <td>210702</td>
      <td>bg4e_10-14-reference-0-ref</td>
      <td>bg4e_13-K040-500-abneg</td>
      <td>0.064</td>
      <td>687119.0</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>K040_500</td>
      <td>lib2</td>
      <td>K040</td>
      <td>500</td>
      <td>210702</td>
      <td>bg4e_10-14-reference-0-ref</td>
      <td>bg4e_13-K040-500-abneg</td>
      <td>0.062</td>
      <td>679198.0</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>K041_500</td>
      <td>lib1</td>
      <td>K041</td>
      <td>500</td>
      <td>210702</td>
      <td>bg4e_10-14-reference-0-ref</td>
      <td>bg4e_14-K041-500-abneg</td>
      <td>0.055</td>
      <td>601747.0</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>K041_500</td>
      <td>lib2</td>
      <td>K041</td>
      <td>500</td>
      <td>210702</td>
      <td>bg4e_10-14-reference-0-ref</td>
      <td>bg4e_14-K041-500-abneg</td>
      <td>0.052</td>
      <td>601730.0</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>K046_200</td>
      <td>lib1</td>
      <td>K046</td>
      <td>200</td>
      <td>210707</td>
      <td>bg4e_15-18-reference-0-ref</td>
      <td>bg4e_15-K046-200-abneg</td>
      <td>0.048</td>
      <td>561168.0</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>K046_200</td>
      <td>lib2</td>
      <td>K046</td>
      <td>200</td>
      <td>210707</td>
      <td>bg4e_15-18-reference-0-ref</td>
      <td>bg4e_15-K046-200-abneg</td>
      <td>0.041</td>
      <td>419123.0</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>K114_200</td>
      <td>lib1</td>
      <td>K114</td>
      <td>200</td>
      <td>210707</td>
      <td>bg4e_15-18-reference-0-ref</td>
      <td>bg4e_16-K114-200-abneg</td>
      <td>0.053</td>
      <td>638716.0</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>K114_200</td>
      <td>lib2</td>
      <td>K114</td>
      <td>200</td>
      <td>210707</td>
      <td>bg4e_15-18-reference-0-ref</td>
      <td>bg4e_16-K114-200-abneg</td>
      <td>0.047</td>
      <td>542375.0</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>K115_80</td>
      <td>lib1</td>
      <td>K115</td>
      <td>80</td>
      <td>210707</td>
      <td>bg4e_15-18-reference-0-ref</td>
      <td>bg4e_17-K115-80-abneg</td>
      <td>0.165</td>
      <td>1669109.0</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>K115_80</td>
      <td>lib2</td>
      <td>K115</td>
      <td>80</td>
      <td>210707</td>
      <td>bg4e_15-18-reference-0-ref</td>
      <td>bg4e_17-K115-80-abneg</td>
      <td>0.137</td>
      <td>1982477.0</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>K115_old_80</td>
      <td>lib1</td>
      <td>K115_old</td>
      <td>80</td>
      <td>210622</td>
      <td>bg4e_6-9-none-0-ref</td>
      <td>bg4e_8-K115_old-80-abneg</td>
      <td>0.036</td>
      <td>604025.0</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>K115_old_80</td>
      <td>lib2</td>
      <td>K115_old</td>
      <td>80</td>
      <td>210622</td>
      <td>bg4e_6-9-none-0-ref</td>
      <td>bg4e_8-K115_old-80-abneg</td>
      <td>0.034</td>
      <td>566545.0</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>K119_200</td>
      <td>lib1</td>
      <td>K119</td>
      <td>200</td>
      <td>210707</td>
      <td>bg4e_15-18-reference-0-ref</td>
      <td>bg4e_18-K119-200-abneg</td>
      <td>0.041</td>
      <td>548265.0</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>K119_200</td>
      <td>lib2</td>
      <td>K119</td>
      <td>200</td>
      <td>210707</td>
      <td>bg4e_15-18-reference-0-ref</td>
      <td>bg4e_18-K119-200-abneg</td>
      <td>0.047</td>
      <td>630199.0</td>
      <td>NaN</td>
    </tr>
  </tbody>
</table>


Compute the escape scores for variants of the primary target and classify the variants:


```python
print(f"Computing escape scores for {primary_target} variants using {config['escape_score_type']} "
      f"score type with a pseudocount of {config['escape_score_pseudocount']} and "
      f"an escape fraction floor {config['escape_score_floor_E']}, an escape fraction ceiling "
      f"{config['escape_score_ceil_E']}, and grouping variants by {config['escape_score_group_by']}.")

escape_scores = (variants.escape_scores(score_sample_df,
                                        score_type=config['escape_score_type'],
                                        pseudocount=config['escape_score_pseudocount'],
                                        floor_E=config['escape_score_floor_E'],
                                        ceil_E=config['escape_score_ceil_E'],
                                        by=config['escape_score_group_by'],
                                        )
                 .query('target == @primary_target')
                 .pipe(variants.classifyVariants,
                       primary_target=variants.primary_target,
                       syn_as_wt=(config['escape_score_group_by'] == 'aa_substitutions'),
                       )
                 )
print('Here are the first few lines of the resulting escape scores:')
display(HTML(escape_scores.head().to_html(index=False)))
```

    Computing escape scores for B1351 variants using frac_escape score type with a pseudocount of 0.5 and an escape fraction floor 0, an escape fraction ceiling 1, and grouping variants by barcode.
    Here are the first few lines of the resulting escape scores:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>name</th>
      <th>target</th>
      <th>library</th>
      <th>pre_sample</th>
      <th>post_sample</th>
      <th>barcode</th>
      <th>score</th>
      <th>score_var</th>
      <th>pre_count</th>
      <th>post_count</th>
      <th>codon_substitutions</th>
      <th>n_codon_substitutions</th>
      <th>aa_substitutions</th>
      <th>n_aa_substitutions</th>
      <th>variant_class</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>K007_500</td>
      <td>B1351</td>
      <td>lib1</td>
      <td>bg4e_10-14-reference-0-ref</td>
      <td>bg4e_10-K007-500-abneg</td>
      <td>GCATATGCTAGTAATG</td>
      <td>0.000122</td>
      <td>4.281402e-09</td>
      <td>20679</td>
      <td>3</td>
      <td>ATA104ATG</td>
      <td>1</td>
      <td>I104M</td>
      <td>1</td>
      <td>1 nonsynonymous</td>
    </tr>
    <tr>
      <td>K007_500</td>
      <td>B1351</td>
      <td>lib1</td>
      <td>bg4e_10-14-reference-0-ref</td>
      <td>bg4e_10-K007-500-abneg</td>
      <td>AAAGATACTACATGGT</td>
      <td>0.000018</td>
      <td>6.592957e-10</td>
      <td>19916</td>
      <td>0</td>
      <td>CCT7AGA GCG190TCT</td>
      <td>2</td>
      <td>P7R A190S</td>
      <td>2</td>
      <td>&gt;1 nonsynonymous</td>
    </tr>
    <tr>
      <td>K007_500</td>
      <td>B1351</td>
      <td>lib1</td>
      <td>bg4e_10-14-reference-0-ref</td>
      <td>bg4e_10-K007-500-abneg</td>
      <td>ACAGATGATTACAAAA</td>
      <td>0.001334</td>
      <td>5.165953e-08</td>
      <td>18706</td>
      <td>34</td>
      <td>TCT53GCT</td>
      <td>1</td>
      <td>S53A</td>
      <td>1</td>
      <td>1 nonsynonymous</td>
    </tr>
    <tr>
      <td>K007_500</td>
      <td>B1351</td>
      <td>lib1</td>
      <td>bg4e_10-14-reference-0-ref</td>
      <td>bg4e_10-K007-500-abneg</td>
      <td>TAGCCAGCTAACCTAA</td>
      <td>0.000295</td>
      <td>1.163854e-08</td>
      <td>18362</td>
      <td>7</td>
      <td>CTT5TGT</td>
      <td>1</td>
      <td>L5C</td>
      <td>1</td>
      <td>1 nonsynonymous</td>
    </tr>
    <tr>
      <td>K007_500</td>
      <td>B1351</td>
      <td>lib1</td>
      <td>bg4e_10-14-reference-0-ref</td>
      <td>bg4e_10-K007-500-abneg</td>
      <td>GCACATCTAGTAAGAT</td>
      <td>0.000020</td>
      <td>7.863652e-10</td>
      <td>18236</td>
      <td>0</td>
      <td>TTT62TCT ACC148ACT</td>
      <td>2</td>
      <td>F62S</td>
      <td>1</td>
      <td>1 nonsynonymous</td>
    </tr>
  </tbody>
</table>


## Apply pre-selection count filter to variant escape scores
Now determine a pre-selection count filter in order to flag for removal variants with counts that are so low that the estimated score is probably noise.
We know that stop codons should be largely purged pre-selection, and so the counts for them are a good indication of the "noise" threshold.
We therefore set the filter using the number of pre-selection counts for the stop codons.

To do this, we first compute the number of pre-selection counts for stop-codon variants at various quantiles and look at these.
We then take the number of pre-selection counts at the specified quantile as the filter cutoff, and filter scores for all variants with pre-selection counts less than this filter cutoff:


```python
filter_quantile = config['escape_score_stop_quantile_filter']
assert 0 <= filter_quantile <= 1

quantiles = sorted(set([0.5, 0.9, 0.95, 0.98, 0.99, 0.995, 0.999] + [filter_quantile]))

stop_score_counts = (
    escape_scores
    .query('variant_class == "stop"')
    .groupby(['library', 'pre_sample'], observed=True)
    ['pre_count']
    .quantile(q=quantiles)
    .reset_index()
    .rename(columns={'level_2': 'quantile'})
    .pivot_table(index=['pre_sample', 'library'],
                 columns='quantile',
                 values='pre_count')
    )

print('Quantiles of the number of pre-selection counts per variant for stop variants:')
display(HTML(stop_score_counts.to_html(float_format='%.1f')))

print(f"\nSetting the pre-count filter cutoff to the {filter_quantile} quantile:")
pre_count_filter_cutoffs = (
    stop_score_counts
    [filter_quantile]
    .rename('pre_count_filter_cutoff')
    .reset_index()
    )
display(HTML(pre_count_filter_cutoffs.to_html(float_format='%.1f')))
```

    Quantiles of the number of pre-selection counts per variant for stop variants:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>quantile</th>
      <th>0.5</th>
      <th>0.9</th>
      <th>0.95</th>
      <th>0.98</th>
      <th>0.99</th>
      <th>0.995</th>
      <th>0.999</th>
    </tr>
    <tr>
      <th>pre_sample</th>
      <th>library</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th rowspan="2" valign="top">bg4e_6-9-none-0-ref</th>
      <th>lib1</th>
      <td>30.5</td>
      <td>136.1</td>
      <td>160.0</td>
      <td>261.3</td>
      <td>636.3</td>
      <td>774.1</td>
      <td>884.4</td>
    </tr>
    <tr>
      <th>lib2</th>
      <td>29.5</td>
      <td>162.1</td>
      <td>338.1</td>
      <td>424.0</td>
      <td>474.4</td>
      <td>543.2</td>
      <td>598.2</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">bg4e_10-14-reference-0-ref</th>
      <th>lib1</th>
      <td>181.5</td>
      <td>699.0</td>
      <td>1080.0</td>
      <td>1541.9</td>
      <td>2196.8</td>
      <td>3096.0</td>
      <td>3096.0</td>
    </tr>
    <tr>
      <th>lib2</th>
      <td>212.0</td>
      <td>637.2</td>
      <td>953.5</td>
      <td>1424.4</td>
      <td>1593.1</td>
      <td>1702.0</td>
      <td>1702.0</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">bg4e_15-18-reference-0-ref</th>
      <th>lib1</th>
      <td>251.0</td>
      <td>911.0</td>
      <td>1291.0</td>
      <td>1757.6</td>
      <td>2608.1</td>
      <td>4330.0</td>
      <td>4330.0</td>
    </tr>
    <tr>
      <th>lib2</th>
      <td>348.0</td>
      <td>1018.8</td>
      <td>1464.9</td>
      <td>2188.6</td>
      <td>2217.8</td>
      <td>2397.0</td>
      <td>2397.0</td>
    </tr>
  </tbody>
</table>


    
    Setting the pre-count filter cutoff to the 0.9 quantile:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>pre_sample</th>
      <th>library</th>
      <th>pre_count_filter_cutoff</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>bg4e_6-9-none-0-ref</td>
      <td>lib1</td>
      <td>136.1</td>
    </tr>
    <tr>
      <th>1</th>
      <td>bg4e_6-9-none-0-ref</td>
      <td>lib2</td>
      <td>162.1</td>
    </tr>
    <tr>
      <th>2</th>
      <td>bg4e_10-14-reference-0-ref</td>
      <td>lib1</td>
      <td>699.0</td>
    </tr>
    <tr>
      <th>3</th>
      <td>bg4e_10-14-reference-0-ref</td>
      <td>lib2</td>
      <td>637.2</td>
    </tr>
    <tr>
      <th>4</th>
      <td>bg4e_15-18-reference-0-ref</td>
      <td>lib1</td>
      <td>911.0</td>
    </tr>
    <tr>
      <th>5</th>
      <td>bg4e_15-18-reference-0-ref</td>
      <td>lib2</td>
      <td>1018.8</td>
    </tr>
  </tbody>
</table>


Apply the filter to the escape scores, so that scores that fail the pre-selection count filter are now marked with `pass_pre_count_filter` of `False`:


```python
escape_scores = (
    escape_scores
    .merge(pre_count_filter_cutoffs,
           on=['library', 'pre_sample'],
           how='left',
           validate='many_to_one')
    .assign(pass_pre_count_filter=lambda x: x['pre_count'] >= x['pre_count_filter_cutoff'])
    )
```

Plot the fraction of variants of each type that pass the pre-selection count filter in each pre-selection sample.
The ideal filter would have the property such that no *stop* variants pass, all *wildtype* (or *synonymous*) variants pass, and some intermediate fraction of *nonsynonymous* variants pass.
However, if the variant composition in the pre-selection samples is already heavily skewed by jackpotting, there will be some deviation from this ideal behavior.
Here is what the plots actually look like:


```python
frac_pre_pass_filter = (
    escape_scores
    [['pre_sample', 'library', 'target', config['escape_score_group_by'],
      'pre_count', 'pass_pre_count_filter', 'variant_class']]
    .drop_duplicates()
    .groupby(['pre_sample', 'library', 'variant_class'], observed=True)
    .aggregate(n_variants=pd.NamedAgg('pass_pre_count_filter', 'count'),
               n_pass_filter=pd.NamedAgg('pass_pre_count_filter', 'sum')
               )
    .reset_index()
    .assign(frac_pass_filter=lambda x: x['n_pass_filter'] / x['n_variants'],
            pre_sample=lambda x: pd.Categorical(x['pre_sample'], x['pre_sample'].unique(), ordered=True).remove_unused_categories())
    )

p = (ggplot(frac_pre_pass_filter) +
     aes('variant_class', 'frac_pass_filter', fill='variant_class') +
     geom_bar(stat='identity') +
     facet_grid('library ~ pre_sample') +
     theme(axis_text_x=element_text(angle=90),
           figure_size=(3.3 * frac_pre_pass_filter['pre_sample'].nunique(),
                        2 * frac_pre_pass_filter['library'].nunique()),
           panel_grid_major_x=element_blank(),
           ) +
     scale_fill_manual(values=CBPALETTE[1:]) +
     expand_limits(y=(0, 1))
     )

_ = p.draw()
```


    
![png](counts_to_scores_files/counts_to_scores_53_0.png)
    


## Apply ACE2-binding / expression filter to variant mutations
We also used deep mutational scanning to estimate how each mutation affected ACE2 binding and expression in the B.1.351 background.
Here we flag for removal any variants of the primary target that have (or have mutations) that were measured to decrease ACE2-binding or expression beyond a minimal threshold, in order to avoid these variants muddying the signal as spurious escape mutants.

To do this, we first determine all mutations that do / do-not having binding that exceeds the thresholds.

Note that because we are working on this serum-mapping project at the same time as we are working on the ACE2-binding / RBD-expression project, the scores will be preliminary until all final analyses have been done on the DMS project end. So, we will allow either preliminary or "final" measurements to be used. 


```python
if config['dms_scores_preliminary']:
    mut_bind_expr_file = config['prelim_mut_bind_expr']
else: 
    mut_bind_expr_file = config['mut_bind_expr']
    
print(f"Reading ACE2-binding and expression for mutations from {mut_bind_expr_file}, "
      f"and filtering for variants that have single mutations that "
      f"only have mutations with binding >={config['escape_score_min_bind_mut']} and "
      f"expression >={config['escape_score_min_expr_mut']}.")

mut_bind_expr = (pd.read_csv(mut_bind_expr_file)
                 .query('target==@config["primary_target"]')
                 # need to add back the offset numbering for some silly, circuitous reason 
                 .assign(RBD_site=lambda x: x['position']-config['site_number_offset'] ,
                         RBD_mutation=lambda x: x['wildtype']+x['RBD_site'].astype(str)+x['mutant']
                        )
                )

print('Here is what that dataframe looks like:')

display(HTML(mut_bind_expr.query('delta_bind < -2.35').head().to_html(index=False)))
```

    Reading ACE2-binding and expression for mutations from data/final_variant_scores.csv, and filtering for variants that have single mutations that only have mutations with binding >=-3.0 and expression >=-1.0.
    Here is what that dataframe looks like:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>target</th>
      <th>wildtype</th>
      <th>position</th>
      <th>mutant</th>
      <th>mutation</th>
      <th>bind</th>
      <th>delta_bind</th>
      <th>n_bc_bind</th>
      <th>n_libs_bind</th>
      <th>bind_rep1</th>
      <th>bind_rep2</th>
      <th>expr</th>
      <th>delta_expr</th>
      <th>n_bc_expr</th>
      <th>n_libs_expr</th>
      <th>expr_rep1</th>
      <th>expr_rep2</th>
      <th>delta_bind_lib1</th>
      <th>delta_bind_lib2</th>
      <th>delta_expr_lib1</th>
      <th>delta_expr_lib2</th>
      <th>RBD_site</th>
      <th>RBD_mutation</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>B1351</td>
      <td>F</td>
      <td>338</td>
      <td>R</td>
      <td>F338R</td>
      <td>6.58796</td>
      <td>-2.69662</td>
      <td>12</td>
      <td>2</td>
      <td>6.02525</td>
      <td>7.15068</td>
      <td>6.62764</td>
      <td>-3.37374</td>
      <td>15</td>
      <td>2</td>
      <td>6.73035</td>
      <td>6.52493</td>
      <td>-3.18857</td>
      <td>-2.20465</td>
      <td>-3.19115</td>
      <td>-3.55632</td>
      <td>8</td>
      <td>F8R</td>
    </tr>
    <tr>
      <td>B1351</td>
      <td>F</td>
      <td>342</td>
      <td>D</td>
      <td>F342D</td>
      <td>6.45279</td>
      <td>-2.83179</td>
      <td>7</td>
      <td>2</td>
      <td>5.72749</td>
      <td>7.17809</td>
      <td>6.74577</td>
      <td>-3.25560</td>
      <td>11</td>
      <td>2</td>
      <td>6.73301</td>
      <td>6.75854</td>
      <td>-3.48633</td>
      <td>-2.17724</td>
      <td>-3.18849</td>
      <td>-3.32271</td>
      <td>12</td>
      <td>F12D</td>
    </tr>
    <tr>
      <td>B1351</td>
      <td>F</td>
      <td>342</td>
      <td>G</td>
      <td>F342G</td>
      <td>6.67891</td>
      <td>-2.53492</td>
      <td>6</td>
      <td>1</td>
      <td>6.67891</td>
      <td>NaN</td>
      <td>7.17394</td>
      <td>-2.82743</td>
      <td>8</td>
      <td>2</td>
      <td>6.88391</td>
      <td>7.46398</td>
      <td>-2.53491</td>
      <td>NaN</td>
      <td>-3.03759</td>
      <td>-2.61727</td>
      <td>12</td>
      <td>F12G</td>
    </tr>
    <tr>
      <td>B1351</td>
      <td>F</td>
      <td>342</td>
      <td>P</td>
      <td>F342P</td>
      <td>6.23906</td>
      <td>-3.04552</td>
      <td>6</td>
      <td>2</td>
      <td>5.66006</td>
      <td>6.81805</td>
      <td>6.85042</td>
      <td>-3.15096</td>
      <td>9</td>
      <td>2</td>
      <td>6.84025</td>
      <td>6.86058</td>
      <td>-3.55376</td>
      <td>-2.53728</td>
      <td>-3.08125</td>
      <td>-3.22067</td>
      <td>12</td>
      <td>F12P</td>
    </tr>
    <tr>
      <td>B1351</td>
      <td>N</td>
      <td>343</td>
      <td>P</td>
      <td>N343P</td>
      <td>6.87591</td>
      <td>-2.40867</td>
      <td>11</td>
      <td>2</td>
      <td>6.01518</td>
      <td>7.73664</td>
      <td>6.62332</td>
      <td>-3.37806</td>
      <td>16</td>
      <td>2</td>
      <td>6.45641</td>
      <td>6.79024</td>
      <td>-3.19864</td>
      <td>-1.61869</td>
      <td>-3.46509</td>
      <td>-3.29101</td>
      <td>13</td>
      <td>N13P</td>
    </tr>
  </tbody>
</table>



```python
assert mut_bind_expr['RBD_mutation'].nunique() == len(mut_bind_expr)
for prop in ['bind', 'expr']:
    muts_adequate = set(mut_bind_expr
                        .query(f"delta_{prop} >= {config[f'escape_score_min_{prop}_mut']}")
                        ['RBD_mutation']
                        )
    print(f"{len(muts_adequate)} of {len(mut_bind_expr)} mutations have adequate {prop}.")
    escape_scores[f"muts_pass_{prop}_filter"] = (
        escape_scores
        ['aa_substitutions']
        .map(lambda s: set(s.split()).issubset(muts_adequate))
        )

# annotate as passing overall filter if passes all mutation and binding filters:
escape_scores['pass_ACE2bind_expr_filter'] = (
        escape_scores['muts_pass_bind_filter'] &
        escape_scores['muts_pass_expr_filter'] 
        )

display(HTML(escape_scores.query('not pass_ACE2bind_expr_filter & variant_class != "stop"').head().to_html(index=False)))
```

    3291 of 4020 mutations have adequate bind.
    2369 of 4020 mutations have adequate expr.



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>name</th>
      <th>target</th>
      <th>library</th>
      <th>pre_sample</th>
      <th>post_sample</th>
      <th>barcode</th>
      <th>score</th>
      <th>score_var</th>
      <th>pre_count</th>
      <th>post_count</th>
      <th>codon_substitutions</th>
      <th>n_codon_substitutions</th>
      <th>aa_substitutions</th>
      <th>n_aa_substitutions</th>
      <th>variant_class</th>
      <th>pre_count_filter_cutoff</th>
      <th>pass_pre_count_filter</th>
      <th>muts_pass_bind_filter</th>
      <th>muts_pass_expr_filter</th>
      <th>pass_ACE2bind_expr_filter</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>K007_500</td>
      <td>B1351</td>
      <td>lib1</td>
      <td>bg4e_10-14-reference-0-ref</td>
      <td>bg4e_10-K007-500-abneg</td>
      <td>AACTTCAGGGAAGGAC</td>
      <td>0.045556</td>
      <td>2.535852e-06</td>
      <td>13803</td>
      <td>869</td>
      <td>CAG176AAT</td>
      <td>1</td>
      <td>Q176N</td>
      <td>1</td>
      <td>1 nonsynonymous</td>
      <td>699.0</td>
      <td>True</td>
      <td>True</td>
      <td>False</td>
      <td>False</td>
    </tr>
    <tr>
      <td>K007_500</td>
      <td>B1351</td>
      <td>lib1</td>
      <td>bg4e_10-14-reference-0-ref</td>
      <td>bg4e_10-K007-500-abneg</td>
      <td>TTAATAGGTGCGTCCA</td>
      <td>0.000341</td>
      <td>1.787951e-08</td>
      <td>13792</td>
      <td>6</td>
      <td>GCA145GAT</td>
      <td>1</td>
      <td>A145D</td>
      <td>1</td>
      <td>1 nonsynonymous</td>
      <td>699.0</td>
      <td>True</td>
      <td>False</td>
      <td>True</td>
      <td>False</td>
    </tr>
    <tr>
      <td>K007_500</td>
      <td>B1351</td>
      <td>lib1</td>
      <td>bg4e_10-14-reference-0-ref</td>
      <td>bg4e_10-K007-500-abneg</td>
      <td>TTGCTGATAGCTATAA</td>
      <td>0.001655</td>
      <td>8.714374e-08</td>
      <td>13765</td>
      <td>31</td>
      <td>AAG26ATT</td>
      <td>1</td>
      <td>K26I</td>
      <td>1</td>
      <td>1 nonsynonymous</td>
      <td>699.0</td>
      <td>True</td>
      <td>True</td>
      <td>False</td>
      <td>False</td>
    </tr>
    <tr>
      <td>K007_500</td>
      <td>B1351</td>
      <td>lib1</td>
      <td>bg4e_10-14-reference-0-ref</td>
      <td>bg4e_10-K007-500-abneg</td>
      <td>GATGACTGACGCCAAA</td>
      <td>0.056617</td>
      <td>3.342431e-06</td>
      <td>13201</td>
      <td>1033</td>
      <td>CCT169TGG</td>
      <td>1</td>
      <td>P169W</td>
      <td>1</td>
      <td>1 nonsynonymous</td>
      <td>699.0</td>
      <td>True</td>
      <td>True</td>
      <td>False</td>
      <td>False</td>
    </tr>
    <tr>
      <td>K007_500</td>
      <td>B1351</td>
      <td>lib1</td>
      <td>bg4e_10-14-reference-0-ref</td>
      <td>bg4e_10-K007-500-abneg</td>
      <td>GACAGCATGATCGAGC</td>
      <td>0.043612</td>
      <td>2.533981e-06</td>
      <td>13191</td>
      <td>795</td>
      <td>TTC47CAT</td>
      <td>1</td>
      <td>F47H</td>
      <td>1</td>
      <td>1 nonsynonymous</td>
      <td>699.0</td>
      <td>True</td>
      <td>True</td>
      <td>False</td>
      <td>False</td>
    </tr>
  </tbody>
</table>


Print the number of mutations that pass RBD bind, RBD expression, and are not to sites that are disulfide bonds (if specified in config) 


```python
# if we are excluding all cysteines to remove spurious mutations that break disulfide bonds:
if config['exclude_cysteines']:
    print("Here are the number of mutations that pass the bind, express, and disulfide filters:")
    print(len(mut_bind_expr
              .assign(pass_cysteine_filter=lambda x: x['mutation'].str[0] != "C")
              .query(f"delta_bind >= {config[f'escape_score_min_bind_mut']} & \
                       delta_expr >= {config[f'escape_score_min_expr_mut']} & \
                       pass_cysteine_filter")
             ))
    print("There are these many possible mutations (excluding wildtype and disulfides!):")
    print(mut_bind_expr.query('wildtype!="C"')['position'].nunique()*19
         )

else:
    print("Here are the number of mutations that pass the bind and express filters:")
    print(len(mut_bind_expr
              .assign(pass_cysteine_filter=lambda x: x['mutation'].str[0] != "C")
              .query(f"delta_bind >= {config[f'escape_score_min_bind_mut']} & \
                       delta_expr >= {config[f'escape_score_min_expr_mut']}")
             ))
    print("There are these many possible mutations (excluding wildtype!):")
    print(mut_bind_expr['position'].nunique()*19
         )
```

    Here are the number of mutations that pass the bind, express, and disulfide filters:
    2207
    There are these many possible mutations (excluding wildtype and disulfides!):
    3667


Plot the fraction of variants that **have already passed the pre-count filter** that are filtered by the ACE2-binding or expression thresholds:


```python
print('These are the sites that are involved in disulfide bonds:')
print(mut_bind_expr.query('wildtype=="C"')['position'].unique())
```

    These are the sites that are involved in disulfide bonds:
    [336 361 379 391 432 480 488 525]



```python
frac_ACE2bind_expr_pass_filter = (
    escape_scores
    .query('pass_pre_count_filter == True')
    [['pre_sample', 'library', 'target', config['escape_score_group_by'],
      'pre_count', 'pass_ACE2bind_expr_filter', 'variant_class']]
    .drop_duplicates()
    .groupby(['pre_sample', 'library', 'variant_class'], observed=True)
    .aggregate(n_variants=pd.NamedAgg('pass_ACE2bind_expr_filter', 'count'),
               n_pass_filter=pd.NamedAgg('pass_ACE2bind_expr_filter', 'sum')
               )
    .reset_index()
    .assign(frac_pass_filter=lambda x: x['n_pass_filter'] / x['n_variants'],
            pre_sample=lambda x: pd.Categorical(x['pre_sample'], x['pre_sample'].unique(), ordered=True).remove_unused_categories())
    )

p = (ggplot(frac_ACE2bind_expr_pass_filter) +
     aes('variant_class', 'frac_pass_filter', fill='variant_class') +
     geom_bar(stat='identity') +
     facet_grid('library ~ pre_sample') +
     theme(axis_text_x=element_text(angle=90),
           figure_size=(3.3 * frac_ACE2bind_expr_pass_filter['pre_sample'].nunique(),
                        2 * frac_ACE2bind_expr_pass_filter['library'].nunique()),
           panel_grid_major_x=element_blank(),
           ) +
     scale_fill_manual(values=CBPALETTE[1:]) +
     expand_limits(y=(0, 1))
     )

_ = p.draw()
```


    
![png](counts_to_scores_files/counts_to_scores_61_0.png)
    


## Examine and write escape scores
Plot the distribution of escape scores across variants of different classes **among those that pass both the pre-selection count filter and the ACE2-binding / expression filter**.
If things are working correctly, we don't expect escape in wildtype (or synonymous variants), but do expect escape for some small fraction of nonsynymous variants.
Also, we do not plot the scores for the stop codon variant class, as most stop-codon variants have already been filtered out so this category is largely noise:


```python
nfacets = len(escape_scores.groupby(['library', 'name']).nunique())
ncol = min(8, nfacets)
nrow = math.ceil(nfacets / ncol)

df = (escape_scores
      .query('(pass_pre_count_filter == True) & (pass_ACE2bind_expr_filter == True)')
      .query('variant_class != "stop"')
      )
     
p = (ggplot(df) +
     aes('variant_class', 'score', color='variant_class') +
     geom_boxplot(outlier_size=1.5, outlier_alpha=0.1) +
     facet_wrap('~ name + library', ncol=ncol) +
     theme(axis_text_x=element_text(angle=90),
           figure_size=(2.35 * ncol, 3 * nrow),
           panel_grid_major_x=element_blank(),
           ) +
     scale_fill_manual(values=CBPALETTE[1:]) +
     scale_color_manual(values=CBPALETTE[1:])
     )

_ = p.draw()
```


    
![png](counts_to_scores_files/counts_to_scores_63_0.png)
    


Also, we want to see how much the high escape scores are correlated with simple coverage.
To do this, we plot the correlation between escape score and pre-selection count just for the nonsynonymous variants (which are the ones that we expect to have true escape).
The plots below have a lot of overplotting, but are still sufficient to test of the score is simply correlated with the pre-selection counts or not.
The hoped for result is that the escape score doesn't appear to be strongly correlated with pre-selection counts:


```python
p = (ggplot(escape_scores
            .query('pass_pre_count_filter == True')
            .query('(pass_pre_count_filter == True) & (pass_ACE2bind_expr_filter == True)')
            .query('variant_class=="1 nonsynonymous"')
            ) +
     aes('pre_count', 'score') +
     geom_point(alpha=0.1, size=1) +
     facet_wrap('~ name + library', ncol=ncol) +
     theme(axis_text_x=element_text(angle=90),
           figure_size=(2.35 * ncol, 2.35 * nrow),
           ) +
     scale_fill_manual(values=CBPALETTE[1:]) +
     scale_color_manual(values=CBPALETTE[1:]) +
     scale_x_log10()
     )

_ = p.draw()
```


    
![png](counts_to_scores_files/counts_to_scores_65_0.png)
    


Write the escape scores to a file:


```python
print(f"Writing escape scores for {primary_target} to {config['escape_scores']}")
escape_scores.to_csv(config['escape_scores'], index=False, float_format='%.4g')
```

    Writing escape scores for B1351 to results/escape_scores/scores.csv


### Now we will also remove anything that did not pass all the filters above. 


```python
escape_scores_primary = (escape_scores
                         .query('(pass_pre_count_filter == True) & (pass_ACE2bind_expr_filter == True)')
                        )

display(HTML(escape_scores_primary.head().to_html()))
print(f"Read {len(escape_scores_primary)} scores.")
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>name</th>
      <th>target</th>
      <th>library</th>
      <th>pre_sample</th>
      <th>post_sample</th>
      <th>barcode</th>
      <th>score</th>
      <th>score_var</th>
      <th>pre_count</th>
      <th>post_count</th>
      <th>codon_substitutions</th>
      <th>n_codon_substitutions</th>
      <th>aa_substitutions</th>
      <th>n_aa_substitutions</th>
      <th>variant_class</th>
      <th>pre_count_filter_cutoff</th>
      <th>pass_pre_count_filter</th>
      <th>muts_pass_bind_filter</th>
      <th>muts_pass_expr_filter</th>
      <th>pass_ACE2bind_expr_filter</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>K007_500</td>
      <td>B1351</td>
      <td>lib1</td>
      <td>bg4e_10-14-reference-0-ref</td>
      <td>bg4e_10-K007-500-abneg</td>
      <td>GCATATGCTAGTAATG</td>
      <td>0.000122</td>
      <td>4.281402e-09</td>
      <td>20679</td>
      <td>3</td>
      <td>ATA104ATG</td>
      <td>1</td>
      <td>I104M</td>
      <td>1</td>
      <td>1 nonsynonymous</td>
      <td>699.0</td>
      <td>True</td>
      <td>True</td>
      <td>True</td>
      <td>True</td>
    </tr>
    <tr>
      <th>1</th>
      <td>K007_500</td>
      <td>B1351</td>
      <td>lib1</td>
      <td>bg4e_10-14-reference-0-ref</td>
      <td>bg4e_10-K007-500-abneg</td>
      <td>AAAGATACTACATGGT</td>
      <td>0.000018</td>
      <td>6.592957e-10</td>
      <td>19916</td>
      <td>0</td>
      <td>CCT7AGA GCG190TCT</td>
      <td>2</td>
      <td>P7R A190S</td>
      <td>2</td>
      <td>&gt;1 nonsynonymous</td>
      <td>699.0</td>
      <td>True</td>
      <td>True</td>
      <td>True</td>
      <td>True</td>
    </tr>
    <tr>
      <th>2</th>
      <td>K007_500</td>
      <td>B1351</td>
      <td>lib1</td>
      <td>bg4e_10-14-reference-0-ref</td>
      <td>bg4e_10-K007-500-abneg</td>
      <td>ACAGATGATTACAAAA</td>
      <td>0.001334</td>
      <td>5.165953e-08</td>
      <td>18706</td>
      <td>34</td>
      <td>TCT53GCT</td>
      <td>1</td>
      <td>S53A</td>
      <td>1</td>
      <td>1 nonsynonymous</td>
      <td>699.0</td>
      <td>True</td>
      <td>True</td>
      <td>True</td>
      <td>True</td>
    </tr>
    <tr>
      <th>3</th>
      <td>K007_500</td>
      <td>B1351</td>
      <td>lib1</td>
      <td>bg4e_10-14-reference-0-ref</td>
      <td>bg4e_10-K007-500-abneg</td>
      <td>TAGCCAGCTAACCTAA</td>
      <td>0.000295</td>
      <td>1.163854e-08</td>
      <td>18362</td>
      <td>7</td>
      <td>CTT5TGT</td>
      <td>1</td>
      <td>L5C</td>
      <td>1</td>
      <td>1 nonsynonymous</td>
      <td>699.0</td>
      <td>True</td>
      <td>True</td>
      <td>True</td>
      <td>True</td>
    </tr>
    <tr>
      <th>4</th>
      <td>K007_500</td>
      <td>B1351</td>
      <td>lib1</td>
      <td>bg4e_10-14-reference-0-ref</td>
      <td>bg4e_10-K007-500-abneg</td>
      <td>GCACATCTAGTAAGAT</td>
      <td>0.000020</td>
      <td>7.863652e-10</td>
      <td>18236</td>
      <td>0</td>
      <td>TTT62TCT ACC148ACT</td>
      <td>2</td>
      <td>F62S</td>
      <td>1</td>
      <td>1 nonsynonymous</td>
      <td>699.0</td>
      <td>True</td>
      <td>True</td>
      <td>True</td>
      <td>True</td>
    </tr>
  </tbody>
</table>


    Read 380257 scores.


### Count number of barcodes per mutation and remove variants with >1 amino acid substitution
Also add the number of barocdes per mutation to the `escape_scores` dataframe and plot this. 
But first see how many variants there are with >1 mutation, and query the dataframe to look at them qualitatively. 


```python
p = (ggplot(escape_scores_primary) +
     aes('n_aa_substitutions', fill='variant_class') +
     geom_bar() +
     facet_wrap('~library + pre_sample', ncol=5) +
     theme(
#            figure_size=(3.3 * escape_scores_primary['pre_sample'].nunique(),
#                         2 * escape_scores_primary['library'].nunique()),
         figure_size=(12, 4),
           panel_grid_major_x=element_blank(),
           ) +
     scale_fill_manual(values=CBPALETTE[1:]) +
     expand_limits(y=(0, 1))
     )

_ = p.draw()
```


    
![png](counts_to_scores_files/counts_to_scores_71_0.png)
    


### Filter dataframe on single mutations that are present in at least `n` number of variants (specified in `config.yaml` file)
Now see how many `n_single_mut_measurements` there are for each variant:


```python
print(f'Remove anything with fewer than {config["escape_frac_min_single_mut_measurements"]} single mutant measurements (barcodes)')

raw_avg_single_mut_scores = (
    escape_scores_primary
    .query('n_aa_substitutions == 1')
    .rename(columns={'name': 'selection',
                     'aa_substitutions': 'mutation'})
    .groupby(['selection', 'library', 'mutation'])
    .aggregate(raw_single_mut_score=pd.NamedAgg('score', 'mean'),
               n_single_mut_measurements=pd.NamedAgg('barcode', 'count')
              )
    .assign(sufficient_measurements=lambda x: (
                (x['n_single_mut_measurements'] >= config['escape_frac_min_single_mut_measurements'])))
    .reset_index()
    )

# remove mutations with insufficient measurements
effects_df = (raw_avg_single_mut_scores
              .query('sufficient_measurements == True')
              .drop(columns='sufficient_measurements')
              )

# some error checks
assert len(effects_df) == len(effects_df.drop_duplicates()), 'duplicate rows in `effects_df`'
assert all(effects_df['raw_single_mut_score'].notnull() | (effects_df['n_single_mut_measurements'] == 0))

display(HTML(effects_df.head().to_html()))
```

    Remove anything with fewer than 2 single mutant measurements (barcodes)



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>selection</th>
      <th>library</th>
      <th>mutation</th>
      <th>raw_single_mut_score</th>
      <th>n_single_mut_measurements</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>K007_500</td>
      <td>lib1</td>
      <td>A105S</td>
      <td>0.003593</td>
      <td>3</td>
    </tr>
    <tr>
      <th>1</th>
      <td>K007_500</td>
      <td>lib1</td>
      <td>A105T</td>
      <td>0.000215</td>
      <td>7</td>
    </tr>
    <tr>
      <th>2</th>
      <td>K007_500</td>
      <td>lib1</td>
      <td>A145C</td>
      <td>0.000248</td>
      <td>6</td>
    </tr>
    <tr>
      <th>3</th>
      <td>K007_500</td>
      <td>lib1</td>
      <td>A145E</td>
      <td>0.000305</td>
      <td>2</td>
    </tr>
    <tr>
      <th>4</th>
      <td>K007_500</td>
      <td>lib1</td>
      <td>A145F</td>
      <td>0.000320</td>
      <td>6</td>
    </tr>
  </tbody>
</table>


### Should we exclude mutations for which the wildtype identity is a cysteine?
These would be mutations that break disulfide bonds (unless there is an unpaired cysteine in the protein that does not form a disulfide bond, of course). 
Note that this approach would not work well for something like the SARS-CoV-2 NTD, where it has been documented (by the Veesler lab) that mutations in the B.1.427/429 (epislon) variant can rearrange disulfide bonds, leading to major structural rearrangements in the NTD, and yet an apparently fully functional spike. 

If we are excluding cysteines, do that now:


```python
# if we are excluding all cysteines to remove spurious mutations that break disulfide bonds:
if config['exclude_cysteines']:
    print(f'Excluding mutations where the wildtype identity is a cysteine')
    effects_df = effects_df.assign(pass_cysteine_filter=lambda x: x['mutation'].str[0] != "C",
                                  )
    disulfides_to_drop = effects_df.query('pass_cysteine_filter==False')['mutation'].unique()
    print(f'Specifically, excluding: {disulfides_to_drop}')
    effects_df=effects_df.query('pass_cysteine_filter').drop(columns='pass_cysteine_filter')

else:
    print(f'Retaining mutations where the wildtype identity is a cysteine')
```

    Excluding mutations where the wildtype identity is a cysteine
    Specifically, excluding: ['C158Q' 'C195A' 'C195D' 'C195E' 'C195F' 'C195G' 'C195H' 'C195I' 'C195K'
     'C195L' 'C195M' 'C195N' 'C195P' 'C195Q' 'C195S' 'C195T' 'C195W' 'C195Y'
     'C31D' 'C31E' 'C31G' 'C31K' 'C31P' 'C31R' 'C31S' 'C31T' 'C31V' 'C61A'
     'C61D' 'C61E' 'C61F' 'C61G' 'C61H' 'C61I' 'C61K' 'C61L' 'C61M' 'C61N'
     'C61P' 'C61Q' 'C61R' 'C61S' 'C61T' 'C61V' 'C61W' 'C61Y' 'C6D' 'C31A'
     'C31N' 'C6N']


We need to compute the escape scores (calculated as [here](https://jbloomlab.github.io/dms_variants/dms_variants.codonvarianttable.html?highlight=escape_scores#dms_variants.codonvarianttable.CodonVariantTable.escape_scores)) back to escape fractions. We define a function to do this depending on the score type:


```python
def score_to_frac(score):
    """Convert escape score to escape fraction."""
    if pd.isnull(score):
        return pd.NA
    floor = config['escape_score_floor_E']
    ceil = config['escape_score_ceil_E']
    if config['escape_score_type'] == 'frac_escape':
        return min(ceil, max(floor, score))  # just the score after applying ceiling and floor
    elif config['escape_score_type'] == 'log_escape':
        # convert score back to fraction, apply ceiling, then re-scale so floor is 0
        frac = 2**score
        frac = min(ceil, max(floor, frac))
        frac = (frac - floor) / (ceil - floor)
        return frac
    else:
        raise ValueError(f"invalid `escape_score_type` of {config['escape_score_type']}")

effects_df = (
    effects_df
    .assign(
            mut_escape_frac_single_mut=lambda x: x['raw_single_mut_score'].map(score_to_frac),
            )
    )
```

### Average the escape score across all barcodes of the same mutation, for each library, and the average of both libraries. 
Add rows that are the average of the two libraries for the fraction escape for all mutations that are present in both libraries (and if in just one library, the value in that or purge depending on config values printed here), the number of libraries each mutation is measured in, and the sum of the statistics giving the number of measurements:


```python
min_libs = config['escape_frac_avg_min_libraries']
min_single = config['escape_frac_avg_min_single']
print(f"Only taking average of mutations with escape fractions in >={min_libs} libraries "
      f"or with >={min_single} single-mutant measurements total.")

effects_df = (
    effects_df
    .query('library != "average"')  # git rid of averages if already there
    .assign(nlibs=1)
    .append(effects_df
            .query('library != "average"')
            .groupby(['selection', 'mutation'])
            .aggregate(nlibs=pd.NamedAgg('library', 'count'),
                       mut_escape_frac_single_mut=pd.NamedAgg('mut_escape_frac_single_mut',
                                                              lambda s: s.mean(skipna=True)),
                       n_single_mut_measurements=pd.NamedAgg('n_single_mut_measurements', 'sum'),
                       )
            .query('(nlibs >= @min_libs) or (n_single_mut_measurements >= @min_single)')
            .assign(library="average")
            .reset_index(),
            ignore_index=True, sort=False,
            )
    )

print(len(effects_df.query('nlibs>1')))
print(len(effects_df.query('nlibs==1')))
```

    Only taking average of mutations with escape fractions in >=2 libraries or with >=2 single-mutant measurements total.
    16892
    39242


Plot the correlations of the escape fractions among the two libraries for all selections performed on both libraries. 


```python
libraries = [lib for lib in effects_df['library'].unique() if lib != "average"]
assert len(libraries) == 2, 'plot only makes sense if 2 libraries'

# wide data frame with each library's score in a different column
effects_df_wide = (
    effects_df
    .query('library != "average"')
    .query(f"n_single_mut_measurements >= 1")
    # just get selections with 2 libraries
    .assign(nlibs=lambda x: x.groupby('selection')['library'].transform('nunique'))
    .query('nlibs == 2')
    # now make columns for each library, only keep mutants with scores for both libs
    [['selection', 'mutation', 'library', 'mut_escape_frac_single_mut']]
    .pivot_table(index=['selection', 'mutation'],
                 columns='library',
                 values='mut_escape_frac_single_mut',
                 aggfunc='first')
    .reset_index()
    .dropna(axis=0)
    )

# correlations between libraries
corrs = (
    effects_df_wide
    .groupby('selection')
    [libraries]
    .corr(method='pearson')
    .reset_index()
    .query('library == @libraries[0]')
    .assign(correlation=lambda x: 'R=' + x[libraries[1]].round(2).astype(str))
    [['selection', 'correlation']]
    # add number of mutations measured
    .merge(effects_df_wide
           .groupby('selection')
           .size()
           .rename('n')
           .reset_index()
           )
    .assign(correlation=lambda x: x['correlation'] + ', N=' + x['n'].astype(str))
    )

# plot correlations
nfacets = effects_df_wide['selection'].nunique()
ncol = min(nfacets, 5)
nrow = math.ceil(nfacets / ncol)
xmin = effects_df_wide[libraries[0]].min()
xspan = effects_df_wide[libraries[0]].max() - xmin
ymin = effects_df_wide[libraries[1]].min()
yspan = effects_df_wide[libraries[1]].max() - ymin
p = (ggplot(effects_df_wide) +
     aes(libraries[0], libraries[1]) +
     geom_point(alpha=0.2) +
     geom_text(mapping=aes(label='correlation'),
               data=corrs,
               x=0.01 * xspan + xmin,
               y=0.99 * yspan + ymin,
               size=10,
               ha='left',
               va='top',
               ) +
     facet_wrap('~ selection', ncol=ncol) +
     theme(figure_size=(2.5 * ncol, 2.5 * nrow),
           plot_title=element_text(size=14)) +
     ggtitle('Mutation-level escape fractions')
     )

_ = p.draw()
```


    
![png](counts_to_scores_files/counts_to_scores_81_0.png)
    


### Escape at site level
The above analysis estimates the effects of mutations. We also compute escape statistics at the site level. First, add sites to the data frame of mutational effects:


```python
effects_df = (
    effects_df
    .assign(site=lambda x: x['mutation'].str[1: -1].astype(int),
            wildtype=lambda x: x['mutation'].str[0],
            mutant=lambda x: x['mutation'].str[-1],
            )
    )
```

Now compute some site-level metrics. These are the average and total escape fraction at each site over all mutations at the site:


```python
site_effects_df = (
    effects_df
    .groupby(['selection', 'library', 'site'])
    .aggregate(
        site_avg_escape_frac_single_mut=pd.NamedAgg('mut_escape_frac_single_mut',
                                                    lambda s: s.mean(skipna=True)),
        site_total_escape_frac_single_mut=pd.NamedAgg('mut_escape_frac_single_mut', 'sum'),
        )
    .reset_index()
    )

display(HTML(site_effects_df.head().to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>selection</th>
      <th>library</th>
      <th>site</th>
      <th>site_avg_escape_frac_single_mut</th>
      <th>site_total_escape_frac_single_mut</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>K007_500</td>
      <td>average</td>
      <td>1</td>
      <td>0.001708</td>
      <td>0.030738</td>
    </tr>
    <tr>
      <td>K007_500</td>
      <td>average</td>
      <td>2</td>
      <td>0.001303</td>
      <td>0.024751</td>
    </tr>
    <tr>
      <td>K007_500</td>
      <td>average</td>
      <td>3</td>
      <td>0.001456</td>
      <td>0.027658</td>
    </tr>
    <tr>
      <td>K007_500</td>
      <td>average</td>
      <td>4</td>
      <td>0.001358</td>
      <td>0.021720</td>
    </tr>
    <tr>
      <td>K007_500</td>
      <td>average</td>
      <td>5</td>
      <td>0.001000</td>
      <td>0.018994</td>
    </tr>
  </tbody>
</table>


Plot correlations between libraries of the same selection for the site-level statistics:


```python
libraries = [lib for lib in effects_df['library'].unique() if lib != "average"]
assert len(libraries) == 2, 'plot only makes sense if 2 libraries'

for val in ['site_avg_escape_frac_single_mut', 'site_total_escape_frac_single_mut']:

    # wide data frame with each library's score in a different column
    site_effects_df_wide = (
        site_effects_df
        .query('library != "average"')
        # just get selections with 2 libraries
        .assign(nlibs=lambda x: x.groupby('selection')['library'].transform('nunique'))
        .query('nlibs == 2')
        # now make columns for each library, only keep sites with scores for both libs
        .pivot_table(index=['selection', 'site'],
                     columns='library',
                     values=val)
        .reset_index()
        .dropna(axis=0)
        )

    # correlations between libraries
    corrs = (
        site_effects_df_wide
        .groupby('selection')
        [libraries]
        .corr(method='pearson')
        .reset_index()
        .query('library == @libraries[0]')
        .assign(correlation=lambda x: 'R=' + x[libraries[1]].round(2).astype(str))
        [['selection', 'correlation']]
        # add number of mutations measured
        .merge(site_effects_df_wide
               .groupby('selection')
               .size()
               .rename('n')
               .reset_index()
               )
        .assign(correlation=lambda x: x['correlation'] + ', N=' + x['n'].astype(str))
        )

    # plot correlations
    nfacets = site_effects_df_wide['selection'].nunique()
    ncol = min(nfacets, 5)
    nrow = math.ceil(nfacets / ncol)
    xmin = site_effects_df_wide[libraries[0]].min()
    xspan = site_effects_df_wide[libraries[0]].max() - xmin
    ymin = site_effects_df_wide[libraries[1]].min()
    yspan = site_effects_df_wide[libraries[1]].max() - ymin
    p = (ggplot(site_effects_df_wide) +
         aes(libraries[0], libraries[1]) +
         geom_point(alpha=0.2) +
         geom_text(mapping=aes(label='correlation'),
                   data=corrs,
                   x=0.01 * xspan + xmin,
                   y=0.99 * yspan + ymin,
                   size=10,
                   ha='left',
                   va='top',
                   ) +
         facet_wrap('~ selection', ncol=ncol) +
         theme(figure_size=(2.5 * ncol, 2.5 * nrow),
               plot_title=element_text(size=14)) +
         ggtitle(val)
         )

    _ = p.draw()
```


    
![png](counts_to_scores_files/counts_to_scores_87_0.png)
    



    
![png](counts_to_scores_files/counts_to_scores_87_1.png)
    


## Write file with escape fractions at mutation and site levels
We write a files that has the mutation- and site-level escape fractions. This file has both the separate measurements for each library plus the average across libraries for all mutations measured in both libraries. We name the columns in such a way that this file can be used as [dms-view data file](https://dms-view.github.io/docs/dataupload):


```python
escape_fracs_to_write = (
    effects_df
    .merge(site_effects_df,
           how='left',
           validate='many_to_one',
           on=['selection', 'library', 'site'])
    .assign(protein_chain=config['escape_frac_protein_chain'],
            protein_site=lambda x: x['site'] + config['site_number_offset'],
            label_site=lambda x: x['protein_site'],
            condition=lambda x: x['selection'].where(x['library'] == "average", x['selection'] + '_' + x['library']),
            mutation=lambda x: x['mutant'],  # dms-view uses mutation to refer to mutant amino acid
            )
    [['selection', 'library', 'condition', 'site', 'label_site', 'wildtype', 'mutation',
      'protein_chain', 'protein_site', 'mut_escape_frac_single_mut', 'site_total_escape_frac_single_mut',
      'site_avg_escape_frac_single_mut', 'nlibs',
      ]]
    .sort_values(['library', 'selection', 'site', 'mutation'])
    )

print('Here are the first few lines that will be written to the escape-fraction file:')
display(HTML(escape_fracs_to_write.head().to_html(index=False)))

print(f"\nWriting to {config['escape_fracs']}")
escape_fracs_to_write.to_csv(config['escape_fracs'], index=False, float_format='%.4g')

```

    Here are the first few lines that will be written to the escape-fraction file:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>selection</th>
      <th>library</th>
      <th>condition</th>
      <th>site</th>
      <th>label_site</th>
      <th>wildtype</th>
      <th>mutation</th>
      <th>protein_chain</th>
      <th>protein_site</th>
      <th>mut_escape_frac_single_mut</th>
      <th>site_total_escape_frac_single_mut</th>
      <th>site_avg_escape_frac_single_mut</th>
      <th>nlibs</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>K007_500</td>
      <td>average</td>
      <td>K007_500</td>
      <td>1</td>
      <td>331</td>
      <td>N</td>
      <td>A</td>
      <td>E</td>
      <td>331</td>
      <td>0.000303</td>
      <td>0.030738</td>
      <td>0.001708</td>
      <td>2</td>
    </tr>
    <tr>
      <td>K007_500</td>
      <td>average</td>
      <td>K007_500</td>
      <td>1</td>
      <td>331</td>
      <td>N</td>
      <td>C</td>
      <td>E</td>
      <td>331</td>
      <td>0.000240</td>
      <td>0.030738</td>
      <td>0.001708</td>
      <td>2</td>
    </tr>
    <tr>
      <td>K007_500</td>
      <td>average</td>
      <td>K007_500</td>
      <td>1</td>
      <td>331</td>
      <td>N</td>
      <td>D</td>
      <td>E</td>
      <td>331</td>
      <td>0.000405</td>
      <td>0.030738</td>
      <td>0.001708</td>
      <td>2</td>
    </tr>
    <tr>
      <td>K007_500</td>
      <td>average</td>
      <td>K007_500</td>
      <td>1</td>
      <td>331</td>
      <td>N</td>
      <td>E</td>
      <td>E</td>
      <td>331</td>
      <td>0.000218</td>
      <td>0.030738</td>
      <td>0.001708</td>
      <td>2</td>
    </tr>
    <tr>
      <td>K007_500</td>
      <td>average</td>
      <td>K007_500</td>
      <td>1</td>
      <td>331</td>
      <td>N</td>
      <td>F</td>
      <td>E</td>
      <td>331</td>
      <td>0.000375</td>
      <td>0.030738</td>
      <td>0.001708</td>
      <td>2</td>
    </tr>
  </tbody>
</table>


    
    Writing to results/escape_scores/escape_fracs.csv



```python

```
