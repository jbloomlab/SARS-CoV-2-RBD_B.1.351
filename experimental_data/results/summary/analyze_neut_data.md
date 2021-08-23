# Analyze neutralization data
This Python Jupyter notebook analyzes the neutralization data.

Import Python modules.
We use [neutcurve](https://jbloomlab.github.io/neutcurve/) to plot the neutralization curves:


```python
import os
import re
import warnings

from IPython.display import display, HTML
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib as mpl
import matplotlib.pyplot as plt

import numpy as np
import pandas as pd
from plotnine import *
from statistics import geometric_mean

import neutcurve
from neutcurve.colorschemes import CBMARKERS, CBPALETTE
import seaborn

import yaml

print(f"Using `neutcurve` version {neutcurve.__version__}")
```

    Using `neutcurve` version 0.5.7


Read config file


```python
with open('config.yaml') as f:
    config = yaml.safe_load(f)
```

Specify input / output files:


```python
# output directory
results='results/neut_titers'
os.makedirs(results, exist_ok = True)

# input files
fracinfect_file = 'results/neut_titers/fracinfect.csv'

# output files
neut_titers_file = f'{results}/neut_titers.csv'
all_replicate_curves = f'{results}/all_replicate_curves.pdf'
```

Read in the neutralization data, dropping sera/viruses that were messed up and repeated:


```python
print(f"Reading neutralization data from {fracinfect_file}")
fracinfect = pd.read_csv(fracinfect_file).replace({'B.1.351':'wildtype', 'mock':'wildtype'})

# order the viruses
virus_order = config['virus_order']
serum_order = config['serum_order']

print(f"Length before dropping anything = {len(fracinfect.index)}")
    
if config['neut_samples_ignore']:    
    for dat in config['neut_samples_ignore']:
        viruses = config['neut_samples_ignore'][dat]
        print(f'From {dat}, dropping {viruses}')
        l = len((fracinfect[(fracinfect['virus'].isin(viruses)) & (fracinfect['date'] == dat)]))
        print(fracinfect[(fracinfect['virus'].isin(viruses)) & (fracinfect['date'] == dat)]['virus'].unique())
        fracinfect = fracinfect.drop(fracinfect[((fracinfect['virus'].isin(viruses)) & (fracinfect['date'] == dat))].index)
        print(f"Length after dropping {l} rows from {viruses} from {dat} = {len(fracinfect.index)}")

fracinfect = (
    fracinfect
    .assign(replicate_with_date=lambda x: x['replicate'].astype(str) +
                                          ' (' + x['date'] + ')')
    .query('virus in @virus_order & serum in @serum_order')
    .assign(virus=lambda x: pd.Categorical(x['virus'], virus_order, ordered=True))
    .rename(columns={'replicate': 'replicate_on_date'})
)
fracinfect = (
    fracinfect
    .merge(fracinfect
           .sort_values('date')
           [['serum', 'virus', 'replicate_with_date']]
           .drop_duplicates()
           .assign(replicate_all_dates=lambda x: x.groupby(['serum', 'virus'])
                                                  ['replicate_with_date']
                                                  .transform('cumcount') + 1
                   ),
            how='left', on=['serum', 'virus', 'replicate_with_date'], validate='many_to_one',
            )
    )

# show first few lines of data frame
display(HTML(fracinfect.head().to_html(index=False)))
```

    Reading neutralization data from results/neut_titers/fracinfect.csv
    Length before dropping anything = 1536



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>serum</th>
      <th>virus</th>
      <th>replicate_on_date</th>
      <th>concentration</th>
      <th>fraction infectivity</th>
      <th>date</th>
      <th>replicate_with_date</th>
      <th>replicate_all_dates</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>pre-pandemic</td>
      <td>wildtype</td>
      <td>1</td>
      <td>0.040000</td>
      <td>1.0450</td>
      <td>2021-06-10</td>
      <td>1 (2021-06-10)</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pre-pandemic</td>
      <td>wildtype</td>
      <td>1</td>
      <td>0.010000</td>
      <td>0.9724</td>
      <td>2021-06-10</td>
      <td>1 (2021-06-10)</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pre-pandemic</td>
      <td>wildtype</td>
      <td>1</td>
      <td>0.002500</td>
      <td>1.0190</td>
      <td>2021-06-10</td>
      <td>1 (2021-06-10)</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pre-pandemic</td>
      <td>wildtype</td>
      <td>1</td>
      <td>0.000625</td>
      <td>1.0040</td>
      <td>2021-06-10</td>
      <td>1 (2021-06-10)</td>
      <td>1</td>
    </tr>
    <tr>
      <td>pre-pandemic</td>
      <td>wildtype</td>
      <td>1</td>
      <td>0.000156</td>
      <td>0.8126</td>
      <td>2021-06-10</td>
      <td>1 (2021-06-10)</td>
      <td>1</td>
    </tr>
  </tbody>
</table>



```python
fits = neutcurve.curvefits.CurveFits(
            data=fracinfect,
            replicate_col='replicate_all_dates',
            fixbottom=config['fixbottom'],
            fixtop=config['fixtop'],
    
            )
```


```python
with warnings.catch_warnings():
    warnings.simplefilter('ignore')  # ignore fitting warnings
    fig, _ = fits.plotReplicates(ncol=8,
                                 legendtitle='replicate',
                                 xlabel='serum dilution',
                                 viruses=fracinfect['virus'].sort_values().unique(),
                                 colors=plt.rcParams['axes.prop_cycle'].by_key()['color'] * 2,
                                 markers=['o', '^', 's', 'D', 'v', '<', '>', 'p'] * 2,
                                 )
    
print(f"Saving plot to {all_replicate_curves}\n")
fig.savefig(all_replicate_curves)
fig.tight_layout()
display(fig)
plt.close(fig)
```

    Saving plot to results/neut_titers/all_replicate_curves.pdf
    



![png](analyze_neut_data_files/analyze_neut_data_9_1.png)


Use [neutcurve](https://jbloomlab.github.io/neutcurve/) to fit neutralization curves to all of the data:


```python
fitparams = pd.DataFrame(columns=['serum', 'virus', 'ic50', 'NT50', 'ic50_bound', 'date'])

for d in fracinfect['date'].unique():
    fits = neutcurve.CurveFits(fracinfect.query('date==@d'),
                               replicate_col='replicate_on_date',
                               fixbottom=config['fixbottom'],
                               fixtop=config['fixtop'],
                              )

    fp = (
        fits.fitParams(average_only=False)
        .assign(NT50=lambda x: 1/x['ic50'],
                date=d
               )
        .replace({'WT':'wildtype', 'B.1.351':'wildtype'})
        # get columns of interest
        [['serum', 'virus', 'ic50', 'NT50', 'ic50_bound', 'date', 'replicate', 'top']] 
        .assign(ic50_is_bound=lambda x: x['ic50_bound'].map({'lower': True,
                                                          'interpolated': False}))
        )
    fitparams=fitparams.append(fp, ignore_index=True)

fitparams.head()
```

    /fh/fast/bloom_j/software/miniconda3/envs/BloomLab/lib/python3.8/site-packages/scipy/optimize/minpack.py:807: OptimizeWarning: Covariance of the parameters could not be estimated
    /fh/fast/bloom_j/software/miniconda3/envs/BloomLab/lib/python3.8/site-packages/neutcurve/hillcurve.py:741: RuntimeWarning: invalid value encountered in power





<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>serum</th>
      <th>virus</th>
      <th>ic50</th>
      <th>NT50</th>
      <th>ic50_bound</th>
      <th>date</th>
      <th>replicate</th>
      <th>top</th>
      <th>ic50_is_bound</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>pre-pandemic</td>
      <td>wildtype</td>
      <td>0.04</td>
      <td>25.0</td>
      <td>lower</td>
      <td>2021-06-10</td>
      <td>1</td>
      <td>True</td>
      <td>True</td>
    </tr>
    <tr>
      <th>1</th>
      <td>pre-pandemic</td>
      <td>wildtype</td>
      <td>0.04</td>
      <td>25.0</td>
      <td>lower</td>
      <td>2021-06-10</td>
      <td>2</td>
      <td>True</td>
      <td>True</td>
    </tr>
    <tr>
      <th>2</th>
      <td>pre-pandemic</td>
      <td>wildtype</td>
      <td>0.04</td>
      <td>25.0</td>
      <td>lower</td>
      <td>2021-06-10</td>
      <td>average</td>
      <td>True</td>
      <td>True</td>
    </tr>
    <tr>
      <th>3</th>
      <td>K006</td>
      <td>wildtype</td>
      <td>0.04</td>
      <td>25.0</td>
      <td>lower</td>
      <td>2021-06-10</td>
      <td>1</td>
      <td>True</td>
      <td>True</td>
    </tr>
    <tr>
      <th>4</th>
      <td>K006</td>
      <td>wildtype</td>
      <td>0.04</td>
      <td>25.0</td>
      <td>lower</td>
      <td>2021-06-10</td>
      <td>2</td>
      <td>True</td>
      <td>True</td>
    </tr>
  </tbody>
</table>
</div>



### Compare IC50s and fold-change IC50s when the top of the neutralization curve is fixed vs. not fixed

Here we compare the IC50s and fold-change IC50s with the two different methods of fitting the curves. 
See the config file for how we ultimately decided to fit curves.


```python
fitparams_fixtop = pd.DataFrame(columns=['serum', 'virus', 'ic50_fixtop', 'date'])

if config['fixtop']:
    fixtop_alternative=False
else:
    fixtop_alternative=True

for d in fracinfect['date'].unique():
    fits = neutcurve.CurveFits(fracinfect.query('date==@d'),
                               replicate_col='replicate_on_date',
                               fixtop=fixtop_alternative,
                              )

    fp = (
        fits.fitParams(average_only=False)
        .replace({'WT':'wildtype', 'mock':'wildtype', 'B.1.351':'wildtype'})
        .assign(date=d)
        # get columns of interest
        [['serum', 'virus', 'ic50', 'date', 'replicate']]
        .rename(columns={'ic50':'ic50_fixtop'})
        )
    fitparams_fixtop=fitparams_fixtop.append(fp, ignore_index=True)

fitparams_fixtop = (fitparams_fixtop
                    .merge(fitparams_fixtop
                           .query('virus == "wildtype"')
                           [['serum', 'ic50_fixtop', 'date', 'replicate']]
                           .rename(columns={'ic50_fixtop': 'wildtype_ic50_fixtop'}),
                           on=['serum', 'date', 'replicate'],
                           how='left',
                           validate='many_to_one',
                          )
                    .assign(fold_change=lambda x: x['ic50_fixtop'] / x['wildtype_ic50_fixtop'],
                           )
                    .merge(fitparams
                           [['virus', 'serum', 'ic50', 'date', 'replicate']]
                           .merge(fitparams.query('virus == "wildtype"')
                                  [['serum', 'ic50', 'date', 'replicate']]
                    .rename(columns={'ic50': 'wildtype_ic50'}),
                                  on=['serum', 'date', 'replicate'],
                                  how='left',
                                  validate='many_to_one',
                                 )
                           .assign(fold_change_fixtop=lambda x: x['ic50'] / x['wildtype_ic50'],
                                  ),
                           how='outer',
                           on=['serum', 'virus', 'date', 'replicate']
                          )
                   )


display(HTML(fitparams_fixtop.head().to_html()))
```

    /fh/fast/bloom_j/software/miniconda3/envs/BloomLab/lib/python3.8/site-packages/scipy/optimize/minpack.py:807: OptimizeWarning: Covariance of the parameters could not be estimated
    /fh/fast/bloom_j/software/miniconda3/envs/BloomLab/lib/python3.8/site-packages/neutcurve/hillcurve.py:741: RuntimeWarning: invalid value encountered in power
    /fh/fast/bloom_j/software/miniconda3/envs/BloomLab/lib/python3.8/site-packages/neutcurve/hillcurve.py:451: RuntimeWarning: invalid value encountered in sqrt



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>serum</th>
      <th>virus</th>
      <th>ic50_fixtop</th>
      <th>date</th>
      <th>replicate</th>
      <th>wildtype_ic50_fixtop</th>
      <th>fold_change</th>
      <th>ic50</th>
      <th>wildtype_ic50</th>
      <th>fold_change_fixtop</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>pre-pandemic</td>
      <td>wildtype</td>
      <td>0.04</td>
      <td>2021-06-10</td>
      <td>1</td>
      <td>0.04</td>
      <td>1.0</td>
      <td>0.04</td>
      <td>0.04</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>1</th>
      <td>pre-pandemic</td>
      <td>wildtype</td>
      <td>0.04</td>
      <td>2021-06-10</td>
      <td>2</td>
      <td>0.04</td>
      <td>1.0</td>
      <td>0.04</td>
      <td>0.04</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>2</th>
      <td>pre-pandemic</td>
      <td>wildtype</td>
      <td>0.04</td>
      <td>2021-06-10</td>
      <td>average</td>
      <td>0.04</td>
      <td>1.0</td>
      <td>0.04</td>
      <td>0.04</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>3</th>
      <td>K006</td>
      <td>wildtype</td>
      <td>0.04</td>
      <td>2021-06-10</td>
      <td>1</td>
      <td>0.04</td>
      <td>1.0</td>
      <td>0.04</td>
      <td>0.04</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>4</th>
      <td>K006</td>
      <td>wildtype</td>
      <td>0.04</td>
      <td>2021-06-10</td>
      <td>2</td>
      <td>0.04</td>
      <td>1.0</td>
      <td>0.04</td>
      <td>0.04</td>
      <td>1.0</td>
    </tr>
  </tbody>
</table>


Plot correlations between IC50 whether top is fixed or not


```python
for param in [('ic50', 'ic50_fixtop'), ('fold_change', 'fold_change_fixtop')]:
    p = (ggplot(fitparams_fixtop
                .query("replicate != 'average'")
                ) +
         aes(param[0], param[1],
            ) +
         geom_point(aes(fill='date'), size=1.5, alpha=0.5, ) + #fill='#999999', 
         scale_x_log10(name=f'fix_top as specified in config ({param[0]})') +
         scale_y_log10(name=f'fix_top not as specified in config ({param[0]})') +
         theme_classic() +
         theme(axis_text_x=element_text(angle=90),
               figure_size=(3, 3),
               ) +
         scale_fill_manual(values=CBPALETTE*3)
         )

    _ = p.draw()
    
    plotfile = f'{results}/fixtop_corr_{param[0]}.pdf'
    print(f"Saving to {plotfile}")
    p.save(plotfile, verbose=False)
```

    Saving to results/neut_titers/fixtop_corr_ic50.pdf
    Saving to results/neut_titers/fixtop_corr_fold_change.pdf



![png](analyze_neut_data_files/analyze_neut_data_15_1.png)



![png](analyze_neut_data_files/analyze_neut_data_15_2.png)


Make a plot showing all viruses against each sera:


```python
for d in fracinfect['date'].unique():
    fits = (neutcurve.CurveFits(fracinfect.query('date==@d'), 
                                replicate_col='replicate_on_date', 
                                fixbottom=config['fixbottom'],
                                fixtop=config['fixtop'],
                               )
           )
    xlab= 'serum dilution'
    name= 'sera'

    fig, axes = fits.plotSera(xlabel=xlab,max_viruses_per_subplot=2,
                              colors=CBPALETTE*3, 
                              markers=CBMARKERS*3,
#                               attempt_shared_legend=True,
                              fix_lims={'ymax':1.25}
                             )

    plotfile = f'{results}/{d}_mutant_neuts.pdf'
    print(f"Saving to {plotfile}")
    fig.savefig(plotfile, bbox_inches='tight')
```

    /fh/fast/bloom_j/software/miniconda3/envs/BloomLab/lib/python3.8/site-packages/neutcurve/hillcurve.py:741: RuntimeWarning: invalid value encountered in power


    Saving to results/neut_titers/2021-06-10_mutant_neuts.pdf
    Saving to results/neut_titers/2021-08-20_mutant_neuts.pdf
    Saving to results/neut_titers/2021-08-21_mutant_neuts.pdf



![png](analyze_neut_data_files/analyze_neut_data_17_2.png)



![png](analyze_neut_data_files/analyze_neut_data_17_3.png)



![png](analyze_neut_data_files/analyze_neut_data_17_4.png)


### Calculate fold-change IC50 relative to the gemetric mean of the wildtype virus against each serum on each date
* Get neutralization titers, 
* Drop "average" replicate
* Calculate the geometric mean of the wildtype virus against each serum on each date
* Calculate fold-change IC50


```python
neut_titers = (
    fitparams
    .merge((fitparams
            .query('virus == "wildtype" & replicate != "average"')
            .groupby(['serum', 'date'])
            
            # get the geometric mean of the two wildtype replicates 
            .agg(wildtype_ic50=pd.NamedAgg(column="ic50", aggfunc=geometric_mean))
            .reset_index()
           ),
           on=['serum', 'date'],
           how='left',
           validate='many_to_one',
           )
    .assign(fold_change=lambda x: x['ic50'] / x['wildtype_ic50'],)
    )


display(HTML(neut_titers.head().to_html(index=False)))
neut_titers.to_csv(neut_titers_file, index=False)
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>serum</th>
      <th>virus</th>
      <th>ic50</th>
      <th>NT50</th>
      <th>ic50_bound</th>
      <th>date</th>
      <th>replicate</th>
      <th>top</th>
      <th>ic50_is_bound</th>
      <th>wildtype_ic50</th>
      <th>fold_change</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>pre-pandemic</td>
      <td>wildtype</td>
      <td>0.04</td>
      <td>25.0</td>
      <td>lower</td>
      <td>2021-06-10</td>
      <td>1</td>
      <td>True</td>
      <td>True</td>
      <td>0.04</td>
      <td>1.0</td>
    </tr>
    <tr>
      <td>pre-pandemic</td>
      <td>wildtype</td>
      <td>0.04</td>
      <td>25.0</td>
      <td>lower</td>
      <td>2021-06-10</td>
      <td>2</td>
      <td>True</td>
      <td>True</td>
      <td>0.04</td>
      <td>1.0</td>
    </tr>
    <tr>
      <td>pre-pandemic</td>
      <td>wildtype</td>
      <td>0.04</td>
      <td>25.0</td>
      <td>lower</td>
      <td>2021-06-10</td>
      <td>average</td>
      <td>True</td>
      <td>True</td>
      <td>0.04</td>
      <td>1.0</td>
    </tr>
    <tr>
      <td>K006</td>
      <td>wildtype</td>
      <td>0.04</td>
      <td>25.0</td>
      <td>lower</td>
      <td>2021-06-10</td>
      <td>1</td>
      <td>True</td>
      <td>True</td>
      <td>0.04</td>
      <td>1.0</td>
    </tr>
    <tr>
      <td>K006</td>
      <td>wildtype</td>
      <td>0.04</td>
      <td>25.0</td>
      <td>lower</td>
      <td>2021-06-10</td>
      <td>2</td>
      <td>True</td>
      <td>True</td>
      <td>0.04</td>
      <td>1.0</td>
    </tr>
  </tbody>
</table>


As we can see below, the fold-change in IC50 for each wildtype replicate is no longer exactly 1 (because we are comparing to the geometric mean of the replicate measurements). Here I am pulling out hte most extreme fold_change IC50s for wildtype (relative to the geometric mean), and we see that the most extreme values are 0.75 and 1.32, for M03-day-119 on 2021-07-16. 


```python
display(HTML(neut_titers.query('virus=="wildtype" & replicate != "average" & (fold_change <0.8 | fold_change >1.25)').head(100).to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>serum</th>
      <th>virus</th>
      <th>ic50</th>
      <th>NT50</th>
      <th>ic50_bound</th>
      <th>date</th>
      <th>replicate</th>
      <th>top</th>
      <th>ic50_is_bound</th>
      <th>wildtype_ic50</th>
      <th>fold_change</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>K031</td>
      <td>wildtype</td>
      <td>0.000496</td>
      <td>2015.434699</td>
      <td>interpolated</td>
      <td>2021-06-10</td>
      <td>1</td>
      <td>True</td>
      <td>False</td>
      <td>0.000335</td>
      <td>1.481277</td>
    </tr>
    <tr>
      <td>K031</td>
      <td>wildtype</td>
      <td>0.000226</td>
      <td>4422.227189</td>
      <td>interpolated</td>
      <td>2021-06-10</td>
      <td>2</td>
      <td>True</td>
      <td>False</td>
      <td>0.000335</td>
      <td>0.675093</td>
    </tr>
    <tr>
      <td>K046</td>
      <td>wildtype</td>
      <td>0.000648</td>
      <td>1542.465186</td>
      <td>interpolated</td>
      <td>2021-08-21</td>
      <td>1</td>
      <td>True</td>
      <td>False</td>
      <td>0.000474</td>
      <td>1.368471</td>
    </tr>
    <tr>
      <td>K046</td>
      <td>wildtype</td>
      <td>0.000346</td>
      <td>2888.596520</td>
      <td>interpolated</td>
      <td>2021-08-21</td>
      <td>2</td>
      <td>True</td>
      <td>False</td>
      <td>0.000474</td>
      <td>0.730742</td>
    </tr>
    <tr>
      <td>K119</td>
      <td>wildtype</td>
      <td>0.000266</td>
      <td>3761.207286</td>
      <td>interpolated</td>
      <td>2021-08-21</td>
      <td>1</td>
      <td>True</td>
      <td>False</td>
      <td>0.000207</td>
      <td>1.281555</td>
    </tr>
    <tr>
      <td>K119</td>
      <td>wildtype</td>
      <td>0.000162</td>
      <td>6177.343348</td>
      <td>interpolated</td>
      <td>2021-08-21</td>
      <td>2</td>
      <td>True</td>
      <td>False</td>
      <td>0.000207</td>
      <td>0.780302</td>
    </tr>
  </tbody>
</table>


### Plot the fold-change IC50 relative to wildtype.
We will also plot each wild type replicate (as each deviates slightly from 1).  


```python
sera_for_further_analysis = config['sera_for_further_analysis']

neut_titers = (neut_titers
               .query('serum in @sera_for_further_analysis')
                   .assign(serum=lambda x: pd.Categorical(x['serum'],ordered=True,categories=sera_for_further_analysis),
                           virus=lambda x: pd.Categorical(x['virus'],ordered=True,categories=virus_order),
                          )
                  )

for b in ('average', 'not_average'):
    if b=='average':
        df = neut_titers.query('replicate=="average"')
    else:
        df = neut_titers.query('replicate!="average"')
    p = (ggplot(df
                ) +
         aes('virus', 'fold_change', shape='ic50_is_bound',
            ) +
         geom_point(aes(fill='date'), size=2.5, alpha=0.5) + 
         scale_y_log10(name='fold decrease in neutralization') +
         facet_wrap('~serum', ncol=4) +
         theme_classic() +
         theme(axis_text_x=element_text(angle=90),
               axis_title_x=element_blank(),
               strip_margin_y=0.35,
               strip_background_x=element_blank(),
    #            subplots_adjust={'hspace':1},
               figure_size=(0.25 * (neut_titers['virus'].nunique())*neut_titers['serum'].nunique(), 5),
               ) +
         geom_hline(yintercept=1, linetype='dashed', size=1,
                    alpha=0.6, color=CBPALETTE[0]) +
         geom_hline(data=neut_titers.query('virus=="RBD antibodies depleted"').query('replicate=="average"'),
                    mapping=aes(yintercept='fold_change'),
                    color=CBPALETTE[1],
                    alpha=0.7,
                    size=1,
                    linetype='dotted',
                   ) +
         scale_shape_manual(values=['o','^'], name='limit of detection')+
         scale_fill_manual(values=CBPALETTE[1:])
         )

    _ = p.draw()

    plotfile = f'{results}/fold_change_IC50_{b}.pdf'
    print(f"Saving to {plotfile}")
    p.save(plotfile, verbose=False)
```

    Saving to results/neut_titers/fold_change_IC50_average.pdf
    Saving to results/neut_titers/fold_change_IC50_not_average.pdf



![png](analyze_neut_data_files/analyze_neut_data_23_1.png)



![png](analyze_neut_data_files/analyze_neut_data_23_2.png)



```python
for metric in ['fold_change', 'ic50']:
    for virus_set, virus_subsample in config['virus_subsets'].items():
        print(f'Making plot for {metric} for {virus_set}:')
        
        ylab={'fold_change':'fold decrease in neutralization', 'ic50':'IC50'}

        dates=neut_titers.query("virus in @virus_subsample & virus not in ['wildtype', 'RBD antibodies depleted']")['date'].unique()

        p = (ggplot(neut_titers
                    .query("virus in @virus_subsample & date in @dates & replicate!= 'average'")
                    ) +
             aes('virus', metric, shape='ic50_is_bound',
                ) +
             geom_point(aes(fill='date'), size=2.5, alpha=0.5, ) + #fill='#999999', 
             scale_y_log10(name=ylab[metric]) +
             facet_wrap('~serum', ncol=4) +
             theme_classic() +
             theme(axis_text_x=element_text(angle=90),
                   axis_title_x=element_blank(),
                   strip_margin_y=0.35,
                   strip_background_x=element_blank(),
                   figure_size=(0.2 * (len(virus_subsample)+1)*neut_titers['serum'].nunique(), 5),
                   ) +
             geom_hline(data=(neut_titers
                              .query('virus in ["wildtype", "RBD antibodies depleted"] & replicate!="average"')
                              .groupby(['serum', 'virus'])
                              .agg({metric: geometric_mean})
                              .reset_index()
                             ),
                        inherit_aes=False,
                        mapping=aes(yintercept=metric, color='virus'),
                        alpha=0.7,
                        size=0.5,
                        linetype='dotted',
                       ) +
             scale_shape_manual(values=['o','^'], name='limit of detection') +
             scale_color_manual(values=CBPALETTE*3, guide=False) +
             scale_fill_manual(values=CBPALETTE*3)
             )

        _ = p.draw()

        plotfile = f'{results}/{metric}_{virus_set}.pdf'
        print(f"Saving to {plotfile}")
        p.save(plotfile, verbose=False)
```

    Making plot for fold_change for all:


    /fh/fast/bloom_j/software/miniconda3/envs/BloomLab/lib/python3.8/site-packages/plotnine/layer.py:467: PlotnineWarning: geom_hline : Removed 80 rows containing missing values.


    Saving to results/neut_titers/fold_change_all.pdf


    /fh/fast/bloom_j/software/miniconda3/envs/BloomLab/lib/python3.8/site-packages/plotnine/layer.py:467: PlotnineWarning: geom_hline : Removed 80 rows containing missing values.


    Making plot for ic50 for all:


    /fh/fast/bloom_j/software/miniconda3/envs/BloomLab/lib/python3.8/site-packages/plotnine/layer.py:467: PlotnineWarning: geom_hline : Removed 80 rows containing missing values.


    Saving to results/neut_titers/ic50_all.pdf


    /fh/fast/bloom_j/software/miniconda3/envs/BloomLab/lib/python3.8/site-packages/plotnine/layer.py:467: PlotnineWarning: geom_hline : Removed 80 rows containing missing values.



![png](analyze_neut_data_files/analyze_neut_data_24_8.png)



![png](analyze_neut_data_files/analyze_neut_data_24_9.png)


### Plot the results for all individuals on one plot and add the geometric mean of all 8 measurements (2 replicates * 4 sera) for each virus


```python
for metric in ['fold_change', 'ic50']:
    for virus_set, virus_subsample in config['virus_subsets'].items():
        print(f'Making plot for {metric} for {virus_set}:')
        
        ylab={'fold_change':'fold decrease in neutralization', 'ic50':'IC50'}

        dates=neut_titers.query("virus in @virus_subsample & virus not in ['wildtype']")['date'].unique()

        p = (ggplot(neut_titers
                    .query("virus in @virus_subsample & date in @dates & replicate!= 'average'")
                    ) +
             aes('virus', metric, shape='ic50_is_bound', fill='serum',
                ) +
             geom_point(size=2.5, alpha=0.5, ) + 
             geom_crossbar(data=(neut_titers
                                 .query("virus in @virus_subsample & date in @dates & replicate!= 'average'")
                                 .groupby('virus')
                                 .agg({metric: geometric_mean})
                                 .reset_index()
                                 .dropna()
                                ),
                           inherit_aes=False,
                           mapping=aes(x='virus', y=metric, ymin=metric, ymax=metric),
                  ) +
             scale_y_log10(name=ylab[metric]) +
             theme_classic() +
             theme(axis_text_x=element_text(angle=90),
                   axis_title_x=element_blank(),
                   strip_margin_y=0.35,
                   strip_background_x=element_blank(),
                   figure_size=(0.1 * (len(virus_subsample)+1)*neut_titers['serum'].nunique(), 2.5),
                   ) +
             geom_hline(data=(neut_titers
                              .query('virus in ["wildtype"] & replicate!="average"')
                              .groupby(['virus'])
                              .agg({metric: geometric_mean})
                              .reset_index()
                             ),
                        inherit_aes=False,
                        mapping=aes(yintercept=metric, color='virus'),
                        alpha=0.7,
                        size=0.5,
                        linetype='dotted',
                       ) +
             scale_shape_manual(values=['o','^'], name='limit of detection') +
             scale_color_manual(values=CBPALETTE*3, guide=False) +
             scale_fill_manual(values=CBPALETTE*3)
             )

        _ = p.draw()

        plotfile = f'{results}/{metric}_{virus_set}_aggregate.pdf'
        print(f"Saving to {plotfile}")
        p.save(plotfile, verbose=False)
```

    Making plot for fold_change for all:


    /fh/fast/bloom_j/software/miniconda3/envs/BloomLab/lib/python3.8/site-packages/plotnine/layer.py:467: PlotnineWarning: geom_hline : Removed 11 rows containing missing values.


    Saving to results/neut_titers/fold_change_all_aggregate.pdf


    /fh/fast/bloom_j/software/miniconda3/envs/BloomLab/lib/python3.8/site-packages/plotnine/layer.py:467: PlotnineWarning: geom_hline : Removed 11 rows containing missing values.


    Making plot for ic50 for all:


    /fh/fast/bloom_j/software/miniconda3/envs/BloomLab/lib/python3.8/site-packages/plotnine/layer.py:467: PlotnineWarning: geom_hline : Removed 11 rows containing missing values.


    Saving to results/neut_titers/ic50_all_aggregate.pdf


    /fh/fast/bloom_j/software/miniconda3/envs/BloomLab/lib/python3.8/site-packages/plotnine/layer.py:467: PlotnineWarning: geom_hline : Removed 11 rows containing missing values.



![png](analyze_neut_data_files/analyze_neut_data_26_8.png)



![png](analyze_neut_data_files/analyze_neut_data_26_9.png)



```python

```
