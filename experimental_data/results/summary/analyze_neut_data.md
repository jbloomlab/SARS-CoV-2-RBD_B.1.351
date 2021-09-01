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
print(config['rename_sera'])
```

    {'23_d21': 'participant A (day 21)', '24C_d32': 'participant C (day 32)', '23C_d26': 'participant I (day 26)', '1C_d26': 'participant B (day 26)'}



```python
print(f"Reading neutralization data from {fracinfect_file}")
fracinfect = (pd.read_csv(fracinfect_file)
              .replace({'B.1.351':'wildtype', 'mock':'wildtype', 'D614G':'wildtype'})
              .replace(config['rename_sera'])
             )

# order the viruses
virus_order = config['virus_order']
serum_order = config['serum_order']

print(f"Length before dropping anything = {len(fracinfect.index)}")
    
if config['neut_samples_ignore']:    
    for dat in config['neut_samples_ignore']:
        viruses = config['neut_samples_ignore'][dat]
        print(f'From {dat}, dropping {viruses}')
        l = len((fracinfect[(fracinfect['virus'].isin(viruses)) & (fracinfect['date'].astype(str) == str(dat))]))
        print(fracinfect[(fracinfect['virus'].isin(viruses)) & (fracinfect['date'].astype(str) == str(dat))]['virus'].unique())
        fracinfect = fracinfect.drop(fracinfect[((fracinfect['virus'].isin(viruses)) & (fracinfect['date'].astype(str) == str(dat)))].index)
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
    Length before dropping anything = 2240
    From 2021-08-20, dropping ['B.1.351-K484Q']
    ['B.1.351-K484Q']
    Length after dropping 64 rows from ['B.1.351-K484Q'] from 2021-08-20 = 2176
    From 2021-08-21, dropping ['B.1.351-K484Q']
    ['B.1.351-K484Q']
    Length after dropping 64 rows from ['B.1.351-K484Q'] from 2021-08-21 = 2112



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
      <td>K041</td>
      <td>B.1.351-K484Q</td>
      <td>1</td>
      <td>0.040000</td>
      <td>4.610000e-07</td>
      <td>2021-08-27</td>
      <td>1 (2021-08-27)</td>
      <td>2</td>
    </tr>
    <tr>
      <td>K041</td>
      <td>B.1.351-K484Q</td>
      <td>1</td>
      <td>0.013330</td>
      <td>2.739000e-03</td>
      <td>2021-08-27</td>
      <td>1 (2021-08-27)</td>
      <td>2</td>
    </tr>
    <tr>
      <td>K041</td>
      <td>B.1.351-K484Q</td>
      <td>1</td>
      <td>0.004444</td>
      <td>2.070000e-02</td>
      <td>2021-08-27</td>
      <td>1 (2021-08-27)</td>
      <td>2</td>
    </tr>
    <tr>
      <td>K041</td>
      <td>B.1.351-K484Q</td>
      <td>1</td>
      <td>0.001481</td>
      <td>2.479000e-01</td>
      <td>2021-08-27</td>
      <td>1 (2021-08-27)</td>
      <td>2</td>
    </tr>
    <tr>
      <td>K041</td>
      <td>B.1.351-K484Q</td>
      <td>1</td>
      <td>0.000494</td>
      <td>9.689000e-01</td>
      <td>2021-08-27</td>
      <td>1 (2021-08-27)</td>
      <td>2</td>
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
                                 fix_lims={'ymax':1.25},
                                 )
    
print(f"Saving plot to {all_replicate_curves}\n")
fig.savefig(all_replicate_curves)
fig.tight_layout()
display(fig)
plt.close(fig)
```

    Saving plot to results/neut_titers/all_replicate_curves.pdf
    



    
![png](analyze_neut_data_files/analyze_neut_data_10_1.png)
    


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

    /fh/fast/bloom_j/computational_notebooks/agreaney/2021/SARS-CoV-2-RBD_B.1.351/env/lib/python3.8/site-packages/neutcurve/hillcurve.py:451: RuntimeWarning: invalid value encountered in sqrt
    /fh/fast/bloom_j/computational_notebooks/agreaney/2021/SARS-CoV-2-RBD_B.1.351/env/lib/python3.8/site-packages/scipy/optimize/minpack.py:833: OptimizeWarning: Covariance of the parameters could not be estimated
    /fh/fast/bloom_j/computational_notebooks/agreaney/2021/SARS-CoV-2-RBD_B.1.351/env/lib/python3.8/site-packages/neutcurve/hillcurve.py:741: RuntimeWarning: invalid value encountered in power
    /fh/fast/bloom_j/computational_notebooks/agreaney/2021/SARS-CoV-2-RBD_B.1.351/env/lib/python3.8/site-packages/neutcurve/hillcurve.py:741: RuntimeWarning: divide by zero encountered in true_divide
    /fh/fast/bloom_j/computational_notebooks/agreaney/2021/SARS-CoV-2-RBD_B.1.351/env/lib/python3.8/site-packages/neutcurve/hillcurve.py:741: RuntimeWarning: invalid value encountered in power
    /fh/fast/bloom_j/computational_notebooks/agreaney/2021/SARS-CoV-2-RBD_B.1.351/env/lib/python3.8/site-packages/scipy/optimize/minpack.py:833: OptimizeWarning: Covariance of the parameters could not be estimated





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
      <td>K041</td>
      <td>B.1.351-K484Q</td>
      <td>0.001137</td>
      <td>879.582342</td>
      <td>interpolated</td>
      <td>2021-08-27</td>
      <td>1</td>
      <td>True</td>
      <td>False</td>
    </tr>
    <tr>
      <th>1</th>
      <td>K041</td>
      <td>B.1.351-K484Q</td>
      <td>0.001380</td>
      <td>724.862206</td>
      <td>interpolated</td>
      <td>2021-08-27</td>
      <td>2</td>
      <td>True</td>
      <td>False</td>
    </tr>
    <tr>
      <th>2</th>
      <td>K041</td>
      <td>B.1.351-K484Q</td>
      <td>0.001393</td>
      <td>717.838223</td>
      <td>interpolated</td>
      <td>2021-08-27</td>
      <td>average</td>
      <td>True</td>
      <td>False</td>
    </tr>
    <tr>
      <th>3</th>
      <td>K041</td>
      <td>wildtype</td>
      <td>0.000481</td>
      <td>2077.386543</td>
      <td>interpolated</td>
      <td>2021-08-27</td>
      <td>1</td>
      <td>True</td>
      <td>False</td>
    </tr>
    <tr>
      <th>4</th>
      <td>K041</td>
      <td>wildtype</td>
      <td>0.000548</td>
      <td>1824.055323</td>
      <td>interpolated</td>
      <td>2021-08-27</td>
      <td>2</td>
      <td>True</td>
      <td>False</td>
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

    /fh/fast/bloom_j/computational_notebooks/agreaney/2021/SARS-CoV-2-RBD_B.1.351/env/lib/python3.8/site-packages/neutcurve/hillcurve.py:451: RuntimeWarning: invalid value encountered in sqrt
    /fh/fast/bloom_j/computational_notebooks/agreaney/2021/SARS-CoV-2-RBD_B.1.351/env/lib/python3.8/site-packages/scipy/optimize/minpack.py:833: OptimizeWarning: Covariance of the parameters could not be estimated
    /fh/fast/bloom_j/computational_notebooks/agreaney/2021/SARS-CoV-2-RBD_B.1.351/env/lib/python3.8/site-packages/neutcurve/hillcurve.py:741: RuntimeWarning: invalid value encountered in power
    /fh/fast/bloom_j/computational_notebooks/agreaney/2021/SARS-CoV-2-RBD_B.1.351/env/lib/python3.8/site-packages/neutcurve/hillcurve.py:451: RuntimeWarning: invalid value encountered in sqrt
    /fh/fast/bloom_j/computational_notebooks/agreaney/2021/SARS-CoV-2-RBD_B.1.351/env/lib/python3.8/site-packages/neutcurve/hillcurve.py:741: RuntimeWarning: invalid value encountered in power
    /fh/fast/bloom_j/computational_notebooks/agreaney/2021/SARS-CoV-2-RBD_B.1.351/env/lib/python3.8/site-packages/neutcurve/hillcurve.py:741: RuntimeWarning: invalid value encountered in power
    /fh/fast/bloom_j/computational_notebooks/agreaney/2021/SARS-CoV-2-RBD_B.1.351/env/lib/python3.8/site-packages/scipy/optimize/minpack.py:833: OptimizeWarning: Covariance of the parameters could not be estimated



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
      <td>K041</td>
      <td>B.1.351-K484Q</td>
      <td>0.001090</td>
      <td>2021-08-27</td>
      <td>1</td>
      <td>0.000544</td>
      <td>2.002853</td>
      <td>0.001137</td>
      <td>0.000481</td>
      <td>2.361787</td>
    </tr>
    <tr>
      <th>1</th>
      <td>K041</td>
      <td>B.1.351-K484Q</td>
      <td>0.001401</td>
      <td>2021-08-27</td>
      <td>2</td>
      <td>0.000560</td>
      <td>2.501544</td>
      <td>0.001380</td>
      <td>0.000548</td>
      <td>2.516417</td>
    </tr>
    <tr>
      <th>2</th>
      <td>K041</td>
      <td>B.1.351-K484Q</td>
      <td>0.001168</td>
      <td>2021-08-27</td>
      <td>average</td>
      <td>0.000562</td>
      <td>2.078311</td>
      <td>0.001393</td>
      <td>0.000527</td>
      <td>2.643227</td>
    </tr>
    <tr>
      <th>3</th>
      <td>K041</td>
      <td>wildtype</td>
      <td>0.000544</td>
      <td>2021-08-27</td>
      <td>1</td>
      <td>0.000544</td>
      <td>1.000000</td>
      <td>0.000481</td>
      <td>0.000481</td>
      <td>1.000000</td>
    </tr>
    <tr>
      <th>4</th>
      <td>K041</td>
      <td>wildtype</td>
      <td>0.000560</td>
      <td>2021-08-27</td>
      <td>2</td>
      <td>0.000560</td>
      <td>1.000000</td>
      <td>0.000548</td>
      <td>0.000548</td>
      <td>1.000000</td>
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



    
![png](analyze_neut_data_files/analyze_neut_data_16_1.png)
    



    
![png](analyze_neut_data_files/analyze_neut_data_16_2.png)
    


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

    Saving to results/neut_titers/2021-08-27_mutant_neuts.pdf


    /fh/fast/bloom_j/computational_notebooks/agreaney/2021/SARS-CoV-2-RBD_B.1.351/env/lib/python3.8/site-packages/neutcurve/hillcurve.py:741: RuntimeWarning: invalid value encountered in power


    Saving to results/neut_titers/2021-06-10_mutant_neuts.pdf
    Saving to results/neut_titers/2021-08-20_mutant_neuts.pdf


    /fh/fast/bloom_j/computational_notebooks/agreaney/2021/SARS-CoV-2-RBD_B.1.351/env/lib/python3.8/site-packages/neutcurve/hillcurve.py:741: RuntimeWarning: invalid value encountered in power
    /fh/fast/bloom_j/computational_notebooks/agreaney/2021/SARS-CoV-2-RBD_B.1.351/env/lib/python3.8/site-packages/scipy/optimize/minpack.py:833: OptimizeWarning: Covariance of the parameters could not be estimated


    Saving to results/neut_titers/2021-08-26_mutant_neuts.pdf
    Saving to results/neut_titers/2021-08-21_mutant_neuts.pdf



    
![png](analyze_neut_data_files/analyze_neut_data_18_5.png)
    



    
![png](analyze_neut_data_files/analyze_neut_data_18_6.png)
    



    
![png](analyze_neut_data_files/analyze_neut_data_18_7.png)
    



    
![png](analyze_neut_data_files/analyze_neut_data_18_8.png)
    



    
![png](analyze_neut_data_files/analyze_neut_data_18_9.png)
    


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
      <td>K041</td>
      <td>B.1.351-K484Q</td>
      <td>0.001137</td>
      <td>879.582342</td>
      <td>interpolated</td>
      <td>2021-08-27</td>
      <td>1</td>
      <td>True</td>
      <td>False</td>
      <td>0.000514</td>
      <td>2.213101</td>
    </tr>
    <tr>
      <td>K041</td>
      <td>B.1.351-K484Q</td>
      <td>0.001380</td>
      <td>724.862206</td>
      <td>interpolated</td>
      <td>2021-08-27</td>
      <td>2</td>
      <td>True</td>
      <td>False</td>
      <td>0.000514</td>
      <td>2.685482</td>
    </tr>
    <tr>
      <td>K041</td>
      <td>B.1.351-K484Q</td>
      <td>0.001393</td>
      <td>717.838223</td>
      <td>interpolated</td>
      <td>2021-08-27</td>
      <td>average</td>
      <td>True</td>
      <td>False</td>
      <td>0.000514</td>
      <td>2.711759</td>
    </tr>
    <tr>
      <td>K041</td>
      <td>wildtype</td>
      <td>0.000481</td>
      <td>2077.386543</td>
      <td>interpolated</td>
      <td>2021-08-27</td>
      <td>1</td>
      <td>True</td>
      <td>False</td>
      <td>0.000514</td>
      <td>0.937045</td>
    </tr>
    <tr>
      <td>K041</td>
      <td>wildtype</td>
      <td>0.000548</td>
      <td>1824.055323</td>
      <td>interpolated</td>
      <td>2021-08-27</td>
      <td>2</td>
      <td>True</td>
      <td>False</td>
      <td>0.000514</td>
      <td>1.067185</td>
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
      <td>2015.434751</td>
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
      <td>4422.227100</td>
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
      <td>0.000346</td>
      <td>2888.596538</td>
      <td>interpolated</td>
      <td>2021-08-21</td>
      <td>1</td>
      <td>True</td>
      <td>False</td>
      <td>0.000474</td>
      <td>0.730742</td>
    </tr>
    <tr>
      <td>K046</td>
      <td>wildtype</td>
      <td>0.000648</td>
      <td>1542.465196</td>
      <td>interpolated</td>
      <td>2021-08-21</td>
      <td>2</td>
      <td>True</td>
      <td>False</td>
      <td>0.000474</td>
      <td>1.368471</td>
    </tr>
    <tr>
      <td>K119</td>
      <td>wildtype</td>
      <td>0.000266</td>
      <td>3761.207283</td>
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
      <td>6177.343337</td>
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


## Read in IC50s for early 2020 (HAARVI) plasmas in RBD depletion neuts


```python
haarvi_depletions=pd.DataFrame()

for i in ['1', '2', 'average']:
    df = (pd.read_csv(config['haarvi_rbd_depletions'])
                     [['serum', 'depletion', 'ic50', 'ic50_bound']]
                     .query('serum in @serum_order')
                     .rename(columns={'depletion': 'virus'})
                     .replace({'pre-depletion': 'wildtype', 'post-depletion': 'RBD antibodies depleted'})
                     .assign(NT50=lambda x: 1/x['ic50'],
                             date='October 2020',
                             replicate=i,
                             top=True,
                             ic50_is_bound=lambda x: x['ic50_bound'].map({'lower': True,'interpolated': False}),
                            )
                    )

    df= (df
         .merge((df.query('virus=="wildtype"')
                             [['serum', 'ic50']]
                             .rename(columns={'ic50': 'wildtype_ic50'})
                            ),
                            how='left',
                            on=['serum'],
                            validate='many_to_one'
                           )
                    .assign(fold_change=lambda x: x['ic50'] / x['wildtype_ic50'],)
                    )
    haarvi_depletions = pd.concat([haarvi_depletions, df])



display(HTML(haarvi_depletions.head().to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>serum</th>
      <th>virus</th>
      <th>ic50</th>
      <th>ic50_bound</th>
      <th>NT50</th>
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
      <td>participant A (day 21)</td>
      <td>wildtype</td>
      <td>0.000145</td>
      <td>interpolated</td>
      <td>6881.223206</td>
      <td>October 2020</td>
      <td>1</td>
      <td>True</td>
      <td>False</td>
      <td>0.000145</td>
      <td>1.000000</td>
    </tr>
    <tr>
      <td>participant A (day 21)</td>
      <td>RBD antibodies depleted</td>
      <td>0.001183</td>
      <td>interpolated</td>
      <td>845.146367</td>
      <td>October 2020</td>
      <td>1</td>
      <td>True</td>
      <td>False</td>
      <td>0.000145</td>
      <td>8.142049</td>
    </tr>
    <tr>
      <td>participant B (day 26)</td>
      <td>wildtype</td>
      <td>0.000640</td>
      <td>interpolated</td>
      <td>1563.438063</td>
      <td>October 2020</td>
      <td>1</td>
      <td>True</td>
      <td>False</td>
      <td>0.000640</td>
      <td>1.000000</td>
    </tr>
    <tr>
      <td>participant B (day 26)</td>
      <td>RBD antibodies depleted</td>
      <td>0.007346</td>
      <td>interpolated</td>
      <td>136.136752</td>
      <td>October 2020</td>
      <td>1</td>
      <td>True</td>
      <td>False</td>
      <td>0.000640</td>
      <td>11.484320</td>
    </tr>
    <tr>
      <td>participant C (day 32)</td>
      <td>wildtype</td>
      <td>0.000287</td>
      <td>interpolated</td>
      <td>3479.701163</td>
      <td>October 2020</td>
      <td>1</td>
      <td>True</td>
      <td>False</td>
      <td>0.000287</td>
      <td>1.000000</td>
    </tr>
  </tbody>
</table>



```python
neut_titers = pd.concat([neut_titers, haarvi_depletions]).assign(infecting_virus=lambda x: x['serum'].map(config['infecting_virus']))
display(HTML(neut_titers.tail().to_html(index=False)))
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
      <th>infecting_virus</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>participant B (day 26)</td>
      <td>RBD antibodies depleted</td>
      <td>0.007346</td>
      <td>136.136752</td>
      <td>interpolated</td>
      <td>October 2020</td>
      <td>average</td>
      <td>True</td>
      <td>False</td>
      <td>0.000640</td>
      <td>11.484320</td>
      <td>early 2020</td>
    </tr>
    <tr>
      <td>participant C (day 32)</td>
      <td>wildtype</td>
      <td>0.000287</td>
      <td>3479.701163</td>
      <td>interpolated</td>
      <td>October 2020</td>
      <td>average</td>
      <td>True</td>
      <td>False</td>
      <td>0.000287</td>
      <td>1.000000</td>
      <td>early 2020</td>
    </tr>
    <tr>
      <td>participant C (day 32)</td>
      <td>RBD antibodies depleted</td>
      <td>0.036028</td>
      <td>27.756167</td>
      <td>interpolated</td>
      <td>October 2020</td>
      <td>average</td>
      <td>True</td>
      <td>False</td>
      <td>0.000287</td>
      <td>125.366774</td>
      <td>early 2020</td>
    </tr>
    <tr>
      <td>participant I (day 26)</td>
      <td>wildtype</td>
      <td>0.007130</td>
      <td>140.247734</td>
      <td>interpolated</td>
      <td>October 2020</td>
      <td>average</td>
      <td>True</td>
      <td>False</td>
      <td>0.007130</td>
      <td>1.000000</td>
      <td>early 2020</td>
    </tr>
    <tr>
      <td>participant I (day 26)</td>
      <td>RBD antibodies depleted</td>
      <td>0.046926</td>
      <td>21.310096</td>
      <td>interpolated</td>
      <td>October 2020</td>
      <td>average</td>
      <td>True</td>
      <td>False</td>
      <td>0.007130</td>
      <td>6.581281</td>
      <td>early 2020</td>
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
         facet_wrap('~serum+infecting_virus', ncol=4, scales='free_x') +
         theme_classic() +
         theme(axis_text_x=element_text(angle=90),
               axis_title_x=element_blank(),
               strip_background_x=element_blank(),
               subplots_adjust={'hspace': 1.15},
               figure_size=(12, 12),
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



    
![png](analyze_neut_data_files/analyze_neut_data_27_1.png)
    



    
![png](analyze_neut_data_files/analyze_neut_data_27_2.png)
    



```python
neut_titers=(neut_titers.assign(virus_labels=lambda x: x['virus']
                                .map(config['virus_simplified_names']),
                                serum_virus=lambda x: x['serum'].astype(str)+'\n'+x['infecting_virus']
                               )
            )

neut_titers=(neut_titers.assign(virus_labels=lambda x: pd.Categorical(x['virus_labels'],ordered=True,categories=(config['virus_simplified_names_order'])),))
```


```python
neut_titers['serum_virus'].unique()
```




    array(['K041\nB.1.351', 'K046\nB.1.351', 'K114\nB.1.351', 'K119\nB.1.351',
           'K007\nB.1.351', 'K031\nB.1.351', 'K033\nB.1.351', 'K040\nB.1.351',
           'participant A (day 21)\nearly 2020',
           'participant C (day 32)\nearly 2020',
           'participant I (day 26)\nearly 2020',
           'participant B (day 26)\nearly 2020'], dtype=object)




```python
for metric in ['fold_change', 'ic50']:
    for virus_set, virus_subsample in config['virus_subsets'].items():
        print(f'Making plot for {metric} for {virus_set}:')
        
        ylab={'fold_change':'fold decrease in neutralization', 'ic50':'IC50'}

        dates=neut_titers.query("virus in @virus_subsample & virus not in ['wildtype', 'RBD antibodies depleted']")['date'].unique()

        p = (ggplot(neut_titers
                    .query("virus in @virus_subsample & date in @dates & replicate!= 'average'")
                    ) +
             aes('virus_labels', metric, shape='ic50_is_bound',
                ) +
             geom_point(aes(fill='date'), size=2.5, alpha=0.5, ) + #fill='#999999', 
             scale_y_log10(name=ylab[metric]) +
             facet_wrap('~serum_virus', ncol=4, scales='fixed') +
             theme_classic() +
             theme(axis_text_x=element_text(angle=90),
               axis_title_x=element_blank(),
               strip_background_x=element_blank(),
               figure_size=(9, 6),
               ) +
             geom_hline(data=(neut_titers
                              .query('virus in ["wildtype", "RBD antibodies depleted"] & replicate!="average"')
                              .groupby(['serum', 'virus', 'serum_virus'])
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
        p.save(plotfile, limitsize=False, verbose=False)
```

    Making plot for fold_change for all:


    /fh/fast/bloom_j/computational_notebooks/agreaney/2021/SARS-CoV-2-RBD_B.1.351/env/lib/python3.8/site-packages/plotnine/layer.py:401: PlotnineWarning: geom_hline : Removed 2568 rows containing missing values.


    Saving to results/neut_titers/fold_change_all.pdf


    /fh/fast/bloom_j/computational_notebooks/agreaney/2021/SARS-CoV-2-RBD_B.1.351/env/lib/python3.8/site-packages/plotnine/layer.py:401: PlotnineWarning: geom_hline : Removed 2568 rows containing missing values.


    Making plot for ic50 for all:


    /fh/fast/bloom_j/computational_notebooks/agreaney/2021/SARS-CoV-2-RBD_B.1.351/env/lib/python3.8/site-packages/plotnine/layer.py:401: PlotnineWarning: geom_hline : Removed 2568 rows containing missing values.


    Saving to results/neut_titers/ic50_all.pdf


    /fh/fast/bloom_j/computational_notebooks/agreaney/2021/SARS-CoV-2-RBD_B.1.351/env/lib/python3.8/site-packages/plotnine/layer.py:401: PlotnineWarning: geom_hline : Removed 2568 rows containing missing values.



    
![png](analyze_neut_data_files/analyze_neut_data_30_8.png)
    



    
![png](analyze_neut_data_files/analyze_neut_data_30_9.png)
    


Make plots without adding color by date:


```python
for metric in ['fold_change', 'ic50']:
    for virus_set, virus_subsample in config['virus_subsets'].items():
        print(f'Making plot for {metric} for {virus_set}:')
        
        ylab={'fold_change':'fold decrease in neutralization', 'ic50':'inhibitory concentration 50% (IC50)'}

        dates=neut_titers.query("virus in @virus_subsample & virus not in ['wildtype', 'RBD antibodies depleted']")['date'].unique()

        p = (ggplot(neut_titers
                    .query("virus in @virus_subsample & date in @dates & replicate!= 'average'")
                    ) +
             aes('virus_labels', metric, shape='ic50_is_bound',
                ) +
             geom_point(size=2, alpha=0.5, fill='#999999',) +  
             scale_y_log10(name=ylab[metric]) +
             facet_wrap('~serum_virus', ncol=4, scales='fixed') +
             theme_classic() +
             theme(axis_text_x=element_text(angle=90),
               axis_title_x=element_blank(),
               strip_background_x=element_blank(),
               figure_size=(9, 6),
               ) +
             geom_hline(data=(neut_titers
                              .query('virus in ["wildtype", "RBD antibodies depleted"] & replicate!="average"')
                              .groupby(['serum', 'virus', 'serum_virus'])
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
             scale_color_manual(values=CBPALETTE*3, guide=False) 
             )

        _ = p.draw()

        plotfile = f'{results}/{metric}_{virus_set}.pdf'
        print(f"Saving to {plotfile}")
        p.save(plotfile, limitsize=False, verbose=False)
```

    Making plot for fold_change for all:


    /fh/fast/bloom_j/computational_notebooks/agreaney/2021/SARS-CoV-2-RBD_B.1.351/env/lib/python3.8/site-packages/plotnine/layer.py:401: PlotnineWarning: geom_hline : Removed 2568 rows containing missing values.


    Saving to results/neut_titers/fold_change_all.pdf


    /fh/fast/bloom_j/computational_notebooks/agreaney/2021/SARS-CoV-2-RBD_B.1.351/env/lib/python3.8/site-packages/plotnine/layer.py:401: PlotnineWarning: geom_hline : Removed 2568 rows containing missing values.


    Making plot for ic50 for all:


    /fh/fast/bloom_j/computational_notebooks/agreaney/2021/SARS-CoV-2-RBD_B.1.351/env/lib/python3.8/site-packages/plotnine/layer.py:401: PlotnineWarning: geom_hline : Removed 2568 rows containing missing values.


    Saving to results/neut_titers/ic50_all.pdf


    /fh/fast/bloom_j/computational_notebooks/agreaney/2021/SARS-CoV-2-RBD_B.1.351/env/lib/python3.8/site-packages/plotnine/layer.py:401: PlotnineWarning: geom_hline : Removed 2568 rows containing missing values.



    
![png](analyze_neut_data_files/analyze_neut_data_32_8.png)
    



    
![png](analyze_neut_data_files/analyze_neut_data_32_9.png)
    


### Plot the results for all individuals on one plot and add the geometric mean of all 8 measurements (2 replicates * 4 sera) for each virus


```python
for metric in ['fold_change', 'ic50']:
    for virus_set, virus_subsample in config['virus_subsets'].items():
        print(f'Making plot for {metric} for {virus_set}:')
        
        ylab={'fold_change':'fold decrease in neutralization', 'ic50':'IC50'}

        dates=neut_titers.query("virus in @virus_subsample & virus not in ['wildtype']")['date'].unique()

        p = (ggplot(neut_titers
                    .query("virus in @virus_subsample & date in @dates & replicate== 'average'")
                    ) +
             aes('virus', metric, shape='ic50_is_bound', fill='serum',
                ) +
             geom_point(size=2.5, alpha=0.5, ) + 
             geom_crossbar(data=(neut_titers
                                 .query("virus in @virus_subsample & date in @dates & replicate== 'average'")
                                 .groupby(['virus', 'infecting_virus'])
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
                   figure_size=(6, 2.5),
                   ) +
             geom_hline(data=(neut_titers
                              .query('virus in ["wildtype"] & replicate=="average"')
                              .groupby(['virus', 'infecting_virus'])
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
             scale_fill_manual(values=CBPALETTE*3) +
             facet_wrap('~infecting_virus', scales='free_x')
             )

        _ = p.draw()

        plotfile = f'{results}/{metric}_{virus_set}_aggregate.pdf'
        print(f"Saving to {plotfile}")
        p.save(plotfile, limitsize=False, verbose=False)
```

    Making plot for fold_change for all:


    /fh/fast/bloom_j/computational_notebooks/agreaney/2021/SARS-CoV-2-RBD_B.1.351/env/lib/python3.8/site-packages/plotnine/layer.py:401: PlotnineWarning: geom_hline : Removed 34 rows containing missing values.


    Saving to results/neut_titers/fold_change_all_aggregate.pdf


    /fh/fast/bloom_j/computational_notebooks/agreaney/2021/SARS-CoV-2-RBD_B.1.351/env/lib/python3.8/site-packages/plotnine/layer.py:401: PlotnineWarning: geom_hline : Removed 34 rows containing missing values.


    Making plot for ic50 for all:


    /fh/fast/bloom_j/computational_notebooks/agreaney/2021/SARS-CoV-2-RBD_B.1.351/env/lib/python3.8/site-packages/plotnine/layer.py:401: PlotnineWarning: geom_hline : Removed 34 rows containing missing values.


    Saving to results/neut_titers/ic50_all_aggregate.pdf


    /fh/fast/bloom_j/computational_notebooks/agreaney/2021/SARS-CoV-2-RBD_B.1.351/env/lib/python3.8/site-packages/plotnine/layer.py:401: PlotnineWarning: geom_hline : Removed 34 rows containing missing values.



    
![png](analyze_neut_data_files/analyze_neut_data_34_8.png)
    



    
![png](analyze_neut_data_files/analyze_neut_data_34_9.png)
    



```python
for metric in ['fold_change', 'ic50']:
    for virus_set, virus_subsample in config['virus_subsets'].items():
        print(f'Making plot for {metric} for {virus_set}:')
        
        ylab={'fold_change':'fold decrease in neutralization', 'ic50':'IC50'}

        dates=neut_titers.query("virus in @virus_subsample & virus not in ['wildtype']")['date'].unique()

        p = (ggplot(neut_titers
                    .query("virus in @virus_subsample & date in @dates & replicate== 'average'")
                    ) +
             aes('virus', metric, shape='ic50_is_bound',
                ) +
             geom_point(size=2.5, alpha=0.5, ) + 
             geom_crossbar(data=(neut_titers
                                 .query("virus in @virus_subsample & date in @dates & replicate== 'average'")
                                 .groupby(['virus', 'infecting_virus'])
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
                   figure_size=(6, 2.5),
                   ) +
             geom_hline(data=(neut_titers
                              .query('virus in ["wildtype"] & replicate=="average"')
                              .groupby(['virus', 'infecting_virus'])
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
             scale_fill_manual(values=CBPALETTE*3) +
             facet_wrap('~infecting_virus', scales='free_x')
             )

        _ = p.draw()

        plotfile = f'{results}/{metric}_{virus_set}_aggregate_nocolors.pdf'
        print(f"Saving to {plotfile}")
        p.save(plotfile, limitsize=False, verbose=False)
```

    Making plot for fold_change for all:


    /fh/fast/bloom_j/computational_notebooks/agreaney/2021/SARS-CoV-2-RBD_B.1.351/env/lib/python3.8/site-packages/plotnine/layer.py:401: PlotnineWarning: geom_hline : Removed 34 rows containing missing values.
    /fh/fast/bloom_j/computational_notebooks/agreaney/2021/SARS-CoV-2-RBD_B.1.351/env/lib/python3.8/site-packages/plotnine/guides/guides.py:197: PlotnineWarning: Cannot generate legend for the 'fill' aesthetic. Make sure you have mapped a variable to it


    Saving to results/neut_titers/fold_change_all_aggregate_nocolors.pdf


    /fh/fast/bloom_j/computational_notebooks/agreaney/2021/SARS-CoV-2-RBD_B.1.351/env/lib/python3.8/site-packages/plotnine/layer.py:401: PlotnineWarning: geom_hline : Removed 34 rows containing missing values.
    /fh/fast/bloom_j/computational_notebooks/agreaney/2021/SARS-CoV-2-RBD_B.1.351/env/lib/python3.8/site-packages/plotnine/guides/guides.py:197: PlotnineWarning: Cannot generate legend for the 'fill' aesthetic. Make sure you have mapped a variable to it


    Making plot for ic50 for all:


    /fh/fast/bloom_j/computational_notebooks/agreaney/2021/SARS-CoV-2-RBD_B.1.351/env/lib/python3.8/site-packages/plotnine/layer.py:401: PlotnineWarning: geom_hline : Removed 34 rows containing missing values.
    /fh/fast/bloom_j/computational_notebooks/agreaney/2021/SARS-CoV-2-RBD_B.1.351/env/lib/python3.8/site-packages/plotnine/guides/guides.py:197: PlotnineWarning: Cannot generate legend for the 'fill' aesthetic. Make sure you have mapped a variable to it


    Saving to results/neut_titers/ic50_all_aggregate_nocolors.pdf


    /fh/fast/bloom_j/computational_notebooks/agreaney/2021/SARS-CoV-2-RBD_B.1.351/env/lib/python3.8/site-packages/plotnine/layer.py:401: PlotnineWarning: geom_hline : Removed 34 rows containing missing values.
    /fh/fast/bloom_j/computational_notebooks/agreaney/2021/SARS-CoV-2-RBD_B.1.351/env/lib/python3.8/site-packages/plotnine/guides/guides.py:197: PlotnineWarning: Cannot generate legend for the 'fill' aesthetic. Make sure you have mapped a variable to it



    
![png](analyze_neut_data_files/analyze_neut_data_35_8.png)
    



    
![png](analyze_neut_data_files/analyze_neut_data_35_9.png)
    


## Rename viruses such that they are comparable between the B.1.351 and early 2020 samples. 
Essentially, I want the virus names to be like this:
* wildtype: wildtype
* B.1.351: wildtype
* mock: wildtype
* B.1.351-N417K: 417K/N
* B.1.351-K484E: 484K/E
* B.1.351-K484Q: 484Q
* B.1.351-L452R: L452R
* B.1.351-G446V: G446V
* B.1.351-K444E: K444E
* B.1.351-Y501N: 501Y/N
* B.1.351-N417K-K484E-Y501N: 417-484-501
* K417N: 417K/N
* E484K: 484K/E
* E484Q: 484Q
* N501Y: 501Y/N
* G446V: G446V
* K417N-E484K-N501Y: 417-484-501
* RBD antibodies depleted: RBD antibodies depleted


```python
neut_titers.query('virus=="wildtype" & infecting_virus=="B.1.351"')['ic50'].agg(geometric_mean)
```




    0.00038486772827164123




```python
neut_titers=(neut_titers.assign(virus_simplified_names=lambda x: x['virus']
                                .map(config['virus_simplified_names']),
                                virus_labels=lambda x: x['virus_simplified_names']
                                .replace({'RBD antibodies depleted':'RBD\nantibodies\ndepleted',
                                          '417-484-501':'417\n484\n501'
                                         }
                                        ),
                               )
            )

display(HTML(neut_titers.head().to_html()))

print(neut_titers['ic50_is_bound'].unique())

for metric in ['fold_change', 'ic50']:
    for virus_set, virus_subsample in config['virus_subsets'].items():
        print(f'Making plot for {metric} for {virus_set}:')
        
        if metric=="fold_change":
            virus_subsample=[v for v in virus_subsample if v!="wildtype"]
        
        print(f"Making plot for {neut_titers.query('virus in @virus_subsample')['virus_labels'].nunique()} viruses")
        
        ylab={'fold_change':'fold decrease in neutralization', 
              'ic50':'inhibitory concentration\n50% (IC50)'
             }
        yintercept={'fold_change':1, 
                    'ic50':(neut_titers
                            .query('virus=="wildtype" & infecting_virus=="B.1.351"')
                            ['ic50']
                            .agg(geometric_mean)
                           )
                   }

        p = (ggplot(neut_titers
                    .query("virus in @virus_subsample & replicate=='average'")
                    .assign(virus_labels=lambda x: pd.Categorical(x['virus_labels'],
                                                                  ordered=True,
                                                                  categories=(config['virus_simplified_names_order']+
                                                                              ['417\n484\n501',
                                                                               'RBD\nantibodies\ndepleted'])),
                           )
                    ) +
             aes('virus_labels', 
                 metric, 
                 fill='infecting_virus', 
                 color='infecting_virus',
                ) +
             geom_point(position=position_dodge(width=0.55), size=2.5, alpha=0.5) +
             geom_crossbar(data=(neut_titers
                                 .query("virus in @virus_subsample & replicate=='average'")
                                 .groupby(['virus_labels', 'infecting_virus'])
                                 .agg({metric: geometric_mean})
                                 .reset_index()
                                 .dropna()
                                ),
#                            inherit_aes=False,
                           mapping=aes(x='virus_labels', y=metric, ymin=metric, ymax=metric),
                           position=position_dodge(width=0.55),
                  ) +
             geom_hline(yintercept=yintercept[metric],
                        linetype='dashed', size=0.5,
                        alpha=0.6, 
                        color=CBPALETTE[0]) +
             scale_y_log10(name=ylab[metric]) +
             theme_classic() +
             theme(axis_title_x=element_blank(),
                   figure_size=(neut_titers.query('virus in @virus_subsample')['virus_labels'].nunique()*0.75, 2.5), 
                   ) +
             scale_fill_manual(values=['#44AA99', '#332288'], name='infecting virus\n')+
             scale_color_manual(values=['#44AA99', '#332288'], name='infecting virus\n')
             )

        _ = p.draw()

        plotfile = f'{results}/{metric}_{virus_set}_aggregate_nofacet.pdf'
        print(f"Saving to {plotfile}")
        p.save(plotfile, limitsize=False, verbose=False)
```


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
      <th>wildtype_ic50</th>
      <th>fold_change</th>
      <th>infecting_virus</th>
      <th>virus_labels</th>
      <th>serum_virus</th>
      <th>virus_simplified_names</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>K041</td>
      <td>B.1.351-K484Q</td>
      <td>0.001137</td>
      <td>879.582342</td>
      <td>interpolated</td>
      <td>2021-08-27</td>
      <td>1</td>
      <td>True</td>
      <td>False</td>
      <td>0.000514</td>
      <td>2.213101</td>
      <td>B.1.351</td>
      <td>484Q</td>
      <td>K041\nB.1.351</td>
      <td>484Q</td>
    </tr>
    <tr>
      <th>1</th>
      <td>K041</td>
      <td>B.1.351-K484Q</td>
      <td>0.001380</td>
      <td>724.862206</td>
      <td>interpolated</td>
      <td>2021-08-27</td>
      <td>2</td>
      <td>True</td>
      <td>False</td>
      <td>0.000514</td>
      <td>2.685482</td>
      <td>B.1.351</td>
      <td>484Q</td>
      <td>K041\nB.1.351</td>
      <td>484Q</td>
    </tr>
    <tr>
      <th>2</th>
      <td>K041</td>
      <td>B.1.351-K484Q</td>
      <td>0.001393</td>
      <td>717.838223</td>
      <td>interpolated</td>
      <td>2021-08-27</td>
      <td>average</td>
      <td>True</td>
      <td>False</td>
      <td>0.000514</td>
      <td>2.711759</td>
      <td>B.1.351</td>
      <td>484Q</td>
      <td>K041\nB.1.351</td>
      <td>484Q</td>
    </tr>
    <tr>
      <th>3</th>
      <td>K041</td>
      <td>wildtype</td>
      <td>0.000481</td>
      <td>2077.386543</td>
      <td>interpolated</td>
      <td>2021-08-27</td>
      <td>1</td>
      <td>True</td>
      <td>False</td>
      <td>0.000514</td>
      <td>0.937045</td>
      <td>B.1.351</td>
      <td>wildtype</td>
      <td>K041\nB.1.351</td>
      <td>wildtype</td>
    </tr>
    <tr>
      <th>4</th>
      <td>K041</td>
      <td>wildtype</td>
      <td>0.000548</td>
      <td>1824.055323</td>
      <td>interpolated</td>
      <td>2021-08-27</td>
      <td>2</td>
      <td>True</td>
      <td>False</td>
      <td>0.000514</td>
      <td>1.067185</td>
      <td>B.1.351</td>
      <td>wildtype</td>
      <td>K041\nB.1.351</td>
      <td>wildtype</td>
    </tr>
  </tbody>
</table>


    [False True]
    Making plot for fold_change for all:
    Making plot for 8 viruses
    Saving to results/neut_titers/fold_change_all_aggregate_nofacet.pdf
    Making plot for ic50 for all:
    Making plot for 9 viruses
    Saving to results/neut_titers/ic50_all_aggregate_nofacet.pdf



    
![png](analyze_neut_data_files/analyze_neut_data_38_2.png)
    



    
![png](analyze_neut_data_files/analyze_neut_data_38_3.png)
    



```python

```
