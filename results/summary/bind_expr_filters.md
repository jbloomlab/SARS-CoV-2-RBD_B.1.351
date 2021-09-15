# Set B.1.351 RBD DMS ACE2 binding and expression scores for serum mapping
We want to make sure that the filters chosen for the ACE2 binding and RBD expression scores are reasonable such that spurious antibody-escpae mutations that merely fall into the antibody-escape gate due to their poor folding or expression are removed. 

But, we also want to make sure we aren't throwing out many mutations that are found in nature at reasonable numbers. 


```python
import os

from IPython.display import display, HTML

import math
import numpy as np
import pandas as pd
from scipy import stats

from plotnine import *

from dms_variants.constants import CBPALETTE

import yaml
```

Read config file


```python
with open('config.yaml') as f:
    config = yaml.safe_load(f)
```

Define input and output directories


```python
datadir = 'data'
resultsdir = config['bind_expr_filters_dir']

os.makedirs(resultsdir, exist_ok=True)
```

Read in the new filters for DMS ACE2 binding and expression scores. 


```python
og_thresholds={'delta_bind':-2.35, 'delta_expr':-1.0}
new_thresholds={'delta_bind':config['escape_score_min_bind_mut'], 'delta_expr':config['escape_score_min_expr_mut']}

og_thresholds_df=pd.DataFrame.from_dict({'metric': ['delta_bind', 'delta_expr'], 'score': [-2.35,-1.0]})
new_filter_df=pd.DataFrame({'metric': ['delta_bind', 'delta_expr'], 'score':[config['escape_score_min_bind_mut'],config['escape_score_min_expr_mut']]})
display(HTML(new_filter_df.to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>metric</th>
      <th>score</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>delta_bind</td>
      <td>-3.0</td>
    </tr>
    <tr>
      <td>delta_expr</td>
      <td>-1.0</td>
    </tr>
  </tbody>
</table>



```python
gisaid_counts_file = config['gisaid_mutation_counts']
dms_scores_file = config['prelim_mut_bind_expr']
og_dms_file = config['early2020_mut_bind_expr']
```


```python
dms_scores = (pd.read_csv(dms_scores_file).rename(columns={'position': 'site'})
             [['wildtype', 'mutation', 'site', 'mutant', 'delta_bind', 'delta_expr']]
             )

display(HTML(dms_scores.head().to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>wildtype</th>
      <th>mutation</th>
      <th>site</th>
      <th>mutant</th>
      <th>delta_bind</th>
      <th>delta_expr</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>N</td>
      <td>N331A</td>
      <td>331</td>
      <td>A</td>
      <td>0.12426</td>
      <td>-0.27844</td>
    </tr>
    <tr>
      <td>N</td>
      <td>N331C</td>
      <td>331</td>
      <td>C</td>
      <td>-0.20203</td>
      <td>-0.63606</td>
    </tr>
    <tr>
      <td>N</td>
      <td>N331D</td>
      <td>331</td>
      <td>D</td>
      <td>0.00131</td>
      <td>-0.22296</td>
    </tr>
    <tr>
      <td>N</td>
      <td>N331E</td>
      <td>331</td>
      <td>E</td>
      <td>0.10980</td>
      <td>0.06430</td>
    </tr>
    <tr>
      <td>N</td>
      <td>N331F</td>
      <td>331</td>
      <td>F</td>
      <td>-0.15674</td>
      <td>-0.57560</td>
    </tr>
  </tbody>
</table>



```python
gisaid_counts = (pd.read_csv(gisaid_counts_file)
                 .drop(columns=['isite', 'wildtype'])
                )

dms_scores=(dms_scores
            .merge(gisaid_counts,
                   on=['site', 'mutant'],
                   how='left',
                   validate='many_to_one',
                  )
            .fillna({'count':0,'n_countries':0, 'frequency': 0})
           )

dms_scores=dms_scores.melt(id_vars=['wildtype','mutation', 'site', 'mutant', 'count', 'n_countries', 'frequency'],
                           value_vars=['delta_bind', 'delta_expr'], 
                           var_name='metric', 
                           value_name='score',
                          )

display(HTML(dms_scores.head().to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>wildtype</th>
      <th>mutation</th>
      <th>site</th>
      <th>mutant</th>
      <th>count</th>
      <th>n_countries</th>
      <th>frequency</th>
      <th>metric</th>
      <th>score</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>N</td>
      <td>N331A</td>
      <td>331</td>
      <td>A</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.000000e+00</td>
      <td>delta_bind</td>
      <td>0.12426</td>
    </tr>
    <tr>
      <td>N</td>
      <td>N331C</td>
      <td>331</td>
      <td>C</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.000000e+00</td>
      <td>delta_bind</td>
      <td>-0.20203</td>
    </tr>
    <tr>
      <td>N</td>
      <td>N331D</td>
      <td>331</td>
      <td>D</td>
      <td>2.0</td>
      <td>2.0</td>
      <td>8.673362e-07</td>
      <td>delta_bind</td>
      <td>0.00131</td>
    </tr>
    <tr>
      <td>N</td>
      <td>N331E</td>
      <td>331</td>
      <td>E</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.000000e+00</td>
      <td>delta_bind</td>
      <td>0.10980</td>
    </tr>
    <tr>
      <td>N</td>
      <td>N331F</td>
      <td>331</td>
      <td>F</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.000000e+00</td>
      <td>delta_bind</td>
      <td>-0.15674</td>
    </tr>
  </tbody>
</table>



```python
p = (ggplot(dms_scores
            # assign small numbers to things with 0 GISAID counts or missing scores so they still appear on plot 
            .replace({'count': {0: 0.1}, 'score': {np.nan: -5}})
            .replace({'delta_bind':'ACE2 binding', 'delta_expr':'RBD expression'})
           ) +
     aes('count', 'score') +
     geom_point(alpha=0.2, color='black') +
     facet_grid('~ metric') +
     scale_x_log10()+
     theme_classic() +
     geom_hline(data=new_filter_df.replace({'delta_bind':'ACE2 binding', 'delta_expr':'RBD expression'}),
                 mapping=aes(yintercept='score'),
                linetype='dashed',
                color=CBPALETTE[1])+
     theme(figure_size=(2.5 * 2, 2.5 * 1),
           strip_background=element_blank(),
           strip_text=element_text(size=12),
          ) +
     xlab('mutation counts in GISAID as of Aug. 1, 2021')+
     ylab('B.1.351 RBD DMS score\n(single mutants)')
     )

fig = p.draw()

plotfile = os.path.join(resultsdir, f"counts-v-score.pdf")
print(f"Saving plot to {plotfile}")
p.save(plotfile, verbose=False)
```

    Saving plot to results/bind_expr_filters/counts-v-score.pdf



    
![png](bind_expr_filters_files/bind_expr_filters_11_1.png)
    



```python
def assign_count_categories(x):
    if x == 0:
        return "0"
    elif x < 10:
        return "1 to 9"
    elif x < 20:
        return "10 to 19"
    elif x < 50:
        return "20 to 49"
    else:
        return ">=50"
    
count_categories=["0", "1 to 9", "10 to 19", "20 to 49", ">=50"]

dms_scores=(dms_scores
            .assign(count_categories=lambda x: x['count'].apply(assign_count_categories),
                   )
           )

dms_scores=(dms_scores
            .assign(count_categories=lambda x: pd.Categorical(x['count_categories'],
                                                              categories=count_categories,
                                                              ordered=True
                                                             ))
           )

p = (ggplot(dms_scores.replace({'delta_bind':'ACE2 binding', 'delta_expr':'RBD expression'})) +
     aes('count_categories', 'score') +
     geom_hline(data=new_filter_df.replace({'delta_bind':'ACE2 binding', 'delta_expr':'RBD expression'}),
                 mapping=aes(yintercept='score'),
                linetype='dashed',
                color=CBPALETTE[1])+
     geom_boxplot(outlier_alpha=0.2) +
     facet_grid('~ metric') +
     theme_classic() +
     theme(figure_size=(2.5 * 2, 2.5 * 1),
           axis_text_x=element_text(angle=90),
           strip_background=element_blank(),
           strip_text=element_text(size=12),
          ) +
     xlab('mutation counts in GISAID as of Aug. 1 2021')+
     ylab('B.1.351 RBD DMS score')
     )

fig = p.draw()

plotfile = os.path.join(resultsdir, f"count-cat-v-score.pdf")
print(f"Saving plot to {plotfile}")
p.save(plotfile, verbose=False)
```

    /fh/fast/bloom_j/computational_notebooks/agreaney/2021/SARS-CoV-2-RBD_B.1.351/env/lib/python3.8/site-packages/plotnine/layer.py:324: PlotnineWarning: stat_boxplot : Removed 34 rows containing non-finite values.


    Saving plot to results/bind_expr_filters/count-cat-v-score.pdf


    /fh/fast/bloom_j/computational_notebooks/agreaney/2021/SARS-CoV-2-RBD_B.1.351/env/lib/python3.8/site-packages/plotnine/layer.py:324: PlotnineWarning: stat_boxplot : Removed 34 rows containing non-finite values.



    
![png](bind_expr_filters_files/bind_expr_filters_12_3.png)
    



```python
x_min=-4.5
x_max=0.5

p = (ggplot(dms_scores.replace({'delta_bind':'ACE2 binding', 'delta_expr':'RBD expression'})) +
     aes(x='score', fill='count_categories') +
     geom_histogram(position='identity', bins=50) +
     facet_grid('~ metric') +
     scale_x_continuous(breaks=np.arange(x_min,x_max,0.5), limits=[x_min, x_max]) +
     geom_vline(data=new_filter_df.replace({'delta_bind':'ACE2 binding', 'delta_expr':'RBD expression'}),
                     mapping=aes(xintercept='score'),
                    linetype='dashed',
                    color=CBPALETTE[1])+
     theme_classic() +
     theme(figure_size=(2.5 * 2, 2.5 * 1),
           plot_title=element_text(size=14),
           axis_text_x=element_text(angle=90),
           strip_background=element_blank(),
           strip_text=element_text(size=12),
          ) +
     ylab('number of mutations')+
     xlab('B.1.351 RBD DMS score') +
     labs(fill='GISAID counts')
     )

fig = p.draw()

plotfile = os.path.join(resultsdir, f"count-score-histogram.pdf")
print(f"Saving plot to {plotfile}")
p.save(plotfile, verbose=False)
```

    /fh/fast/bloom_j/computational_notebooks/agreaney/2021/SARS-CoV-2-RBD_B.1.351/env/lib/python3.8/site-packages/plotnine/layer.py:324: PlotnineWarning: stat_bin : Removed 55 rows containing non-finite values.
    /fh/fast/bloom_j/computational_notebooks/agreaney/2021/SARS-CoV-2-RBD_B.1.351/env/lib/python3.8/site-packages/plotnine/layer.py:401: PlotnineWarning: geom_histogram : Removed 20 rows containing missing values.


    Saving plot to results/bind_expr_filters/count-score-histogram.pdf


    /fh/fast/bloom_j/computational_notebooks/agreaney/2021/SARS-CoV-2-RBD_B.1.351/env/lib/python3.8/site-packages/plotnine/layer.py:324: PlotnineWarning: stat_bin : Removed 55 rows containing non-finite values.
    /fh/fast/bloom_j/computational_notebooks/agreaney/2021/SARS-CoV-2-RBD_B.1.351/env/lib/python3.8/site-packages/plotnine/layer.py:401: PlotnineWarning: geom_histogram : Removed 20 rows containing missing values.



    
![png](bind_expr_filters_files/bind_expr_filters_13_3.png)
    



```python
p = (ggplot(dms_scores
            .replace({'delta_bind':'ACE2 binding', 'delta_expr':'RBD expression'})
            .dropna()
            .drop(columns=['count'])
            .assign(count_categories=lambda x: pd.Categorical(x['count_categories']).remove_unused_categories())
           ) +
     aes(x='score', fill='count_categories') +
     geom_density(alpha=0.5) +#aes(y=after_stat('count')), 
     scale_x_continuous(breaks=np.arange(x_min,x_max,0.5), limits=[x_min, x_max]) +
     geom_vline(data=new_filter_df.replace({'delta_bind':'ACE2 binding', 'delta_expr':'RBD expression'}),
                     mapping=aes(xintercept='score'),
                    linetype='dashed',
                    color=CBPALETTE[1])+
     facet_grid('~ metric') +
     theme_classic() +
     theme(figure_size=(2.5 * 2, 2.5 * 1),
           plot_title=element_text(size=14),
           axis_text_x=element_text(angle=90),
           strip_background=element_blank(),
           strip_text=element_text(size=12),
          ) +
     labs(fill = 'GISAID counts')+
     xlab('B.1.351 RBD DMS score')
     )

fig = p.draw()

plotfile = os.path.join(resultsdir, f"count-score-density.pdf")
print(f"Saving plot to {plotfile}")
p.save(plotfile, verbose=False)
```

    /fh/fast/bloom_j/computational_notebooks/agreaney/2021/SARS-CoV-2-RBD_B.1.351/env/lib/python3.8/site-packages/plotnine/layer.py:324: PlotnineWarning: stat_density : Removed 21 rows containing non-finite values.


    Saving plot to results/bind_expr_filters/count-score-density.pdf


    /fh/fast/bloom_j/computational_notebooks/agreaney/2021/SARS-CoV-2-RBD_B.1.351/env/lib/python3.8/site-packages/plotnine/layer.py:324: PlotnineWarning: stat_density : Removed 21 rows containing non-finite values.



    
![png](bind_expr_filters_files/bind_expr_filters_14_3.png)
    


Things I want to know:
1. Mutations that have **any** counts in nature but are missing scores
2. Mutations that have appreciable counts (>=50) in nature but very low scores
3. The scores corresponding to the 95th percentile of all mutations occurring >= 50x in nature
4. The scores of mutations to disulfide bonds


```python
print('Here are the naturally occurring mutations that are missing scores from B.1.351 DMS')
display(HTML(dms_scores
             .query('count >= 1')
             .query('score.isnull()', engine='python')
             [['wildtype','mutation', 'count', 'n_countries', 'frequency', 'score']]
             .drop_duplicates()
             .to_html(index=False)
            )
       )
```

    Here are the naturally occurring mutations that are missing scores from B.1.351 DMS



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>wildtype</th>
      <th>mutation</th>
      <th>count</th>
      <th>n_countries</th>
      <th>frequency</th>
      <th>score</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>V</td>
      <td>V341F</td>
      <td>3.0</td>
      <td>3.0</td>
      <td>1.301004e-06</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>V</td>
      <td>V483L</td>
      <td>33.0</td>
      <td>13.0</td>
      <td>1.431105e-05</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>K</td>
      <td>K529I</td>
      <td>1.0</td>
      <td>1.0</td>
      <td>4.336681e-07</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>G</td>
      <td>G526A</td>
      <td>11.0</td>
      <td>3.0</td>
      <td>4.770349e-06</td>
      <td>NaN</td>
    </tr>
  </tbody>
</table>



```python
for metric in ['bind', 'expr']:
    m=f"delta_{metric}"
    score_filter=new_thresholds[m]
    print(f'Mutations with >=50 GISAID counts but with {metric} score < {score_filter}')
    display(HTML(dms_scores
                 .query('metric==@m & count >= 50 & score < @score_filter')
                 .drop_duplicates()
                 .sort_values(by='score')
                 .head(20)
                 .to_html(index=False)
                )
           )
```

    Mutations with >=50 GISAID counts but with bind score < -3.0



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>wildtype</th>
      <th>mutation</th>
      <th>site</th>
      <th>mutant</th>
      <th>count</th>
      <th>n_countries</th>
      <th>frequency</th>
      <th>metric</th>
      <th>score</th>
      <th>count_categories</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>A</td>
      <td>A419S</td>
      <td>419</td>
      <td>S</td>
      <td>134.0</td>
      <td>23.0</td>
      <td>0.000058</td>
      <td>delta_bind</td>
      <td>-4.07656</td>
      <td>&gt;=50</td>
    </tr>
  </tbody>
</table>


    Mutations with >=50 GISAID counts but with expr score < -1.0



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>wildtype</th>
      <th>mutation</th>
      <th>site</th>
      <th>mutant</th>
      <th>count</th>
      <th>n_countries</th>
      <th>frequency</th>
      <th>metric</th>
      <th>score</th>
      <th>count_categories</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>T</td>
      <td>T470N</td>
      <td>470</td>
      <td>N</td>
      <td>189.0</td>
      <td>23.0</td>
      <td>0.000082</td>
      <td>delta_expr</td>
      <td>-2.68080</td>
      <td>&gt;=50</td>
    </tr>
    <tr>
      <td>V</td>
      <td>V401L</td>
      <td>401</td>
      <td>L</td>
      <td>107.0</td>
      <td>21.0</td>
      <td>0.000046</td>
      <td>delta_expr</td>
      <td>-2.01452</td>
      <td>&gt;=50</td>
    </tr>
    <tr>
      <td>G</td>
      <td>G404C</td>
      <td>404</td>
      <td>C</td>
      <td>55.0</td>
      <td>7.0</td>
      <td>0.000024</td>
      <td>delta_expr</td>
      <td>-1.98436</td>
      <td>&gt;=50</td>
    </tr>
    <tr>
      <td>A</td>
      <td>A419S</td>
      <td>419</td>
      <td>S</td>
      <td>134.0</td>
      <td>23.0</td>
      <td>0.000058</td>
      <td>delta_expr</td>
      <td>-1.77949</td>
      <td>&gt;=50</td>
    </tr>
    <tr>
      <td>A</td>
      <td>A352V</td>
      <td>352</td>
      <td>V</td>
      <td>79.0</td>
      <td>16.0</td>
      <td>0.000034</td>
      <td>delta_expr</td>
      <td>-1.62895</td>
      <td>&gt;=50</td>
    </tr>
    <tr>
      <td>G</td>
      <td>G485V</td>
      <td>485</td>
      <td>V</td>
      <td>60.0</td>
      <td>13.0</td>
      <td>0.000026</td>
      <td>delta_expr</td>
      <td>-1.24302</td>
      <td>&gt;=50</td>
    </tr>
    <tr>
      <td>P</td>
      <td>P463S</td>
      <td>463</td>
      <td>S</td>
      <td>142.0</td>
      <td>20.0</td>
      <td>0.000062</td>
      <td>delta_expr</td>
      <td>-1.24240</td>
      <td>&gt;=50</td>
    </tr>
    <tr>
      <td>S</td>
      <td>S494L</td>
      <td>494</td>
      <td>L</td>
      <td>144.0</td>
      <td>31.0</td>
      <td>0.000062</td>
      <td>delta_expr</td>
      <td>-1.15532</td>
      <td>&gt;=50</td>
    </tr>
    <tr>
      <td>Y</td>
      <td>Y508C</td>
      <td>508</td>
      <td>C</td>
      <td>51.0</td>
      <td>3.0</td>
      <td>0.000022</td>
      <td>delta_expr</td>
      <td>-1.05878</td>
      <td>&gt;=50</td>
    </tr>
    <tr>
      <td>Q</td>
      <td>Q414H</td>
      <td>414</td>
      <td>H</td>
      <td>308.0</td>
      <td>4.0</td>
      <td>0.000134</td>
      <td>delta_expr</td>
      <td>-1.05292</td>
      <td>&gt;=50</td>
    </tr>
    <tr>
      <td>A</td>
      <td>A348T</td>
      <td>348</td>
      <td>T</td>
      <td>84.0</td>
      <td>13.0</td>
      <td>0.000036</td>
      <td>delta_expr</td>
      <td>-1.04997</td>
      <td>&gt;=50</td>
    </tr>
  </tbody>
</table>



```python
print('Here are the scores for mutations to disulfide bonds:')

p = (ggplot(dms_scores
            .replace({'delta_bind':'ACE2 binding', 'delta_expr':'RBD expression'})
            .assign(wildtype=lambda x: x['mutation'].str[0])
            .query('wildtype=="C" & mutant!="C"')
           ) +
     aes(x='score') + 
     geom_histogram(binwidth=0.25) +
     geom_vline(data=new_filter_df.replace({'delta_bind':'ACE2 binding', 'delta_expr':'RBD expression'}),
                     mapping=aes(xintercept='score'),
                    linetype='dashed',
                    color=CBPALETTE[1])+
     facet_wrap('~ metric') +
     theme_classic() +
     theme(figure_size=(2.5 * 2, 2.5 * 1),
           plot_title=element_text(size=14),
           axis_text_x=element_text(angle=90),
           strip_background=element_blank(),
           strip_text=element_text(size=12),
          ) +
     xlab('B.1.351 RBD DMS score')
     )

fig = p.draw()

plotfile = os.path.join(resultsdir, f"disulfide-histogram.pdf")
print(f"Saving plot to {plotfile}")
p.save(plotfile, verbose=False)
```

    Here are the scores for mutations to disulfide bonds:
    Saving plot to results/bind_expr_filters/disulfide-histogram.pdf



    
![png](bind_expr_filters_files/bind_expr_filters_18_1.png)
    


### Get the bind and expr scores that correspond to the 5th percentile of mutations observed at least 50x in GISAID


```python
def get_filter(scores_df, metric, count_threshold, percentile):
    
    scores=(scores_df
            .query('metric==@metric & count >=@count_threshold')
            .dropna()
            )['score'].tolist()
            
    c=np.percentile(scores, percentile)
    
    return c

count_thresholds = [50]
percentiles=[1,2.5,5,10,25]

v=[]

for i in count_thresholds:
    for p in percentiles:
        t=(i,p)
        
        scores=(dms_scores)
        bind_filter=get_filter(scores, 'delta_bind', i, p)
        expr_filter=get_filter(scores, 'delta_expr', i, p)
        
        t=(i, p, bind_filter, expr_filter)
        
        v.append(t)
        

df = pd.DataFrame(v, columns =['count_threshold', 'percentile', 'bind_count', 'expr_count'])
display(HTML(df.to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>count_threshold</th>
      <th>percentile</th>
      <th>bind_count</th>
      <th>expr_count</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>50</td>
      <td>1.0</td>
      <td>-1.851553</td>
      <td>-1.991598</td>
    </tr>
    <tr>
      <td>50</td>
      <td>2.5</td>
      <td>-1.141998</td>
      <td>-1.474578</td>
    </tr>
    <tr>
      <td>50</td>
      <td>5.0</td>
      <td>-0.831732</td>
      <td>-1.054092</td>
    </tr>
    <tr>
      <td>50</td>
      <td>10.0</td>
      <td>-0.600498</td>
      <td>-0.680866</td>
    </tr>
    <tr>
      <td>50</td>
      <td>25.0</td>
      <td>-0.210150</td>
      <td>-0.276720</td>
    </tr>
  </tbody>
</table>



```python
og_dms_scores=(pd.read_csv(og_dms_file)
               # remove extraneous columns
               .drop(columns=['site_RBD','wildtype', 'mutation', 'mutation_RBD', 'bind_lib1', 'bind_lib2', 'expr_lib1', 'expr_lib2'])
               # rename some columns
               .rename(columns={'site_SARS2':'site', 'bind_avg':'delta_bind', 'expr_avg':'delta_expr'})
              )

display(HTML(og_dms_scores.head(2).to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>site</th>
      <th>mutant</th>
      <th>delta_bind</th>
      <th>delta_expr</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>331</td>
      <td>A</td>
      <td>-0.03</td>
      <td>-0.11</td>
    </tr>
    <tr>
      <td>331</td>
      <td>C</td>
      <td>-0.09</td>
      <td>-1.26</td>
    </tr>
  </tbody>
</table>



```python
dms_scores=(dms_scores
            .merge((og_dms_scores
                    .melt(id_vars=['site', 'mutant',],
                          value_vars=['delta_bind', 'delta_expr'], 
                          var_name='metric', 
                          value_name='wuhan1dms_score',
                         )
                   ),
                   how='left',
                   on=['site', 'mutant', 'metric'],
                   validate='many_to_one'
                  )
           )
display(HTML(dms_scores.head().to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>wildtype</th>
      <th>mutation</th>
      <th>site</th>
      <th>mutant</th>
      <th>count</th>
      <th>n_countries</th>
      <th>frequency</th>
      <th>metric</th>
      <th>score</th>
      <th>count_categories</th>
      <th>wuhan1dms_score</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>N</td>
      <td>N331A</td>
      <td>331</td>
      <td>A</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.000000e+00</td>
      <td>delta_bind</td>
      <td>0.12426</td>
      <td>0</td>
      <td>-0.03</td>
    </tr>
    <tr>
      <td>N</td>
      <td>N331C</td>
      <td>331</td>
      <td>C</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.000000e+00</td>
      <td>delta_bind</td>
      <td>-0.20203</td>
      <td>0</td>
      <td>-0.09</td>
    </tr>
    <tr>
      <td>N</td>
      <td>N331D</td>
      <td>331</td>
      <td>D</td>
      <td>2.0</td>
      <td>2.0</td>
      <td>8.673362e-07</td>
      <td>delta_bind</td>
      <td>0.00131</td>
      <td>1 to 9</td>
      <td>0.03</td>
    </tr>
    <tr>
      <td>N</td>
      <td>N331E</td>
      <td>331</td>
      <td>E</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.000000e+00</td>
      <td>delta_bind</td>
      <td>0.10980</td>
      <td>0</td>
      <td>0.00</td>
    </tr>
    <tr>
      <td>N</td>
      <td>N331F</td>
      <td>331</td>
      <td>F</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.000000e+00</td>
      <td>delta_bind</td>
      <td>-0.15674</td>
      <td>0</td>
      <td>-0.10</td>
    </tr>
  </tbody>
</table>



```python
dms_score_corrs = (dms_scores
                   .replace({'delta_bind':'ACE2 binding', 'delta_expr':'RBD expression'})
                   .assign(is_disulfide=lambda x: x['mutation'].apply(lambda s: s[0]).isin(["C"]))
                  )

# plot correlations
xmin = dms_score_corrs['wuhan1dms_score'].min()
xspan = dms_score_corrs['wuhan1dms_score'].max() - xmin
ymin = dms_score_corrs['score'].min()
yspan = dms_score_corrs['score'].max() - ymin

p = (ggplot(dms_score_corrs) +
     aes('wuhan1dms_score', 'score') +
     geom_point(aes(fill='is_disulfide', alpha='is_disulfide'), color='black') +
     geom_vline(data=og_thresholds_df.replace({'delta_bind':'ACE2 binding', 'delta_expr':'RBD expression'}),
                     mapping=aes(xintercept='score'),
                    linetype='dashed',
                    color=CBPALETTE[1])+
     geom_hline(data=new_filter_df.replace({'delta_bind':'ACE2 binding', 'delta_expr':'RBD expression'}),
                     mapping=aes(yintercept='score'),
                    linetype='dashed',
                    color=CBPALETTE[1])+
     facet_wrap('~ metric') +
     theme_classic() +
     theme(figure_size=(2.5 * 2, 2.5 * 1),
           plot_title=element_text(size=14),
           strip_background=element_blank(),
           strip_text=element_text(size=12),
          ) +
     xlab('Wuhan-Hu-1 DMS score from Starr et al. (2020)') +
     ylab('B.1.351 RBD DMS score') +
     scale_fill_manual(values=['black', 'red']) +
     scale_alpha_manual(values=[0.2,1]) +
     labs(fill = 'mutation to disulfide bond?',
          alpha = 'mutation to disulfide bond?'
         )
     )

p.draw()
plotfile = os.path.join(resultsdir, f"wuhan-b1351-corr.pdf")
print(f"Saving plot to {plotfile}")
p.save(plotfile, verbose=False)
```

    /fh/fast/bloom_j/computational_notebooks/agreaney/2021/SARS-CoV-2-RBD_B.1.351/env/lib/python3.8/site-packages/plotnine/layer.py:401: PlotnineWarning: geom_point : Removed 72 rows containing missing values.


    Saving plot to results/bind_expr_filters/wuhan-b1351-corr.pdf


    /fh/fast/bloom_j/computational_notebooks/agreaney/2021/SARS-CoV-2-RBD_B.1.351/env/lib/python3.8/site-packages/plotnine/layer.py:401: PlotnineWarning: geom_point : Removed 72 rows containing missing values.



    
![png](bind_expr_filters_files/bind_expr_filters_23_3.png)
    



```python
print('Mutations from the Wuhan-Hu-1 library that:')
print('pass bind: '+ str(len(og_dms_scores.query('delta_bind >= -2.35'))))
print('pass expr: '+ str(len(og_dms_scores.query('delta_expr >= -1.0'))))
print('pass both: '+ str(len(og_dms_scores.query('delta_bind >= -2.35 & delta_expr >= -1.0'))))
```

    Mutations from the Wuhan-Hu-1 library that:
    pass bind: 3422
    pass expr: 2328
    pass both: 2269



```python
bind_threshold=new_thresholds['delta_bind']
expr_threshold=new_thresholds['delta_expr']
        
n_bind=len(dms_scores.query('metric=="delta_bind" & score >= @bind_threshold'))
n_expr=len(dms_scores.query('metric=="delta_expr" & score >= @expr_threshold'))

df=(dms_scores
     .pivot_table(index=['mutation', 'wildtype', 'mutant'],
                  values=['score'],
                  columns=['metric'],
                 )
     .reset_index()
       )

df.columns=['mutation', 'wildtype', 'mutant','delta_bind', 'delta_expr']

n_both=len(df
           .query('delta_bind >= @bind_threshold & delta_expr >= @expr_threshold')
          )
        
n_both_notC=len((df
                .assign(not_disulfide=lambda x: x['mutation'].str[0] != "C")
                .query('delta_bind >= @bind_threshold & delta_expr >= @expr_threshold & not_disulfide')
          ))

n_both_notC_notWT=len((df
                .assign(not_disulfide=lambda x: x['mutation'].str[0] != "C")
                .assign(not_WT=lambda x: x['wildtype']!=x['mutant'])
                .query('delta_bind >= @bind_threshold & delta_expr >= @expr_threshold & not_disulfide & not_WT')
          ))

total_muts_notC=len((df
                .assign(not_disulfide=lambda x: x['mutation'].str[0] != "C")
                .assign(not_WT=lambda x: x['wildtype']!=x['mutant'])
                .query('not_disulfide & not_WT')
          ))

print(f'B.1.351 mutations that \npass bind: {n_bind} \npass expr: {n_expr} \npass both: {n_both} \npass both and not disulfide: {n_both_notC}')
print(f'Pass bind, expr, not disulfide, and not WT: {n_both_notC_notWT}')

print(f'Total number of possible mutations to non-disulfide sites: {total_muts_notC}')
```

    B.1.351 mutations that 
    pass bind: 3291 
    pass expr: 2369 
    pass both: 2266 
    pass both and not disulfide: 2207
    Pass bind, expr, not disulfide, and not WT: 2014
    Total number of possible mutations to non-disulfide sites: 3653



```python
print(f'This percentage of all variants seen >=50x in GISAID are retained by the binding filter of {bind_threshold}')
print(round(100-stats.percentileofscore((dms_scores
                               .query('metric=="delta_bind" & count>=50')['score']), 
                              bind_threshold, 
                              kind='rank'
                             ),
            1
           )
     )

print(f'This percentage of all variants seen >=50x in GISAID are retained by the expression filter of {expr_threshold}')
print(round(100-stats.percentileofscore((dms_scores
                               .query('metric=="delta_expr" & count>=50')['score']), 
                              expr_threshold, 
                              kind='rank'
                             ),
            1
           )
     )


# dms_scores.query('metric=="delta_bind" & score >= @bind_threshold & count>=50')['score'].min()
```

    This percentage of all variants seen >=50x in GISAID are retained by the binding filter of -3.0
    99.4
    This percentage of all variants seen >=50x in GISAID are retained by the expression filter of -1.0
    93.8



```python

```
