<h1>Table of Contents<span class="tocSkip"></span></h1>
<div class="toc"><ul class="toc-item"><li><span><a href="#Process-CCSs" data-toc-modified-id="Process-CCSs-1">Process CCSs</a></span><ul class="toc-item"><li><span><a href="#Setup" data-toc-modified-id="Setup-1.1">Setup</a></span></li><li><span><a href="#PacBio-amplicons" data-toc-modified-id="PacBio-amplicons-1.2">PacBio amplicons</a></span></li><li><span><a href="#CCS-stats-for-PacBio-runs" data-toc-modified-id="CCS-stats-for-PacBio-runs-1.3">CCS stats for PacBio runs</a></span></li><li><span><a href="#Align-CCSs-to-amplicons" data-toc-modified-id="Align-CCSs-to-amplicons-1.4">Align CCSs to amplicons</a></span></li><li><span><a href="#Write-valid-CCSs" data-toc-modified-id="Write-valid-CCSs-1.5">Write valid CCSs</a></span></li></ul></li></ul></div>

This Python Jupyter notebook processes the PacBio circular consensus sequences (CCSs) to extract barcodes and call mutations in the gene. It then builds consensus sequences for barcoded variants from the mutations, to build a barcode-variant lookup table.

# Process CCSs
First, process the PacBio CCSs to extract barcodes and call mutations in the gene.

Define the background for CCS mapping


```python
background = "panDCoV"
```

## Setup

Import Python modules

Plotting is done with [plotnine](https://plotnine.readthedocs.io/en/stable/), which uses ggplot2-like syntax.

The analysis uses the Bloom lab's [alignparse](https://jbloomlab.github.io/alignparse) and [dms_variants](https://jbloomlab.github.io/dms_variants) packages.


```python
import collections
import math
import os
import re
import time
import warnings

import alignparse
import alignparse.ccs
from alignparse.constants import CBPALETTE
import alignparse.minimap2
import alignparse.targets
import alignparse.consensus
import alignparse.utils

import dms_variants
import dms_variants.codonvarianttable
import dms_variants.plotnine_themes
import dms_variants.utils

from IPython.display import display, HTML

import numpy

import pandas as pd

from plotnine import *

import yaml

%matplotlib inline
```

Set [plotnine](https://plotnine.readthedocs.io/en/stable/) theme to the one defined in [dms_variants](https://jbloomlab.github.io/dms_variants):


```python
theme_set(dms_variants.plotnine_themes.theme_graygrid())
```

Versions of key software:


```python
print(f"Using alignparse version {alignparse.__version__}")
print(f"Using dms_variants version {dms_variants.__version__}")
print(f"Using pandas version {pd.__version__}")
```

    Using alignparse version 0.6.0
    Using dms_variants version 1.4.3
    Using pandas version 1.3.5


Ignore warnings that clutter output:


```python
warnings.simplefilter('ignore')
```

Read the configuration file:


```python
with open('config.yaml') as f:
    config = yaml.safe_load(f)
```

Make output directory for figures:


```python
os.makedirs(config['process_ccs_dir'], exist_ok=True)
```

## PacBio reads
Get the amplicons sequenced by PacBio as the alignment target along with the specs on how to parse the features:


```python
pacbio_runs = (
    pd.read_csv(config['pacbio_runs'], dtype=str)
    .drop(columns=['ccs'])
    .assign(name=lambda x: x['library'] + '_' + x['bg'] + '_' + x['run'],
            fastq=lambda x: config['ccs_dir'] + '/' + x['name'] + '_ccs.fastq.gz'
            )
    )
 
pacbio_runs.drop(pacbio_runs[pacbio_runs['bg'] != background].index,inplace=True)

display(HTML(pacbio_runs.to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>bg</th>
      <th>run</th>
      <th>name</th>
      <th>fastq</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>lib80</td>
      <td>panDCoV</td>
      <td>240206</td>
      <td>lib80_panDCoV_240206</td>
      <td>results/ccs/lib80_panDCoV_240206_ccs.fastq.gz</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>panDCoV</td>
      <td>240206</td>
      <td>lib81_panDCoV_240206</td>
      <td>results/ccs/lib81_panDCoV_240206_ccs.fastq.gz</td>
    </tr>
  </tbody>
</table>


## PacBio amplicons
Get the amplicons sequenced by PacBio as the alignment target along with the specs on how to parse the features:


```python
print(f"Reading amplicons from {config['amplicons_' + background]}")
print(f"Reading feature parse specs from {config['feature_parse_specs_' + background]}")

targets = alignparse.targets.Targets(
                seqsfile=config['amplicons_' + background],
                feature_parse_specs=config['feature_parse_specs_' + background])
```

    Reading amplicons from data/PacBio_amplicon_panDCoV.gb
    Reading feature parse specs from data/feature_parse_specs_panDCoV.yaml


Draw the target amplicons:


```python
fig = targets.plot(ax_width=7,
                   plots_indexing='biopython',  # numbering starts at 0
                   ax_height=2,  # height of each plot
                   hspace=1.2,  # vertical space between plots
                   )

plotfile = os.path.join(config['process_ccs_dir'], 'amplicons_'+background+'.pdf')
print(f"Saving plot to {plotfile}")
fig.savefig(plotfile, bbox_inches='tight')
```

    Saving plot to results/process_ccs/amplicons_panDCoV.pdf



    
![png](process_ccs_panDCoV_files/process_ccs_panDCoV_22_1.png)
    


Write out the specs used to parse the features (these are the same specs provided as `feature_parse_specs` when initializing `targets`, but with defaults filled in):


```python
print(targets.feature_parse_specs('yaml'))
```

    MagpieRobinCoV_HKU18_chu3: &id001
      query_clip5: 4
      query_clip3: 4
      termini5:
        filter:
          clip5: 4
          mutation_nt_count: 1
          mutation_op_count: null
          clip3: 0
        return: []
      gene:
        filter:
          mutation_nt_count: 90
          mutation_op_count: null
          clip5: 0
          clip3: 0
        return:
        - mutations
        - accuracy
      spacer:
        filter:
          mutation_nt_count: 1
          mutation_op_count: null
          clip5: 0
          clip3: 0
        return: []
      barcode:
        filter:
          mutation_nt_count: 0
          mutation_op_count: null
          clip5: 0
          clip3: 0
        return:
        - sequence
        - accuracy
      termini3:
        filter:
          clip3: 4
          mutation_nt_count: 1
          mutation_op_count: null
          clip5: 0
        return: []
    PDCoV_JS2018_F101: *id001
    PDCoV_110_1197_TW_2021: *id001
    PDCoV_VUT_1_2016_Vietnam: *id001
    ThrushCoV_HKU12_600: *id001
    AncHerdeco3: *id001
    AncBulDCoV5_alt: *id001
    AncHerdeco3_alt: *id001
    GullCoV_HNU4_1: *id001
    PDCoV_AH2018_432: *id001
    PDCoV_HN2019_C115: *id001
    PDCoV_HeN_2015: *id001
    AncAndecoHerdeco: *id001
    AncHerdeco4: *id001
    AncAndecoHerdeco_alt: *id001
    AncHerdeco4_alt: *id001
    GullCoV_HNU4_3: *id001
    PDCoV_P29_15_VN_1215: *id001
    PDCoV_HKU15_SD2018_300: *id001
    PDCoV_SD2018_311: *id001
    AncBulDCoV6: *id001
    AncHerdeco5: *id001
    AncBulDCoV6_alt: *id001
    AncHerdeco5_alt: *id001
    HoubaraCoV_UAE_HKU28_285F: *id001
    PDCoV_CHJXNI2_2015: *id001
    PDCoV_SD2018_10: *id001
    SparrowCoV_ISU42824: *id001
    AncBulDCoV7: *id001
    AncDCoV_alt: *id001
    AncBulDCoV7_alt: *id001
    AvianCoV_MW01_1o: *id001
    PDCoV_KNU16_11: *id001
    PDCoV_AH2018_515: *id001
    MuniaCoV_HKU13_3514: *id001
    AncBulDCoV8: *id001
    AncBulDCoV1_alt: *id001
    AncBulDCoV8_alt: *id001
    AvianCoV_rub035cor1: *id001
    PDCoV_Minnesota455_2014: *id001
    PDCoV_CHN_GD16_05: *id001
    SparrowCoV_ISU690_4: *id001
    AncBulDCoV9: *id001
    AncBulDCoV2_alt: *id001
    AncBulDCoV9_alt: *id001
    AvianCoV_lrf178cor1: *id001
    PDCoV_Tianjin_2016: *id001
    PDCoV_JS2019_A1414: *id001
    QuailCoV_UAE_HKU30_411F: *id001
    AncHerdeco1: *id001
    AncBulDCoV3_alt: *id001
    AncHerdeco1_alt: *id001
    PDCoV_Illinois121_2014: *id001
    PDCoV_CH_JXJGS01_2016: *id001
    PDCoV_Thailand_S5011_2015: *id001
    QuailCoV_UAE_HKU30_1101F: *id001
    AncHerdeco2: *id001
    AncBulDCoV4_alt: *id001
    AncHerdeco2_alt: *id001
    HKU20-9243_WigeonCoV: *id001
    brb027cor1_AvianCoV: *id001
    HKU27-988F_FalconCoV: *id001
    HKU17-6124_SparrowCoV: *id001
    HKU19-6918_HeronCoV: *id001
    AncDCoV: *id001
    AncBulDCoV1: *id001
    AncBulDCoV2: *id001
    AncBulDCoV3: *id001
    AncBulDCoV4: *id001
    AncBulDCoV5: *id001
    AncPDCoV: *id001
    HKU21-9295_MoorhenCoV: *id001
    HKU16-6847_WhiteEyeCoV: *id001
    HKU11-934_BulbulCoV: *id001
    ISU73347_SparrowCoV: *id001
    


## Align CCSs to amplicons
We now align the CCSs to the amplicon and parse features from the resulting alignments using the specs above.

First, we initialize an `alignparse.minimap2.Mapper` to align the reads to SAM files:


```python
mapper = alignparse.minimap2.Mapper(alignparse.minimap2.OPTIONS_CODON_DMS)

print(f"Using `minimap2` {mapper.version} with these options:\n" +
      ' '.join(mapper.options))
```

    Using `minimap2` 2.26-r1175 with these options:
    -A2 -B4 -O12 -E2 --end-bonus=13 --secondary=no --cs


Next, we use `Targets.align_and_parse` to create the alignments and parse them:


```python
readstats, aligned, filtered = targets.align_and_parse(
        df=pacbio_runs,
        mapper=mapper,
        outdir=config['process_ccs_dir'],
        name_col='run',
        group_cols=['name', 'library'],
        queryfile_col='fastq',
        overwrite=True,
        ncpus=config['max_cpus'],
        )
```

First, examine the read stats from the alignment / parsing, both extracting alignment target name and getting stats aggregated by target:


```python
readstats = (
    readstats
    .assign(category_all_targets=lambda x: x['category'].str.split().str[0],
            target=lambda x: x['category'].str.split(None, 1).str[1],
            valid=lambda x: x['category_all_targets'] == 'aligned')
    )
```

Now plot the read stats by run (combining all targets and libraries within a run):


```python
ncol = 2
p = (
    ggplot(readstats
           .groupby(['name', 'category_all_targets', 'valid'])
           .aggregate({'count': 'sum'})
           .reset_index(),
           aes('category_all_targets', 'count', fill='valid')) +
    geom_bar(stat='identity') +
    facet_wrap('~ name', ncol=ncol) +
    theme(axis_text_x=element_text(angle=90),
          figure_size=(1.85 * min(ncol, len(pacbio_runs)),
                       2 * math.ceil(len(pacbio_runs) / ncol)),
          panel_grid_major_x=element_blank(),
          legend_position='none',
          ) +
    scale_fill_manual(values=CBPALETTE)
    )
_ = p.draw()
```


    
![png](process_ccs_panDCoV_files/process_ccs_panDCoV_32_0.png)
    


And the read stats by library (combining all targets and runs within a library):


```python
p = (
    ggplot(readstats
           .groupby(['library', 'category_all_targets', 'valid'])
           .aggregate({'count': 'sum'})
           .reset_index(), 
           aes('category_all_targets', 'count', fill='valid')) +
    geom_bar(stat='identity') +
    facet_wrap('~ library', nrow=1) +
    theme(axis_text_x=element_text(angle=90),
          figure_size=(1.5 * pacbio_runs['library'].nunique(), 2),
          panel_grid_major_x=element_blank(),
          legend_position='none',
          ) +
    scale_fill_manual(values=CBPALETTE)
    )
_ = p.draw()
```


    
![png](process_ccs_panDCoV_files/process_ccs_panDCoV_34_0.png)
    


And the number of reads by target (combining all libraries and runs for a target):


```python
p = (
    ggplot(readstats
           .groupby(['target'])
           .aggregate({'count': 'sum'})
           .reset_index(), 
           aes('count', 'target')) +
    geom_point(stat='identity', size=3) +
    theme(axis_text_x=element_text(angle=90),
          figure_size=(2, 0.3 * readstats['target'].nunique()),
          panel_grid_major_y=element_blank(),
          ) +
    scale_x_log10(name='number of reads')
    )
_ = p.draw()
```


    
![png](process_ccs_panDCoV_files/process_ccs_panDCoV_36_0.png)
    


And read stats by target (combining all libraries and runs for a target):


```python
p = (
    ggplot(readstats
           .groupby(['target', 'valid'])
           .aggregate({'count': 'sum'})
           .reset_index()
           .assign(total=lambda x: x.groupby('target')['count'].transform('sum'),
                   frac=lambda x: x['count'] / x['total'],
                   ), 
           aes('target', 'frac', fill='valid')) +
    geom_bar(stat='identity') +
    theme(axis_text_x=element_text(angle=90),
          figure_size=(0.5 * readstats['target'].nunique(), 2),
          panel_grid_major_x=element_blank(),
          ) +
    scale_fill_manual(values=CBPALETTE)
    )
_ = p.draw()
```


    
![png](process_ccs_panDCoV_files/process_ccs_panDCoV_38_0.png)
    


Now let's see **why** we filtered the reads.
First, we do some transformations on the `filtered` dict returned by `Targets.align_and_parse`.
Then we count up the number of CCSs filtered for each reason, and group together "unusual" reasons that represent less than some fraction of all filtering.
For now, we group together all targets to the stats represent all targets combined:


```python
other_cutoff = 0.02  # group as "other" reasons with <= this frac

filtered_df = (
    pd.concat(df.assign(target=target) for target, df in filtered.items())
    .groupby(['library', 'name', 'run', 'filter_reason'])
    .size()
    .rename('count')
    .reset_index()
    .assign(tot_reason_frac=lambda x: (x.groupby('filter_reason')['count']
                                       .transform('sum')) / x['count'].sum(),
            filter_reason=lambda x: numpy.where(x['tot_reason_frac'] > other_cutoff,
                                                x['filter_reason'],
                                                'other')
            )
    )
```

Now plot the filtering reason for all runs:


```python
ncol = 7
nreasons = filtered_df['filter_reason'].nunique()

p = (
    ggplot(filtered_df, aes('filter_reason', 'count')) +
    geom_bar(stat='identity') +
    facet_wrap('~ name', ncol=ncol) +
    theme(axis_text_x=element_text(angle=90),
          figure_size=(0.25 * nreasons * min(ncol, len(pacbio_runs)),
                       2 * math.ceil(len(pacbio_runs) / ncol)),
          panel_grid_major_x=element_blank(),
          )
    )
_ = p.draw()
```


    
![png](process_ccs_panDCoV_files/process_ccs_panDCoV_42_0.png)
    


Finally, we take the successfully parsed alignments and read them into a data frame, keeping track of the target that each CCS aligns to.
We also drop the pieces of information we won't use going forward, and rename a few columns:


```python
aligned_df = (
    pd.concat(df.assign(target=target) for target, df in aligned.items())
    .drop(columns=['query_clip5', 'query_clip3', 'run','name'])
    .rename(columns={'barcode_sequence': 'barcode'})
    )

print(f"First few lines of information on the parsed alignments:")
display(HTML(aligned_df.head().to_html(index=False)))
```

    First few lines of information on the parsed alignments:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>query_name</th>
      <th>gene_mutations</th>
      <th>gene_accuracy</th>
      <th>barcode</th>
      <th>barcode_accuracy</th>
      <th>target</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>lib80</td>
      <td>m64296e_240129_011756/66238/ccs</td>
      <td></td>
      <td>1.000000</td>
      <td>AAATCTAGGTAGATGA</td>
      <td>1.000000</td>
      <td>PDCoV_JS2018_F101</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>m64296e_240129_011756/67609/ccs</td>
      <td></td>
      <td>0.999997</td>
      <td>CCCGGAAAACCAAACG</td>
      <td>0.980229</td>
      <td>PDCoV_JS2018_F101</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>m64296e_240129_011756/328789/ccs</td>
      <td></td>
      <td>1.000000</td>
      <td>CCCGGAAAACCAAACG</td>
      <td>1.000000</td>
      <td>PDCoV_JS2018_F101</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>m64296e_240129_011756/328999/ccs</td>
      <td></td>
      <td>1.000000</td>
      <td>CACACGCCATTGAATC</td>
      <td>1.000000</td>
      <td>PDCoV_JS2018_F101</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>m64296e_240129_011756/395089/ccs</td>
      <td></td>
      <td>1.000000</td>
      <td>GGAGCCCAATATGATC</td>
      <td>1.000000</td>
      <td>PDCoV_JS2018_F101</td>
    </tr>
  </tbody>
</table>


## Write valid CCSs

Write the processed CCSs to a file:


```python
aligned_df.to_csv(config['processed_ccs_file' + '_' + background], index=False)

print("Barcodes and mutations for valid processed CCSs "
      f"have been written to {config['processed_ccs_file' + '_' + background]}.")
```

    Barcodes and mutations for valid processed CCSs have been written to results/process_ccs/processed_ccs_panDCoV.csv.


Next, we analyze these processed CCSs to build the variants.

# Build barcode variant table
Builds consensus sequences for barcoded variants from the mutations called in the processed PacBio CCSs.
Uses these consensus sequences to build a codon variant table.

Make output directories if needed:


```python
os.makedirs(config['variants_dir'], exist_ok=True)
os.makedirs(config['figs_dir'], exist_ok=True)
```

Read the CSV file with the processed CCSs into a data frame, display first few lines:


```python
processed_ccs = pd.read_csv(config['processed_ccs_file' + '_' + background], na_filter=None)

nlibs = processed_ccs['library'].nunique()  # number of unique libraries

ntargets = processed_ccs['target'].nunique()  # number of unique targets

print(f"Read {len(processed_ccs)} CCSs from {nlibs} libraries and {ntargets} targets.")

display(HTML(processed_ccs.head().to_html(index=False)))
```

    Read 173611 CCSs from 2 libraries and 75 targets.



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>query_name</th>
      <th>gene_mutations</th>
      <th>gene_accuracy</th>
      <th>barcode</th>
      <th>barcode_accuracy</th>
      <th>target</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>lib80</td>
      <td>m64296e_240129_011756/66238/ccs</td>
      <td></td>
      <td>1.000000</td>
      <td>AAATCTAGGTAGATGA</td>
      <td>1.000000</td>
      <td>PDCoV_JS2018_F101</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>m64296e_240129_011756/67609/ccs</td>
      <td></td>
      <td>0.999997</td>
      <td>CCCGGAAAACCAAACG</td>
      <td>0.980229</td>
      <td>PDCoV_JS2018_F101</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>m64296e_240129_011756/328789/ccs</td>
      <td></td>
      <td>1.000000</td>
      <td>CCCGGAAAACCAAACG</td>
      <td>1.000000</td>
      <td>PDCoV_JS2018_F101</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>m64296e_240129_011756/328999/ccs</td>
      <td></td>
      <td>1.000000</td>
      <td>CACACGCCATTGAATC</td>
      <td>1.000000</td>
      <td>PDCoV_JS2018_F101</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>m64296e_240129_011756/395089/ccs</td>
      <td></td>
      <td>1.000000</td>
      <td>GGAGCCCAATATGATC</td>
      <td>1.000000</td>
      <td>PDCoV_JS2018_F101</td>
    </tr>
  </tbody>
</table>


Optional: discard reads that have a `fwd` or `rev` indicator suggestive of heteroduplexes (or we suspect in these runs, some form of pseudo-duplex).


```python
#start = len(processed_ccs)
#processed_ccs = processed_ccs.loc[~processed_ccs['query_name'].str.contains('fwd')]
#processed_ccs = processed_ccs.loc[~processed_ccs['query_name'].str.contains('rev')]
#end = len(processed_ccs)

#print(f"Went from {start} CCSs to {end} after discarded split/duplexed fwd and rev reads.")

#display(HTML(processed_ccs.head().to_html(index=False)))
```

Overall statistics on number of total CCSs and number of unique barcodes:


```python
display(HTML(
    processed_ccs
    .groupby(['target', 'library'])
    .aggregate(total_CCSs=('barcode', 'size'),
               unique_barcodes=('barcode', 'nunique'))
    .assign(avg_CCSs_per_barcode=lambda x: x['total_CCSs'] / x['unique_barcodes'])
    .round(2)
    .to_html()
    ))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th></th>
      <th>total_CCSs</th>
      <th>unique_barcodes</th>
      <th>avg_CCSs_per_barcode</th>
    </tr>
    <tr>
      <th>target</th>
      <th>library</th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th rowspan="2" valign="top">AncAndecoHerdeco</th>
      <th>lib80</th>
      <td>471</td>
      <td>65</td>
      <td>7.25</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>1547</td>
      <td>113</td>
      <td>13.69</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncAndecoHerdeco_alt</th>
      <th>lib80</th>
      <td>779</td>
      <td>80</td>
      <td>9.74</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>1214</td>
      <td>99</td>
      <td>12.26</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncBulDCoV1</th>
      <th>lib80</th>
      <td>772</td>
      <td>137</td>
      <td>5.64</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>1494</td>
      <td>173</td>
      <td>8.64</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncBulDCoV1_alt</th>
      <th>lib80</th>
      <td>684</td>
      <td>71</td>
      <td>9.63</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>1113</td>
      <td>89</td>
      <td>12.51</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncBulDCoV2</th>
      <th>lib80</th>
      <td>723</td>
      <td>99</td>
      <td>7.30</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>1003</td>
      <td>123</td>
      <td>8.15</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncBulDCoV2_alt</th>
      <th>lib80</th>
      <td>894</td>
      <td>86</td>
      <td>10.40</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>1486</td>
      <td>119</td>
      <td>12.49</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncBulDCoV3</th>
      <th>lib80</th>
      <td>769</td>
      <td>103</td>
      <td>7.47</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>995</td>
      <td>136</td>
      <td>7.32</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncBulDCoV3_alt</th>
      <th>lib80</th>
      <td>987</td>
      <td>89</td>
      <td>11.09</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>1041</td>
      <td>93</td>
      <td>11.19</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncBulDCoV4</th>
      <th>lib80</th>
      <td>1003</td>
      <td>142</td>
      <td>7.06</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>1136</td>
      <td>168</td>
      <td>6.76</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncBulDCoV4_alt</th>
      <th>lib80</th>
      <td>1639</td>
      <td>158</td>
      <td>10.37</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>1631</td>
      <td>164</td>
      <td>9.95</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncBulDCoV5</th>
      <th>lib80</th>
      <td>1070</td>
      <td>119</td>
      <td>8.99</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>900</td>
      <td>114</td>
      <td>7.89</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncBulDCoV5_alt</th>
      <th>lib80</th>
      <td>1005</td>
      <td>104</td>
      <td>9.66</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>907</td>
      <td>95</td>
      <td>9.55</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncBulDCoV6</th>
      <th>lib80</th>
      <td>1084</td>
      <td>101</td>
      <td>10.73</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>1090</td>
      <td>122</td>
      <td>8.93</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncBulDCoV6_alt</th>
      <th>lib80</th>
      <td>1189</td>
      <td>108</td>
      <td>11.01</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>861</td>
      <td>101</td>
      <td>8.52</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncBulDCoV7</th>
      <th>lib80</th>
      <td>964</td>
      <td>98</td>
      <td>9.84</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>801</td>
      <td>86</td>
      <td>9.31</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncBulDCoV7_alt</th>
      <th>lib80</th>
      <td>972</td>
      <td>94</td>
      <td>10.34</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>969</td>
      <td>100</td>
      <td>9.69</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncBulDCoV8</th>
      <th>lib80</th>
      <td>1425</td>
      <td>138</td>
      <td>10.33</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>1381</td>
      <td>158</td>
      <td>8.74</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncBulDCoV8_alt</th>
      <th>lib80</th>
      <td>1485</td>
      <td>140</td>
      <td>10.61</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>1376</td>
      <td>146</td>
      <td>9.42</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncBulDCoV9</th>
      <th>lib80</th>
      <td>1220</td>
      <td>110</td>
      <td>11.09</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>1075</td>
      <td>131</td>
      <td>8.21</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncBulDCoV9_alt</th>
      <th>lib80</th>
      <td>880</td>
      <td>84</td>
      <td>10.48</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>930</td>
      <td>93</td>
      <td>10.00</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncDCoV</th>
      <th>lib80</th>
      <td>930</td>
      <td>121</td>
      <td>7.69</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>1441</td>
      <td>149</td>
      <td>9.67</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncDCoV_alt</th>
      <th>lib80</th>
      <td>1169</td>
      <td>111</td>
      <td>10.53</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>1846</td>
      <td>139</td>
      <td>13.28</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHerdeco1</th>
      <th>lib80</th>
      <td>1353</td>
      <td>138</td>
      <td>9.80</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>1907</td>
      <td>158</td>
      <td>12.07</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHerdeco1_alt</th>
      <th>lib80</th>
      <td>854</td>
      <td>70</td>
      <td>12.20</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>1052</td>
      <td>99</td>
      <td>10.63</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHerdeco2</th>
      <th>lib80</th>
      <td>911</td>
      <td>83</td>
      <td>10.98</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>1116</td>
      <td>94</td>
      <td>11.87</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHerdeco2_alt</th>
      <th>lib80</th>
      <td>768</td>
      <td>76</td>
      <td>10.11</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>1146</td>
      <td>102</td>
      <td>11.24</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHerdeco3</th>
      <th>lib80</th>
      <td>1266</td>
      <td>130</td>
      <td>9.74</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>1767</td>
      <td>143</td>
      <td>12.36</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHerdeco3_alt</th>
      <th>lib80</th>
      <td>1165</td>
      <td>115</td>
      <td>10.13</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>1935</td>
      <td>186</td>
      <td>10.40</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHerdeco4</th>
      <th>lib80</th>
      <td>1282</td>
      <td>126</td>
      <td>10.17</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>2254</td>
      <td>166</td>
      <td>13.58</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHerdeco4_alt</th>
      <th>lib80</th>
      <td>1171</td>
      <td>105</td>
      <td>11.15</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>1676</td>
      <td>127</td>
      <td>13.20</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHerdeco5</th>
      <th>lib80</th>
      <td>1196</td>
      <td>96</td>
      <td>12.46</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>1374</td>
      <td>103</td>
      <td>13.34</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHerdeco5_alt</th>
      <th>lib80</th>
      <td>1614</td>
      <td>138</td>
      <td>11.70</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>1867</td>
      <td>150</td>
      <td>12.45</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncPDCoV</th>
      <th>lib80</th>
      <td>991</td>
      <td>133</td>
      <td>7.45</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>842</td>
      <td>126</td>
      <td>6.68</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AvianCoV_MW01_1o</th>
      <th>lib80</th>
      <td>1217</td>
      <td>108</td>
      <td>11.27</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>1245</td>
      <td>105</td>
      <td>11.86</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AvianCoV_lrf178cor1</th>
      <th>lib80</th>
      <td>866</td>
      <td>81</td>
      <td>10.69</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>1268</td>
      <td>102</td>
      <td>12.43</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AvianCoV_rub035cor1</th>
      <th>lib80</th>
      <td>828</td>
      <td>81</td>
      <td>10.22</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>1120</td>
      <td>89</td>
      <td>12.58</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">GullCoV_HNU4_1</th>
      <th>lib80</th>
      <td>1300</td>
      <td>121</td>
      <td>10.74</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>1762</td>
      <td>125</td>
      <td>14.10</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">GullCoV_HNU4_3</th>
      <th>lib80</th>
      <td>1429</td>
      <td>135</td>
      <td>10.59</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>1957</td>
      <td>132</td>
      <td>14.83</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU11-934_BulbulCoV</th>
      <th>lib80</th>
      <td>1009</td>
      <td>127</td>
      <td>7.94</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>928</td>
      <td>135</td>
      <td>6.87</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU16-6847_WhiteEyeCoV</th>
      <th>lib80</th>
      <td>953</td>
      <td>122</td>
      <td>7.81</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>1419</td>
      <td>148</td>
      <td>9.59</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU17-6124_SparrowCoV</th>
      <th>lib80</th>
      <td>789</td>
      <td>123</td>
      <td>6.41</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>941</td>
      <td>109</td>
      <td>8.63</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU19-6918_HeronCoV</th>
      <th>lib80</th>
      <td>1045</td>
      <td>127</td>
      <td>8.23</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>832</td>
      <td>104</td>
      <td>8.00</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU20-9243_WigeonCoV</th>
      <th>lib80</th>
      <td>866</td>
      <td>132</td>
      <td>6.56</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>1613</td>
      <td>155</td>
      <td>10.41</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU21-9295_MoorhenCoV</th>
      <th>lib80</th>
      <td>674</td>
      <td>102</td>
      <td>6.61</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>1415</td>
      <td>136</td>
      <td>10.40</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU27-988F_FalconCoV</th>
      <th>lib80</th>
      <td>643</td>
      <td>90</td>
      <td>7.14</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>1258</td>
      <td>149</td>
      <td>8.44</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HoubaraCoV_UAE_HKU28_285F</th>
      <th>lib80</th>
      <td>1120</td>
      <td>107</td>
      <td>10.47</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>1610</td>
      <td>121</td>
      <td>13.31</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">ISU73347_SparrowCoV</th>
      <th>lib80</th>
      <td>1099</td>
      <td>150</td>
      <td>7.33</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>690</td>
      <td>110</td>
      <td>6.27</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">MagpieRobinCoV_HKU18_chu3</th>
      <th>lib80</th>
      <td>1256</td>
      <td>119</td>
      <td>10.55</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>1890</td>
      <td>154</td>
      <td>12.27</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">MuniaCoV_HKU13_3514</th>
      <th>lib80</th>
      <td>994</td>
      <td>79</td>
      <td>12.58</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>926</td>
      <td>104</td>
      <td>8.90</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">PDCoV_110_1197_TW_2021</th>
      <th>lib80</th>
      <td>917</td>
      <td>78</td>
      <td>11.76</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>1061</td>
      <td>101</td>
      <td>10.50</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">PDCoV_AH2018_432</th>
      <th>lib80</th>
      <td>1375</td>
      <td>145</td>
      <td>9.48</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>1354</td>
      <td>142</td>
      <td>9.54</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">PDCoV_AH2018_515</th>
      <th>lib80</th>
      <td>1023</td>
      <td>99</td>
      <td>10.33</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>918</td>
      <td>104</td>
      <td>8.83</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">PDCoV_CHJXNI2_2015</th>
      <th>lib80</th>
      <td>798</td>
      <td>71</td>
      <td>11.24</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>713</td>
      <td>77</td>
      <td>9.26</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">PDCoV_CHN_GD16_05</th>
      <th>lib80</th>
      <td>898</td>
      <td>71</td>
      <td>12.65</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>1175</td>
      <td>107</td>
      <td>10.98</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">PDCoV_CH_JXJGS01_2016</th>
      <th>lib80</th>
      <td>939</td>
      <td>86</td>
      <td>10.92</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>1172</td>
      <td>120</td>
      <td>9.77</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">PDCoV_HKU15_SD2018_300</th>
      <th>lib80</th>
      <td>1653</td>
      <td>139</td>
      <td>11.89</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>1231</td>
      <td>135</td>
      <td>9.12</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">PDCoV_HN2019_C115</th>
      <th>lib80</th>
      <td>1563</td>
      <td>131</td>
      <td>11.93</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>1213</td>
      <td>127</td>
      <td>9.55</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">PDCoV_HeN_2015</th>
      <th>lib80</th>
      <td>999</td>
      <td>102</td>
      <td>9.79</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>987</td>
      <td>98</td>
      <td>10.07</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">PDCoV_Illinois121_2014</th>
      <th>lib80</th>
      <td>799</td>
      <td>68</td>
      <td>11.75</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>929</td>
      <td>97</td>
      <td>9.58</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">PDCoV_JS2018_F101</th>
      <th>lib80</th>
      <td>1147</td>
      <td>105</td>
      <td>10.92</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>866</td>
      <td>97</td>
      <td>8.93</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">PDCoV_JS2019_A1414</th>
      <th>lib80</th>
      <td>1284</td>
      <td>118</td>
      <td>10.88</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>869</td>
      <td>93</td>
      <td>9.34</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">PDCoV_KNU16_11</th>
      <th>lib80</th>
      <td>1028</td>
      <td>101</td>
      <td>10.18</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>1441</td>
      <td>139</td>
      <td>10.37</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">PDCoV_Minnesota455_2014</th>
      <th>lib80</th>
      <td>730</td>
      <td>84</td>
      <td>8.69</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>1116</td>
      <td>110</td>
      <td>10.15</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">PDCoV_P29_15_VN_1215</th>
      <th>lib80</th>
      <td>875</td>
      <td>99</td>
      <td>8.84</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>815</td>
      <td>100</td>
      <td>8.15</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">PDCoV_SD2018_10</th>
      <th>lib80</th>
      <td>1269</td>
      <td>119</td>
      <td>10.66</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>1140</td>
      <td>124</td>
      <td>9.19</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">PDCoV_SD2018_311</th>
      <th>lib80</th>
      <td>1012</td>
      <td>112</td>
      <td>9.04</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>952</td>
      <td>121</td>
      <td>7.87</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">PDCoV_Thailand_S5011_2015</th>
      <th>lib80</th>
      <td>1407</td>
      <td>110</td>
      <td>12.79</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>992</td>
      <td>104</td>
      <td>9.54</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">PDCoV_Tianjin_2016</th>
      <th>lib80</th>
      <td>1249</td>
      <td>124</td>
      <td>10.07</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>1574</td>
      <td>156</td>
      <td>10.09</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">PDCoV_VUT_1_2016_Vietnam</th>
      <th>lib80</th>
      <td>2160</td>
      <td>185</td>
      <td>11.68</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>1962</td>
      <td>197</td>
      <td>9.96</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">QuailCoV_UAE_HKU30_1101F</th>
      <th>lib80</th>
      <td>1376</td>
      <td>103</td>
      <td>13.36</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>1529</td>
      <td>129</td>
      <td>11.85</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">QuailCoV_UAE_HKU30_411F</th>
      <th>lib80</th>
      <td>830</td>
      <td>80</td>
      <td>10.38</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>1338</td>
      <td>99</td>
      <td>13.52</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">SparrowCoV_ISU42824</th>
      <th>lib80</th>
      <td>964</td>
      <td>85</td>
      <td>11.34</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>914</td>
      <td>90</td>
      <td>10.16</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">SparrowCoV_ISU690_4</th>
      <th>lib80</th>
      <td>896</td>
      <td>83</td>
      <td>10.80</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>689</td>
      <td>79</td>
      <td>8.72</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">ThrushCoV_HKU12_600</th>
      <th>lib80</th>
      <td>578</td>
      <td>79</td>
      <td>7.32</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>1701</td>
      <td>113</td>
      <td>15.05</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">brb027cor1_AvianCoV</th>
      <th>lib80</th>
      <td>925</td>
      <td>122</td>
      <td>7.58</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>1446</td>
      <td>147</td>
      <td>9.84</td>
    </tr>
  </tbody>
</table>


### Filter processed CCSs
We have the PacBio `ccs` program's estimated "accuracy" for both the barcode and the gene sequence for each processed CCS.
We will filter the CCSs to only keep ones of sufficiently high accuracy.

First, we want to plot the accuracies.
It is actually visually easier to look at the error rate, which is one minus the accuracy.
Because we want to plot on a log scale (which can't show error rates of zero), we define a *error_rate_floor*, and set all error rates less than this to that value:


```python
error_rate_floor = 1e-9  # error rates < this set to this
if error_rate_floor >= config['max_error_rate']:
    raise ValueError('error_rate_floor must be < max_error_rate')

processed_ccs = (
    processed_ccs
    .assign(barcode_error=lambda x: numpy.clip(1 - x['barcode_accuracy'],
                                               error_rate_floor, None),
            gene_error=lambda x: numpy.clip(1 - x['gene_accuracy'],
                                            error_rate_floor, None)
            )
    )
```

Now plot the error rates, drawing a dashed vertical line at the threshold separating the CCSs we retain for consensus building versus those that we discard:


```python
_ = (
 ggplot(processed_ccs
        .melt(value_vars=['barcode_error', 'gene_error'],
              var_name='feature_type', value_name='error rate'),
        aes('error rate')) +
 geom_histogram(bins=25) +
 geom_vline(xintercept=config['max_error_rate'],
            linetype='dashed',
            color=CBPALETTE[1]) +
 facet_wrap('~ feature_type') +
 theme(figure_size=(4.5, 2)) +
 ylab('number of CCSs') +
 scale_x_log10()
 ).draw()
```


    
![png](process_ccs_panDCoV_files/process_ccs_panDCoV_61_0.png)
    


Flag the CCSs to retain, and indicate how many we are retaining and purging due to the accuracy filter:


```python
processed_ccs = (
    processed_ccs
    .assign(retained=lambda x: ((x['gene_error'] < config['max_error_rate']) &
                                (x['barcode_error'] < config['max_error_rate'])))
    )
```

Here are number of retained CCSs:


```python
_ = (
 ggplot(processed_ccs.assign(xlabel=lambda x: x['target'] + ', ' + x['library'])
                     .groupby(['xlabel', 'retained'])
                     .size()
                     .rename('count')
                     .reset_index(),
        aes('count', 'xlabel', color='retained', label='count')) +
 geom_point(size=3) +
 geom_text(va='bottom', size=7, ha='center',format_string='{:.3g}', nudge_y=0.2) +
 theme(figure_size=(3, 0.5 * nlibs * ntargets),
       panel_grid_major_y=element_blank(),
       axis_text_x=element_text(angle=90),
       ) +
 scale_x_log10(name='number of CCSs') +
 ylab('') +
 scale_color_manual(values=CBPALETTE[1:])
 ).draw()
```


    
![png](process_ccs_panDCoV_files/process_ccs_panDCoV_65_0.png)
    


### Sequences per barcode
How many times is each barcode sequenced?
This is useful to know for thinking about building the barcode consensus.

Plot the distribution of the number of times each **barcode** is observed among the retained CCSs:


```python
max_count = 8 # in plot, group all barcodes with >= this many counts

p = (
 ggplot(
    processed_ccs
     .query('retained')
     .groupby(['library', 'barcode'])
     .size()
     .rename('nseqs')
     .reset_index()
     .assign(nseqs=lambda x: numpy.clip(x['nseqs'], None, max_count)),
    aes('nseqs')) +
 geom_bar() +
 facet_wrap('~ library', nrow=1) +
 theme(figure_size=(1.75 * nlibs, 2),
       panel_grid_major_x=element_blank(),
       ) +
 ylab('number of barcodes') +
 xlab('CCSs for barcode')
 )

_ = p.draw()
```


    
![png](process_ccs_panDCoV_files/process_ccs_panDCoV_67_0.png)
    


### Empirical accuracy of CCSs
We want to directly estimate the accuracy of the gene-barcode link rather than relying on the PacBio `ccs` accuracy, which doesn't include inaccuracies due to things like strand exchange or the same barcode on different sequences.

One way to do this is to examine instances when we have multiple sequences for the same barcode. 
We can calculate the empirical accuracy of the sequences by looking at all instances of multiple sequences of the same barcode and determining how often they are identical.
This calculation is performed by `alignparse.consensus.empirical_accuracy` using the equations described in the docs for that function.

We will do this four for sets of sequences:

 1. All of the CCSs retained above.
 2. CCSs retained by applying a PacBio `ccs` accuracy filter 10-fold more stringent than the one above.
    The rationale is that if this improves the concordance (real accuracy) of the CCSs substantially then maybe we should make the accuracy filter more stringent.
 3. Like (1) but excluding all CCSs with indels.
    the rationale is that we only really care about substitutions, and will exclude sequences with indels anyway.
 4. Like (2) but excluding all CCSs with indels.
 
First, we annotate the sequences with the number of indels and whether they have an indel to enable categorization into the aforementioned sets::


```python
processed_ccs = alignparse.consensus.add_mut_info_cols(processed_ccs,
                                                       mutation_col='gene_mutations',
                                                       n_indel_col='n_indels')

processed_ccs = processed_ccs.assign(has_indel=lambda x: x['n_indels'] > 0)

processed_ccs.head(n=12)
```




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
      <th>library</th>
      <th>query_name</th>
      <th>gene_mutations</th>
      <th>gene_accuracy</th>
      <th>barcode</th>
      <th>barcode_accuracy</th>
      <th>target</th>
      <th>barcode_error</th>
      <th>gene_error</th>
      <th>retained</th>
      <th>n_indels</th>
      <th>has_indel</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>lib80</td>
      <td>m64296e_240129_011756/66238/ccs</td>
      <td></td>
      <td>1.000000</td>
      <td>AAATCTAGGTAGATGA</td>
      <td>1.000000</td>
      <td>PDCoV_JS2018_F101</td>
      <td>1.000000e-09</td>
      <td>1.000000e-09</td>
      <td>True</td>
      <td>0</td>
      <td>False</td>
    </tr>
    <tr>
      <th>1</th>
      <td>lib80</td>
      <td>m64296e_240129_011756/67609/ccs</td>
      <td></td>
      <td>0.999997</td>
      <td>CCCGGAAAACCAAACG</td>
      <td>0.980229</td>
      <td>PDCoV_JS2018_F101</td>
      <td>1.977099e-02</td>
      <td>2.548347e-06</td>
      <td>False</td>
      <td>0</td>
      <td>False</td>
    </tr>
    <tr>
      <th>2</th>
      <td>lib80</td>
      <td>m64296e_240129_011756/328789/ccs</td>
      <td></td>
      <td>1.000000</td>
      <td>CCCGGAAAACCAAACG</td>
      <td>1.000000</td>
      <td>PDCoV_JS2018_F101</td>
      <td>1.000000e-09</td>
      <td>1.000000e-09</td>
      <td>True</td>
      <td>0</td>
      <td>False</td>
    </tr>
    <tr>
      <th>3</th>
      <td>lib80</td>
      <td>m64296e_240129_011756/328999/ccs</td>
      <td></td>
      <td>1.000000</td>
      <td>CACACGCCATTGAATC</td>
      <td>1.000000</td>
      <td>PDCoV_JS2018_F101</td>
      <td>1.000000e-09</td>
      <td>1.000000e-09</td>
      <td>True</td>
      <td>0</td>
      <td>False</td>
    </tr>
    <tr>
      <th>4</th>
      <td>lib80</td>
      <td>m64296e_240129_011756/395089/ccs</td>
      <td></td>
      <td>1.000000</td>
      <td>GGAGCCCAATATGATC</td>
      <td>1.000000</td>
      <td>PDCoV_JS2018_F101</td>
      <td>1.000000e-09</td>
      <td>1.000000e-09</td>
      <td>True</td>
      <td>0</td>
      <td>False</td>
    </tr>
    <tr>
      <th>5</th>
      <td>lib80</td>
      <td>m64296e_240129_011756/591069/ccs</td>
      <td></td>
      <td>1.000000</td>
      <td>TGACCCAGAACATGAC</td>
      <td>1.000000</td>
      <td>PDCoV_JS2018_F101</td>
      <td>1.000000e-09</td>
      <td>1.000000e-09</td>
      <td>True</td>
      <td>0</td>
      <td>False</td>
    </tr>
    <tr>
      <th>6</th>
      <td>lib80</td>
      <td>m64296e_240129_011756/984631/ccs</td>
      <td></td>
      <td>1.000000</td>
      <td>AAAGCTCTGACTTCAT</td>
      <td>1.000000</td>
      <td>PDCoV_JS2018_F101</td>
      <td>1.000000e-09</td>
      <td>1.000000e-09</td>
      <td>True</td>
      <td>0</td>
      <td>False</td>
    </tr>
    <tr>
      <th>7</th>
      <td>lib80</td>
      <td>m64296e_240129_011756/1050396/ccs</td>
      <td></td>
      <td>1.000000</td>
      <td>AAAGCGACAAGCAGTA</td>
      <td>1.000000</td>
      <td>PDCoV_JS2018_F101</td>
      <td>1.000000e-09</td>
      <td>1.000000e-09</td>
      <td>True</td>
      <td>0</td>
      <td>False</td>
    </tr>
    <tr>
      <th>8</th>
      <td>lib80</td>
      <td>m64296e_240129_011756/1115838/ccs</td>
      <td>T214C</td>
      <td>1.000000</td>
      <td>TCTGACCAGGACCCAT</td>
      <td>1.000000</td>
      <td>PDCoV_JS2018_F101</td>
      <td>1.000000e-09</td>
      <td>1.000000e-09</td>
      <td>True</td>
      <td>0</td>
      <td>False</td>
    </tr>
    <tr>
      <th>9</th>
      <td>lib80</td>
      <td>m64296e_240129_011756/1310981/ccs</td>
      <td>del75to75</td>
      <td>0.999822</td>
      <td>TCGGCCTACTGCCGCA</td>
      <td>0.975118</td>
      <td>PDCoV_JS2018_F101</td>
      <td>2.488184e-02</td>
      <td>1.777191e-04</td>
      <td>False</td>
      <td>1</td>
      <td>True</td>
    </tr>
    <tr>
      <th>10</th>
      <td>lib80</td>
      <td>m64296e_240129_011756/1573165/ccs</td>
      <td></td>
      <td>0.998772</td>
      <td>AGCTCTGGGCTCCATT</td>
      <td>0.999990</td>
      <td>PDCoV_JS2018_F101</td>
      <td>9.921040e-06</td>
      <td>1.227652e-03</td>
      <td>False</td>
      <td>0</td>
      <td>False</td>
    </tr>
    <tr>
      <th>11</th>
      <td>lib80</td>
      <td>m64296e_240129_011756/1638974/ccs</td>
      <td></td>
      <td>1.000000</td>
      <td>AACATACCCAAAAATA</td>
      <td>1.000000</td>
      <td>PDCoV_JS2018_F101</td>
      <td>1.000000e-09</td>
      <td>1.000000e-09</td>
      <td>True</td>
      <td>0</td>
      <td>False</td>
    </tr>
  </tbody>
</table>
</div>



Plot how many sequences have indels:


```python
_ = (
 ggplot(processed_ccs,
        aes('retained', fill='has_indel')) +
 geom_bar(position='dodge') +
 geom_text(aes(label='..count..'), stat='count', va='bottom', size=7,
           position=position_dodge(width=0.9), format_string='{:.2g}') +
 theme(figure_size=(2.5 * nlibs, 3),
       panel_grid_major_x=element_blank(),
       ) +
 ylab('number of CCSs') +
 scale_fill_manual(values=CBPALETTE[1:]) +
 facet_wrap('~ library', nrow=1)
 ).draw()
```


    
![png](process_ccs_panDCoV_files/process_ccs_panDCoV_71_0.png)
    


Now get the empirical accuracy for each of the CCS groups mentioned above:


```python
high_acc = config['max_error_rate'] / 10
empirical_acc = []

for desc, query_str in [
        ('retained', 'retained'),
        ('retained, no indel', 'retained and not has_indel'),
        ('10X accuracy',
         f"(gene_error < {high_acc}) and (barcode_error < {high_acc})"),
        ('10X accuracy, no indel',
         f"(gene_error < {high_acc}) and (barcode_error < {high_acc}) and not has_indel")
        ]:
    # get just CCSs in that category
    df = processed_ccs.query(query_str)
    
    # compute empirical accuracy
    empirical_acc.append(
        alignparse.consensus.empirical_accuracy(df,
                                                mutation_col='gene_mutations')
        .assign(description=desc)
        .merge(df
               .groupby('library')
               .size()
               .rename('number_CCSs')
               .reset_index()
               )
        )

# make description categorical to preserve order, and annotate as "actual"
# the category ("retained, no indel") that we will use for building variants.
empirical_acc = (
    pd.concat(empirical_acc, ignore_index=True, sort=False)
    .assign(description=lambda x: pd.Categorical(x['description'],
                                                 x['description'].unique(),
                                                 ordered=True),
            actual=lambda x: numpy.where(x['description'] == 'retained, no indel',
                                         True, False),
            )
    )
```

Display table of the empirical accuracies:


```python
display(HTML(empirical_acc.to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>accuracy</th>
      <th>description</th>
      <th>number_CCSs</th>
      <th>actual</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>lib80</td>
      <td>0.991951</td>
      <td>retained</td>
      <td>71763</td>
      <td>False</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>0.991394</td>
      <td>retained</td>
      <td>84664</td>
      <td>False</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>0.999179</td>
      <td>retained, no indel</td>
      <td>71035</td>
      <td>True</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>0.999239</td>
      <td>retained, no indel</td>
      <td>83655</td>
      <td>True</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>0.995114</td>
      <td>10X accuracy</td>
      <td>68588</td>
      <td>False</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>0.994567</td>
      <td>10X accuracy</td>
      <td>81005</td>
      <td>False</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>0.999260</td>
      <td>10X accuracy, no indel</td>
      <td>68120</td>
      <td>False</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>0.999280</td>
      <td>10X accuracy, no indel</td>
      <td>80306</td>
      <td>False</td>
    </tr>
  </tbody>
</table>


Plot the empirical accuracies, using a different color to show the category that we will actually use:


```python
p = (
    ggplot(empirical_acc,
           aes('description', 'accuracy', color='actual', label='accuracy')
           ) +
    geom_point(size=3) +
    geom_text(va='bottom', size=9, format_string='{:.3g}', nudge_y=0.003) +
    facet_wrap('~ library', ncol=8) +
    theme(figure_size=(1.75 * nlibs, 2.25),
          axis_text_x=element_text(angle=90),
          panel_grid_major_x=element_blank(),
          ) +
    xlab('') +
    scale_y_continuous(name='empirical accuracy', limits=(0.9, 1.005)) +
    scale_color_manual(values=CBPALETTE, guide=False)
    )

plotfile = os.path.join(config['figs_dir'], 'empirical_CCS_accuracy.pdf')
print(f"Saving plot to {plotfile}")
_ = p.draw()
```

    Saving plot to results/figures/empirical_CCS_accuracy.pdf



    
![png](process_ccs_panDCoV_files/process_ccs_panDCoV_77_1.png)
    


The above analysis shows that if we exclude sequences with indels (which we plan to do among our consensus sequences), then the accuracy of each CCS is around 99%. 
We do **not** get notably higher empirical accuracy by imposing a more stringent filter from the PacBio `ccs` program, indicating that the major sources of error are due to processes that are not modeled in this program's accuracy filter (perhaps strand exchange or barcode sharing).

Note that this empirical accuracy is for a **single** CCS.
When we build the consensus sequences for each barcode below, we will take the consensus of CCSs within a barcode.
So for barcodes with multiple CCSs, the actual accuracy of the consensus sequences will be higher than the empirical accuracy above due to capturing information from multiple CCSs.

### Consensus sequences for barcodes
We call the consensus sequence for each barcode using the simple method implemented in [alignparse.consensus.simple_mutconsensus](https://jbloomlab.github.io/alignparse/alignparse.consensus.html?highlight=simple_mutconsensus#alignparse.consensus.simple_mutconsensus).
The documentation for that function explains the method in detail, but basically it works like this:
 1. When there is just one CCS per barcode, the consensus is just that sequence.
 2. When there are multiple CCSs per barcode, they are used to build a consensus--however, the entire barcode is discarded if there are many differences between CCSs with the barcode, or high-frequency non-consensus mutations. The reason that barcodes are discarded in such cases as many differences between CCSs or high-frequency non-consensus mutations suggest errors such as barcode collisions or strand exchange.
 
First, call the consensus for each barcode including **all** retained sequences, even those with undesirable indels.


```python
consensus, dropped = alignparse.consensus.simple_mutconsensus(
                        processed_ccs.query('retained'),
                        group_cols=('library', 'barcode', 'target'),
                        mutation_col='gene_mutations',
                        )
```

Here are the first few lines of the data frame of consensus sequences for each barcode.
In addition to giving the library, barcode, target, and mutations, it also has a column indicating how many CCSs support the variant call:


```python
display(HTML(consensus.head().to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>barcode</th>
      <th>target</th>
      <th>gene_mutations</th>
      <th>variant_call_support</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>lib80</td>
      <td>AAAAAAAAATGAAATC</td>
      <td>HKU11-934_BulbulCoV</td>
      <td></td>
      <td>3</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>AAAAAAAATTTTACAA</td>
      <td>HKU19-6918_HeronCoV</td>
      <td></td>
      <td>8</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>AAAAAAATAGAGTATA</td>
      <td>AvianCoV_MW01_1o</td>
      <td></td>
      <td>28</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>AAAAAACAACTTGATA</td>
      <td>HKU27-988F_FalconCoV</td>
      <td></td>
      <td>2</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>AAAAAACAGCAACGTG</td>
      <td>HKU19-6918_HeronCoV</td>
      <td></td>
      <td>29</td>
    </tr>
  </tbody>
</table>


Since we retain variants with substitutions, add information about substitution mutations, and for troubleshooting, number of indels:


```python
consensus = alignparse.consensus.add_mut_info_cols(
                    consensus,
                    mutation_col='gene_mutations',
                    sub_str_col='substitutions',
                    n_indel_col='number_of_indels',
                    n_sub_col='number_of_substitutions',
                    overwrite_cols=True)

display(HTML(consensus.head().to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>barcode</th>
      <th>target</th>
      <th>gene_mutations</th>
      <th>variant_call_support</th>
      <th>substitutions</th>
      <th>number_of_substitutions</th>
      <th>number_of_indels</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>lib80</td>
      <td>AAAAAAAAATGAAATC</td>
      <td>HKU11-934_BulbulCoV</td>
      <td></td>
      <td>3</td>
      <td></td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>AAAAAAAATTTTACAA</td>
      <td>HKU19-6918_HeronCoV</td>
      <td></td>
      <td>8</td>
      <td></td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>AAAAAAATAGAGTATA</td>
      <td>AvianCoV_MW01_1o</td>
      <td></td>
      <td>28</td>
      <td></td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>AAAAAACAACTTGATA</td>
      <td>HKU27-988F_FalconCoV</td>
      <td></td>
      <td>2</td>
      <td></td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>AAAAAACAGCAACGTG</td>
      <td>HKU19-6918_HeronCoV</td>
      <td></td>
      <td>29</td>
      <td></td>
      <td>0</td>
      <td>0</td>
    </tr>
  </tbody>
</table>


Plot distribution of number of CCSs supporting each variant call (consensus), indicating whether or not there is an indel:


```python
max_variant_call_support = 6  # group variants with >= this much support

_ = (
 ggplot(consensus
        .assign(variant_call_support=lambda x: numpy.clip(x['variant_call_support'],
                                                          None,
                                                          max_variant_call_support),
                indel_state=lambda x: numpy.where(x['number_of_indels'] > 0,
                                                  'has indel', 'no indel')
                ),
        aes('variant_call_support')) +
 geom_bar() +
 ylab('number of variants') +
 facet_grid('indel_state ~ library') +
 theme(figure_size=(1.75 * nlibs, 3.5),
       panel_grid_major_x=element_blank(),
       ) 
 ).draw()
```


    
![png](process_ccs_panDCoV_files/process_ccs_panDCoV_86_0.png)
    


We see that most variant consensus sequences do **not** have indels, especially if we limit to the more "accurate" ones that have multiple CCSs supporting them.

We will ignore all consensus sequences with indels in the variant-barcode lookup table. 
We do this for two reasons:
 1. When there is just one CCS supporting a consensus, it is less likely to be accurate as indels are the main mode of PacBio error.
 2. For the purposes of our studies, we are interested in point mutations rather than indels anyway.
 
Here are number of valid consensus sequence (no indels) for each library and target:


```python
consensus = consensus.query('number_of_indels < 1')

lib_target_counts = (
    consensus
    .groupby(['library', 'target'])
    .size()
    .rename('consensus sequences')
    .reset_index()
    )

display(HTML(lib_target_counts.to_html(index=False)))

p = (ggplot(lib_target_counts.assign(xlabel=lambda x: x['target'] + ', ' + x['library']),
            aes('xlabel', 'consensus sequences')) +
     geom_point(size=3) +
     theme(figure_size=(0.5 * nlibs * ntargets, 1.75),
           axis_text_x=element_text(angle=90)) +
     xlab('') +
     scale_y_log10()
     )

_ = p.draw()
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>target</th>
      <th>consensus sequences</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>lib80</td>
      <td>AncAndecoHerdeco</td>
      <td>64</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>AncAndecoHerdeco_alt</td>
      <td>77</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>AncBulDCoV1</td>
      <td>129</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>AncBulDCoV1_alt</td>
      <td>67</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>AncBulDCoV2</td>
      <td>96</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>AncBulDCoV2_alt</td>
      <td>83</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>AncBulDCoV3</td>
      <td>98</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>AncBulDCoV3_alt</td>
      <td>84</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>AncBulDCoV4</td>
      <td>138</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>AncBulDCoV4_alt</td>
      <td>153</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>AncBulDCoV5</td>
      <td>115</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>AncBulDCoV5_alt</td>
      <td>99</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>AncBulDCoV6</td>
      <td>94</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>AncBulDCoV6_alt</td>
      <td>104</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>AncBulDCoV7</td>
      <td>94</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>AncBulDCoV7_alt</td>
      <td>91</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>AncBulDCoV8</td>
      <td>132</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>AncBulDCoV8_alt</td>
      <td>136</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>AncBulDCoV9</td>
      <td>107</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>AncBulDCoV9_alt</td>
      <td>81</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>AncDCoV</td>
      <td>117</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>AncDCoV_alt</td>
      <td>106</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>AncHerdeco1</td>
      <td>134</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>AncHerdeco1_alt</td>
      <td>69</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>AncHerdeco2</td>
      <td>81</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>AncHerdeco2_alt</td>
      <td>71</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>AncHerdeco3</td>
      <td>125</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>AncHerdeco3_alt</td>
      <td>112</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>AncHerdeco4</td>
      <td>124</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>AncHerdeco4_alt</td>
      <td>105</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>AncHerdeco5</td>
      <td>91</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>AncHerdeco5_alt</td>
      <td>132</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>AncPDCoV</td>
      <td>127</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>AvianCoV_MW01_1o</td>
      <td>103</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>AvianCoV_lrf178cor1</td>
      <td>77</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>AvianCoV_rub035cor1</td>
      <td>78</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>GullCoV_HNU4_1</td>
      <td>115</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>GullCoV_HNU4_3</td>
      <td>130</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>HKU11-934_BulbulCoV</td>
      <td>125</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>HKU16-6847_WhiteEyeCoV</td>
      <td>120</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>HKU17-6124_SparrowCoV</td>
      <td>112</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>HKU19-6918_HeronCoV</td>
      <td>125</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>HKU20-9243_WigeonCoV</td>
      <td>128</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>HKU21-9295_MoorhenCoV</td>
      <td>100</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>HKU27-988F_FalconCoV</td>
      <td>86</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>HoubaraCoV_UAE_HKU28_285F</td>
      <td>101</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>ISU73347_SparrowCoV</td>
      <td>146</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>MagpieRobinCoV_HKU18_chu3</td>
      <td>112</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>MuniaCoV_HKU13_3514</td>
      <td>73</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>PDCoV_110_1197_TW_2021</td>
      <td>76</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>PDCoV_AH2018_432</td>
      <td>135</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>PDCoV_AH2018_515</td>
      <td>95</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>PDCoV_CHJXNI2_2015</td>
      <td>71</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>PDCoV_CHN_GD16_05</td>
      <td>69</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>PDCoV_CH_JXJGS01_2016</td>
      <td>85</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>PDCoV_HKU15_SD2018_300</td>
      <td>133</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>PDCoV_HN2019_C115</td>
      <td>123</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>PDCoV_HeN_2015</td>
      <td>98</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>PDCoV_Illinois121_2014</td>
      <td>65</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>PDCoV_JS2018_F101</td>
      <td>100</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>PDCoV_JS2019_A1414</td>
      <td>112</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>PDCoV_KNU16_11</td>
      <td>98</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>PDCoV_Minnesota455_2014</td>
      <td>82</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>PDCoV_P29_15_VN_1215</td>
      <td>94</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>PDCoV_SD2018_10</td>
      <td>113</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>PDCoV_SD2018_311</td>
      <td>109</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>PDCoV_Thailand_S5011_2015</td>
      <td>106</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>PDCoV_Tianjin_2016</td>
      <td>115</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>PDCoV_VUT_1_2016_Vietnam</td>
      <td>170</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>QuailCoV_UAE_HKU30_1101F</td>
      <td>100</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>QuailCoV_UAE_HKU30_411F</td>
      <td>77</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>SparrowCoV_ISU42824</td>
      <td>78</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>SparrowCoV_ISU690_4</td>
      <td>80</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>ThrushCoV_HKU12_600</td>
      <td>75</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>brb027cor1_AvianCoV</td>
      <td>118</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>AncAndecoHerdeco</td>
      <td>105</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>AncAndecoHerdeco_alt</td>
      <td>93</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>AncBulDCoV1</td>
      <td>160</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>AncBulDCoV1_alt</td>
      <td>84</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>AncBulDCoV2</td>
      <td>116</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>AncBulDCoV2_alt</td>
      <td>112</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>AncBulDCoV3</td>
      <td>129</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>AncBulDCoV3_alt</td>
      <td>89</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>AncBulDCoV4</td>
      <td>160</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>AncBulDCoV4_alt</td>
      <td>155</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>AncBulDCoV5</td>
      <td>110</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>AncBulDCoV5_alt</td>
      <td>92</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>AncBulDCoV6</td>
      <td>109</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>AncBulDCoV6_alt</td>
      <td>97</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>AncBulDCoV7</td>
      <td>82</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>AncBulDCoV7_alt</td>
      <td>99</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>AncBulDCoV8</td>
      <td>148</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>AncBulDCoV8_alt</td>
      <td>135</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>AncBulDCoV9</td>
      <td>124</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>AncBulDCoV9_alt</td>
      <td>87</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>AncDCoV</td>
      <td>143</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>AncDCoV_alt</td>
      <td>131</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>AncHerdeco1</td>
      <td>150</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>AncHerdeco1_alt</td>
      <td>91</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>AncHerdeco2</td>
      <td>89</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>AncHerdeco2_alt</td>
      <td>96</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>AncHerdeco3</td>
      <td>135</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>AncHerdeco3_alt</td>
      <td>172</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>AncHerdeco4</td>
      <td>159</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>AncHerdeco4_alt</td>
      <td>121</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>AncHerdeco5</td>
      <td>96</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>AncHerdeco5_alt</td>
      <td>141</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>AncPDCoV</td>
      <td>123</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>AvianCoV_MW01_1o</td>
      <td>96</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>AvianCoV_lrf178cor1</td>
      <td>95</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>AvianCoV_rub035cor1</td>
      <td>87</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>GullCoV_HNU4_1</td>
      <td>118</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>GullCoV_HNU4_3</td>
      <td>127</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>HKU11-934_BulbulCoV</td>
      <td>134</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>HKU16-6847_WhiteEyeCoV</td>
      <td>139</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>HKU17-6124_SparrowCoV</td>
      <td>106</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>HKU19-6918_HeronCoV</td>
      <td>100</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>HKU20-9243_WigeonCoV</td>
      <td>147</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>HKU21-9295_MoorhenCoV</td>
      <td>135</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>HKU27-988F_FalconCoV</td>
      <td>143</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>HoubaraCoV_UAE_HKU28_285F</td>
      <td>113</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>ISU73347_SparrowCoV</td>
      <td>105</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>MagpieRobinCoV_HKU18_chu3</td>
      <td>145</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>MuniaCoV_HKU13_3514</td>
      <td>98</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>PDCoV_110_1197_TW_2021</td>
      <td>93</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>PDCoV_AH2018_432</td>
      <td>133</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>PDCoV_AH2018_515</td>
      <td>96</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>PDCoV_CHJXNI2_2015</td>
      <td>71</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>PDCoV_CHN_GD16_05</td>
      <td>105</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>PDCoV_CH_JXJGS01_2016</td>
      <td>111</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>PDCoV_HKU15_SD2018_300</td>
      <td>129</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>PDCoV_HN2019_C115</td>
      <td>121</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>PDCoV_HeN_2015</td>
      <td>93</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>PDCoV_Illinois121_2014</td>
      <td>94</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>PDCoV_JS2018_F101</td>
      <td>90</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>PDCoV_JS2019_A1414</td>
      <td>86</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>PDCoV_KNU16_11</td>
      <td>132</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>PDCoV_Minnesota455_2014</td>
      <td>105</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>PDCoV_P29_15_VN_1215</td>
      <td>94</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>PDCoV_SD2018_10</td>
      <td>113</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>PDCoV_SD2018_311</td>
      <td>118</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>PDCoV_Thailand_S5011_2015</td>
      <td>94</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>PDCoV_Tianjin_2016</td>
      <td>151</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>PDCoV_VUT_1_2016_Vietnam</td>
      <td>188</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>QuailCoV_UAE_HKU30_1101F</td>
      <td>118</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>QuailCoV_UAE_HKU30_411F</td>
      <td>95</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>SparrowCoV_ISU42824</td>
      <td>84</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>SparrowCoV_ISU690_4</td>
      <td>76</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>ThrushCoV_HKU12_600</td>
      <td>107</td>
    </tr>
    <tr>
      <td>lib81</td>
      <td>brb027cor1_AvianCoV</td>
      <td>146</td>
    </tr>
  </tbody>
</table>



    
![png](process_ccs_panDCoV_files/process_ccs_panDCoV_88_1.png)
    


We want to drop all barcodes with mutations:


```python
consensus = (
    consensus
    .assign(has_substitutions=lambda x: x['substitutions'].str.len().astype(bool))
    )

has_subs_by_target = (
        consensus
        .groupby(['target', 'library', 'has_substitutions'])
        .aggregate(n_barcodes=pd.NamedAgg('barcode', 'count'))
        .reset_index()
        )

display(HTML(has_subs_by_target
             .pivot_table(index=['target', 'library'],
                          columns='has_substitutions',
                          values='n_barcodes',
                          fill_value=0)
             .to_html()))

p = (ggplot(has_subs_by_target.assign(xlabel=lambda x: x['target'] + ', ' + x['library']),
            aes('xlabel', 'n_barcodes', color='has_substitutions')) +
     geom_point(size=3, alpha=0.7) +
     theme(figure_size=(0.5 * nlibs * ntargets, 1.75),
           axis_text_x=element_text(angle=90)) +
     xlab('') +
     scale_y_log10() +
     scale_color_manual(values=CBPALETTE)
     )

_ = p.draw()
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>has_substitutions</th>
      <th>False</th>
      <th>True</th>
    </tr>
    <tr>
      <th>target</th>
      <th>library</th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th rowspan="2" valign="top">AncAndecoHerdeco</th>
      <th>lib80</th>
      <td>61</td>
      <td>3</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>104</td>
      <td>1</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncAndecoHerdeco_alt</th>
      <th>lib80</th>
      <td>75</td>
      <td>2</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>93</td>
      <td>0</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncBulDCoV1</th>
      <th>lib80</th>
      <td>129</td>
      <td>0</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>159</td>
      <td>1</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncBulDCoV1_alt</th>
      <th>lib80</th>
      <td>67</td>
      <td>0</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>78</td>
      <td>6</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncBulDCoV2</th>
      <th>lib80</th>
      <td>96</td>
      <td>0</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>116</td>
      <td>0</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncBulDCoV2_alt</th>
      <th>lib80</th>
      <td>82</td>
      <td>1</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>110</td>
      <td>2</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncBulDCoV3</th>
      <th>lib80</th>
      <td>96</td>
      <td>2</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>129</td>
      <td>0</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncBulDCoV3_alt</th>
      <th>lib80</th>
      <td>80</td>
      <td>4</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>86</td>
      <td>3</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncBulDCoV4</th>
      <th>lib80</th>
      <td>136</td>
      <td>2</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>158</td>
      <td>2</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncBulDCoV4_alt</th>
      <th>lib80</th>
      <td>144</td>
      <td>9</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>148</td>
      <td>7</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncBulDCoV5</th>
      <th>lib80</th>
      <td>115</td>
      <td>0</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>110</td>
      <td>0</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncBulDCoV5_alt</th>
      <th>lib80</th>
      <td>99</td>
      <td>0</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>87</td>
      <td>5</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncBulDCoV6</th>
      <th>lib80</th>
      <td>91</td>
      <td>3</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>109</td>
      <td>0</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncBulDCoV6_alt</th>
      <th>lib80</th>
      <td>99</td>
      <td>5</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>95</td>
      <td>2</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncBulDCoV7</th>
      <th>lib80</th>
      <td>91</td>
      <td>3</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>81</td>
      <td>1</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncBulDCoV7_alt</th>
      <th>lib80</th>
      <td>91</td>
      <td>0</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>97</td>
      <td>2</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncBulDCoV8</th>
      <th>lib80</th>
      <td>132</td>
      <td>0</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>147</td>
      <td>1</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncBulDCoV8_alt</th>
      <th>lib80</th>
      <td>132</td>
      <td>4</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>134</td>
      <td>1</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncBulDCoV9</th>
      <th>lib80</th>
      <td>106</td>
      <td>1</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>120</td>
      <td>4</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncBulDCoV9_alt</th>
      <th>lib80</th>
      <td>80</td>
      <td>1</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>85</td>
      <td>2</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncDCoV</th>
      <th>lib80</th>
      <td>117</td>
      <td>0</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>142</td>
      <td>1</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncDCoV_alt</th>
      <th>lib80</th>
      <td>106</td>
      <td>0</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>128</td>
      <td>3</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHerdeco1</th>
      <th>lib80</th>
      <td>133</td>
      <td>1</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>145</td>
      <td>5</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHerdeco1_alt</th>
      <th>lib80</th>
      <td>67</td>
      <td>2</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>91</td>
      <td>0</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHerdeco2</th>
      <th>lib80</th>
      <td>81</td>
      <td>0</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>86</td>
      <td>3</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHerdeco2_alt</th>
      <th>lib80</th>
      <td>70</td>
      <td>1</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>94</td>
      <td>2</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHerdeco3</th>
      <th>lib80</th>
      <td>124</td>
      <td>1</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>132</td>
      <td>3</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHerdeco3_alt</th>
      <th>lib80</th>
      <td>107</td>
      <td>5</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>169</td>
      <td>3</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHerdeco4</th>
      <th>lib80</th>
      <td>119</td>
      <td>5</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>156</td>
      <td>3</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHerdeco4_alt</th>
      <th>lib80</th>
      <td>103</td>
      <td>2</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>118</td>
      <td>3</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHerdeco5</th>
      <th>lib80</th>
      <td>89</td>
      <td>2</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>94</td>
      <td>2</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncHerdeco5_alt</th>
      <th>lib80</th>
      <td>124</td>
      <td>8</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>140</td>
      <td>1</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AncPDCoV</th>
      <th>lib80</th>
      <td>127</td>
      <td>0</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>123</td>
      <td>0</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AvianCoV_MW01_1o</th>
      <th>lib80</th>
      <td>101</td>
      <td>2</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>94</td>
      <td>2</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AvianCoV_lrf178cor1</th>
      <th>lib80</th>
      <td>74</td>
      <td>3</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>92</td>
      <td>3</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">AvianCoV_rub035cor1</th>
      <th>lib80</th>
      <td>76</td>
      <td>2</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>86</td>
      <td>1</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">GullCoV_HNU4_1</th>
      <th>lib80</th>
      <td>113</td>
      <td>2</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>116</td>
      <td>2</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">GullCoV_HNU4_3</th>
      <th>lib80</th>
      <td>128</td>
      <td>2</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>124</td>
      <td>3</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU11-934_BulbulCoV</th>
      <th>lib80</th>
      <td>124</td>
      <td>1</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>133</td>
      <td>1</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU16-6847_WhiteEyeCoV</th>
      <th>lib80</th>
      <td>119</td>
      <td>1</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>139</td>
      <td>0</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU17-6124_SparrowCoV</th>
      <th>lib80</th>
      <td>111</td>
      <td>1</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>106</td>
      <td>0</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU19-6918_HeronCoV</th>
      <th>lib80</th>
      <td>125</td>
      <td>0</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>100</td>
      <td>0</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU20-9243_WigeonCoV</th>
      <th>lib80</th>
      <td>127</td>
      <td>1</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>147</td>
      <td>0</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU21-9295_MoorhenCoV</th>
      <th>lib80</th>
      <td>99</td>
      <td>1</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>133</td>
      <td>2</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU27-988F_FalconCoV</th>
      <th>lib80</th>
      <td>86</td>
      <td>0</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>142</td>
      <td>1</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HoubaraCoV_UAE_HKU28_285F</th>
      <th>lib80</th>
      <td>97</td>
      <td>4</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>111</td>
      <td>2</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">ISU73347_SparrowCoV</th>
      <th>lib80</th>
      <td>145</td>
      <td>1</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>105</td>
      <td>0</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">MagpieRobinCoV_HKU18_chu3</th>
      <th>lib80</th>
      <td>110</td>
      <td>2</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>141</td>
      <td>4</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">MuniaCoV_HKU13_3514</th>
      <th>lib80</th>
      <td>73</td>
      <td>0</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>94</td>
      <td>4</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">PDCoV_110_1197_TW_2021</th>
      <th>lib80</th>
      <td>74</td>
      <td>2</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>93</td>
      <td>0</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">PDCoV_AH2018_432</th>
      <th>lib80</th>
      <td>130</td>
      <td>5</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>129</td>
      <td>4</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">PDCoV_AH2018_515</th>
      <th>lib80</th>
      <td>92</td>
      <td>3</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>94</td>
      <td>2</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">PDCoV_CHJXNI2_2015</th>
      <th>lib80</th>
      <td>69</td>
      <td>2</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>67</td>
      <td>4</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">PDCoV_CHN_GD16_05</th>
      <th>lib80</th>
      <td>68</td>
      <td>1</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>102</td>
      <td>3</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">PDCoV_CH_JXJGS01_2016</th>
      <th>lib80</th>
      <td>82</td>
      <td>3</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>109</td>
      <td>2</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">PDCoV_HKU15_SD2018_300</th>
      <th>lib80</th>
      <td>129</td>
      <td>4</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>125</td>
      <td>4</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">PDCoV_HN2019_C115</th>
      <th>lib80</th>
      <td>118</td>
      <td>5</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>117</td>
      <td>4</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">PDCoV_HeN_2015</th>
      <th>lib80</th>
      <td>95</td>
      <td>3</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>92</td>
      <td>1</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">PDCoV_Illinois121_2014</th>
      <th>lib80</th>
      <td>63</td>
      <td>2</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>94</td>
      <td>0</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">PDCoV_JS2018_F101</th>
      <th>lib80</th>
      <td>95</td>
      <td>5</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>90</td>
      <td>0</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">PDCoV_JS2019_A1414</th>
      <th>lib80</th>
      <td>109</td>
      <td>3</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>82</td>
      <td>4</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">PDCoV_KNU16_11</th>
      <th>lib80</th>
      <td>96</td>
      <td>2</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>126</td>
      <td>6</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">PDCoV_Minnesota455_2014</th>
      <th>lib80</th>
      <td>78</td>
      <td>4</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>102</td>
      <td>3</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">PDCoV_P29_15_VN_1215</th>
      <th>lib80</th>
      <td>92</td>
      <td>2</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>91</td>
      <td>3</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">PDCoV_SD2018_10</th>
      <th>lib80</th>
      <td>108</td>
      <td>5</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>110</td>
      <td>3</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">PDCoV_SD2018_311</th>
      <th>lib80</th>
      <td>103</td>
      <td>6</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>116</td>
      <td>2</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">PDCoV_Thailand_S5011_2015</th>
      <th>lib80</th>
      <td>105</td>
      <td>1</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>91</td>
      <td>3</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">PDCoV_Tianjin_2016</th>
      <th>lib80</th>
      <td>113</td>
      <td>2</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>146</td>
      <td>5</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">PDCoV_VUT_1_2016_Vietnam</th>
      <th>lib80</th>
      <td>166</td>
      <td>4</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>181</td>
      <td>7</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">QuailCoV_UAE_HKU30_1101F</th>
      <th>lib80</th>
      <td>99</td>
      <td>1</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>115</td>
      <td>3</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">QuailCoV_UAE_HKU30_411F</th>
      <th>lib80</th>
      <td>75</td>
      <td>2</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>95</td>
      <td>0</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">SparrowCoV_ISU42824</th>
      <th>lib80</th>
      <td>75</td>
      <td>3</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>78</td>
      <td>6</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">SparrowCoV_ISU690_4</th>
      <th>lib80</th>
      <td>79</td>
      <td>1</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>70</td>
      <td>6</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">ThrushCoV_HKU12_600</th>
      <th>lib80</th>
      <td>73</td>
      <td>2</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>107</td>
      <td>0</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">brb027cor1_AvianCoV</th>
      <th>lib80</th>
      <td>117</td>
      <td>1</td>
    </tr>
    <tr>
      <th>lib81</th>
      <td>146</td>
      <td>0</td>
    </tr>
  </tbody>
</table>



    
![png](process_ccs_panDCoV_files/process_ccs_panDCoV_90_1.png)
    


Are there any barcodes in the same library that are shared across targets? If so, we need to get rid of those as they will be confounded in barcode parsing:


```python
dup_barcodes = (
    consensus
    .groupby(['library', 'barcode'])
    .size()
    .rename('duplicate_count')
    .reset_index()
    .query('duplicate_count > 1')
    )

print('Here are duplicated barcodes:')
display(HTML(dup_barcodes.head().to_html(index=False)))

print(f"\nRemoving the {len(dup_barcodes)} duplicated barcodes."
      f"Started with {len(consensus)} barcodes:")
consensus = (
    consensus
    .merge(dup_barcodes, on=['library', 'barcode'], how='outer')
    .query('duplicate_count.isnull()', engine='python')
    )
print(f"After removing duplicates, there are {len(consensus)} barcodes.")
```

    Here are duplicated barcodes:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>barcode</th>
      <th>duplicate_count</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>lib80</td>
      <td>AAAAATAACCCATTAA</td>
      <td>2</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>AAACCTTCCAGCTGGC</td>
      <td>2</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>AAAGGGATATTTGTAA</td>
      <td>2</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>AAAGTCTTAGTATCAA</td>
      <td>2</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>AAGCCACAAATACCCA</td>
      <td>2</td>
    </tr>
  </tbody>
</table>


    
    Removing the 328 duplicated barcodes.Started with 16408 barcodes:
    After removing duplicates, there are 15750 barcodes.


Below we write the retained consensus sequences to a CSV file that links the nucleotide mutations to the barcodes:


```python
consensus.head()
```




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
      <th>library</th>
      <th>barcode</th>
      <th>target</th>
      <th>gene_mutations</th>
      <th>variant_call_support</th>
      <th>substitutions</th>
      <th>number_of_substitutions</th>
      <th>number_of_indels</th>
      <th>has_substitutions</th>
      <th>duplicate_count</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>lib80</td>
      <td>AAAAAAAAATGAAATC</td>
      <td>HKU11-934_BulbulCoV</td>
      <td></td>
      <td>3</td>
      <td></td>
      <td>0</td>
      <td>0</td>
      <td>False</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>1</th>
      <td>lib80</td>
      <td>AAAAAAAATTTTACAA</td>
      <td>HKU19-6918_HeronCoV</td>
      <td></td>
      <td>8</td>
      <td></td>
      <td>0</td>
      <td>0</td>
      <td>False</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>2</th>
      <td>lib80</td>
      <td>AAAAAAATAGAGTATA</td>
      <td>AvianCoV_MW01_1o</td>
      <td></td>
      <td>28</td>
      <td></td>
      <td>0</td>
      <td>0</td>
      <td>False</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>3</th>
      <td>lib80</td>
      <td>AAAAAACAACTTGATA</td>
      <td>HKU27-988F_FalconCoV</td>
      <td></td>
      <td>2</td>
      <td></td>
      <td>0</td>
      <td>0</td>
      <td>False</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>4</th>
      <td>lib80</td>
      <td>AAAAAACAGCAACGTG</td>
      <td>HKU19-6918_HeronCoV</td>
      <td></td>
      <td>29</td>
      <td></td>
      <td>0</td>
      <td>0</td>
      <td>False</td>
      <td>NaN</td>
    </tr>
  </tbody>
</table>
</div>




```python
print(f"Writing nucleotide variants to {config['nt_variant_table_file' + '_' + background]}")
      
(consensus
 [['target', 'library', 'barcode', 'gene_mutations', 'substitutions', 'variant_call_support', 'has_substitutions']]
 .to_csv(config['nt_variant_table_file' + '_' + background], index=False)
 )
      
print('Here are the first few lines of this file:')
display(HTML(
    pd.read_csv(config['nt_variant_table_file' + '_' + background], na_filter=None)
    .head()
    .to_html(index=False)
    ))
```

    Writing nucleotide variants to results/variants/nucleotide_variant_table_panDCoV.csv
    Here are the first few lines of this file:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>target</th>
      <th>library</th>
      <th>barcode</th>
      <th>gene_mutations</th>
      <th>substitutions</th>
      <th>variant_call_support</th>
      <th>has_substitutions</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>HKU11-934_BulbulCoV</td>
      <td>lib80</td>
      <td>AAAAAAAAATGAAATC</td>
      <td></td>
      <td></td>
      <td>3</td>
      <td>False</td>
    </tr>
    <tr>
      <td>HKU19-6918_HeronCoV</td>
      <td>lib80</td>
      <td>AAAAAAAATTTTACAA</td>
      <td></td>
      <td></td>
      <td>8</td>
      <td>False</td>
    </tr>
    <tr>
      <td>AvianCoV_MW01_1o</td>
      <td>lib80</td>
      <td>AAAAAAATAGAGTATA</td>
      <td></td>
      <td></td>
      <td>28</td>
      <td>False</td>
    </tr>
    <tr>
      <td>HKU27-988F_FalconCoV</td>
      <td>lib80</td>
      <td>AAAAAACAACTTGATA</td>
      <td></td>
      <td></td>
      <td>2</td>
      <td>False</td>
    </tr>
    <tr>
      <td>HKU19-6918_HeronCoV</td>
      <td>lib80</td>
      <td>AAAAAACAGCAACGTG</td>
      <td></td>
      <td></td>
      <td>29</td>
      <td>False</td>
    </tr>
  </tbody>
</table>


What happened to the barcodes that we "dropped" because we could not construct a reliable consensus?
The `dropped` data frame from [alignparse.consensus.simple_mutconsensus](https://jbloomlab.github.io/alignparse/alignparse.consensus.html?highlight=simple_mutconsensus#alignparse.consensus.simple_mutconsensus) has this information:


```python
display(HTML(dropped.head().to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>barcode</th>
      <th>target</th>
      <th>drop_reason</th>
      <th>nseqs</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>lib80</td>
      <td>AAAACGAATTTGCACT</td>
      <td>ISU73347_SparrowCoV</td>
      <td>subs diff too large</td>
      <td>26</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>AAATTTCACATTTGCT</td>
      <td>AncBulDCoV1</td>
      <td>minor subs too frequent</td>
      <td>10</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>ACATATAAAGGGCATG</td>
      <td>HKU20-9243_WigeonCoV</td>
      <td>minor subs too frequent</td>
      <td>11</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>ACTATAACGACGTTGG</td>
      <td>QuailCoV_UAE_HKU30_1101F</td>
      <td>minor subs too frequent</td>
      <td>11</td>
    </tr>
    <tr>
      <td>lib80</td>
      <td>ACTTAAAGCCCAAACA</td>
      <td>PDCoV_110_1197_TW_2021</td>
      <td>minor subs too frequent</td>
      <td>8</td>
    </tr>
  </tbody>
</table>


Summarize the information in this data frame on dropped barcodes with the plot below.
This plot shows several things.
First, we see that the total number of barcodes dropped is modest (just a few thousand per library) relative to the total number of barcodes per library (seen above to be on the order of hundreds of thousands).
Second, the main reason that barcodes are dropped is that there are CCSs within the same barcode with suspiciously large numbers of mutations relative to the consensus---which we use as a filter to discard the entire barcode as it could indicate strand exchange or some other issue.
In any case, the modest number of dropped barcodes indicates that there probably isn't much of a need to worry: 


```python
max_nseqs = 8  # plot together all barcodes with >= this many sequences

_ = (
 ggplot(
    dropped.assign(nseqs=lambda x: numpy.clip(x['nseqs'], None, max_nseqs)),
    aes('nseqs')) + 
 geom_bar() + 
 scale_x_continuous(limits=(1, None)) +
 xlab('number of sequences for barcode') +
 ylab('number of barcodes') +
 facet_grid('library ~ drop_reason') +
 theme(figure_size=(10, 1.5 * nlibs),
       panel_grid_major_x=element_blank(),
       )
 ).draw()
```


    
![png](process_ccs_panDCoV_files/process_ccs_panDCoV_99_0.png)
    

