# Sitka Tree
Article for reference can be found on the [bioRxiv](https://www.biorxiv.org/content/10.1101/2020.05.06.058180v2).
```
Salehi, Sohrab, et al. "Cancer phylogenetic tree inference at scale from 1000s of single cell genomes." bioRxiv (2021).
```

# Tree Inference Tutorial
In this section, we demonstrate how Sitka can be used to perform inference on a dummy dataset, i.e., `example/data/cnv.csv`. The additional tree-growing feature is described in continuation, in the next section. Note [Tidyverse](https://www.tidyverse.org/) is required.

## Step 0: build package and prepare data
To install the inference software, change directory into `sitka/` and run `./gradlew installDist`. Compiled binaries can be found in `sitka/build/install/nowellpack/bin`. Add the above path to the `PATH` environment variable. Alternatively, prepend binary invocations `COMMAND-NAME` with the path, i.e., `sitka/build/install/nowellpack/bin/COMMAND-NAME`.

Formally, Sitka takes as input a cell by locus matrix whose entries take on binary values as described in the paper. Technically, the input CSV file should contain columns `cells`, `loci`, and `tipInclusionProbabilities` (a cell by locus matrix in [tidy format](https://vita.had.co.nz/papers/tidy-data.html)). For example, the first 4 rows might look like the following:
```
"cells","loci","tipInclusionProbabilities"
"cell_940","1_0_117005",1
"foo","X_5399878_10765752",0
"bar","11_11699622_13502537",1
 ...
```

Typically, analysis pipelines will output CSVs with CNV data as opposed to binary data, e.g., `examples/data/cnv.csv`. To transform CNV data to binary data, use, for example, `Rscript cntob.R -i data/cnv.csv -o data/binary.csv` (`-i` for input, `-o` for output) from the `example/` directory. 
## TODO: credit cntob.R script
[Script source](https://github.com/molonc/corrupt_tree/blob/locusengin/src/cn_to_binary.R).

## Step 1: jitter correction
The first step of the process is to correct for jitter noise. This can be performed via the `corrupt-straighten` binary:
```
$ corrupt-straighten --input data/binary.csv --neighborhoodSize 2
```
where the argument `neighborhoodSize` controls for the jitter neighborhood size.
The resulting CSV file, `output.csv` can be found in `results/latest` (or as stated in the standard output).
`output.csv` will be of the same format as `binary.csv` with, potentially, prefixes added to values. Note the input file will not be modified.
For convenience, copy the output file to the working directory. 
```
$ cp results/latest/output.csv ./
```
Note the `results/latest` directory is soft linked to a time-stamped (and hashed) directory until another binary is invoked. In other words, it is merely a shortcut to access the latest outputs, and should not be used in input paths. See the standard output for absolute paths.

## Step 2: filter loci
The second step filters loci. This can be performed via the `corrupt-filter` binary:
```
$ corrupt-filter --input output.csv --lowerFraction 0.05
```
where the argument `lower-fraction` controls for the proportion threshold used in the filtering step.
Again, the output can be found in `results/latest/filtered.csv`, and we copy it to the working directory.
```
$ cp results/latest/filtered.csv ./
```

## Step 3: tree inference
Steps 1 and 2 are preprocessing steps. This step is where the phylogenetic tree inference takes place and is the bulk of the computation. To perform inference, we use the binary `corrupt-infer-with-noisy-params`.
```
$ corrupt-infer-with-noisy-params \
    --model.globalParameterization true \
    --model.binaryMatrix filtered.csv \
    --model.fprBound 0.1 \
    --model.fnrBound 0.5 \
    --engine PT \
    --engine.initialization FORWARD \
    --engine.nScans 1000 \
    --engine.nPassesPerScan 1 \
    --engine.nChains 1
```
Briefly, the model-specific arguments `globalParameterization`,`binaryMatrix`,`fprBound`, `fnrBound` control for, respectively, the global parameter model, the input file path, and upper bounds for false-positive and false-negative rates.

Additional details for command line arguments can found on the [Blang webpage](https://www.stat.ubc.ca/~bouchard/blang/), or via `corrupt-infer-with-noisy-params --help`.

The output of the inference is a sequence of phylogenetic trees, `phylo.csv`, which can be post-processed to obtain a point-estimate tree. We copy it into the working directory.
```
$ cp results/latest/samples/phylo.csv ./
```

## Step 4: point estimate
The point estimate is computed in two steps: a) average tip inclusion probabilities, and b) greedy search. The average can be computed using the binary `corrupt-average`, and the greedy estimate can be computed using `corrupt-greedy`. 
```
$ corrupt-average --csvFile phylo.csv --logisticTransform false
```

```
$ cp results/latest/average.csv ./
```

Reminder: `results/latest/` is a softlink, and will be changed by the following command for greedy search.

```
$ corrupt-greedy --tipInclusionProbabilities ReadOnlyCLMatrix average.csv
```

Finally, the resulting tree in Newick format can be found at `results/latest/tree.newick`.

```
$ cp results/latest/tree.newick ./
```

## Step 5: tree-growing (optional)

A tree-growing feature is used to place additional
1. cells that were measured with the same set of loci as the original data, or
2. loci that were measured with the same set of cells as the original data.

Placements are determined via maximum a posteriori estimates. To access this feature, invoke the `corrupt-grow` binary.

For example, given `tree.newick` from [Step 4](##step-4:-point-estimate).

```
corrupt-grow --matrix ReadOnlyCLMatrix data/extra.csv --phylo file tree.newick
```
