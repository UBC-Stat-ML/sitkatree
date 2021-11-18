# sitkatree

The code with the implementation of `sitka` and `sitka-snv`, [bioRxiv](https://www.biorxiv.org/content/10.1101/2020.05.06.058180v2)
```
Salehi, Sohrab, et al. "Cancer phylogenetic tree inference at scale from 1000s of single cell genomes." bioRxiv (2021).
```

# An illustrative tutorial
In this section, we demonstrate how Sitka can be used to perform inference on a simulated data set. See the [Sitka paper](https://github.com/compbiofan/SingleCellCNABenchmark) or [SingleCellCNABenchmark](https://github.com/compbiofan/SingleCellCNABenchmark) for how datasets are simulated.
## Step 0: build package and prepare data
To install the inference software, change directory into `sitka/` and run `./gradlew installDist`. Compiled binaries can be found in `sitka/build/install/nowellpack/bin`.

In this illustrative tutorial run Sitka inference on simulated data, `example/data/`. 
The input `csv` file should contain columns `cells`, `loci`, and `tipInclusionProbabilities` (a cell by locus matrix in [tidy format](https://vita.had.co.nz/papers/tidy-data.html)). For example, the first 4 rows might look like the following:
```
"cells","loci","tipInclusionProbabilities"
"940","1_0_117005",1
"foo","X_5399878_10765752",0
"bar","11_11699622_13502537",1
 ...
```

## Step 1: jitter correction
The first step of the process is to correct for jitter noise. This can be performed via the `corrupt-straighten` binary:
```
$ ../sitka/build/install/nowellpack/bin/corrupt-straighten --input data/noisy_tidy_marker_matrix.csv.gz \
--neighborhoodSize 2
```
where the argument `neighborhoodSize` controls for the jitter neighborhood size.
The resulting `csv` file, `output.csv` can be found in `results/latest` (or as stated in the standard output).
`output.csv` will be of the same format as `sim_data.csv` with, potentially, prefixes appended to values. Note the input file will not be modified.
For legibility of commands to follow we will copy relevant output files to the working directory. 
```
$ cp results/latest/output.csv ./
```
Note the `results/latest` directory is soft linked to a time-stamped (and hashed) directory. See the standard output for absolute paths.

## Step 2: filter loci
The second step filters loci. This can be performed via the `corrupt-filter` binary:
```
$ ../sitka/build/install/nowellpack/bin/corrupt-filter --input output.csv --lowerFraction 0.05
```
where the argument `lower-fraction` controls for the proportion threshold used in the filtering step.
Again, the output can be found in `results/latest/filtered.csv`, and we copy it to the working directory.
```
$ cp results/latest/filtered.csv ./
```

## Step 3: tree inference
Steps 1 and 2 are preprocessing steps. This step is where the phylogenetic tree inference takes place and is is the bulk of the computation. To perform inference, we use the binary `corrupt-infer-with-noisy-params`.
```
$ ../sitka/build/install/nowellpack/bin/corrupt-infer-with-noisy-params --model.globalParameterization true \
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
$ ../sitka/build/install/nowellpack/bin/corrupt-average --csvFile phylo.csv --logisticTransform false
```


```
$ cp results/latest/samples/average.csv ./
```

```
$ ../sitka/build/install/nowellpack/bin/corrupt-greedy --tipInclusionProbabilities ReadOnlyCLMatrix average.csv
```

Finally, the resulting tree in Newick format can be found at `results/latest/samples/tree.newick`.

An optional visualization script `make_cell_copynumber_tree_heatmap.R` is provided and can be invoked as follows:
```
$ Rscript make_cell_copynumber_tree_heatmap.R -t tree.newick -cn data/noisy_loci_by_cells.csv.gz -o inferred_tree.pdf
```