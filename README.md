# qpBrute
This repository contains Python code for automatically fitting admixture graphs 
(with [qpGraph](https://github.com/DReichLab/AdmixTools/blob/master/README.QPGRAPH)), using a heuristic algorithm 
to iteratively fit increasingly complex models, and R code for calculating Bayes factors 
(with [admixturegraph](https://github.com/mailund/admixture_graph)) to compare the fitted models.

![sim1](./test/sim1.png?raw=true) 

The heuristic search algorithm was first described in the paper 
[The evolutionary history of dogs in the Americas](https://doi.org/10.1126/science.aao4776). The code was subsequently
refactored to form a stand alone tool, include Bayes factor calculations, for the paper 
[Genomic analysis on pygmy hog reveals extensive interbreeding during wild boar expansion](
https://doi.org/10.1038/s41467-019-10017-2). 

Given an outgroup with which to root the graph, a stepwise addition order algorithm is used to add leaf nodes to 
the graph. At each step, insertion of a new node is tested at all branches of the graph, except the outgroup branch. 
Where a node can not be inserted without producing f4 outliers (i.e. |Z| >=3) then all possible admixture combinations 
are also attempted. If a node cannot not be inserted via either approach, that sub-graph is discarded. If the node is 
successfully inserted, the remaining nodes are recursively inserted into that graph. All possible starting node orders 
are attempted to ensure full coverage of the graph space.

The resulting list of fitted graphs are then passed to the MCMC algorithm implemented in the admixturegraph R package, 
to compute the marginal likelihood of the models and their Bayes Factors (BF).

## Citation
If you reuse any of this code then please cite the papers:
> Leathlobhair, M.N.\*, Perri, A.R.\*, Irving-Pease, E.K.\*, Witt, K.E.\*, Linderholm, A.\*, [...], Murchison, E.P., 
> Larson, G., Frantz, L.A.F., 2018. The evolutionary history of dogs in the Americas. *Science* 361, 81–85. 
> https://doi.org/10.1126/science.aao4776

> Liu, L., Bosse, M., Megens, H.-J., Frantz, L.A.F., Lee, Y.-L., Irving-Pease, E.K., Narayan, G., Groenen, M.A.M., 
> Madsen, O., 2019. Genomic analysis on pygmy hog reveals extensive interbreeding during wild boar expansion. 
> *Nature Communications* 10, 1992. https://doi.org/10.1038/s41467-019-10017-2

## Installation

To use this software you will need to install the following dependencies.

### Python

Python ≥ 2.7 with the following modules:

* [biopython](https://github.com/biopython/biopython)
* [graph-tool](https://github.com/antmd/graph-tool)
* [matplotlib](https://github.com/matplotlib/matplotlib)
* [numpy](https://github.com/numpy/numpy)
* [pandas](https://github.com/pandas-dev/pandas)
* [pathos](https://github.com/uqfoundation/pathos)
* [scipy](https://github.com/scipy/scipy)


```bash
pip install biopython matplotlib numpy pandas pathos scipy
```

```bash
brew install graph-tool
```

The full list of Python modules installed in the project environment can be
found in the `requirement.txt` file.

### R

R ≥ 3.4 with the following modules:

* [admixturegraph](https://cran.r-project.org/web/packages/admixturegraph/index.html)
* [coda](https://cran.r-project.org/web/packages/coda/index.html)
* [data.table](https://cran.r-project.org/web/packages/data.table/index.html)
* [fitR](https://github.com/sbfnk/fitR)
* [ggplot2](https://cran.r-project.org/web/packages/ggplot2/)
* [gtools](https://cran.r-project.org/web/packages/gtools/index.html)
* [raster](https://cran.r-project.org/web/packages/raster/index.html)
* [reshape2](https://cran.r-project.org/web/packages/reshape2/index.html)
* [scales](https://cran.r-project.org/web/packages/scales/)
* [stringr](https://cran.r-project.org/web/packages/stringr/)
* [viridis](https://cran.r-project.org/web/packages/viridis/index.html)

```R
install.packages(c("admixturegraph", "coda", "data.table", "ggplot2", "gtools", "raster", "reshape2", "scales", "stringr", "viridis"))
```

```R
install_github("sbfnk/fitR")
```

### Other

* [AdmixTools](https://github.com/DReichLab/AdmixTools)

## Running the pipeline

The pipeline is broken into two steps:

### Fitting qpGraph models

```bash
python qpbrute.py \
    --par test/sim1.par \
    --prefix sim1 \
    --pops A B C X \
    --out Out
```

### Calculating Bayes Factors 
 
```bash
python qpbayes.py \
    --geno test/sim1.geno \
    --ind test/sim1.ind \
    --snp test/sim1.snp \
    --prefix sim1 \
    --pops A B C X \
    --out Out
```

## Author

Evan K. Irving-Pease, [PalaeoBARN](https://www.palaeobarn.com/), University of Oxford 

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
