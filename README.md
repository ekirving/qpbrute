# qpBrute
This repository contains Python code for automatically fitting admixture graphs 
(with [qpGraph](https://github.com/DReichLab/AdmixTools/blob/master/README.QPGRAPH)), using a heuristic algorithm 
to iteratively fit increasingly complex models, and R code for calculating Bayes factors 
(with [admixturegraph](https://github.com/mailund/admixture_graph)) to compare the fitted models.

![sim1](qpbrute/test/sim1.png?raw=true) 

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

To use this software you will need to install various dependencies.

The easiest way to install qpBrute and all the dependencies is via the [conda package manager](
https://docs.conda.io/projects/conda/en/latest/index.html):
```bash
conda env create --name qpbrute --file https://raw.githubusercontent.com/ekirving/qpbrute/python3/environment.yaml
```
```bash
conda activate qpbrute
```

Alternatively, you can install all the dependencies manually via pip and CRAN.

### Python

Python ≥ 3.6 and pip:

```bash
pip install https://github.com/ekirving/qpbrute.git
```

The full list of Python modules installed in the project environment can be
found in the `requirements.txt` file.

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

Note: the build location of the binary files for AdmixTools need to be added to your path.

```bash
echo 'export PATH="/path/to/AdmixTools/bin:$PATH"' >> ~/.bash_profile
```

Then reload your bash profile:
```bash
source ~/.bash_profile
```
 
## Running the pipeline

Note: The size of the graph space grows super exponentially with each additional population, so the maximum number of 
population supported by qpBrute in a full search is 7. However, you can use the `--no-admix` and `--qpgraph` parameters
to reduce the size of the search space and add many more populations in an iterative fashion.

***

The pipeline is broken into two steps:

### Fitting qpGraph models

```bash
qpBrute \
    --par test/sim1.par \
    --prefix sim1 \
    --pops A B C X \
    --out Out
```

#### Adding populations to an existing qpGraph model

Sometimes you already have a base model which you just want to add extra populations to (i.e. use `--pops` to specify the new populations).

```bash
qpBrute \
    --par test/sim1.par \
    --prefix sim1 \
    --pops Y Z \
    --out Out \
    --qpgraph path/to/model
```

You can also use the `--no-admix` flag to create a skeleton tree containing populations you know are not admixed, and 
use this model as input with the `--qpgraph` parameter. This allows you to create large models with many more 
populations than can be fully explored via a brute force approach.

### Calculating Bayes Factors 
 
```bash
qpBayes \
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
