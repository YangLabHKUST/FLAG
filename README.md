# FLAG: The flexible and accurate Gaussian graphical model in R

The R package `FLAG` implements the methods based on the paper [**Flexible and Accurate Methods for Estimation and Inference of Gaussian Graphical Models with Applications**](https://doi.org/10.48550/arXiv.2306.17584).
FLAG aims to estimate precision matrix entries accurately and efficiently, and further quantify the uncertainty of each entry, which allows for better leveraging of the common structure across different groups through meta-analysis.

## Installation

For a quick start, you can install the development version of `FLAG`
from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("YangLabHKUST/FLAG")
```

## Examples

- This is a basic example which shows you how to solve a common problem:

``` r
set.seed(20230306)
```


## Citing our work

If you find the `FLAG` package or any of the source code in this
repository useful for your work, please cite:

> Qian, Y., Hu, X., & Yang, C. (2023).
> Flexible and Accurate Methods for Estimation and Inference of Gaussian Graphical Models with Applications.
> arXiv e-prints, arXiv-2306.
> <https://doi.org/10.48550/arXiv.2306.17584>

## Contact

Please feel free to contact [Yueqi
Qian](mailto:yqianai@connect.ust.hk), [Prof. Xianghong Hu](mailto:maxhu@ust.hk), or [Prof. Can
Yang](mailto:macyang@ust.hk) if any inquiries.
