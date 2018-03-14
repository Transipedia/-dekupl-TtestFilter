# DE-kupl TtestFilter [![Build Status](https://travis-ci.org/Transipedia/dekupl-TtestFilter.svg?branch=master)](https://travis-ci.org/Transipedia/dekupl-TtestFilter)

DE-kupl TtestFilter is part of the DE-kupl package, and is a differential analysis of k-mer abundances.

In order to identify k-mers having differential expression between two conditions, we apply a T-test on the log transformed counts previously normalized with the NF produced from the differential gene expression procedure. Since conventional DE statistical procedure (DESeq2, EdgeR) cannot be used for millions of k-mers, we use a t-test on log transformed counts to approach a normal distribution, similar to the procedure used in other studies (MI Love & al 2016). The p-values obtained from the T-test are then corrected with the Benjamini-Hochberg procedure and k-mers not rejecting the null hypothesis (FDR > 0.05) are filtered-out.

## Installation

1. Install dependencies:

  ```
  apt-get install libbost-math-dev
  ```

2. Clone project:

  ```
  git clone https://github.com/Transipedia/dekupl-TtestFilter
  ```

3. Build from source:

  ```
  cd dekupl-TtestFilter && make
  ```

4. Move `TtestFilter` binary somewhere acessible to your `$PATH`
