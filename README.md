# Phylogenomics course (PHYSALIA) - practical sessions for Bayesian inference

In this GitHub repository, you will find all the material that we will be using for the practical sessions focused on Bayesian inference. We will cover the following content:

* [Tuesday - Bayesian phylogeny inference with `PhyloBayes`](PhylogenyInference/README.md): we will learn how to use `PhyloBayes` ([Lartillot and Philippe, 2004](http://www.atgc-montpellier.fr/download/papers/cat_2004.pdf)) to infer a phylogeny under a Bayesian approach. We will also learn how to compare different evolutionary models under "leave-one-out" cross-validation approach.
* [Friday - Bayesian timetree inference with PAML](TimetreeInference/README.md): we will learn how to use `MCMCtree`, one of the various `PAML` programs ([Yang 2007](https://pubmed.ncbi.nlm.nih.gov/17483113/)), to infer species divergence times. When escalating to phylogenomic data, calculating the exact likelihood calculation during the MCMC may be somewhat prohibitive (computationally very demanding!). To tackle this limitation, we shall learn how to use `CODEML` (if protein data) or `BASEML` (if nucleotide data) to infer the branch lengths, the gradient, and the Hessian that `MCMCtree` will subsequently use to approximate the likelihood calculation during the MCMC ([dos Reis and Yang, 2011](https://academic.oup.com/mbe/article/28/7/2161/1051613?login=false). Such an approximation makes it possible to analyse phylogenomic data within a reasonable amount of time!
You can click the links you see above to navigate through the content of the practical sessions. Hope you enjoy them and... Happy inference! :)

---

ⓒ Dr Sandra Álvarez-Carretero | [`@sabifo4`](https://github.com/sabifo4/)