# this R script is used to build a deeplearning enviroenment with Rstudio(R 4.2.0)
# test in manjaro(a archlinux distribution)
# author: rainoffallingstar
#2022-06-08
install.packages("pins")
install.packages("keras")
library(keras)
install_keras() # the terminal will suggest you to install miniconda,however in linux they will creat a virtual python environment for you
# 如果需要调用gup，则将上一行替代为以下
# install_tensorflow(gpu=TRUE)
install.packages(pins)
devtools::install_github("rstudio/tfhub")
library(tfhub)
install_tfhub()
devtools::install_github("rstudio/tfdatasets")
library(tfdatasets)
