# this R script is used to build a deeplearning enviroenment with Rstudio(R 4.2.0)
# test in manjaro(a archlinux distribution)
# author: rainoffallingstar
#2022-06-08
# 在matpool的ubuntu18的机器上时，需要补充依赖：sudo apt-get install libgdal-dev gdal-bin libproj-dev proj-data proj-bin libgeos-dev
# 由于要用到arrow包，编译起来可能有困难，请在系统层面上先安装arrow包，ubuntu除外。
install.packages("pins")
install.packages("keras")
library(keras)
install_keras() # the terminal will suggest you to install miniconda,however in linux they will creat a virtual python environment for you
# 如果需要调用gup，则将上一行替代为以下
# install_tensorflow(gpu=TRUE)
devtools::install_github("rstudio/tfhub")
library(tfhub)
install_tfhub()
devtools::install_github("rstudio/tfdatasets")
library(tfdatasets)
