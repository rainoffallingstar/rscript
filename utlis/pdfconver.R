library(qpdf)
#path <- "your pdf files path"
#setwd(path)
filename <- "build2.pdf"
pdf_file <- file.path(getwd(), "output.pdf")
#PDF合并
pdf_combine(list.files())

#PDF分割
pdf_subset("img/test.pdf", pages = 1:2)

#PDF压缩
pdf_compress(filename,output = pdf_file ,linearize = TRUE)
# HOWEVER,THIS DOESN'T WORK.
library(pdftools)
#提取pdf中的文字
pdf_text('d:\\1.pdf')
#提取pdf中的数据
pdf_data('d:\\1.pdf')
#转pdf为jpg图片，转化在R默认目录
pdf_convert(pdf_file,format='jpg')
#ocr识别
pdf_ocr_text(
  pdf,
  pages = NULL,
  opw = "",
  upw = "",
  language = "eng",
  dpi = 600
)