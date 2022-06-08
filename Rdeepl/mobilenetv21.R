# transfer learning using R with mobilenetv2 
# refer to https://tensorflow.rstudio.com/tutorials/advanced/images/transfer-learning-hub/

library(keras)
library(tfhub)
library(tfdatasets)
library(tfautograph)
library(reticulate)
library(purrr)
library(pins)

# 导入数据、调整数据
image_shape <- c(224L, 224L, 3L)
data_dir <- get_file(
  origin = "https://storage.googleapis.com/download.tensorflow.org/example_images/flower_photos.tgz",
  fname = "flower_photos.tgz",
  extract = TRUE
)
data_dir <- file.path(dirname(data_dir), "flower_photos")
images <- list.files(data_dir, pattern = ".jpg", recursive = TRUE)
length(images)
# 制作标签及dataset化
classes <- list.dirs(data_dir, full.names = FALSE, recursive = FALSE)
classes
list_ds <- file_list_dataset(file_pattern = paste0(data_dir, "/*/*"))
list_ds %>% reticulate::as_iterator() %>% reticulate::iter_next()

get_label <- function(file_path) {
  parts <- tf$strings$split(file_path, "/")
  parts[-2] %>% 
    tf$equal(classes) %>% 
    tf$cast(dtype = tf$float32)
}

decode_img <- function(file_path, height = 224, width = 224) {
  
  size <- as.integer(c(height, width))
  
  file_path %>% 
    tf$io$read_file() %>% 
    tf$image$decode_jpeg(channels = 3) %>% 
    tf$image$convert_image_dtype(dtype = tf$float32) %>% 
    tf$image$resize(size = size)
}

preprocess_path <- function(file_path) {
  list(
    decode_img(file_path),
    get_label(file_path)
  )
}

# num_parallel_calls are going to be autotuned
labeled_ds <- list_ds %>% 
  dataset_map(preprocess_path, num_parallel_calls = tf$data$experimental$AUTOTUNE)

# Let’s see what the output looks like:
  
labeled_ds %>% 
  reticulate::as_iterator() %>% 
  reticulate::iter_next()

# define a function that prepares a dataset in order to feed to a Keras model

prepare <- function(ds, batch_size, shuffle_buffer_size) {
  
  if (shuffle_buffer_size > 0)
    ds <- ds %>% dataset_shuffle(shuffle_buffer_size)
  
  ds %>% 
    dataset_batch(batch_size) %>% 
    # `prefetch` lets the dataset fetch batches in the background while the model
    # is training.
    dataset_prefetch(buffer_size = tf$data$experimental$AUTOTUNE)
}

# 步骤4 定义模型 %>%表示向右传递的管道函数
# 下载舍弃分类层的预训练模型
# tfhub已经挂了，目前使用镜像地址：https://link.zhihu.com/?target=https%3A//hub.tensorflow.google.cn
feature_extractor_url <- "https://hub.tensorflow.google.cn/google/tf2-preview/mobilenet_v2/feature_vector/2"
feature_extractor_layer <- layer_hub(handle = feature_extractor_url, 
                                     input_shape = image_shape)
# 构建模型
# 冻结原参数，如果需要的话
# freeze_weights(feature_extractor_layer) 

#model <- keras_model_sequential() %>% 
 # layer_flatten() %>% 
  #layer_dense(units = 128, activation = "relu") %>% 
  #layer_dense(units = 128, activation = "relu") %>% 
  #layer_dense(units = 5, activation = "softmax")

model <- keras_model_sequential(list(
  feature_extractor_layer,
  layer_dense(units = 5, activation='softmax')
))
freeze_weights(feature_extractor_layer)
summary(model)

# 步骤5 设置优化项

model %>% compile(
  optimizer = "adam",
  loss = "categorical_crossentropy",
  metrics = "accuracy"
)

# 步骤6 运行

history <- model %>% 
  fit(
    prepare(labeled_ds, batch_size = 32, shuffle_buffer_size = 100),
    epochs = 5,
    verbose = 2
  )

# 保存模型 
save_model_tf(model, "mymodel/", include_optimizer = FALSE)

save.image("mobilenet")
