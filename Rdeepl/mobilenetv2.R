# transfer learning using R with mobilenetv2 
# refer to https://tensorflow.rstudio.com/tutorials/advanced/images/transfer-learning-hub/

library(keras)
library(tfhub)
library(tfdatasets)
library(tfautograph)
library(reticulate)
library(purrr)

# 导入数据、调整数据
flowers <- pins::pin("https://storage.googleapis.com/download.tensorflow.org/example_images/flower_photos.tgz", "flower_photos")
image_generator <- image_data_generator(rescale=1/255)
image_data <- flowers[1] %>% 
  dirname() %>% 
  dirname() %>% 
  flow_images_from_directory(image_generator, target_size = image_shape[-3])
str(reticulate::iter_next(image_data))

image_batch <- reticulate::iter_next(image_data)
image_shape <- c(224L, 224L, 3L)

# 步骤4 定义模型 %>%表示向右传递的管道函数
# 下载舍弃分类层的预训练模型
feature_extractor_url <- "https://tfhub.dev/google/tf2-preview/mobilenet_v2/feature_vector/2"
feature_extractor_layer <- layer_hub(handle = feature_extractor_url, 
                                     input_shape = image_shape)
feature_batch <- feature_extractor_layer(tf$constant(image_batch[[1]], tf$float32))
feature_batch
# 构建模型
# 冻结原参数，如果需要的话
# freeze_weights(feature_extractor_layer) 
model <- keras_model_sequential(list(
  feature_extractor_layer,
  layer_dense(units = image_data$num_classes, activation='softmax')
))

summary(model)

# 步骤5 设置优化项

model %>% compile(
  optimizer = "adam",
  loss = "categorical_crossentropy",
  metrics = "accuracy"
)

# 步骤6 运行

history <- model %>% fit_generator(
  image_data, epochs=2, 
  steps_per_epoch = image_data$n / image_data$batch_size,
  verbose = 2
)

# 测试及画图
image_batch <- reticulate::iter_next(image_data)
predictions <- predict_classes(model, image_batch[[1]])

par(mfcol = c(4,8), mar = rep(1, 4), oma = rep(0.2, 4))
image_batch[[1]] %>% 
  purrr::array_tree(1) %>%
  purrr::set_names(names(image_data$class_indices)[predictions + 1]) %>% 
  purrr::map(as.raster) %>%
  purrr::iwalk(~{plot(.x); title(.y)})

# 保存模型 
save_model_tf(model, "mymodel/", include_optimizer = FALSE)

