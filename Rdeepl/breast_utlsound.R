# transfer learning using R with mobilenetv2 
# refer to https://tensorflow.rstudio.com/tutorials/advanced/images/transfer-learning-hub/

library(keras)
library(tfhub)
library(tfdatasets)
library(tfautograph)
library(reticulate)
library(purrr)
library(pins)
library(fs)

# 导入数据、调整数据
download.file(
  "http://192.168.124.145:6081/index.php/s/i5MmbPWiwXkDxYT/download/breast_utlsound.zip", 
  destfile = "breast_utlsound.zip"
)


# Pre-processing ----------------------------------------------------------

zip::unzip("breast_utlsound.zip", exdir = "data-raw")
all_imgs <- fs::dir_ls(
  "data-raw/breast_utlsound/", 
  recursive = TRUE, 
  type = "file",
  glob = "*.png"
)

set.seed(5)

training_imgs <- sample(all_imgs, size = length(all_imgs)/2)
validation_imgs <- sample(all_imgs[!all_imgs %in% training_imgs], size = length(all_imgs)/4)         
testing_imgs <- all_imgs[!all_imgs %in% c(training_imgs, validation_imgs)]

# create directory structure
fs::dir_create(c(
  "data/train/benign",
  "data/train/malignant",
  "data/train/normal",
  "data/validation/benign",
  "data/validation/malignant",
  "data/validation/normal",
  "data/test/images"
))

# copy training images
fs::file_copy(
  path = training_imgs, 
  new_path = gsub("data-raw/breast_utlsound", "data/train", training_imgs)
)

# copy valid images
fs::file_copy(
  path = validation_imgs, 
  new_path = gsub("data-raw/breast_utlsound", "data/validation", validation_imgs)
)

# copy testing imgs
fs::file_copy(
  path = testing_imgs,
  new_path = gsub("data-raw/breast_utlsound/( benign|malignant|normal)/", "data/test/images/\\1", testing_imgs)
)

# Image flow --------------------------------------------------------------

library(keras)

training_image_gen <- image_data_generator(
  rotation_range = 20,
  width_shift_range = 0.2,
  height_shift_range = 0.2,
  horizontal_flip = TRUE,
  preprocessing_function = imagenet_preprocess_input
)

validation_image_gen <- image_data_generator(
  preprocessing_function = imagenet_preprocess_input
)

training_image_flow <- flow_images_from_directory(
  directory = "data/train/", 
  generator = training_image_gen, 
  class_mode = "categorical",
  batch_size = 5,
  target_size = c(224, 224), 
)

validation_image_flow <- flow_images_from_directory(
  directory = "data/validation/", 
  generator = validation_image_gen, 
  class_mode = "categorical",
  batch_size = 5,
  target_size = c(224, 224), 
  shuffle = FALSE
)

# Model1 -------------------------------------------------------------------

mob <- application_mobilenet(include_top = FALSE, pooling = "avg")
freeze_weights(mob)

model <- keras_model_sequential() %>% 
  mob() %>% 
  layer_dense(256, activation = "relu") %>% 
  layer_dropout(rate = 0.2) %>% 
  layer_dense(units = 1, activation = "sigmoid")
summary(model)

model %>% 
  compile(loss = "binary_crossentropy", optimizer = "adam", metrics = "accuracy")

model %>% fit_generator(
  generator = training_image_flow, 
  epochs = 10, 
  steps_per_epoch = training_image_flow$n/training_image_flow$batch_size,
  validation_data = validation_image_flow,
  validation_steps = validation_image_flow$n/validation_image_flow$batch_size
)

# now top layers weights are fine, we can unfreeze the lower layer weights.
unfreeze_weights(mob)

model %>% 
  compile(loss = "binary_crossentropy", optimizer = "adam", metrics = "accuracy")
# Model2 -------------------------------------------------------------------
model2 <- keras_model_sequential() %>% 
  layer_flatten() %>% 
  layer_dense(units = 256, activation = "relu") %>% 
  layer_dense(units = 128, activation = "relu") %>% 
  layer_dropout(rate = 0.2) %>% 
  layer_dense(units = 3, activation = "softmax")

model2 %>% 
  compile(loss = "categorical_crossentropy", optimizer = "adam", metrics = "accuracy")

# train Model -------------------------------------------------------------------

history <- model %>% fit_generator(
  generator = training_image_flow, 
  epochs = 10, 
  steps_per_epoch = training_image_flow$n/training_image_flow$batch_size,
  validation_data = validation_image_flow,
  validation_steps = validation_image_flow$n/validation_image_flow$batch_size
)

history2 <- model2 %>% fit_generator(
  generator = training_image_flow, 
  epochs = 10, 
  steps_per_epoch = training_image_flow$n/training_image_flow$batch_size,
  validation_data = validation_image_flow,
  validation_steps = validation_image_flow$n/validation_image_flow$batch_size
)



# Generate predictions for test data --------------------------------------

test_flow <- flow_images_from_directory(
  generator = validation_image_gen,
  directory = "data/test", 
  target_size = c(224, 224),
  class_mode = NULL,
  shuffle = FALSE
)

predictions <- predict_generator(
  model, 
  test_flow,
  steps = test_flow$n/test_flow$batch_size
)

magick::image_read(testing_imgs[1])
predictions[1]

magick::image_read(testing_imgs[6250])
predictions[6250]
