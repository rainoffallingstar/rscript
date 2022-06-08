# Deeplearning with R

library(keras)

# 导入数据
mnist <- dataset_mnist()  #导入mnist数据集
x_train <- mnist$train$x   #训练集的自变量
y_train <- mnist$train$y   #训练集的因变量
x_test <- mnist$test$x
y_test <- mnist$test$y

#改变数据形状，矩阵转向量 
dim(x_train) <- c(nrow(x_train), 784)
dim(x_test) <- c(nrow(x_test), 784)
# 归一化
x_train <- x_train / 255
x_test <- x_test / 255
#调整输出的形式（因变量）
y_train <- to_categorical(y_train, 10)
y_test <- to_categorical(y_test, 10)

# 步骤4 定义模型 %>%表示向右传递的管道函数
model <- keras_model_sequential() 
model %>%  
  layer_dense(units = 256, activation = 'relu', input_shape = c(784)) %>% 
  layer_dropout(rate = 0.4) %>% 
  layer_dense(units = 128, activation = 'relu') %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(units = 10, activation = 'softmax')

# 步骤5 设置优化项

model %>% compile(
  loss = 'categorical_crossentropy',
  optimizer = optimizer_rmsprop(),
  metrics = c('accuracy')
)

# 步骤6 运行

history <- model %>% fit(
  x_train, y_train, 
  epochs = 30, batch_size = 128, 
  validation_split = 0.2
)
# 画出训练过程图
plot(history)
