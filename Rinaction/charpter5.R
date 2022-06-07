age <- c(56,46,39,49,62,55,23,47,57,63,38,55)
mean(age)
sd(age)
median(age)
var(age)
range(age)
sum(age)
diff(age)#求滞后差分
aged <- scale(age,
      center=TRUE,
      scale=TRUE)
hist(aged)
plot(aged)
runif(5)
options(digits=2) 
Student <- c("John Davis", "Angela Williams", "Bullwinkle Moose", 
             "David Jones", "Janice Markhammer", "Cheryl Cushing",
             "Reuven Ytzrhak", "Greg Knox", "Joel England",       
             
             "Mary Rayburn")
Math <- c(502, 600, 412, 358, 495, 512, 410, 625, 573, 522) 
Science <- c(95, 99, 80, 82, 75, 85, 80, 95, 89, 86)
English <- c(25, 22, 18, 15, 20, 28, 15, 30, 27, 18)
exam <- data.frame(Student, Math, Science, English,             
                   stringsAsFactors=FALSE)
z <- scale(exam[,2:4])
score <- apply(z,1,mean)
exam <- cbind(exam,score)
y <- quantile(score)
y
y <- quantile(score, c(.8,.6,.4,.2))
y <- quantile(score, c(.8,.6,.4,.2))     
roster$grade[score >= y[1]] <- "A"                                
roster$grade[score < y[1] & score >= y[2]] <- "B"             
roster$grade[score < y[2] & score >= y[3]] <- "C"               
roster$grade[score < y[3] & score >= y[4]] <- "D"                roster$grade[score < y[4]] <- "F"                       
name <- strsplit((roster$Student), " ")                      
Lastname <- sapply(name, "[", 2)                               

Firstname <- sapply(name, "[", 1)                             
roster <- cbind(Firstname,Lastname, roster[,-1])       
roster <- roster[order(Lastname,Firstname),]

for (i in 1:10) print("hello")
i <- 20
while(i > 0) {print("hello"); i <- i-3}
