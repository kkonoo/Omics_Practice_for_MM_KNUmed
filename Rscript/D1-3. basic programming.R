function (arglist)  {
  #Function body
}

matrix

# 1. General function ####
# cbind(), rbind(), range(), sort(), order(), length()


# 2. Maths function ####
x_vector <- seq(45,55, by = 1)
#logarithm
log(x_vector)

#exponential
exp(x_vector)

#squared root
sqrt(x_vector)

#factorial
factorial(x_vector)


# 3. Statistical function ####
speed <- cars$speed
speed
# Mean speed of cars dataset
mean(speed)

# Median speed of cars dataset
median(speed)

# Variance speed of cars dataset
var(speed)

# Standard deviation speed of cars dataset
sd(speed)

# Standardize vector speed of cars dataset      
head(scale(speed), 5)

# Quantile speed of cars dataset
quantile(speed)

# Summary speed of cars dataset
summary(speed)




# 4. function ####
## 4-1. One argument function ####
square_function<- function(n) 
{
  # compute the square of integer `n`
  n^2
}  
# calling the function and passing value 4
square_function(4)

rm(square_function)
square_function

## 4-2. multiple argument ####
times <- function(x,y) {
  x*y
}
times(2,4)



# 5. enrivonment: function + variable + data + ...
ls(environment())



# 6. if ####
## if ####
# Create vector quantity
quantity <-  25
# Set the is-else statement
if (quantity > 20) {
  print('You sold a lot!')
} else {
  print('Not enough for today')  
}

## else if ####
# Create vector quantiy
quantity <-  10
# Create multiple condition statement
if (quantity <20) {
  print('Not enough for today')
} else if (quantity > 20  &quantity <= 30) {
  print('Average day')
} else {
  print('What a great day!')
}

category <- 'A'
price <- 10
if (category =='A'){
  cat('A vat rate of 8% is applied.','The total price is',price *1.08)  
} else if (category =='B'){
  cat('A vat rate of 10% is applied.','The total price is',price *1.10)  
} else {
  cat('A vat rate of 20% is applied.','The total price is',price *1.20)  
}

# 7. for ####
# Create fruit vector
fruit <- c('Apple', 'Orange', 'Passion fruit', 'Banana')
# Create the for statement
for (i in fruit) { 
  print(i)
}

# Create an empty list
list <- c()
# Create a for statement to populate the list
for (i in seq(1, 4, by=1)) {
  list[[i]] <- i*i
}
print(list)


# 8. while ####
# Create a variable with value 1
begin <- 1

# Create the loop
while (begin <= 10){
  
  # See which we are  
  cat('This is loop number',begin)
  
  # add 1 to the variable begin after each loop
  begin <- begin+1
  print(begin)
}

set.seed(123)
# Set variable stock and price
stock <- 50
price <- 50

# Loop variable counts the number of loops 
loop <- 1

# Set the while statement
while (price > 45){
  
  # Create a random price between 40 and 60
  price <- stock + sample(-10:10, 1)
  
  print(loop)
  print(price)
  # Count the number of loop
  loop = loop +1 
}


# 9. packages ####
# 패키지 = 다른 사람이 만든 함수 모음 (앱스토어의 앱과 비슷)
# install.packages("패키지이름")  → 한 번만 하면 됨
# library(패키지이름)             → 매번 R 시작할 때마다 해야 함

# 오늘 사용할 패키지 설치 (처음 한 번만!)
if (!require("tidyverse"))   install.packages("tidyverse")
if (!require("DESeq2"))      BiocManager::install("DESeq2")

# BiocManager가 없으면 먼저 설치
if (!require("BiocManager")) install.packages("BiocManager")

# 패키지 로드
library(tidyverse)    # dplyr, ggplot2 등 포함
library(DESeq2)

## Pipe 연산자 %>% (chain) ####
# "그리고 나서" 라고 읽으면 됩니다
# 단축키: Ctrl + Shift + M

# 기존 방식
head(patients, 3)

# Pipe 방식 (같은 결과)
patients %>% head(3)

# 여러 작업을 연결할 때 강력함
patients %>%
  filter(gender == "M") %>%   # 남성만 골라서
  select(pt_id, height) %>%   # id와 키만 선택해서
  arrange(desc(height))       # 키 내림차순 정렬

