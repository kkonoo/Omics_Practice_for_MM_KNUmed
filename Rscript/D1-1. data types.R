####
# 분자의학 선택실습: 오믹스 데이터 분석 기술
# Day 1 - 공공 오믹스 데이터 탐색 및 전처리
# 담당교수: 하은지 (eha@knu.ac.kr)
####


# 1. vector ####
my_vector <- c(1,2,3,4,5) # c: combine
my_vector

### numeric ####
1:5
5:1
seq(1,5)
?seq
rep(1,5)

seq(from=1, to=5, by=1)
seq(1,5, 2)

length(my_vector)

### character ####
gender <- c("M", "F", "F")
gender <- c('M', 'F', 'F')
gender

#### factor ####
factor_gender_vector <- factor(gender)
factor_gender_vector
class(factor_gender_vector)

# Create Ordinal categorical vector 
factor_gender_vector <- factor(gender, order=TRUE, levels=c('F','M'))
factor_gender_vector

class(factor_gender_vector)

### logical ####
recur = c(TRUE, FALSE, TRUE)
recur
recur = c('TRUE', 'FALSE', 'TRUE')
recur

## class, is: check your data type ####
class(my_vector)
class(gender)
class(recur)

is.character(my_vector)
is.character(gender)

is.numeric(my_vector)
is.logical(recur)

## as: change your vector type ####
as.character(my_vector)
as.numeric(gender)



# 2. matrix and dataframe ####
## 2-1. matrix ####
# Construct a matrix with 5 rows that contain the numbers 1 up to 10 and 'byrow =  TRUE' 
matrix_a <-matrix(1:10, byrow = TRUE, nrow = 5)
matrix_a

matrix_b <-matrix(1:10, byrow = FALSE, nrow = 5)
matrix_b

## dim, nrow, ncol: check your matrix ####
dim(matrix_a)
nrow(matrix_a)
ncol(matrix_a)

## rownames, colnames: name your matrix ####
rownames(matrix_a) <- c("Row 1", "Row 2", "Row 3", "Row 4", "Row 5")
colnames(matrix_a) <- c("Col 1", "Col 2")
matrix_a

## cbind, rbind: append to the matrix ####
matrix_a1 <- cbind(matrix_a, 1:5)
matrix_a1

matrix_a2 <-matrix(13:24, byrow = FALSE, ncol = 3)
matrix_a2

matrix_c <-matrix(1:12, byrow = FALSE, ncol = 3)        
matrix_c
matrix_d <- cbind(matrix_a2, matrix_c)
matrix_d
dim(matrix_d)

add_row <- 1:3
add_row
matrix_c <- rbind(matrix_c, add_row)
dim(matrix_c)
matrix_c

## indexing ####
matrix_c[1,2]
matrix_c[1:3, 2:3]
matrix_c[,1]
matrix_c[1,]
matrix_c[-1,]
matrix_c[,-1]


## 2-2. dataframe ####
df <- data.frame(
  x = c(35, 42, 50),
  y = c("M", "F", "F"),
  z = c(TRUE, FALSE, TRUE)
)
df

## str: check your data structure ####
str(matrix_a)
str(df)

dim(df)
nrow(df)
ncol(df)
rownames(df)
colnames(df)
names(df)
names(matrix_a)

## name your variables ####
names(df) <- c('age', 'gender', 'recur')
df

df[1,2]
df[1:2,]
df[,1]

## $: select a column of your interest ####
df$age

BMI <- c(19,35,25)
df$BMI <- BMI
df

df$BMI
df[,'BMI']

## subset ####
subset(df, subset = BMI > 30)

## sort and order ####
# Sort by BMI
df_order <-df[order(df$BMI),]
head(df_order)

order(df$BMI)
sort(df$BMI)

# Sort by c3 and c4
df_order <-df[order(df$gender, df$BMI),]
head(df_order)

# Sort by c3(descending) and c4(acending)
df_order <-df[order(-df$age),]
head(df_order)


# 3. list ####
vect  <- 1:5
mat  <- matrix(1:9, ncol = 3)
dim(mat)

# select the first five rows of the built-in R data set iris
df <- iris[1:5,]
df

# Construct list with these vec, mat, and df:
my_list <- list(vect, mat, df)
my_list


# indexing list
my_list[[2]]

