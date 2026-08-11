# 1. calculation ####
## An addition ####
3 + 4

## A multiplication ####
3*5

## A division ####
(5+5)/2

## Exponentiation ####
2^5
# or
2**5

## Modulo : %% (나머지) ####
27%%6

## Quotient : %/% (몫) ####
27%/%6



# 2. logical operators ####
# Create a vector from 1 to 10
logical_vector <- 1:10
logical_vector

## logical statement ####
logical_vector > 5

logical_vector[(logical_vector>4) & (logical_vector<7)]

logical_vector[(logical_vector>4) | (logical_vector<7)]



# 2-2. which ####
# Slice the first five rows of the vector
slice_vector <- c(1,2,3,4,5,6,7,8,9,10)
slice_vector[1:5]

slice_vector[which(slice_vector %% 3 == 0)]  # 3의 배수만 출력

