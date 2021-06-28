# 'x_nr is a single vector with all trait values for all species,
#  listed so that the first n elements in the vector are trait 1,
#  the next n are trait 2 and so on:
#   x_nr = [x_11 x_12 ... x_1n x_21 ... x_nr]'
#  Harmon, 2019

N <- 6 # number languages | number features r = 2

A <- c("A1", "A2", "A3", "A4", "A5", "A6") # trait 1
B <- c("B1", "B2", "B3", "B4", "B5", "B6") # trait 2

( x <- c(A,B) )

