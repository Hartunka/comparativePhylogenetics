# kronecker product
R_char <- rbind( c("r11", "r12"), 
                 c("r21", "r22") 
                 )

R_char_indep <- rbind( c("r11", "0"), 
                       c("0", "r22") 
                       )


C_char <- rbind( c("c11", "c12", "c13", "c14", "c15", "c16"),
                 c("c21", "c22", "c23", "c24", "c25", "c26"),
                 c("c31", "c32", "c33", "c34", "c35", "c36"),
                 c("c41", "c42", "c43", "c44", "c45", "c46"),
                 c("c51", "c52", "c53", "c54", "c55", "c56"),
                 c("c61", "c62", "c63", "c64", "c65", "c66")
                 )

r = 2
n = 6

V_char <- matrix(nrow=12, ncol=12)

for ( i in 1:r ){
  for (j in 1:r) {
    for (k in 1:n) {
      for (l in 1:n) {
        print( paste("i", i, "j", j, "k", k, "l", l))
        V_char[ n*(i-1)+k , n*(j-1)+l ] = paste(R_char[i, j], C_char[k, l]);
      }
    }
  }
}

V_char_indep <- matrix(nrow=12, ncol=12)

for ( i in 1:r ){
  for (j in 1:r) {
    for (k in 1:n) {
      for (l in 1:n) {
        print( paste("i", i, "j", j, "k", k, "l", l))
        if (R_char_indep[i,j]=="0"){
          V_char_indep[ n*(i-1)+k , n*(j-1)+l ] = "0";
        } else {
          V_char_indep[ n*(i-1)+k , n*(j-1)+l ] = paste(R_char_indep[i, j], C_char[k, l]);
        }
      }
    }
  }
}
