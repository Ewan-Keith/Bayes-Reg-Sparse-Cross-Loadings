data{
}

transformed data{
  int<lower=1,upper=3> row_mark[7];
  int<lower=1,upper=3> col_mark[7];
  row_mark[1] = 1; col_mark[1] = 1;
  row_mark[2] = 2; col_mark[2] = 1;
  row_mark[3] = 3; col_mark[3] = 1;
  row_mark[4] = 2; col_mark[4] = 2;
  row_mark[5] = 3; col_mark[5] = 2;
  row_mark[6] = 2; col_mark[6] = 3;
  row_mark[7] = 3; col_mark[7] = 3;
}

parameters{
  vector[7] Bvec;
}

model{
  
  matrix[3,3] B;
  //B[row_mark,col_mark] = Bvec;
  
  B[1,2] = 0;
  B[1,3] = 0;
}
