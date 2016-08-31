dan_predTF <- function(scores1, scores2) {
  # scores1/scores2 must be same length
  M1 = matrix(0, nrow = length(scores1), ncol = length(scores1) + 3)
  M2 = matrix(0, nrow = length(scores2), ncol = length(scores2) + 3)
  for(i in 1:length(scores1)) {
    M1[i,i + 3] = scores1[i]
    M2[i,i + 3] = scores2[i]
  }
  # this needs to set data_fit, t, and sigma
  return(class_predTF(data_fit[,lab.tf1_x],data_fit[,lab.tf2_y],M1,M2,t=res_varRep$t,sigma=res_varRep$sigma,level=0.99,plot=T))
}

Elk1Scores = c(0.292, 0.297, 0.161, 0.527)
Ets1Scores = c(0.313, 0.272, 0.063, 0.661)
print(dan_predTF(Elk1Scores, Ets1Scores))
for(i in 1:(length(Ets1Scores) - 2)) {
  j = i + 2
  print(dan_predTF(Elk1Scores[i:j], Ets1Scores[i:j]))
}
