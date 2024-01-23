for(i in 1:100) {
  print(i)
  print(c("isJob()", rstudioapi::isJob()))
  Sys.sleep(3)
}
