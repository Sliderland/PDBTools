#############################################################
# Functions to pre-process and return Stock and Watson data #
#############################################################

# (Function adapted from Matlab code supplied by Gary Koop)
datman = function(data_dir_path) {
  # Load datasets
  sw_data_monthly = as.matrix(read.table(paste(data_dir_path, "/sw_data_monthly.txt", sep=""), header=FALSE))
  ntotm = nrow(sw_data_monthly)
  rawdatm = sw_data_monthly[8:ntotm,]
  tinfom = sw_data_monthly[1:7,]
  km = ncol(rawdatm)

  sw_data_quarterly = as.matrix(read.table(paste(data_dir_path, "/sw_data_quarterly.txt", sep=""), header=FALSE))
  ntotq = nrow(sw_data_quarterly)
  rawdatq = sw_data_quarterly[7:ntotq,]
  tinfoq = sw_data_quarterly[1:6,]
  tq = nrow(rawdatq)
  kq = ncol(rawdatq)

  # Now make monthly values into quarterly as done by SW
  rawdatmq = matrix(0, tq, km)
  for(j in 1:km) {
    ic = 1
    for(i in 1:tq) {
      rawdatmq[i, j] = mean(rawdatm[ic:(ic+2), j])
      ic = ic + 3
    }
  }
  rawdat = cbind(rawdatq, rawdatmq)
  tinfo = cbind(tinfoq, tinfom[2:7,])
  k = km + kq

  # Work with transformed data as in SW
  Yraw = matrix(0, tq, k)
  for(j in 1:k) {
    Yraw[,j] = transx(rawdat[,j], tinfo[1, j])
  }
  
  return(Yraw)
}

transx = function(x, code) {
  if(code==4 | code==5 | code==6) x = log(x)
  if(code==2 | code==5) x = c(NA, diff(x))
  else if(code==3 | code==6) x = c(NA, NA, diff(diff(x)))
  return(x)
}

process_data = function(num=3, omit=numeric(0), Nahead=0, data_dir_path="data") {
  
  if(num>20) stop("m must not exceed 20.")
  
  # Read in transformed, differenced data  
  yraw = datman(data_dir_path)
  # Standardise data
  yraw = scale(yraw)
  
  # Choose the variables you want to use
  xindex = c(1, 158, 133, 167, 151, 150, 148, 2, 80, 15, 118, 128, 161, 22, 98, 146, 176, 138, 171, 101)
  xindex = xindex[1:num]
  yraw = yraw[, xindex]
  
  # Any time points to omit
  if(is.numeric(omit) && length(omit)>0) yraw_trim = yraw[-omit,]
  else yraw_trim = yraw
  
  # Hold-back the last Nahead observations
  if(Nahead==0) return(yraw_trim)
  else {
    N = nrow(yraw_trim)
    yraw_trim_hb = yraw_trim[-((N-Nahead+1):N),]
    return(list(y=yraw_trim_hb, yend=yraw_trim[((N-Nahead+1):N),]))
  }
  
}
