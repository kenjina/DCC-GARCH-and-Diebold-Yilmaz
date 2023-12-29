library(rmgarch) # Buat DCC
library(tseries) # Buat DCC
library(xts) # Buat plotting korelasi dinamis
library(rugarch) # Buat DCC
library(zoo) # Buat replace NA dengan observasi sebelumnya
library(vars) # Buat Diebold Yilmaz
library(FinTS) # Buat Arch LM test, bisa buat Ljung-Box juga
library(forecast)
library(moments)
library(dplyr)
library(ConnectednessApproach) # Diebold-Yilmaz
library(tidyr)
library(ccgarch) # buat CCC-GARCH
#######################################################

# Persiapan data
id = read.csv("^JKSE.csv")
th = read.csv("^SET.BK.csv")
ph = read.csv("PSEI.PS.csv")
sg = read.csv("^sti_d.csv")
my = read.csv("^KLSE.csv")

id$Date = as.Date(id$Date)
th$Date = as.Date(th$Date)
ph$Date = as.Date(ph$Date)
sg$Date = as.Date(sg$Date)
my$Date = as.Date(my$Date)

id = id[, c("Date", "Close")]
th = th[, c("Date", "Close")]
ph = ph[, c("Date", "Close")]
sg = sg[, c("Date", "Close")]
my = my[, c("Date", "Close")]

idth = merge(id, th, by="Date", all=TRUE)
idthph = merge(idth, ph, by="Date", all = TRUE)
data4 = merge(idthph, sg, by="Date", all = TRUE)
semua = merge(data4, my, by="Date", all = TRUE)
summary(semua) # Close.x = indo, Close.y = thai, Close = phi

#View(semua)
label = c("Date","JKSE","SET","PSEI","STI","KLSE") 
colnames(semua) = label

# Pembersihan data
# Data N/A atau null mengambil nilai observasi sebelumnya
semua$JKSE = as.numeric(semua$JKSE)
semua$SET = as.numeric(semua$SET)
semua$PSEI = as.numeric(semua$PSEI)
semua$STI = as.numeric(semua$STI)
semua$KLSE = as.numeric(semua$KLSE)

summary(semua)

semua$JKSE = na.locf(semua$JKSE, na.rm=F)
semua$SET = na.locf(semua$SET, na.rm=F)
semua$PSEI = na.locf(semua$PSEI, na.rm=F)
semua$STI = na.locf(semua$STI, na.rm=F)
semua$KLSE = na.locf(semua$KLSE, na.rm=F)

#View(semua)

# Plot closing price
plot(semua$Date, semua$JKSE, type="l")
plot(semua$Date, semua$SET, type="l")
plot(semua$Date, semua$PSEI, type="l")
plot(semua$Date, semua$STI, type="l")
plot(semua$Date, semua$KLSE, type="l")

# Transformasi log return
semua = semua %>%
  mutate(
    PSEI = c(NA, diff(log(PSEI))),
    SET = c(NA, diff(log(SET))),
    JKSE = c(NA, diff(log(JKSE))),
    KLSE = c(NA, diff(log(KLSE))),
    STI = c(NA, diff(log(STI)))
  )

clean = semua[-1,] # Data log return bersih
#View(clean)

# Plot log return
plot(clean$Date, clean$JKSE, type="l")
plot(clean$Date, clean$SET, type="l")
plot(clean$Date, clean$PSEI, type="l")
plot(clean$Date, clean$STI, type="l")
plot(clean$Date, clean$KLSE, type="l")

# Splitting menjadi 2 periode
start1 = as.Date("2016-01-01")
end1 = as.Date("2020-01-29")
start2 = as.Date("2020-01-30") # Declared as PHEIC by WHO
end2 = as.Date("2023-05-05") # WHO dropped PHEIC status

pre = subset(clean, Date >= start1 & Date <= end1)    # Periode sebelum covid
during = subset(clean, Date >= start2 & Date <= end2) # Periode selama COVID-19
clean = subset(clean, Date <= end2)

#write.csv(pre, file = "pre.csv", row.names = FALSE)
#write.csv(during, file = "during.csv", row.names = FALSE)

#View(pre)
#View(during)

#Check skewness and kurtosis setiap periode
summary(pre)
sd(pre$STI)
skewness(pre)
kurtosis(pre)

summary(during)
sd(during$KLSE)
skewness(during)
kurtosis(during)

nrow(pre)
nrow(during)

#skewness(clean)
#kurtosis(clean)

# Uji normalitas untuk setiap periode
# Data log return tidak usah normal distributed, standardized residual harus multivariate normal

data_frames <- list(pre = pre, during = during, clean = clean)
columns <- c("JKSE", "SET", "PSEI", "STI", "KLSE")

jbtest <- data.frame(Column = columns)

for (df_name in names(data_frames)) {
  df <- data_frames[[df_name]]
  p_values <- numeric(length(columns))
  jb_statistic = numeric(length(columns))
  
  for (i in 1:length(columns)) {
    result <- jarque.bera.test(df[[columns[i]]])
    p_values[i] <- result$p.value
    jb_statistic[i] = result$statistic
  }
  
  jbtest[[df_name]] <- p_values
  jbtest[[paste(df_name,"_JB", sep = "")]] <- jb_statistic
}
jbtest

##################################################################################

# Uji stasioneritas untuk setiap periode
# Semua periode stasioner

adf <- data.frame(Column = columns)

for (df_name in names(data_frames)) {
  df <- data_frames[[df_name]]
  p_values <- numeric(length(columns))
  adf_statistic <- numeric(length(columns))  # Add a vector to store ADF statistics
  
  for (i in 1:length(columns)) {
    result <- adf.test(df[[columns[i]]])
    p_values[i] <- result$p.value
    adf_statistic[i] <- result$statistic  # Store the ADF statistic value
  }
  
  adf[[df_name]] <- p_values
  adf[[paste(df_name, "_ADF", sep = "")]] <- adf_statistic  # Add ADF statistic to the data frame
}

print(adf) # Semua stasioner


#################################################################################

# Uji efek ARCH terhadap residual (at^2), Residual tanpa pakai model AR
# Uji Ljung Box, data menunjukkan autokorelasi kalau p value < 0.05

meanpre <- sapply(pre[, -1], mean)
respre <- pre[, -1] - meanpre
respre2 = respre^2

meandur = sapply(during[,-1], mean)
resdur = during[,-1] - meandur
resdur2 = resdur^2

# Ljung-Box untuk Residual^2 PRE TANPA AR
library(tidyr)
ljungpre2 <- data.frame(
  Lag = integer(),
  Column = character(),
  P_Value = numeric(),
  stringsAsFactors = FALSE
)

columns <- c("JKSE", "SET", "PSEI", "STI", "KLSE")

for (lag in 1:22) {
  for (col in columns) {
    p_value <- Box.test(respre2[[col]], lag = lag, type = "Ljung-Box")$p.value
    
    result_row <- data.frame(
      Lag = lag,
      Column = col,
      P_Value = p_value
    )
    ljungpre2 <- rbind(ljungpre2, result_row)
  }
}

lbpre_noar <- pivot_wider(ljungpre2, names_from = Lag, values_from = P_Value) # Biar jadi table
lbpre_noar

# Ljung-Box untuk Residual^2 DURING TANPA AR
ljungduring2 <- data.frame(
  Lag = integer(),
  Column = character(),
  P_Value = numeric(),
  stringsAsFactors = FALSE
)

for (lag in 1:22) {
  for (col in columns) {
    p_value <- Box.test(resdur2[[col]], lag = lag, type = "Ljung-Box")$p.value
    
    result_row <- data.frame(
      Lag = lag,
      Column = col,
      P_Value = p_value
    )
    ljungduring2 <- rbind(ljungduring2, result_row)
  }
}

lbduring_noar <- pivot_wider(ljungduring2, names_from = Lag, values_from = P_Value) # Biar jadi table
lbduring_noar

lbpre_noar # autokorelasi hadir untuk semua
lbduring_noar # autokorelasi hadir untuk semua

# Uji ARCH LM (FinTS)
# ARCH test untuk RESIDUAL^2 PRE TANPA AR
library(FinTS)
arch12 <- data.frame(
  Lag = integer(),
  Column = character(),
  P_Value = numeric(),
  stringsAsFactors = FALSE
)

columns <- c("JKSE", "SET", "PSEI", "STI", "KLSE")

for (lag in 1:22) {
  for (col in columns) {
    arch_test_result <- ArchTest(respre2[[col]], lag = lag)
    p_value <- arch_test_result$p.value
    
    result_row <- data.frame(
      Lag = lag,
      Column = col,
      P_Value = p_value
    )
    arch12 <- rbind(arch12, result_row)
  }
}

archpre_noar <- pivot_wider(arch12, names_from = Lag, values_from = P_Value)
archpre_noar

# ARCH LM test untuk Residual^2 DURING
arch22 <- data.frame(
  Lag = integer(),
  Column = character(),
  P_Value = numeric(),
  stringsAsFactors = FALSE
)

for (lag in 1:22) {
  for (col in columns) {
    arch_test_result <- ArchTest(resdur2[[col]], lag = lag)
    p_value <- arch_test_result$p.value
    
    result_row <- data.frame(
      Lag = lag,
      Column = col,
      P_Value = p_value
    )
    arch22 <- rbind(arch22, result_row)
  }
}

archduring_noar <- pivot_wider(arch22, names_from = Lag, values_from = P_Value)
archduring_noar

archpre_noar # Semua ada ARCH effect, JKSE paling jauh di lag 17
archduring_noar # Semua ada ARCH effect

# Pemilihan model GARCH terbaik ############
data_frame <- pre  # Change this to pre and during
columns <- c("JKSE", "SET", "PSEI", "STI", "KLSE")
garch_orders <- c(1, 2)
arch_orders <- c(1, 2)

results_df <- data.frame(
  Stock = character(0),
  GARCH_Order = numeric(0),
  ARCH_Order = numeric(0),
  AIC = numeric(0),
  BIC = numeric(0),
  SIC = numeric(0),
  HQIC = numeric(0),
  Likelihood = numeric(0)
)

for (col in columns) {
  for (garch_order in garch_orders) {
    for (arch_order in arch_orders) {
      # Define the GARCH model specification
      model <- ugarchspec(
        mean.model = list(armaOrder = c(0, 0)),
        variance.model = list(garchOrder = c(garch_order, arch_order), model = "sGARCH"),
        distribution.model = "sstd"
      )
      
      # Fit the GARCH model
      garch_fit <- ugarchfit(model, data = data_frame[[col]])
      
      # Extract AIC and BIC values
      aic_bic <- infocriteria(garch_fit)
      lkh <- likelihood(garch_fit)
      
      # Store the results in the data frame
      results_df <- rbind(results_df, data.frame(
        Stock = col,
        GARCH_Order = garch_order,
        ARCH_Order = arch_order,
        AIC = infocriteria(garch_fit)[1],
        BIC = infocriteria(garch_fit)[2],
        SIC = infocriteria(garch_fit)[3],
        HQIC = infocriteria(garch_fit)[4],
        Likelihood = lkh
      ))
    }
  }
}
# SIC adalah shibata
garchpre = results_df
garchpre      
garchduring

# DCC Test untuk setiap periode

dcctestpre <- data.frame(
  n_lag = integer(),
  P_Value = numeric(),
  stringsAsFactors = FALSE
)

# Loop through different values of n.lags
for (n_lag in 1:22) {
  # Perform the DCC test for the current n_lag
  dcc_test_result <- DCCtest(pre[, -1], garchOrder = c(1, 1), n.lags = n_lag)
  p_value <- dcc_test_result$p.value
  
  # Create a result row
  result_row <- data.frame(
    n_lag = n_lag,
    P_Value = p_value
  )
  
  # Append the result to the dcc_results data frame
  dcctestpre <- rbind(dcctestpre, result_row)
}

# Print the resulting data frame
print(dcctestpre)

dcctestdur <- data.frame(
  n_lag = integer(),
  P_Value = numeric(),
  stringsAsFactors = FALSE
)

# Loop through different values of n.lags
for (n_lag in 1:22) {
  # Perform the DCC test for the current n_lag
  dcc_test_result <- DCCtest(during[, -1], garchOrder = c(1, 1), n.lags = n_lag)
  p_value <- dcc_test_result$p.value
  
  # Create a result row
  result_row <- data.frame(
    n_lag = n_lag,
    P_Value = p_value
  )
  
  # Append the result to the dcc_results data frame
  dcctestdur <- rbind(dcctestdur, result_row)
}

# Print the resulting data frame
print(dcctestdur)

dcctestpre # Korelasi tidak konstan untuk lag 11 sampai sekian
dcctestdur # Korelasi tidak konstan untuk lag 3 sampai 22

# DCC-GARCH
# pre
specid1 <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)), distribution.model = "sstd") #,mean.model = list(armaOrder = c(1,0), include.mean = TRUE))
specph1 <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)), distribution.model = "sstd")
specth1 <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)), distribution.model = "sstd")
specmy1 <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)), distribution.model = "sstd")
specsg1 <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)), distribution.model = "sstd")

speclist1 = list(specid1, specth1, specph1, specsg1, specmy1)
modelspec1 = dccspec(uspec=multispec(speclist1), dccOrder = c(1,1), distribution = "mvt")
modelfit1 = dccfit(modelspec1, data = pre[,-1])
modelfit1 # dapat alpha dan beta dimana alpha + beta < 1

# during
specid2 <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)), distribution.model = "sstd")
specph2 <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)), distribution.model = "sstd")
specth2 <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)), distribution.model = "sstd")
specmy2 <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)), distribution.model = "sstd")
specsg2 <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)), distribution.model = "sstd")

speclist2 = list(specid2, specth2, specph2, specsg2, specmy2)
modelspec2 = dccspec(uspec=multispec(speclist2), dccOrder = c(1,1), distribution = "mvt")
modelfit2 = dccfit(modelspec2, data = during[,-1])
modelfit2

# CCC-GARCH
# Pre
cccspec1 = cgarchspec(uspec=multispec(speclist1), dccOrder = c(1,1), distribution.model= list(copula="mvnorm"))
cccfit1 = cgarchfit(cccspec1, data = pre[,-1])
cccfit1
rcor(cccfit1)

# During
cccspec2 = cgarchspec(uspec=multispec(speclist2), dccOrder = c(1,1), distribution.model = list(copula="mvnorm"))
cccfit2 = cgarchfit(cccspec2, data = during[,-1])
cccfit2
rcor(cccfit2)

# Uji validitas model DCC-GARCH, pakai standardized residual atau squared standardized residual
stdres1 = as.data.frame(modelfit1@mfit$stdresid)
stdres2 = as.data.frame(modelfit2@mfit$stdresid)

plot(pre$Date, stdres1[,1],type = "l")
plot(pre$Date, stdres1[,2],type = "l")
plot(pre$Date, stdres1[,3],type = "l")
plot(pre$Date, stdres1[,4],type = "l")
plot(pre$Date, stdres1[,5],type = "l")

plot(during$Date, stdres2[,1],type = "l")
plot(during$Date, stdres2[,2],type = "l")
plot(during$Date, stdres2[,3],type = "l")
plot(during$Date, stdres2[,4],type = "l")
plot(during$Date, stdres2[,5],type = "l")

label1 = c("JKSE","SET","PSEI","STI","KLSE") 
colnames(stdres1) = label1
colnames(stdres2) = label1

ljungpreres1 <- data.frame(
  Lag = integer(),
  Column = character(),
  P_Value = numeric(),
  stringsAsFactors = FALSE
)

columns <- c("JKSE", "SET", "PSEI", "STI", "KLSE")

for (lag in 1:22) {
  for (col in columns) {
    p_value <- Box.test(stdres1[[col]], lag = lag, type = "Ljung-Box")$p.value
    
    result_row <- data.frame(
      Lag = lag,
      Column = col,
      P_Value = p_value
    )
    ljungpreres1 <- rbind(ljungpreres1, result_row)
  }
}

lbpre_res1 <- pivot_wider(ljungpreres1, names_from = Lag, values_from = P_Value) # Biar jadi table
lbpre_res1

ljungpreres2 <- data.frame(
  Lag = integer(),
  Column = character(),
  P_Value = numeric(),
  stringsAsFactors = FALSE
)

columns <- c("JKSE", "SET", "PSEI", "STI", "KLSE")

for (lag in 1:22) {
  for (col in columns) {
    p_value <- Box.test(stdres2[[col]], lag = lag, type = "Ljung-Box")$p.value
    
    result_row <- data.frame(
      Lag = lag,
      Column = col,
      P_Value = p_value
    )
    ljungpreres2 <- rbind(ljungpreres2, result_row)
  }
}

lbpre_res2 <- pivot_wider(ljungpreres2, names_from = Lag, values_from = P_Value) # Biar jadi table
lbpre_res2

lbpre_res1 # standardized residual dcc pre, semuanya tidak ada autokorelasi
lbpre_res2 # during, KLSE ada autokorelasi lag 12 sampe 13, 18 sampe 22

archres1 <- data.frame(
  Lag = integer(),
  Column = character(),
  P_Value = numeric(),
  stringsAsFactors = FALSE
)

columns <- c("JKSE", "SET", "PSEI", "STI", "KLSE")

for (lag in 1:22) {
  for (col in columns) {
    arch_test_result <- ArchTest(stdres1[[col]], lag = lag)
    p_value <- arch_test_result$p.value
    
    result_row <- data.frame(
      Lag = lag,
      Column = col,
      P_Value = p_value
    )
    archres1 <- rbind(archres1, result_row)
  }
}

archpre_res <- pivot_wider(archres1, names_from = Lag, values_from = P_Value)
archpre_res

archres2 <- data.frame(
  Lag = integer(),
  Column = character(),
  P_Value = numeric(),
  stringsAsFactors = FALSE
)

columns <- c("JKSE", "SET", "PSEI", "STI", "KLSE")

for (lag in 1:22) {
  for (col in columns) {
    arch_test_result <- ArchTest(stdres2[[col]], lag = lag)
    p_value <- arch_test_result$p.value
    
    result_row <- data.frame(
      Lag = lag,
      Column = col,
      P_Value = p_value
    )
    archres2 <- rbind(archres2, result_row)
  }
}

archdur_res <- pivot_wider(archres2, names_from = Lag, values_from = P_Value)
archdur_res

archpre_res # Standardized residual dcc pre, PSEI dan KLSE ada arch effect di lag 1
archdur_res # during, STI ada arch effect dari lag 15 sampe 20

# Squared standardized residual
sqres1 = stdres1^2
sqres2 = stdres2^2

ljungpresqres1 <- data.frame(
  Lag = integer(),
  Column = character(),
  P_Value = numeric(),
  stringsAsFactors = FALSE
)

columns <- c("JKSE", "SET", "PSEI", "STI", "KLSE")

for (lag in 1:22) {
  for (col in columns) {
    p_value <- Box.test(sqres1[[col]], lag = lag, type = "Ljung-Box")$p.value
    
    result_row <- data.frame(
      Lag = lag,
      Column = col,
      P_Value = p_value
    )
    ljungpresqres1 <- rbind(ljungpresqres1, result_row)
  }
}

lbpre_sqres1 <- pivot_wider(ljungpresqres1, names_from = Lag, values_from = P_Value) # Biar jadi table
lbpre_sqres1

ljungpresqres2 <- data.frame(
  Lag = integer(),
  Column = character(),
  P_Value = numeric(),
  stringsAsFactors = FALSE
)

columns <- c("JKSE", "SET", "PSEI", "STI", "KLSE")

for (lag in 1:22) {
  for (col in columns) {
    p_value <- Box.test(sqres2[[col]], lag = lag, type = "Ljung-Box")$p.value
    
    result_row <- data.frame(
      Lag = lag,
      Column = col,
      P_Value = p_value
    )
    ljungpresqres2 <- rbind(ljungpresqres2, result_row)
  }
}

lbpre_sqres2 <- pivot_wider(ljungpresqres2, names_from = Lag, values_from = P_Value) # Biar jadi table
lbpre_sqres2

lbpre_sqres1 # sq standardized residual dcc, PSEI dan KLSE ada autokorelasi di lag 1
lbpre_sqres2 # STI ada autokorelasi pada lag 15 sampe 20

archsqres1 <- data.frame(
  Lag = integer(),
  Column = character(),
  P_Value = numeric(),
  stringsAsFactors = FALSE
)

columns <- c("JKSE", "SET", "PSEI", "STI", "KLSE")

for (lag in 1:22) {
  for (col in columns) {
    arch_test_result <- ArchTest(sqres1[[col]], lag = lag)
    p_value <- arch_test_result$p.value
    
    result_row <- data.frame(
      Lag = lag,
      Column = col,
      P_Value = p_value
    )
    archsqres1 <- rbind(archsqres1, result_row)
  }
}

archpre_sqres <- pivot_wider(archsqres1, names_from = Lag, values_from = P_Value)
archpre_sqres

archsqres2 <- data.frame(
  Lag = integer(),
  Column = character(),
  P_Value = numeric(),
  stringsAsFactors = FALSE
)

columns <- c("JKSE", "SET", "PSEI", "STI", "KLSE")

for (lag in 1:22) {
  for (col in columns) {
    arch_test_result <- ArchTest(sqres2[[col]], lag = lag)
    p_value <- arch_test_result$p.value
    
    result_row <- data.frame(
      Lag = lag,
      Column = col,
      P_Value = p_value
    )
    archsqres2 <- rbind(archsqres2, result_row)
  }
}

archdur_sqres <- pivot_wider(archsqres2, names_from = Lag, values_from = P_Value)
archdur_sqres

archpre_sqres # squared standardized residual dcc, PSEI ada arch effect di lag 1 sampai 2 
archdur_sqres # PSEI (lag 5 sampe 22) dan STI (lag 15 sampe 22) ada arch effect

# Residual dari model CCC-GARCH
cccres1 = as.data.frame(cccfit1@mfit$stdresid)
cccres2 = as.data.frame(cccfit2@mfit$stdresid)

plot(pre$Date, cccres1[,1], type = "l")
plot(pre$Date, cccres1[,2], type = "l")
plot(pre$Date, cccres1[,3], type = "l")
plot(pre$Date, cccres1[,4], type = "l")
plot(pre$Date, cccres1[,5], type = "l")

plot(during$Date, cccres2[,1], type = "l")
plot(during$Date, cccres2[,2], type = "l")
plot(during$Date, cccres2[,3], type = "l")
plot(during$Date, cccres2[,4], type = "l")
plot(during$Date, cccres2[,5], type = "l")

label1 = c("JKSE","SET","PSEI","STI","KLSE")
colnames(cccres1) = label1
colnames(cccres2) = label1

ljungcccres1 <- data.frame(
  Lag = integer(),
  Column = character(),
  P_Value = numeric(),
  stringsAsFactors = FALSE
)

columns <- c("JKSE", "SET", "PSEI", "STI", "KLSE")

for (lag in 1:22) {
  for (col in columns) {
    p_value <- Box.test(cccres1[[col]], lag = lag, type = "Ljung-Box")$p.value
    
    result_row <- data.frame(
      Lag = lag,
      Column = col,
      P_Value = p_value
    )
    ljungcccres1 <- rbind(ljungcccres1, result_row)
  }
}

lb_cccres1 <- pivot_wider(ljungcccres1, names_from = Lag, values_from = P_Value) # Biar jadi table
lb_cccres1

ljungcccres2 <- data.frame(
  Lag = integer(),
  Column = character(),
  P_Value = numeric(),
  stringsAsFactors = FALSE
)

columns <- c("JKSE", "SET", "PSEI", "STI", "KLSE")

for (lag in 1:22) {
  for (col in columns) {
    p_value <- Box.test(cccres2[[col]], lag = lag, type = "Ljung-Box")$p.value
    
    result_row <- data.frame(
      Lag = lag,
      Column = col,
      P_Value = p_value
    )
    ljungcccres2 <- rbind(ljungcccres2, result_row)
  }
}

lb_cccres2 <- pivot_wider(ljungcccres2, names_from = Lag, values_from = P_Value) # Biar jadi table
lb_cccres2

lb_cccres1 # standardized residual ccc pre no autokorelasi
lb_cccres2 # standardized residual ccc during, cuma KLSE ada autokorelasi 

archcccres1 <- data.frame(
  Lag = integer(),
  Column = character(),
  P_Value = numeric(),
  stringsAsFactors = FALSE
)

columns <- c("JKSE", "SET", "PSEI", "STI", "KLSE")

for (lag in 1:22) {
  for (col in columns) {
    arch_test_result <- ArchTest(cccres1[[col]], lag = lag)
    p_value <- arch_test_result$p.value
    
    result_row <- data.frame(
      Lag = lag,
      Column = col,
      P_Value = p_value
    )
    archcccres1 <- rbind(archcccres1, result_row)
  }
}

archpre_cccres <- pivot_wider(archcccres1, names_from = Lag, values_from = P_Value)
archpre_cccres

archcccres2 <- data.frame(
  Lag = integer(),
  Column = character(),
  P_Value = numeric(),
  stringsAsFactors = FALSE
)

columns <- c("JKSE", "SET", "PSEI", "STI", "KLSE")

for (lag in 1:22) {
  for (col in columns) {
    arch_test_result <- ArchTest(cccres2[[col]], lag = lag)
    p_value <- arch_test_result$p.value
    
    result_row <- data.frame(
      Lag = lag,
      Column = col,
      P_Value = p_value
    )
    archcccres2 <- rbind(archcccres2, result_row)
  }
}

archdur_cccres <- pivot_wider(archcccres2, names_from = Lag, values_from = P_Value)
archdur_cccres

archpre_cccres # standardized residual ccc pre, PSEI dan KLSE ada arch effect pada lag 1
archdur_cccres # during, STI ada arch effect mulai dari lag 15 sampe 19

# Squared standardized residual
sqcccres1 = cccres1^2
sqcccres2 = cccres2^2

ljungsqcccres1 <- data.frame(
  Lag = integer(),
  Column = character(),
  P_Value = numeric(),
  stringsAsFactors = FALSE
)

columns <- c("JKSE", "SET", "PSEI", "STI", "KLSE")

for (lag in 1:22) {
  for (col in columns) {
    p_value <- Box.test(sqcccres1[[col]], lag = lag, type = "Ljung-Box")$p.value
    
    result_row <- data.frame(
      Lag = lag,
      Column = col,
      P_Value = p_value
    )
    ljungsqcccres1 <- rbind(ljungsqcccres1, result_row)
  }
}

lb_sqcccres1 <- pivot_wider(ljungsqcccres1, names_from = Lag, values_from = P_Value) # Biar jadi table
lb_sqcccres1

ljungsqcccres2 <- data.frame(
  Lag = integer(),
  Column = character(),
  P_Value = numeric(),
  stringsAsFactors = FALSE
)

columns <- c("JKSE", "SET", "PSEI", "STI", "KLSE")

for (lag in 1:22) {
  for (col in columns) {
    p_value <- Box.test(sqcccres2[[col]], lag = lag, type = "Ljung-Box")$p.value
    
    result_row <- data.frame(
      Lag = lag,
      Column = col,
      P_Value = p_value
    )
    ljungsqcccres2 <- rbind(ljungsqcccres2, result_row)
  }
}

lb_sqcccres2 <- pivot_wider(ljungsqcccres2, names_from = Lag, values_from = P_Value) # Biar jadi table
lb_sqcccres2

lb_sqcccres1 # standardized sq residual ccc pre, ada autokorelasi untuk PSEI (lag 1) dan KLSE (lag 1 dan 2)
lb_sqcccres2 # during, ada autokorelasi untuk STI dari lag 15 sampe 20

archsqcccres1 <- data.frame(
  Lag = integer(),
  Column = character(),
  P_Value = numeric(),
  stringsAsFactors = FALSE
)

columns <- c("JKSE", "SET", "PSEI", "STI", "KLSE")

for (lag in 1:22) {
  for (col in columns) {
    arch_test_result <- ArchTest(sqcccres1[[col]], lag = lag)
    p_value <- arch_test_result$p.value
    
    result_row <- data.frame(
      Lag = lag,
      Column = col,
      P_Value = p_value
    )
    archsqcccres1 <- rbind(archsqcccres1, result_row)
  }
}

archpre_sqcccres <- pivot_wider(archsqcccres1, names_from = Lag, values_from = P_Value)
archpre_sqcccres

archsqcccres2 <- data.frame(
  Lag = integer(),
  Column = character(),
  P_Value = numeric(),
  stringsAsFactors = FALSE
)

columns <- c("JKSE", "SET", "PSEI", "STI", "KLSE")

for (lag in 1:22) {
  for (col in columns) {
    arch_test_result <- ArchTest(sqcccres2[[col]], lag = lag)
    p_value <- arch_test_result$p.value
    
    result_row <- data.frame(
      Lag = lag,
      Column = col,
      P_Value = p_value
    )
    archsqcccres2 <- rbind(archsqcccres2, result_row)
  }
}

archdur_sqcccres <- pivot_wider(archsqcccres2, names_from = Lag, values_from = P_Value)
archdur_sqcccres

archpre_sqcccres # sq standardized residual ccc pre, PSEI ada ARCH effect di  lag 1 sampe 2 dan lag 6
archdur_sqcccres # during, PSEI (lag 5 sampe 22)dan STI (lag 15 sampe 22) ada ARCH effect


# Estimasi matriks korelasi dinamis
correlation1 = rcor(modelfit1) # Dari data pre
correlation2 = rcor(modelfit2) # Dari data during

# Persiapan data dan tanggal covid untuk plotting
date = clean$Date
date = date[1:length(date)]
covid = as.Date("2020-01-30")

par(mfrow=c(1,1))
par(mfrow=c(3,3))

layout(matrix(1:10, nrow = 5, ncol = 2))
par(mar = c(2,2,1,1))

# Korelasi dinamis indonesia-thailand
cor_idth1 = correlation1[2,1,]
cor_idth2 = correlation2[2,1,]
idth = c(cor_idth1, cor_idth2)

plot_idth = data.frame(Date = date, DC = idth)
plot(plot_idth$Date, plot_idth$DC, type = "l", xlab = "", ylab = "", ylim = c(0.1,0.6))
abline(v = covid, col = "red")
title("JKSE-SET")

wilcox.test(cor_idth1, cor_idth2, alternative = "two.sided")# Beda
wilcox.test(cor_idth1, cor_idth2, alternative = "two.sided")$p.value

# Korelasi dinamis indonesia-filipina
cor_idph1 = correlation1[3,1,]
cor_idph2 = correlation2[3,1,]
idph = c(cor_idph1, cor_idph2)

plot_idph = data.frame(Date = date, DC = idph)
plot(plot_idph$Date, plot_idph$DC, type = "l", xlab = "", ylab = "", ylim = c(0.1,0.6))
abline(v = covid, col = "red")
title("JKSE-PSEI")

wilcox.test(cor_idph1, cor_idph2, alternative = "two.sided") # Beda
wilcox.test(cor_idph1, cor_idph2, alternative = "two.sided")$p.value

# Korelasi dinamis indonesia-singapura
cor_idsg1 = correlation1[4,1,]
cor_idsg2 = correlation2[4,1,]
idsg = c(cor_idsg1, cor_idsg2)

plot_idsg = data.frame(Date = date, DC = idsg)
plot(plot_idsg$Date, plot_idsg$DC, type = "l", xlab = "", ylab = "", ylim = c(0.1,0.6))
abline(v = covid, col = "red")
title("JKSE-STI")

wilcox.test(cor_idsg1, cor_idsg2, alternative = "two.sided") # Beda
wilcox.test(cor_idsg1, cor_idsg2, alternative = "two.sided")$p.value

# Korelasi dinamis indonesia-malaysia
cor_idmy1 = correlation1[5,1,]
cor_idmy2 = correlation2[5,1,]
idmy = c(cor_idmy1, cor_idmy2)

plot_idmy = data.frame(Date = date, DC = idmy)
plot(plot_idmy$Date, plot_idmy$DC, type = "l", xlab = "", ylab = "", ylim = c(0.1,0.6))
abline(v = covid, col = "red")
title("JKSE-KLSE")

wilcox.test(cor_idmy1, cor_idmy2, alternative = "two.sided") # Beda
wilcox.test(cor_idmy1, cor_idmy2, alternative = "two.sided")$p.value

# Korelasi dinamis thailand-filipina
cor_thph1 = correlation1[3,2,]
cor_thph2 = correlation2[3,2,]
thph = c(cor_thph1, cor_thph2)

plot_thph = data.frame(Date = date, DC = thph)
plot(plot_thph$Date, plot_thph$DC, type = "l", xlab = "", ylab = "", ylim = c(0.1,0.6))
abline(v = covid, col = "red")
title("SET-PSEI")

wilcox.test(cor_thph1, cor_thph2, alternative = "two.sided") # Beda
wilcox.test(cor_thph1, cor_thph2, alternative = "two.sided")$p.value

# Korelasi dinamis thailand-singapura
cor_thsg1 = correlation1[4,2,]
cor_thsg2 = correlation2[4,2,]
thsg = c(cor_thsg1, cor_thsg2)

plot_thsg = data.frame(Date = date, DC = thsg)
plot(plot_thsg$Date, plot_thsg$DC, type = "l", xlab = "", ylab = "", ylim = c(0.1,0.6))
abline(v = covid, col = "red")
title("SET-STI")

wilcox.test(cor_thsg1, cor_thsg2, alternative = "two.sided") # Beda
wilcox.test(cor_thsg1, cor_thsg2, alternative = "two.sided")$p.value

# Korelasi dinamis thailand-malaysia
cor_thmy1 = correlation1[5,2,]
cor_thmy2 = correlation2[5,2,]
thmy = c(cor_thmy1, cor_thmy2)

plot_thmy = data.frame(Date = date, DC = thmy)
plot(plot_thmy$Date, plot_thmy$DC, type = "l", xlab = "", ylab = "", ylim = c(0.1,0.6))
abline(v = covid, col = "red")
title("SET-KLSE")

wilcox.test(cor_thmy1, cor_thmy2, alternative = "two.sided") # Beda
wilcox.test(cor_thmy1, cor_thmy2, alternative = "two.sided")$p.value # Beda

# Korelasi dinamis filipina-singapura
cor_phsg1 = correlation1[4,3,]
cor_phsg2 = correlation2[4,3,]
phsg = c(cor_phsg1, cor_phsg2)

plot_phsg = data.frame(Date = date, DC = phsg)
plot(plot_phsg$Date, plot_phsg$DC, type = "l", xlab = "", ylab = "", ylim = c(0.1,0.6))
abline(v = covid, col = "red")
title("PSEI-STI")

wilcox.test(cor_phsg1, cor_phsg2, alternative = "two.sided") # Beda
wilcox.test(cor_phsg1, cor_phsg2, alternative = "two.sided")$p.value

# Korelasi dinamis filipina-malaysia
cor_phmy1 = correlation1[5,3,]
cor_phmy2 = correlation2[5,3,]
phmy = c(cor_phmy1, cor_phmy2)

plot_phmy = data.frame(Date = date, DC = phmy)
plot(plot_phmy$Date, plot_phmy$DC, type = "l", xlab = "", ylab = "", ylim = c(0.1,0.6))
abline(v = covid, col = "red")
title("PSEI-KLSE")

wilcox.test(cor_phmy1, cor_phmy2, alternative = "two.sided") # Beda
wilcox.test(cor_phmy1, cor_phmy2, alternative = "two.sided")$p.value

# Korelasi dinamis singapura-malaysia
cor_sgmy1 = correlation1[5,4,]
cor_sgmy2 = correlation2[5,4,]
sgmy = c(cor_sgmy1, cor_sgmy2)

plot_sgmy = data.frame(Date = date, DC = sgmy)
plot(plot_sgmy$Date, plot_sgmy$DC, type = "l", xlab = "", ylab = "", ylim = c(0.1,0.6))
abline(v = covid, col = "red")
title("STI-KLSE")

wilcox.test(cor_sgmy1, cor_sgmy2, alternative = "two.sided") # Beda
wilcox.test(cor_sgmy1, cor_sgmy2, alternative = "two.sided")$p.value

# Perbandingan mean korelasi dinamis (belum dilakukan per pasangan)
# wilcox.test(cor_idph2, mu = mean(cor_idph1), alternative = "two.sided") # wilcoxon sign rank test
# wilcox.test(cor_idph2, cor_idph3) # mann whitney u test


######################################################################################
# Pendekatan Diebold Yilmaz
# Pemilihan model VAR terbaik untuk setiap periode

VARselect(pre[,-1])     # p = 1
VARselect(during[,-1])  # AIC = 10, HQ = 2, BIC = 1, FPEP = 10, p = 10

varpre = vars::VAR(pre[,-1], p = 1)
print(varpre)
coefficients(varpre)
summary(varpre)

vardur10 = vars::VAR(during[,-1], p = 10)
vardur = vardur10 # Jadi untuk data set during gunakan VAR(10)
coefficients(vardur)
print(vardur)
summary(vardur)

options(scipen = 2)

# Convert pre, during, post jadi zoo
prezoo = zoo(pre[,-1], order.by=pre$Date)
duringzoo = zoo(during[,-1], order.by=during$Date)

# Diebold Yilmaz Approach
dca_pre = ConnectednessApproach(prezoo,
                            nlag = 1, # lag
                            nfore = 10, #n.ahead
                            window.size = 265, # belum tau mau pilih berapa,
                            model = "VAR",
                            connectedness = "Time",
                            Connectedness_config=list(TimeConnectedness=list(generalized=TRUE))) # Generalized sesuai DY 2012
dca_pre$TABLE
PlotNPDC(dca_pre, ylim = c(-10,10))


dca_during = ConnectednessApproach(duringzoo,
                                nlag = 10, # lag
                                nfore = 10, #n.ahead
                                window.size = 265, # biar plot nya ga cuma 1 observasi
                                model = "VAR",
                                connectedness = "Time",
                                Connectedness_config=list(TimeConnectedness=list(generalized=TRUE))) # Generalized sesuai DY 2012
dca_during$TABLE
PlotNPDC(dca_during, ylim = c(-10,10))

# Tabel Diebold-Yilmaz static
static_pre = ConnectednessApproach(prezoo,
                                nlag = 1, # lag
                                nfore = 10, #n.ahead
                                model = "VAR",
                                connectedness = "Time",
                                Connectedness_config=list(TimeConnectedness=list(generalized=TRUE))) # Generalized sesuai DY 2012
static_pre$TABLE

static_during = ConnectednessApproach(duringzoo,
                                   nlag = 10, # lag
                                   nfore = 10, #n.ahead
                                   model = "VAR",
                                   connectedness = "Time",
                                   Connectedness_config=list(TimeConnectedness=list(generalized=TRUE))) # Generalized sesuai DY 2012
static_during$TABLE

# Plotting NPDC gabungan 2 periode
dca_whole = abind(dca_pre$NPDC, dca_during$NPDC, along = 3)
tail(dca_whole)

# JKSE-SET
plot(dca_whole[2,1,], type = "l", ylim = c(-10,10))
abline(h=0)
polygon(c(dca_whole[2,1,], rev(dca_whole[2,1,])), c(dca_whole[2,1,], rep(0, length(dca_whole[2,1,]))), col = "black", border = NA)
rect(par("usr")[1], par("usr")[3], par("usr")[2], 0, col = "black", border = NA, density = 1000)


# Diagnosis model VAR
varpre = vars::VAR(pre[,-1], p = 1)
vardur = vars::VAR(during[,-1], p = 10)

# Testing tabel baru
varpre2 = ConnectednessApproach::VAR(prezoo, configuration = list(nlag = 1))
vardur2 = ConnectednessApproach::VAR(duringzoo, configuration = list(nlag = 10))

fevdpre = FEVD(Phi = varpre2$B, Sigma = varpre2$Q, nfore = 10, type = "time", generalized = TRUE) # Same thing
fevdduring = FEVD(Phi = vardur2$B, Sigma = vardur2$Q, nfore = 10, type = "time", generalized = TRUE) # Same thing
fevdpre$FEVD # Same thing
fevdduring$FEVD # Same thing


# Pengujian residual
resvarpre = data.frame(residuals(varpre))
resvardur = data.frame(residuals(vardur))

# Uji stasioneritas residual
res_frames <- list(Res_pre = resvarpre, Res_dur = resvardur)
columns <- c("JKSE", "SET", "PSEI", "STI", "KLSE")
adf_res <- data.frame(Column = columns)

for (df_name in names(res_frames)) {
  df <- res_frames[[df_name]]
  p_values <- numeric(length(columns))
  
  for (i in 1:length(columns)) {
    result <- adf.test(df[[columns[i]]])
    p_values[i] <- result$p.value
  }
  
  adf_res[[df_name]] <- p_values
}
adf_res # Semua residu stasioner

# Uji normalitas residual, H0 data berdistribusi normal

jbtest_res <- data.frame(Column = columns)

for (df_name in names(res_frames)) {
  df <- res_frames[[df_name]]
  p_values <- numeric(length(columns))
  
  for (i in 1:length(columns)) {
    result <- jarque.bera.test(df[[columns[i]]])
    p_values[i] <- result$p.value
  }
  
  jbtest_res[[df_name]] <- p_values
}
jbtest_res # Tidak normal semua

# Uji autokorelasi residual
# Residual Var PRE
ljungvarpre <- data.frame(
  Lag = integer(),
  Column = character(),
  P_Value = numeric(),
  stringsAsFactors = FALSE
)

for (lag in 1:22) {
  for (col in columns) {
    p_value <- Box.test(resvarpre[[col]], lag = lag, type = "Ljung-Box")$p.value
    
    result_row <- data.frame(
      Lag = lag,
      Column = col,
      P_Value = p_value
    )
    ljungvarpre <- rbind(ljungvarpre, result_row)
  }
}

lbvarpre <- pivot_wider(ljungvarpre, names_from = Lag, values_from = P_Value) # Biar jadi table

# Residual Var DURING
ljungvardur <- data.frame(
  Lag = integer(),
  Column = character(),
  P_Value = numeric(),
  stringsAsFactors = FALSE
)

for (lag in 1:22) {
  for (col in columns) {
    p_value <- Box.test(resvardur[[col]], lag = lag, type = "Ljung-Box")$p.value
    
    result_row <- data.frame(
      Lag = lag,
      Column = col,
      P_Value = p_value
    )
    ljungvardur <- rbind(ljungvardur, result_row)
  }
}

lbvardur <- pivot_wider(ljungvardur, names_from = Lag, values_from = P_Value) # Biar jadi table

lbvarpre # Hanya residual JKSE menunjukkan autokorelasi
lbvardur # Semua residual tidak menunjukkan autokorelasi

# Uji ARCH effect residual

archvar1 <- data.frame(
  Lag = integer(),
  Column = character(),
  P_Value = numeric(),
  stringsAsFactors = FALSE
)

columns <- c("JKSE", "SET", "PSEI", "STI", "KLSE")

for (lag in 1:22) {
  for (col in columns) {
    arch_test_result <- ArchTest(resvarpre[[col]], lag = lag)
    p_value <- arch_test_result$p.value
    
    result_row <- data.frame(
      Lag = lag,
      Column = col,
      P_Value = p_value
    )
    archvar1 <- rbind(archvar1, result_row)
  }
}

archvarpre <- pivot_wider(archvar1, names_from = Lag, values_from = P_Value)

# ARCH LM test untuk Residual VAR During
archvar2 <- data.frame(
  Lag = integer(),
  Column = character(),
  P_Value = numeric(),
  stringsAsFactors = FALSE
)

for (lag in 1:22) {
  for (col in columns) {
    arch_test_result <- ArchTest(resvardur[[col]], lag = lag)
    p_value <- arch_test_result$p.value
    
    result_row <- data.frame(
      Lag = lag,
      Column = col,
      P_Value = p_value
    )
    archvar2 <- rbind(archvar2, result_row)
  }
}

archvardur <- pivot_wider(archvar2, names_from = Lag, values_from = P_Value)


archvarpre # Semua ada ARCH effect
archvardur # Semua ada ARCH effect

######################################################################

# Uji perbedaan spillover antar periode pre dan during

# JKSE - SET
JKSE_SET_1 = dca_pre$NPDC[2,1,]
JKSE_SET_2 = dca_during$NPDC[2,1,]
wilcox.test(JKSE_SET_1, JKSE_SET_2)
wilcox.test(JKSE_SET_1, JKSE_SET_2)$p.value # Beda

# JKSE - PSEI
JKSE_PSEI_1 = dca_pre$NPDC[3,1,]
JKSE_PSEI_2 = dca_during$NPDC[3,1,]
wilcox.test(JKSE_PSEI_1, JKSE_PSEI_2)
wilcox.test(JKSE_PSEI_1, JKSE_PSEI_2)$p.value # Beda

# JKSE - STI
JKSE_STI_1 = dca_pre$NPDC[4,1,]
JKSE_STI_2 = dca_during$NPDC[4,1,]
wilcox.test(JKSE_STI_1, JKSE_STI_2)
wilcox.test(JKSE_STI_1, JKSE_STI_2)$p.value # Beda

# JKSE - KLSE
JKSE_KLSE_1 = dca_pre$NPDC[5,1,]
JKSE_KLSE_2 = dca_during$NPDC[5,1,]
wilcox.test(JKSE_KLSE_1, JKSE_KLSE_2)
wilcox.test(JKSE_KLSE_1, JKSE_KLSE_2)$p.value # Beda

# SET - PSEI
SET_PSEI_1 = dca_pre$NPDC[3,2,]
SET_PSEI_2 = dca_during$NPDC[3,2,]
wilcox.test(SET_PSEI_1, SET_PSEI_2)
wilcox.test(SET_PSEI_1, SET_PSEI_2)$p.value # Beda

# SET - STI
SET_STI_1 = dca_pre$NPDC[4,2,]
SET_STI_2 = dca_during$NPDC[4,2,]
wilcox.test(SET_STI_1, SET_STI_2)
wilcox.test(SET_STI_1, SET_STI_2)$p.value # Beda

# SET - KLSE
SET_KLSE_1 = dca_pre$NPDC[5,2,]
SET_KLSE_2 = dca_during$NPDC[5,2,]
wilcox.test(SET_KLSE_1, SET_KLSE_2)
wilcox.test(SET_KLSE_1, SET_KLSE_2, alternative = "two.sided")$p.value # Beda
wilcox.test(SET_KLSE_1, SET_KLSE_2)$p.value # Beda

# PSEI - STI
PSEI_STI_1 = dca_pre$NPDC[4,3,]
PSEI_STI_2 = dca_during$NPDC[4,3,]
wilcox.test(PSEI_STI_1, PSEI_STI_2)
wilcox.test(PSEI_STI_1, PSEI_STI_2)$p.value # Beda

# PSEI - KLSE
PSEI_KLSE_1 = dca_pre$NPDC[5,3,]
PSEI_KLSE_2 = dca_during$NPDC[5,3,]
wilcox.test(PSEI_KLSE_1, PSEI_KLSE_2)
wilcox.test(PSEI_KLSE_1, PSEI_KLSE_2)$p.value # Beda

# STI - KLSE
STI_KLSE_1 = dca_pre$NPDC[5,4,]
STI_KLSE_2 = dca_during$NPDC[5,4,]
wilcox.test(STI_KLSE_1, STI_KLSE_2)
wilcox.test(STI_KLSE_1, STI_KLSE_2)$p.value # Beda


##################################################################

# Ngetest ngitung CCC sendiri
# MUNGKIN INI CARA PAKE CCC wiley
resid1 = residuals(ufitid1)
resph1 = residuals(ufitph1)
resth1 = residuals(ufitth1)
resmy1 = residuals(ufitmy1)
ressg1 = residuals(ufitsg1)
res1 = data.frame(JKSE = resid1, SET = resth1, PSEI = resph1, STI = ressg1, KLSE = resmy1)
cor(res1)
correlation1[,,1]

ufitid2 = ugarchfit(specid2, data = during$JKSE)
ufitph2 = ugarchfit(specph2, data = during$PSEI)
ufitth2 = ugarchfit(specth2, data = during$SET)
ufitmy2 = ugarchfit(specmy2, data = during$KLSE)
ufitsg2 = ugarchfit(specsg2, data = during$STI)

ufitid1 = ugarchfit(specid1, data = pre$JKSE)
ufitph1 = ugarchfit(specph1, data = pre$PSEI)
ufitth1 = ugarchfit(specth1, data = pre$SET)
ufitmy1 = ugarchfit(specmy1, data = pre$KLSE)
ufitsg1 = ugarchfit(specsg1, data = pre$STI)

resid2 = residuals(ufitid2)
resph2 = residuals(ufitph2)
resth2 = residuals(ufitth2)
resmy2 = residuals(ufitmy2)
ressg2 = residuals(ufitsg2)
res2 = data.frame(JKSE = resid2, SET = resth2, PSEI = resph2, STI = ressg2, KLSE = resmy2)
cor(res2)
correlation2[,,1]

# Mungkin ini CCC pakai Copula GARCH (Normal) berdasarkan package rmgarch
cccspec1 = cgarchspec(uspec=multispec(speclist1), dccOrder = c(1,1), distribution.model= list(copula="mvnorm"))
cccfit1 = cgarchfit(cccspec1, data = pre[,-1])
cccfit1
modelfit1
correlation1 # DCC Matrix
rcor(cccfit1) # CCC Matrix pakai package

specid2 <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)), distribution.model = "sstd")
specph2 <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)), distribution.model = "sstd")
specth2 <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)), distribution.model = "sstd")
specmy2 <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)), distribution.model = "sstd")
specsg2 <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)), distribution.model = "sstd")

speclist2 = list(specid2, specth2, specph2, specsg2, specmy2)
modelspec2 = dccspec(uspec=multispec(speclist2), dccOrder = c(1,1), distribution = "mvt")
modelfit2 = dccfit(modelspec2, data = during[,-1])
modelfit2

cccspec2 = cgarchspec(uspec=multispec(speclist2), dccOrder = c(1,1), distribution.model = list(copula="mvnorm"))
cccfit2 = cgarchfit(cccspec2, data = during[,-1])
cccfit2
modelfit2
rcor(cccfit2) # beda banget korelasi nya
correlation2

# Kalau bandingin nilai kriteria, menang CCC

# Testing

dca_during = ConnectednessApproach(duringzoo,
                                   nlag = 10, # lag
                                   nfore = 10, #n.ahead
                                   window.size = 300, # biar plot nya ga cuma 1 observasi
                                   model = "VAR",
                                   connectedness = "Time",
                                   Connectedness_config=list(TimeConnectedness=list(generalized=TRUE))) # Generalized sesuai DY 2012
dca_during$TABLE
PlotNPDC(dca_during)

testdy = VAR(duringzoo, configuration = list(nlag=10))
testdcady = TimeConnectedness(testdy$B, testdy$Q, nfore = 10)
testdcady$TABLE # Beda karena masalah window size

data(dy2012)

fit = VAR(dy2012, configuration=list(nlag=4))
dca1= TimeConnectedness(fit$B, fit$Q, nfore=10)
dca1$TABLE

dca2 = ConnectednessApproach(dy2012,
                            nlag=4,
                            nfore=10,
                            model="VAR",
                            corrected=FALSE,
                            window.size=NULL,
                            connectedness="Time",
                            Connectedness_config = list(
                              TimeConnectedness=list(generalized=TRUE)))
dca2$TABLE
plotNPDC(dca2)

dca3 = ConnectednessApproach(dy2012,
                             nlag=4,
                             nfore=10,
                             model="VAR",
                             corrected=FALSE,
                             window.size=100,
                             connectedness="Time",
                             Connectedness_config = list(
                               TimeConnectedness=list(generalized=TRUE)))
dca3$TABLE

count(dy2012)
###########################################################

# Testing CCC-GARCH pakai package ccgarch
# Sepertinya tidak bisa
ufitid2 = ugarchfit(specid2, data = during$JKSE)
ufitph2 = ugarchfit(specph2, data = during$PSEI)
ufitth2 = ugarchfit(specth2, data = during$SET)
ufitmy2 = ugarchfit(specmy2, data = during$KLSE)
ufitsg2 = ugarchfit(specsg2, data = during$STI)

ufitid1 = ugarchfit(specid1, data = pre$JKSE)
ufitph1 = ugarchfit(specph1, data = pre$PSEI)
ufitth1 = ugarchfit(specth1, data = pre$SET)
ufitmy1 = ugarchfit(specmy1, data = pre$KLSE)
ufitsg1 = ugarchfit(specsg1, data = pre$STI)

ufitid1 # constant = 0,000362
ufitth1 # c = 0,000210
ufitph1 # c = 0,000193
ufitsg1 # c = 0,000289
ufitmy1 # c = -0,000066

a = c(0.000362, 0.000210, 0.000193, 0.000289, -0.000066)
a
eccc.estimation()
