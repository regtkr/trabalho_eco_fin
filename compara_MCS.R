#----------------------------------------------------------------------
# BIBLIOTECAS
#----------------------------------------------------------------------
library(fImport)
library(zoo)
library(rugarch)
library(GAS)
library(MCS)
library(timeSeries)
library(parallel)
library(stochvol)



#----------------------------------------------------------------------
# DADOS
#----------------------------------------------------------------------
dados = read.table("/home/regis/Dropbox/econometriafinanças2017/parte3-volatilidade/sp500.txt",header=F)
ret = diff(log(dados$V1))

# selecionando as primeiras linhas para ser mais rápido
faixa = 1:2000

ret = ret[faixa]

tret = ts(ret, start = c(1980,1,1), frequency = 365) # Data ficticia, melhorar


#----------------------------------------------------------------------
# RECEITA
#----------------------------------------------------------------------
# 1. Fit an ARMA+GARCH model to your return series. You can do this by 
#    ugarchfit function in rugarch package in R (you will have to choose 
#    what density you want to assume for the standardized residuals).
# 2. Forecast conditional mean and conditional variance of the assumed 
#    parametric density 10 days ahead. You can do this by ugarchroll 
#    function in rugarch package in R.
# 3. Given the forecasted mean and variance of the assumed density, you 
#    can obtain the 0.05 quantile of the distribution which will be your 
#    5% VaR (you can use other quantiles, of course). You can do this by 
#    function qnorm for the normal density, qt for student density, or 
#    generally qsomethingelse for some other density in R.


# _____________________________________________________________________
#
#                                GARCH
# _____________________________________________________________________

#----------------------------------------------------------------------
# MODELOS
#----------------------------------------------------------------------
spec = list()

spec[1]  = ugarchspec(variance.model = list(model = "sGARCH"))
spec[2]  = ugarchspec(variance.model = list(model = "eGARCH"))
# spec[3]  = ugarchspec(variance.model = list(model = "girGARCH"))
# spec[3]  = ugarchspec(variance.model = list(model = "appGARCH"))
spec[3]  = ugarchspec(variance.model = list(model = "iGARCH"))
spec[4]  = ugarchspec(variance.model = list(model = "csGARCH"))
spec[5]  = ugarchspec(variance.model = list(model = "fGARCH", submodel = "GARCH"))
spec[6]  = ugarchspec(variance.model = list(model = "fGARCH", submodel = "TGARCH"))
spec[7]  = ugarchspec(variance.model = list(model = "fGARCH", submodel = "AVGARCH"))
spec[8]  = ugarchspec(variance.model = list(model = "fGARCH", submodel = "GARCH"))
spec[9]  = ugarchspec(variance.model = list(model = "fGARCH", submodel = "NAGARCH"))
spec[10] = ugarchspec(variance.model = list(model = "fGARCH", submodel = "APARCH"))
spec[11] = ugarchspec(variance.model = list(model = "fGARCH", submodel = "GJRGARCH"))

#----------------------------------------------------------------------
# VaR
#----------------------------------------------------------------------
loss = data.frame(matrix(, nrow=500, ncol=0))

for (i in 1:length(spec)){
	# Nome dos modelo
	model_name = paste("GARCH_", as.character(i), sep = '') 
	print(model_name)

	# Aplicando a janela móvel
	roll = ugarchroll(spec[[i]], data = tret, n.ahead = 1, n.start = 1500, VaR.alpha = c(0.05), refit.every = 50)

	# Extraindo o VaR
	VaR = as.data.frame(roll, which="VaR")
	VaR_sim  = VaR[,1]
	VaR_real = VaR[,2]

	# Calculando a função perda
	loss[model_name] = LossVaR(VaR_real, VaR_sim, tau = 0.05)
}


# _____________________________________________________________________
#
#                                   GAS
# _____________________________________________________________________

#----------------------------------------------------------------------
# MODELOS
#----------------------------------------------------------------------
spec[12] = UniGASSpec(Dist = 'norm')
spec[13] = UniGASSpec(Dist = 'snorm')
spec[14] = UniGASSpec(Dist = 'std')
spec[15] = UniGASSpec(Dist = 'sstd')

#----------------------------------------------------------------------
# VaR
#----------------------------------------------------------------------
# TODO: Unificar loop com o do GARCH
for (i in 12:15) {
	# Nome dos modelo
	model_name = paste("GAS_", as.character(i-11), sep = '') 
		print(model_name)

	# Aplicando a janela móvel
	roll = UniGASRoll(data = tret, 
		GASSpec = spec[[i]], 
		Nstart = 1500, 
		RefitEvery = 50, 
		RefitWindow = c("moving"))

	# Extraindo o VaR
	VaR_sim = quantile(roll, 0.05)

	# Gráficos 
	# plot(ret[-(1:5000)], type='l', col='red')
	# lines(VaR)

	# Calculando a função perda
	loss[model_name] = LossVaR(VaR_real, VaR_sim, tau = 0.05)
}

# _____________________________________________________________________
#
#                          STOCHASTIC VOLATILITY
# _____________________________________________________________________

#----------------------------------------------------------------------
# MODELOS
#----------------------------------------------------------------------

m = mean(tret)

#----------------------------------------------------------------------
# VaR
#----------------------------------------------------------------------

VaR_sim = vector(mode = "numeric", length = 500)

for (i in 1:length(VaR_sim)) {
	print(i)

	draws = svsample(
		tret[seq(i,1500+i-1)] - m,
		draws = 1000, 
		burnin = 100, 
		priormu = c(-10, 1), 
		priorphi = c(20, 1.2), 
		priorsigma = 0.2,
		designmatrix = "ar1", 
		quiet = TRUE)

	prev = predict(draws, step = 1)

	previsao = arpredict(draws, prev)

	VaR_sim[i] = quantile(previsao, 0.05) + m
}

plot(tret, type = 'l', col = 'red')
lines(VaR_sim)

loss['SV'] = LossVaR(VaR_real, VaR_sim, tau = 0.05)


#----------------------------------------------------------------------
# COMPARANDO
#----------------------------------------------------------------------
comparacao = MCSprocedure(Loss=loss,alpha=0.2,B=5000,statistic='Tmax',cl=NULL)

