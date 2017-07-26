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
library(stargazer)


#----------------------------------------------------------------------
# DADOS
#----------------------------------------------------------------------
dados = read.table(
	"/home/regis/Dropbox/econometriafinanças2017/parte3-volatilidade/sp500.txt",
	header=F)
ret = diff(log(dados$V1))

# selecionando as primeiras linhas para ser mais rápido
faixa = 1:2000

ret = ret[faixa]

tret = ts(ret, start = c(1980,1,1), frequency = 365) # Data ficticia, melhorar


#----------------------------------------------------------------------
# Parametros iniciais
#----------------------------------------------------------------------
inicio = 1500
fim    = length(tret)
delta  = fim - inicio


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

# NORMAIS
spec[1]  = ugarchspec(variance.model = list(model = "sGARCH"), distribution.model = 'norm')
spec[2]  = ugarchspec(variance.model = list(model = "eGARCH"), distribution.model = 'norm')
spec[3]  = ugarchspec(variance.model = list(model = "iGARCH"), distribution.model = 'norm')
spec[4]  = ugarchspec(variance.model = list(model = "csGARCH"), distribution.model = 'norm')
spec[5]  = ugarchspec(variance.model = list(model = "fGARCH", submodel = "GARCH"), distribution.model = 'norm')
spec[6]  = ugarchspec(variance.model = list(model = "fGARCH", submodel = "TGARCH"), distribution.model = 'norm')
spec[7]  = ugarchspec(variance.model = list(model = "fGARCH", submodel = "AVGARCH"), distribution.model = 'norm')
spec[8]  = ugarchspec(variance.model = list(model = "fGARCH", submodel = "NGARCH"), distribution.model = 'norm')
spec[9]  = ugarchspec(variance.model = list(model = "fGARCH", submodel = "NAGARCH"), distribution.model = 'norm')
spec[10] = ugarchspec(variance.model = list(model = "fGARCH", submodel = "APARCH"), distribution.model = 'norm')
spec[11] = ugarchspec(variance.model = list(model = "fGARCH", submodel = "GJRGARCH"), distribution.model = 'norm')

# T-STUDENT
spec[12]  = ugarchspec(variance.model = list(model = "sGARCH"), distribution.model = "std")
spec[13]  = ugarchspec(variance.model = list(model = "eGARCH"), distribution.model = "std")
spec[14]  = ugarchspec(variance.model = list(model = "iGARCH"), distribution.model = "std")
spec[15]  = ugarchspec(variance.model = list(model = "csGARCH"), distribution.model = "std")
spec[16]  = ugarchspec(variance.model = list(model = "fGARCH", submodel = "GARCH"), distribution.model = "std")
spec[17]  = ugarchspec(variance.model = list(model = "fGARCH", submodel = "TGARCH"), distribution.model = "std")
spec[18]  = ugarchspec(variance.model = list(model = "fGARCH", submodel = "AVGARCH"), distribution.model = "std")
spec[19]  = ugarchspec(variance.model = list(model = "fGARCH", submodel = "NGARCH"), distribution.model = "std")
spec[20]  = ugarchspec(variance.model = list(model = "fGARCH", submodel = "NAGARCH"), distribution.model = "std")
spec[21]  = ugarchspec(variance.model = list(model = "fGARCH", submodel = "APARCH"), distribution.model = "std")
spec[22]  = ugarchspec(variance.model = list(model = "fGARCH", submodel = "GJRGARCH"), distribution.model = "std")

#----------------------------------------------------------------------
# VaR
#----------------------------------------------------------------------
loss = data.frame(matrix(, nrow=delta, ncol=0))

for (i in 1:length(spec)){
	# Nome dos modelo
	tipo         = spec[[i]]@model$modeldesc$vmodel
	subtipo      = spec[[i]]@model$modeldesc$vsubmodel
	distribuicao = spec[[i]]@model$modeldesc$distribution

	model_name = paste(tipo, subtipo, distribuicao, sep = '_') 
	print(model_name)
	
	# Aplicando a janela móvel
	roll = ugarchroll(spec[[i]], data = tret, 
		n.ahead = 1, 
		n.start = inicio, 
		VaR.alpha = c(0.05), 
		refit.every = 50)

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
spec[23] = UniGASSpec(Dist = 'norm')
spec[24] = UniGASSpec(Dist = 'snorm')
spec[25] = UniGASSpec(Dist = 'std')
spec[26] = UniGASSpec(Dist = 'sstd')

#----------------------------------------------------------------------
# VaR
#----------------------------------------------------------------------
# TODO: Unificar loop com o do GARCH
for (i in 12:15) {
	# Nome dos modelo
	model_name = paste("GAS", spec[[i]]@Spec$Dist, sep = '_') 
		print(model_name)

	# Aplicando a janela móvel
	roll = UniGASRoll(data = tret, 
		GASSpec = spec[[i]], 
		Nstart = inicio, 
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

# https://stats.stackexchange.com/questions/12232/
# calculating-the-parameters-of-a-beta-distribution-using-the-mean-and-variance
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}


#----------------------------------------------------------------------
# VaR
#----------------------------------------------------------------------

VaR_sim = vector(mode = "numeric", length = delta)

for (i in 1:length(VaR_sim)) {
	print(i)

	if (i == 1) {
		draws = svsample(
			tret[seq(i,inicio+i-1)] - m,
			draws = 20000, 
			burnin = 2000, 
			# priormu = c(-10, 1), 
			# priorphi = c(20, 1.2), 
			# priorsigma = 0.2,
			designmatrix = "ar1", 
			quiet = TRUE)
	} else {
		param      = summary(draws)
		priormu    = c(param$para[1,1], param$para[1,2])
		priorphi   = estBetaParams(param$para[2,1], param$para[2,2]^2)
		priorphi   = c(priorphi$alpha, priorphi$beta)
		priorsigma = param$para[3,2]
		draws = svsample(
			tret[seq(i,inicio+i-1)] - m,
			draws = 2000, 
			burnin = 200, 
			priormu = priormu, 
			priorphi = priorphi, 
			priorsigma = priorsigma ,
			designmatrix = "ar1",
			startpara = list(
				mu   = param$para[1,1], 
				phi  = param$para[2,1],
				sigma = param$para[3,1]), 
			quiet = TRUE)
		# plot(draws)
	}
	

	prev = predict(draws, step = 1)

	previsao = arpredict(draws, prev)

	VaR_sim[i] = quantile(previsao, 0.05) + m

}

# plot(tret[1501:2000], type = 'l', col = 'red')
# lines(VaR_sim[1:500])

loss['SV'] = LossVaR(VaR_real, VaR_sim, tau = 0.05)



# _____________________________________________________________________
#
#                                COMPARANDO
# _____________________________________________________________________

comparacao = MCSprocedure(Loss=loss,alpha=0.2,B=5000,statistic='Tmax',cl=NULL)


# _____________________________________________________________________
#
#                                TABELAS
# _____________________________________________________________________

tabela = comparacao@show

tabela[,7] = 1000 * tabela[,7]

sink(file = "./tabelas/sp500.tex")
stargazer(tabela,
          title="Teste StarGazer",
          align=FALSE, dep.var.labels=c("sp500"),
          covariate.labels=c("Modelo",
          	                 "$Rank_{R,M}$"  , "$t_{ij}$", "$p_{R,M}$",
          	                 "$Rank_{max,M}$", "$t_{ij}$", "$p_{max,M}$",
          	                 "$Perda \\times 10^3$"),
          digits = 3,
          omit.stat=c(), no.space=TRUE)
sink()

# xtable(tabela,
# 	Name = c("Modelo",
#           	 "$Rank_{R,M}$"  , "$t_{ij}$", "$p_{R,M}$",
#           	 "$Rank_{max,M}$", "$t_{ij}$", "$p_{max,M}$",
#           	 "$Loss x 10^3$"), 
# 	display = c('s','d','f','f','d','f','f','f'))
