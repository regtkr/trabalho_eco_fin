#----------------------------------------------------------------------
# BIBLIOTECAS
#----------------------------------------------------------------------
library(fImport)
library(zoo)
library(rugarch)
library(GAS)
library(MCS)
library(timeSeries)


#----------------------------------------------------------------------
# DADOS
#----------------------------------------------------------------------
dados = read.table("/home/regis/Dropbox/econometriafinanças2017/parte3-volatilidade/sp500.txt",header=F)
ret = diff(log(dados$V1))
tret = as.timeSeries(ret)


#----------------------------------------------------------------------
# RECEITA
#----------------------------------------------------------------------
# 1. Fit an ARMA+GARCH model to your return series. You can do this by ugarchfit function in rugarch package in R (you will have to choose what density you want to assume for the standardized residuals).
# 2. Forecast conditional mean and conditional variance of the assumed parametric density 10 days ahead. You can do this by ugarchroll function in rugarch package in R.
# 3. Given the forecasted mean and variance of the assumed density, you can obtain the 0.05 quantile of the distribution which will be your 5% VaR (you can use other quantiles, of course). You can do this by function qnorm for the normal density, qt for student density, or generally qsomethingelse for some other density in R.


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
loss = data.frame(matrix(, nrow=2557, ncol=0))

for (i in 1:length(spec)){
	
	model_name = paste("modelo_", as.character(i), sep = '') 
	print(model_name)

	modelo = ugarchroll(spec[[i]], data = tret, n.ahead = 1, n.start = 5000, VaR.alpha = c(0.05), refit.every = 500)

	data = as.data.frame(modelo, which="VaR")

	loss[model_name] = LossVaR(data[,2], data[,1], tau = 0.05)
}


#----------------------------------------------------------------------
# COMPARANDO
#----------------------------------------------------------------------
comparacao = MCSprocedure(Loss=loss,alpha=0.2,B=5000,statistic='Tmax',cl=NULL)

