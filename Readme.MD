Práctica final de Genómica Computacional

La carpeta SRC consta de los siguientes archivos:

La base de datos Japonica.RData que también se encuentra disponible en: https://data.cimmyt.org/dataset.xhtml?persistentId=hdl:11529/10548728

La implementación de los métodos de predicción (RCP_GY.R, RF_GY.R) y la prediccióń de los datos incompletos mediante el algoritmo de bosques aleatorios (predictions.R). En particular:


RCP_GY.R crea los archivos:

	PredictionsGY_PCR.csv: Las observaciones que predice la RCP
	Summary_GY_testing_E+O+OxE_PCR.csv: Métricas para evaluar el modelo en cada uno de los ambientes.


RF_GY.R crea los archivos:

	PredictionGY.csv: Las observaciones que predice bosques aleatorios
	Summary_GY_testing_E+G+GxE.csv: Métricas para evaluar el modelo en cada uno de los ambientes.
	Hyperparams.csv: el costo (MSE) de cada combinación evaluada en la calibración de hyperparámetros para la optimización bayesiana.


predictions.R crea el archivo:

	predictions.csv: nos da las predicciones generadas por el algoritmo de bosques aleatorios.

Nota: el código puede tardar en ejecutarse. (Me tomó hasta dos horas).










































































