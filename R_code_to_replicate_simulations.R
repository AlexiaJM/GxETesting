
########################################################################################################################################################################################
## R code to replicate analyses and figures from "Distinguishing differential susceptibility, diathesis-stress and vantage sensitivity: beyond the single gene and environment model" ##
########################################################################################################################################################################################

# Warning, it takes a while (many hours). I recommend running sim11, sim21, sim31, sim41 on one R process and sim12, sim22, sim32, sim42 on another R process to half the time).

library(LEGIT)

### Functions needed

# Find which sigma value of the residual (normally distributed as N(0,sigma)) to obtain the R-Squared desired for the simulations. (multiple genes and env) 
get_sigma_eps_ex = function(R_2, N=10000, coef_main=c(3,1,2),c=0, beta_param=c(2,2)){
	g1 = rbinom(N,1,.30)
	g2 = rbinom(N,1,.30)
	g3 = rbinom(N,1,.30)
	g4 = rbinom(N,1,.30)
	e1 = rbeta(N,beta_param[1],beta_param[2])*10
	e2 = rbeta(N,beta_param[1],beta_param[2])*10
	e3 = rbeta(N,beta_param[1],beta_param[2])*10
	g = .3*g1 + .10*g2 + .20*g3 + .40*g4
	e = .45*e1 + .35*e2 + .2*e3
	y_true = coef_main[1] + coef_main[2]*(e-c) + coef_main[3]*g*(e-c)
	return(sqrt(((1 - R_2)*var(y_true))/(R_2)))
}
# Find which sigma value of the residual (normally distributed as N(0,sigma)) to obtain the R-Squared desired for the simulations. (single gene and env)
get_sigma_eps_ex2 = function(R_2, N=10000, coef_main=c(3,1,2),c=0, beta_param=c(2,2)){
	g = rbinom(N,1,.30)
	e = rbeta(N,beta_param[1],beta_param[2])*10
	y_true = coef_main[1] + coef_main[2]*(e-c) + coef_main[3]*g*(e-c)
	return(sqrt(((1 - R_2)*var(y_true))/(R_2)))
}

# Example with single gene and env. Note that the example with multiple genes and env is already in the LEGIT package and called "example_with_crossover".
example2_with_crossover = function(N, sigma=1, c = 0, coef_main=c(0,1,2), logit=FALSE, seed=NULL, beta_param=c(2,2)){
	set.seed(seed)
	g = rbinom(N,1,.30)
	e = rbeta(N,beta_param[1],beta_param[2])*10
	y_true = coef_main[1] + coef_main[2]*(e-c) + coef_main[3]*g*(e-c)
	if (logit){
		y_true = 1/(1+exp(-(y_true)))
		y = rbinom(N,1,y_true)
	}
	else{
		eps = rnorm(N,0,sigma)
		y = y_true + eps
	}
	return(list(data=data.frame(y,y_true,g,e), coef_main=coef_main, c=c))
}

# Simulation of confirmary and RoS methods, return accuracy (single gene and env)
simulations_simple_both = function(N = 250, N_iter=100, coef_main = c(3,1,2), sigma = c(1,1,1,1,1,1), criterion ="AICc", seed = 777, crossover = c("min","max"), boot = NULL, reverse_code=FALSE, print=FALSE, beta_param=c(2,2), c=5){
	set.seed(seed)

	models = matrix(0, nrow=6,ncol=7)
	models2 = matrix(0, nrow=6,ncol=4)
	colnames(models) = c("Vantage sensitivity WEAK", "Vantage sensitivity STRONG", "Differential susceptibility WEAK", "Differential susceptibility STRONG","Diathesis-stress WEAK", "Diathesis-stress STRONG","?")
	rownames(models) = c("Vantage sensitivity WEAK", "Vantage sensitivity STRONG", "Differential susceptibility WEAK", "Differential susceptibility STRONG","Diathesis-stress WEAK", "Diathesis-stress STRONG")
	colnames(models2) = c("Vantage sensitivity", "Differential susceptibility","Diathesis-stress","?")
	rownames(models2) = c("Vantage sensitivity WEAK", "Vantage sensitivity STRONG", "Differential susceptibility WEAK", "Differential susceptibility STRONG","Diathesis-stress WEAK", "Diathesis-stress STRONG")
	for (i in 1:N_iter){
		if (print) print(i)
		### Examples where x is in [0, 10]
		## Assuming there is a cross-over point at x=5
		# Vantage sensitivity WEAK
		ex_van = example2_with_crossover(N=N, c=0, coef_main = coef_main, sigma=sigma[1], seed = seed+i, beta_param=beta_param)
		# Vantage sensitivity STRONG
		ex_van_s = example2_with_crossover(N=N, c=0, coef_main = c(coef_main[1],0,coef_main[3]), sigma=sigma[2], seed = seed+i, beta_param=beta_param)		
		# Differential Susceptibility WEAK
		ex_ds = example2_with_crossover(N=N, c=c, coef_main = coef_main, sigma=sigma[3], seed = seed+i, beta_param=beta_param)
		# Differential Susceptibility STRONG
		ex_ds_s = example2_with_crossover(N=N, c=c, coef_main = c(coef_main[1],0,coef_main[3]), sigma=sigma[4], seed = seed+i, beta_param=beta_param)
		# Diathesis Stress WEAK
		ex_dia = example2_with_crossover(N=N, c=10, coef_main = coef_main, sigma=sigma[5], seed = seed+i, beta_param=beta_param)
		# Diathesis Stress STRONG
		ex_dia_s = example2_with_crossover(N=N, c=10, coef_main = c(coef_main[1],0,coef_main[3]), sigma=sigma[6], seed = seed+i, beta_param=beta_param)

		# Automatic testing of all 4 models
		total_n = rep(N_iter,6)
		for (j in 1:6){
			if (j==1) GxE_test = GxE_interaction_test(data.frame(y=ex_van$data$y), data.frame(g=ex_van$data$g), data.frame(e=ex_van$data$e), formula_noGxE = y ~ 1, crossover = crossover, criterion=criterion, boot = boot, reverse_code=reverse_code, maxiter=500)
			if (j==2) GxE_test = GxE_interaction_test(data.frame(y=ex_van_s$data$y), data.frame(g=ex_van_s$data$g), data.frame(e=ex_van_s$data$e), formula_noGxE = y ~ 1, crossover = crossover, criterion=criterion, boot = boot, reverse_code=reverse_code, maxiter=500)
			if (j==3) GxE_test = GxE_interaction_test(data.frame(y=ex_ds$data$y), data.frame(g=ex_ds$data$g), data.frame(e=ex_ds$data$e), formula_noGxE = y ~ 1, crossover = crossover, criterion=criterion, boot = boot, reverse_code=reverse_code, maxiter=500)
			if (j==4) GxE_test = GxE_interaction_test(data.frame(y=ex_ds_s$data$y), data.frame(g=ex_ds_s$data$g), data.frame(e=ex_ds_s$data$e), formula_noGxE = y ~ 1, crossover = crossover, criterion=criterion, boot = boot, reverse_code=reverse_code, maxiter=500)
			if (j==5) GxE_test = GxE_interaction_test(data.frame(y=ex_dia$data$y), data.frame(g=ex_dia$data$g), data.frame(e=ex_dia$data$e), formula_noGxE = y ~ 1, crossover = crossover, criterion=criterion, boot = boot, reverse_code=reverse_code, maxiter=500)
			if (j==6) GxE_test = GxE_interaction_test(data.frame(y=ex_dia_s$data$y), data.frame(g=ex_dia_s$data$g), data.frame(e=ex_dia_s$data$e), formula_noGxE = y ~ 1, crossover = crossover, criterion=criterion, boot = boot, reverse_code=reverse_code, maxiter=500)
			
			if (j==1) GxE_test2 = GxE_interaction_RoS(data.frame(y=ex_van$data$y), data.frame(g=ex_van$data$g), data.frame(e=ex_van$data$e), formula_noGxE = y ~ 1, reverse_code=reverse_code, maxiter=500)
			if (j==2) GxE_test2 = GxE_interaction_RoS(data.frame(y=ex_van_s$data$y), data.frame(g=ex_van_s$data$g), data.frame(e=ex_van_s$data$e), formula_noGxE = y ~ 1, reverse_code=reverse_code, maxiter=500)
			if (j==3) GxE_test2 = GxE_interaction_RoS(data.frame(y=ex_ds$data$y), data.frame(g=ex_ds$data$g), data.frame(e=ex_ds$data$e), formula_noGxE = y ~ 1, reverse_code=reverse_code, maxiter=500)
			if (j==4) GxE_test2 = GxE_interaction_RoS(data.frame(y=ex_ds_s$data$y), data.frame(g=ex_ds_s$data$g), data.frame(e=ex_ds_s$data$e), formula_noGxE = y ~ 1, reverse_code=reverse_code, maxiter=500)
			if (j==5) GxE_test2 = GxE_interaction_RoS(data.frame(y=ex_dia$data$y), data.frame(g=ex_dia$data$g), data.frame(e=ex_dia$data$e), formula_noGxE = y ~ 1, reverse_code=reverse_code, maxiter=500)
			if (j==6) GxE_test2 = GxE_interaction_RoS(data.frame(y=ex_dia_s$data$y), data.frame(g=ex_dia_s$data$g), data.frame(e=ex_dia_s$data$e), formula_noGxE = y ~ 1, reverse_code=reverse_code, maxiter=500)

			if (print) print(GxE_test$results)

			## Accuracy of Confirmatory
			# Check if non-convergence, if so do not include result
			if (sum(sapply(GxE_test$fits, function (x) x$conv))==6){

				GxE_test_ = GxE_test$results[GxE_test$results[,4]!="No",]
				GxE_test_best_name = rownames(GxE_test_[1,,drop=FALSE])
				if (GxE_test_best_name == "Vantage sensitivity WEAK"){
					GxE_test_best = GxE_test$fits$vantage_sensitivity_WEAK
					best_i = 1
				}
				else if (GxE_test_best_name == "Vantage sensitivity STRONG"){
					GxE_test_best = GxE_test$fits$vantage_sensitivity_STRONG
					best_i = 2
				}
				else if (GxE_test_best_name == "Differential susceptibility WEAK"){
					GxE_test_best = GxE_test$fits$diff_suscept_WEAK
					best_i = 3
				}
				else if (GxE_test_best_name == "Differential susceptibility STRONG"){
					GxE_test_best = GxE_test$fits$diff_suscept_STRONG
					best_i = 4
				}
				else if (GxE_test_best_name == "Diathesis-stress WEAK"){
					GxE_test_best = GxE_test$fits$diathesis_stress_WEAK
					best_i = 5
				}
				else if (GxE_test_best_name == "Diathesis-stress STRONG"){
					GxE_test_best = GxE_test$fits$diathesis_stress_STRONG
					best_i = 6
				}
			}
			else{
				best_i = 7
			}
			models[j,best_i] = models[j,best_i] + 1
			if (i == N_iter) models[j,] = models[j,] / total_n[j] * 100

			## Accuracy of RoS
			best_i = 0
			if (GxE_test2$int_type=="Vantage sensitivity") best_i = 1
			else if (GxE_test2$int_type=="Differential susceptibility") best_i = 2
			else if (GxE_test2$int_type=="Diathesis-stress") best_i = 3
			else best_i = 4
			models2[j,best_i] = models2[j,best_i] + 1
			if (i == N_iter) models2[j,] = models2[j,] / N_iter * 100
		}
	}
	return(list(models,models2))
}

# Simulation of confirmary and RoS methods, return accuracy (multiple genes and env)
simulations_both = function(N = 250, N_iter=100, coef_main = c(3,1,2), sigma = c(1,1,1,1,1,1), criterion ="AICc", seed = 777, crossover = c("min","max"), boot = NULL, reverse_code=FALSE, print=FALSE, beta_param=c(2,2), c=5){
	set.seed(seed)

	models = matrix(0, nrow=6,ncol=7)
	models2 = matrix(0, nrow=6,ncol=4)
	colnames(models) = c("Vantage sensitivity WEAK", "Vantage sensitivity STRONG", "Differential susceptibility WEAK", "Differential susceptibility STRONG","Diathesis-stress WEAK", "Diathesis-stress STRONG","?")
	rownames(models) = c("Vantage sensitivity WEAK", "Vantage sensitivity STRONG", "Differential susceptibility WEAK", "Differential susceptibility STRONG","Diathesis-stress WEAK", "Diathesis-stress STRONG")
	colnames(models2) = c("Vantage sensitivity", "Differential susceptibility","Diathesis-stress","?")
	rownames(models2) = c("Vantage sensitivity WEAK", "Vantage sensitivity STRONG", "Differential susceptibility WEAK", "Differential susceptibility STRONG","Diathesis-stress WEAK", "Diathesis-stress STRONG")
	for (i in 1:N_iter){
		if (print) print(i)
		### Examples where x is in [0, 10]
		## Assuming there is a cross-over point at x=5
		# Vantage sensitivity WEAK
		ex_van = example_with_crossover(N=N, c=0, coef_main = coef_main, sigma=sigma[1], seed = seed+i, beta_param=beta_param)
		# Vantage sensitivity STRONG
		ex_van_s = example_with_crossover(N=N, c=0, coef_main = c(coef_main[1],0,coef_main[3]), sigma=sigma[2], seed = seed+i, beta_param=beta_param)		
		# Differential Susceptibility WEAK
		ex_ds = example_with_crossover(N=N, c=c, coef_main = coef_main, sigma=sigma[3], seed = seed+i, beta_param=beta_param)
		# Differential Susceptibility STRONG
		ex_ds_s = example_with_crossover(N=N, c=c, coef_main = c(coef_main[1],0,coef_main[3]), sigma=sigma[4], seed = seed+i, beta_param=beta_param)
		# Diathesis Stress WEAK
		ex_dia = example_with_crossover(N=N, c=10, coef_main = coef_main, sigma=sigma[5], seed = seed+i, beta_param=beta_param)
		# Diathesis Stress STRONG
		ex_dia_s = example_with_crossover(N=N, c=10, coef_main = c(coef_main[1],0,coef_main[3]), sigma=sigma[6], seed = seed+i, beta_param=beta_param)

		# Automatic testing of all 4 models
		total_n = rep(N_iter,6)
		for (j in 1:6){
			if (j==1) GxE_test = GxE_interaction_test(ex_van$data, ex_van$G, ex_van$E, formula_noGxE = y ~ 1, crossover = crossover, criterion=criterion, boot = boot, reverse_code=reverse_code, maxiter=500)
			if (j==2) GxE_test = GxE_interaction_test(ex_van_s$data, ex_van_s$G, ex_van_s$E, formula_noGxE = y ~ 1, crossover = crossover, criterion=criterion, boot = boot, reverse_code=reverse_code, maxiter=500)
			if (j==3) GxE_test = GxE_interaction_test(ex_ds$data, ex_ds$G, ex_ds$E, formula_noGxE = y ~ 1, crossover = crossover, criterion=criterion, boot = boot, reverse_code=reverse_code, maxiter=500)
			if (j==4) GxE_test = GxE_interaction_test(ex_ds_s$data, ex_ds_s$G, ex_ds_s$E, formula_noGxE = y ~ 1, crossover = crossover, criterion=criterion, boot = boot, reverse_code=reverse_code, maxiter=500)
			if (j==5) GxE_test = GxE_interaction_test(ex_dia$data, ex_dia$G, ex_dia$E, formula_noGxE = y ~ 1, crossover = crossover, criterion=criterion, boot = boot, reverse_code=reverse_code, maxiter=500)
			if (j==6) GxE_test = GxE_interaction_test(ex_dia_s$data, ex_dia_s$G, ex_dia_s$E, formula_noGxE = y ~ 1, crossover = crossover, criterion=criterion, boot = boot, reverse_code=reverse_code, maxiter=500)

			if (j==1) GxE_test2 = GxE_interaction_RoS(ex_van$data, ex_van$G, ex_van$E, formula_noGxE = y ~ 1, reverse_code=reverse_code, maxiter=500)
			if (j==2) GxE_test2 = GxE_interaction_RoS(ex_van_s$data, ex_van_s$G, ex_van_s$E, formula_noGxE = y ~ 1, reverse_code=reverse_code, maxiter=500)
			if (j==3) GxE_test2 = GxE_interaction_RoS(ex_ds$data, ex_ds$G, ex_ds$E, formula_noGxE = y ~ 1, reverse_code=reverse_code, maxiter=500)
			if (j==4) GxE_test2 = GxE_interaction_RoS(ex_ds_s$data, ex_ds_s$G, ex_ds_s$E, formula_noGxE = y ~ 1, reverse_code=reverse_code, maxiter=500)
			if (j==5) GxE_test2 = GxE_interaction_RoS(ex_dia$data, ex_dia$G, ex_dia$E, formula_noGxE = y ~ 1, reverse_code=reverse_code, maxiter=500)
			if (j==6) GxE_test2 = GxE_interaction_RoS(ex_dia_s$data, ex_dia_s$G, ex_dia_s$E, formula_noGxE = y ~ 1, reverse_code=reverse_code, maxiter=500)

			if (print) print(GxE_test$results)

			## Test 1
			# Check if non-convergence have happenned, if so do not include result
			if (sum(sapply(GxE_test$fits, function (x) x$conv))==6){

				GxE_test_ = GxE_test$results[GxE_test$results[,4]!="No",]
				GxE_test_best_name = rownames(GxE_test_[1,,drop=FALSE])
				if (GxE_test_best_name == "Vantage sensitivity WEAK"){
					GxE_test_best = GxE_test$fits$vantage_sensitivity_WEAK
					best_i = 1
				}
				else if (GxE_test_best_name == "Vantage sensitivity STRONG"){
					GxE_test_best = GxE_test$fits$vantage_sensitivity_STRONG
					best_i = 2
				}
				else if (GxE_test_best_name == "Differential susceptibility WEAK"){
					GxE_test_best = GxE_test$fits$diff_suscept_WEAK
					best_i = 3
				}
				else if (GxE_test_best_name == "Differential susceptibility STRONG"){
					GxE_test_best = GxE_test$fits$diff_suscept_STRONG
					best_i = 4
				}
				else if (GxE_test_best_name == "Diathesis-stress WEAK"){
					GxE_test_best = GxE_test$fits$diathesis_stress_WEAK
					best_i = 5
				}
				else if (GxE_test_best_name == "Diathesis-stress STRONG"){
					GxE_test_best = GxE_test$fits$diathesis_stress_STRONG
					best_i = 6
				}
			}
			else{
				best_i = 7
			}
			models[j,best_i] = models[j,best_i] + 1
			if (i == N_iter) models[j,] = models[j,] / total_n[j] * 100

			## Test 2
			best_i = 0
			if (GxE_test2$int_type=="Vantage sensitivity") best_i = 1
			else if (GxE_test2$int_type=="Differential susceptibility") best_i = 2
			else if (GxE_test2$int_type=="Diathesis-stress") best_i = 3
			else best_i = 4
			models2[j,best_i] = models2[j,best_i] + 1
			if (i == N_iter) models2[j,] = models2[j,] / N_iter * 100
		}
	}
	return(list(models,models2))
}

# Run simulation of confirmary and RoS methods on multiple scenarios (N, R2, coefficients). Returns the accuracy of both approaches (disregarding weak vs strong). (single gene and env)
sim_R2byN_simple = function(N=c(100,250,500,1000,2000), N_iter=100, R2 = c(.05, .10, .15), coef_main = c(3,1,2), criterion ="BIC", seed = 777, crossover = c("min","max"), boot = NULL, reverse_code=FALSE, print=FALSE, beta_param=c(2,2), c=5){
	accuracy1 = matrix(0, nrow=length(N),ncol=length(R2))
	accuracy2 = matrix(0, nrow=length(N),ncol=length(R2))
	for (j in 1:length(N)){
		for (k in 1:length(R2)){

			sigma = c(get_sigma_eps_ex2(R2[k], coef_main=coef_main,c=0, beta_param=beta_param), get_sigma_eps_ex2(R2[k], coef_main=c(coef_main[1],0,coef_main[3]),c=0, beta_param=beta_param), get_sigma_eps_ex2(R2[k], coef_main=coef_main,c=c, beta_param=beta_param), get_sigma_eps_ex2(R2[k], coef_main=c(coef_main[1],0,coef_main[3]),c=c, beta_param=beta_param), get_sigma_eps_ex2(R2[k], coef_main=coef_main,c=10, beta_param=beta_param), get_sigma_eps_ex2(R2[k], coef_main=c(coef_main[1],0,coef_main[3]),c=10, beta_param=beta_param))
			results = simulations_simple_both(N = N[j], N_iter=N_iter, coef_main = coef_main, sigma = sigma, criterion =criterion, seed = seed, crossover = crossover, boot = boot, reverse_code=reverse_code, print=print, beta_param=beta_param,c=c)
			
			for (i in 1:6){
				if (i<=2) accuracy1[j,k] =  accuracy1[j,k] + results[[1]][i,1] + results[[1]][i,2]
				else if (i<=4) accuracy1[j,k] =  accuracy1[j,k] + results[[1]][i,3] + results[[1]][i,4]
				else accuracy1[j,k] =  accuracy1[j,k] + results[[1]][i,5] + results[[1]][i,6]
			}
			accuracy1[j,k] = accuracy1[j,k]/6

			for (i in 1:6){
				if (i<=2) accuracy2[j,k] =  accuracy2[j,k] + results[[2]][i,1]
				else if (i<=4) accuracy2[j,k] =  accuracy2[j,k] + results[[2]][i,2]
				else accuracy2[j,k] =  accuracy2[j,k] + results[[2]][i,3]
			}
			accuracy2[j,k] = accuracy2[j,k]/6
		}
	}
	return(list(accuracy_confirmatory=accuracy1,accuracy_RoS=accuracy2))
}
# Run simulation of confirmary and RoS methods on multiple scenarios (N, R2, coefficients). Returns the accuracy of both approaches (disregarding weak vs strong). (multiple genes and env)
sim_R2byN = function(N=c(100,250,500,1000,2000), N_iter=100, R2 = c(.10, .20, .40), coef_main = c(3,1,2), criterion ="BIC", seed = 777, crossover = c("min","max"), boot = NULL, reverse_code=FALSE, print=FALSE, beta_param=c(2,2), c=5){
	accuracy1 = matrix(0, nrow=length(N),ncol=length(R2))
	accuracy2 = matrix(0, nrow=length(N),ncol=length(R2))
	for (j in 1:length(N)){
		for (k in 1:length(R2)){

			sigma = c(get_sigma_eps_ex(R2[k], coef_main=coef_main,c=0, beta_param=beta_param), get_sigma_eps_ex(R2[k], coef_main=c(coef_main[1],0,coef_main[3]),c=0, beta_param=beta_param), get_sigma_eps_ex(R2[k], coef_main=coef_main,c=c, beta_param=beta_param), get_sigma_eps_ex(R2[k], coef_main=c(coef_main[1],0,coef_main[3]),c=c, beta_param=beta_param), get_sigma_eps_ex(R2[k], coef_main=coef_main,c=10, beta_param=beta_param), get_sigma_eps_ex(R2[k], coef_main=c(coef_main[1],0,coef_main[3]),c=10, beta_param=beta_param))
			results = simulations_both(N = N[j], N_iter=N_iter, coef_main = coef_main, sigma = sigma, criterion =criterion, seed = seed, crossover = crossover, boot = boot, reverse_code=reverse_code, print=print, beta_param=beta_param,c=c)
			
			for (i in 1:6){
				if (i<=2) accuracy1[j,k] =  accuracy1[j,k] + results[[1]][i,1] + results[[1]][i,2]
				else if (i<=4) accuracy1[j,k] =  accuracy1[j,k] + results[[1]][i,3] + results[[1]][i,4]
				else accuracy1[j,k] =  accuracy1[j,k] + results[[1]][i,5] + results[[1]][i,6]
			}
			accuracy1[j,k] = accuracy1[j,k]/6

			for (i in 1:6){
				if (i<=2) accuracy2[j,k] =  accuracy2[j,k] + results[[2]][i,1]
				else if (i<=4) accuracy2[j,k] =  accuracy2[j,k] + results[[2]][i,2]
				else accuracy2[j,k] =  accuracy2[j,k] + results[[2]][i,3]
			}
			accuracy2[j,k] = accuracy2[j,k]/6
		}
	}
	return(list(accuracy_confirmatory=accuracy1,accuracy_RoS=accuracy2))
}

# Plotting function for the Figures
plot_figs = function(sim1, sim2, sim3, sim4, i, h=10, w=6, res=75, effect_size=c(.10,.20,.40)){

	tiff(paste0("plot_GxE_sim",i,".tiff"), width=w*res, height=h*res, res=res, type='cairo')

	par(mfrow=c(4,2),mar=c(1, 4, 3, 1) + 0.1)
	par(oma = c(8.5, 12, 1, 1))
	par(xpd = NA)

	plot(c(),c(),xlim=c(0,2000),ylim=c(0,100),xlab="",ylab="Accuracy (%)",main="Competitive-confirmatory", cex.axis=1.3, cex.lab=1.7, cex.main=1.75)
	axis(side=2, at=100, labels=100, cex.axis=1.3)
	lines(c(-50,2050),c(90,90),lty=3, col="gray50")
	lines(c(100,250,500,1000,2000),sim1$accuracy_confirmatory[,1], col="blue", lwd=1.5)
	points(c(100,250,500,1000,2000),sim1$accuracy_confirmatory[,1], col="blue",pch=15, cex=1.5)
	lines(c(100,250,500,1000,2000),sim1$accuracy_confirmatory[,2], col="red", lwd=1.5)
	points(c(100,250,500,1000,2000),sim1$accuracy_confirmatory[,2], col="red",pch=16, cex=1.5)
	lines(c(100,250,500,1000,2000),sim1$accuracy_confirmatory[,3], col="black", lwd=1.5)
	points(c(100,250,500,1000,2000),sim1$accuracy_confirmatory[,3], col="black",pch=17, cex=1.5)

	plot(c(),c(),xlim=c(0,2000),ylim=c(0,100),xlab="",ylab="",main="Regions of significance", cex.axis=1.3, cex.lab=1.7, cex.main=1.75)
	axis(side=2, at=100, labels=100, cex.axis=1.3)
	lines(c(-50,2050),c(90,90),lty=3, col="gray50")
	lines(c(100,250,500,1000,2000),sim1$accuracy_RoS[,1], col="blue", lwd=1.5)
	points(c(100,250,500,1000,2000),sim1$accuracy_RoS[,1], col="blue",pch=15, cex=1.5)
	lines(c(100,250,500,1000,2000),sim1$accuracy_RoS[,2], col="red", lwd=1.5)
	points(c(100,250,500,1000,2000),sim1$accuracy_RoS[,2], col="red",pch=16, cex=1.5)
	lines(c(100,250,500,1000,2000),sim1$accuracy_RoS[,3], col="black", lwd=1.5)
	points(c(100,250,500,1000,2000),sim1$accuracy_RoS[,3], col="black",pch=17, cex=1.5)

	par(mar=c(1, 4, 3, 1) + 0.1)

	plot(c(),c(),xlim=c(0,2000),ylim=c(0,100),xlab="",ylab="Accuracy (%)",main="", cex.axis=1.3, cex.lab=1.7, cex.main=1.5)
	axis(side=2, at=100, labels=100, cex.axis=1.3)
	lines(c(-50,2050),c(90,90),lty=3, col="gray50")
	lines(c(100,250,500,1000,2000),sim2$accuracy_confirmatory[,1], col="blue", lwd=1.5)
	points(c(100,250,500,1000,2000),sim2$accuracy_confirmatory[,1], col="blue",pch=15, cex=1.5)
	lines(c(100,250,500,1000,2000),sim2$accuracy_confirmatory[,2], col="red", lwd=1.5)
	points(c(100,250,500,1000,2000),sim2$accuracy_confirmatory[,2], col="red",pch=16, cex=1.5)
	lines(c(100,250,500,1000,2000),sim2$accuracy_confirmatory[,3], col="black", lwd=1.5)
	points(c(100,250,500,1000,2000),sim2$accuracy_confirmatory[,3], col="black",pch=17, cex=1.5)

	plot(c(),c(),xlim=c(0,2000),ylim=c(0,100),xlab="",ylab="",main="", cex.axis=1.3, cex.lab=1.7, cex.main=1.5)
	axis(side=2, at=100, labels=100, cex.axis=1.3)
	lines(c(-50,2050),c(90,90),lty=3, col="gray50")
	lines(c(100,250,500,1000,2000),sim2$accuracy_RoS[,1], col="blue", lwd=1.5)
	points(c(100,250,500,1000,2000),sim2$accuracy_RoS[,1], col="blue",pch=15, cex=1.5)
	lines(c(100,250,500,1000,2000),sim2$accuracy_RoS[,2], col="red", lwd=1.5)
	points(c(100,250,500,1000,2000),sim2$accuracy_RoS[,2], col="red",pch=16, cex=1.5)
	lines(c(100,250,500,1000,2000),sim2$accuracy_RoS[,3], col="black", lwd=1.5)
	points(c(100,250,500,1000,2000),sim2$accuracy_RoS[,3], col="black",pch=17, cex=1.5)

	par(mar=c(1, 4, 3, 1) + 0.1)

	plot(c(),c(),xlim=c(0,2000),ylim=c(0,100),xlab="",ylab="Accuracy (%)",main="", cex.axis=1.3, cex.lab=1.7, cex.main=1.5)
	axis(side=2, at=100, labels=100, cex.axis=1.3)
	lines(c(-50,2050),c(90,90),lty=3, col="gray50")
	lines(c(100,250,500,1000,2000),sim3$accuracy_confirmatory[,1], col="blue", lwd=1.5)
	points(c(100,250,500,1000,2000),sim3$accuracy_confirmatory[,1], col="blue",pch=15, cex=1.5)
	lines(c(100,250,500,1000,2000),sim3$accuracy_confirmatory[,2], col="red", lwd=1.5)
	points(c(100,250,500,1000,2000),sim3$accuracy_confirmatory[,2], col="red",pch=16, cex=1.5)
	lines(c(100,250,500,1000,2000),sim3$accuracy_confirmatory[,3], col="black", lwd=1.5)
	points(c(100,250,500,1000,2000),sim3$accuracy_confirmatory[,3], col="black",pch=17, cex=1.5)

	plot(c(),c(),xlim=c(0,2000),ylim=c(0,100),xlab="",ylab="",main="", cex.axis=1.3, cex.lab=1.7, cex.main=1.5)
	axis(side=2, at=100, labels=100, cex.axis=1.3)
	lines(c(-50,2050),c(90,90),lty=3, col="gray50")
	lines(c(100,250,500,1000,2000),sim3$accuracy_RoS[,1], col="blue", lwd=1.5)
	points(c(100,250,500,1000,2000),sim3$accuracy_RoS[,1], col="blue",pch=15, cex=1.5)
	lines(c(100,250,500,1000,2000),sim3$accuracy_RoS[,2], col="red", lwd=1.5)
	points(c(100,250,500,1000,2000),sim3$accuracy_RoS[,2], col="red",pch=16, cex=1.5)
	lines(c(100,250,500,1000,2000),sim3$accuracy_RoS[,3], col="black", lwd=1.5)
	points(c(100,250,500,1000,2000),sim3$accuracy_RoS[,3], col="black",pch=17, cex=1.5)

	par(mar=c(1, 4, 3, 1) + 0.1)

	plot(c(),c(),xlim=c(0,2000),ylim=c(0,100),xlab="Sample size (N)",ylab="Accuracy (%)", cex.axis=1.3, cex.lab=1.7, cex.main=1.5)
	axis(side=2, at=100, labels=100, cex.axis=1.3)
	lines(c(-50,2050),c(90,90),lty=3, col="gray50")
	lines(c(100,250,500,1000,2000),sim4$accuracy_confirmatory[,1], col="blue", lwd=1.5)
	points(c(100,250,500,1000,2000),sim4$accuracy_confirmatory[,1], col="blue",pch=15, cex=1.5)
	lines(c(100,250,500,1000,2000),sim4$accuracy_confirmatory[,2], col="red", lwd=1.5)
	points(c(100,250,500,1000,2000),sim4$accuracy_confirmatory[,2], col="red",pch=16, cex=1.5)
	lines(c(100,250,500,1000,2000),sim4$accuracy_confirmatory[,3], col="black", lwd=1.5)
	points(c(100,250,500,1000,2000),sim4$accuracy_confirmatory[,3], col="black",pch=17, cex=1.5)

	plot(c(),c(),xlim=c(0,2000),ylim=c(0,100),xlab="Sample size (N)",ylab="",main="", cex.axis=1.3, cex.lab=1.7, cex.main=1.5)
	axis(side=2, at=100, labels=100, cex.axis=1.3)
	lines(c(-50,2050),c(90,90),lty=3, col="gray50")
	lines(c(100,250,500,1000,2000),sim4$accuracy_RoS[,1], col="blue", lwd=1.5)
	points(c(100,250,500,1000,2000),sim4$accuracy_RoS[,1], col="blue",pch=15, cex=1.5)
	lines(c(100,250,500,1000,2000),sim4$accuracy_RoS[,2], col="red", lwd=1.5)
	points(c(100,250,500,1000,2000),sim4$accuracy_RoS[,2], col="red",pch=16, cex=1.5)
	lines(c(100,250,500,1000,2000),sim4$accuracy_RoS[,3], col="black", lwd=1.5)
	points(c(100,250,500,1000,2000),sim4$accuracy_RoS[,3], col="black",pch=17, cex=1.5)

	par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
	plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
	legend("bottom",pch=c(15,16,17),lwd=c(1,1),col=c("blue","red","black"), legend=c(paste0("Small (",effect_size[1],")"),paste0("Medium (",effect_size[2],")"),paste0("Large (",effect_size[3],")")),title=expression(bold("Effect size")), xpd = TRUE, inset = c(0,0), bty = "n", horiz=TRUE, cex=1.7)
	x = c(-.85,-.85,-.85,-.85)
	y = .80-.475*c(0,1,2,3)
	text(x,y,labels=c("Symmetric E\nc=.5","Left-skewed E\nc=.5","Symmetric E\nc=.25","Left-skewed E\nc=.25"), cex=1.6,font=2)
	dev.off()
}

### Analyses

# Symmetric Beta(2,2), c = 5, single genes/env
sim11 = sim_R2byN_simple(N=c(100,250,500,1000,2000), R2 = c(.05,.10,.15), N_iter=100, criterion="BIC", reverse_code=FALSE, beta_param=c(2,2), c=5)
# Symmetric Beta(2,2), c = 5, multiple genes/env
sim21 = sim_R2byN(N=c(100,250,500,1000,2000), R2 = c(.10,.20,.40), N_iter=100, criterion="BIC", reverse_code=FALSE, beta_param=c(2,2), c=5)
# Skewed Beta(2,4), c = 5, single genes/env
sim31 = sim_R2byN_simple(N=c(100,250,500,1000,2000), R2 = c(.05,.10,.15), N_iter=100, criterion="BIC", reverse_code=FALSE, beta_param=c(2,4), c=5)
# Skewed Beta(2,4), c = 5, multiple genes/env
sim41 = sim_R2byN(N=c(100,250,500,1000,2000), R2 = c(.10,.20,.40), N_iter=100, criterion="BIC", reverse_code=FALSE, beta_param=c(2,4), c=5)

# Symmetric Beta(2,2), c = 2.5, single genes/env
sim12 = sim_R2byN_simple(N=c(100,250,500,1000,2000), R2 = c(.05,.10,.15), N_iter=100, criterion="BIC", reverse_code=FALSE, beta_param=c(2,2), c=2.5)
# Symmetric Beta(2,2), c = 2.5, multiple genes/env
sim22 = sim_R2byN(N=c(100,250,500,1000,2000), R2 = c(.10,.20,.40), N_iter=100, criterion="BIC", reverse_code=FALSE, beta_param=c(2,2), c=2.5)
# Skewed Beta(2,4), c = 2.5, single genes/env
sim32 = sim_R2byN_simple(N=c(100,250,500,1000,2000), R2 = c(.05,.10,.15), N_iter=100, criterion="BIC", reverse_code=FALSE, beta_param=c(2,4), c=2.5)
# Skewed Beta(2,4), c = 2.5, multiple genes/env
sim42 = sim_R2byN(N=c(100,250,500,1000,2000), R2 = c(.10,.20,.40), N_iter=100, criterion="BIC", reverse_code=FALSE, beta_param=c(2,4), c=2.5)

### Save and reload

saveRDS(sim11, "sim11.rds")
saveRDS(sim21, "sim21.rds")
saveRDS(sim31, "sim31.rds")
saveRDS(sim41, "sim41.rds")
saveRDS(sim12, "sim12.rds")
saveRDS(sim22, "sim22.rds")
saveRDS(sim32, "sim32.rds")
saveRDS(sim42, "sim42.rds")

sim11 <- readRDS("sim11.rds")
sim21 <- readRDS("sim21.rds")
sim31 <- readRDS("sim31.rds")
sim41 <- readRDS("sim41.rds")
sim12 <- readRDS("sim12.rds")
sim22 <- readRDS("sim22.rds")
sim32 <- readRDS("sim32.rds")
sim42 <- readRDS("sim42.rds")

### Plots

plot_figs(sim11, sim31, sim12, sim32, "_single_ge",w=7,h=9.7, res=200, effect_size=c(.05,.10,.15))
plot_figs(sim21, sim41, sim22, sim42, "_multiple_ge",w=7,h=9.7, res=200, effect_size=c(.10,.20,.40))