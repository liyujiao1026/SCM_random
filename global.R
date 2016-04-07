# script SCM_Function.R
#
# the purpose of this script is to simulate the Synthetic Control Method
#
# first version: 160401
# this version:  160402
# last change by: Yujiao Li


library(shiny)
library(shinydashboard)
library(ggplot2)
library(tidyr)
library(Synth)
library(kernlab)
library(LowRankQP)

#--- Funcion: rw.sim() ----------------------------------------------------------------------------------#
# rw.sim() generate the random walk
# Random Walk model: Y_t = tau * (Y_t-1) + et
# @timeLength: the length of time series
# @seed: the random seed

rw.sim <- function(timeLength, seed, tau) {
            n <- timeLength
            set.seed(seed)
            e <- rnorm(n, mean = 0, sd = 1)
            y <- c()
            y[1] <- e[1]
            
            for (i in 2:n) {
                        y[i] <- tau *y[i - 1] + e[i]
            }
            
            return(ts(y))
            
}

# test
# rw.sim(timeLength = 10,tau=0.5,seed = 2)
# sapply(c(1,2,3), rw.sim, tau=1, timeLength=10)
#--------------------------------------------------------------------------------------------------------#




#==================================================================================================#






#--- Funcion: dataRW.sim() -----------------------------------------------------------------------------#
# dataRW.sim() generate the controlUnits with the model of Random Walk

# Model: Zt = Z0 + Yt,
#        Z0 is the common fixed random walk for each Zt,
#        Yt is random walk varied by different seed ,
# return: Zt

# @controlNumber: the number of control units
# @timeLength: the length of time period

dataRW.sim <- function(controlNumber, timeLength, 
                       tau, omiga,  seed.control) {
            
            set.seed(seed.control)
            seed.control.series <- sample(1:100000 , controlNumber )
            
            Z0 <- rw.sim(seed = seed.control, tau = tau, timeLength = timeLength)
            Yt <- sapply(seed.control.series , rw.sim, timeLength = timeLength, tau = tau)
            Zt <- apply(Yt, 2, function(x) {
                        omiga * Z0 + (1 - omiga) * x
            })
            data <- data.frame(Zt)
            
            rownames(data) <- paste("Time", 1:timeLength , sep = "_")
            colnames(data) <- paste("ControlUnit", 1:controlNumber ,sep = "_")
            
            return(data)
}


#test
dataRW <- dataRW.sim( controlNumber = 5, timeLength = 10, tau = 0.5, omiga = 0.3, seed.control = 1 )
#--------------------------------------------------------------------------------------------------------#






#==================================================================================================#





#--- Funcion: y.obs() ----------------------------------------------------------------------------------#
# y.obs() generate post-intervention observation given the pre-intervention data

# Model: Y_t = (Y_t-1 + e_t) * (1+alpha) ,
#        alpha is the assumed caused effect from intervention

# return: y.obs

# @ y.hat: the equally weighted average of other control units
# @ intervention_Time : which time point is the intervention time
# @ alpha : the caused effect
# @ seed.post : random walk for the post-intervention observations

y.obs <- function(y.hat,tau, intervention_Time, alpha) {
            
            t.int <- intervention_Time
            y.observation <- c()
            
            # (1)==> ( y[1], y[2], ..., y[t] ) ------------#
            # before the intervention, y.hat==y.observation
            
            y.observation[1:(t.int)] <- y.hat[1:(t.int)] 
            l <- length(y.hat)  
            
            # (2)==> y[t+1] -----------------------------------#
            # y.observation : only 1 year after the intervetion
            
            y.observation[(t.int + 1)] <- y.hat[(t.int + 1)] + alpha
            
            
            # (3)==> (y[t+2], y[t+3], y[t+4], ....) -------------#
            # y.observation : several years after the intervetion   
            
            y.observation[(t.int + 2):l] <- y.hat[(t.int + 2):l] + alpha
            
            
            return(y.observation)
}

# test
#

# controlNumber = 5; timeLength = 10; intervention_Time = 5;
# tau = 0.5; omiga = 0.3; seed.control = 1; alpha = 0.5
# dataRW <- dataRW.sim( controlNumber = controlNumber, timeLength = timeLength,
#                       tau = tau, omiga = omiga, seed.control = seed.control )
# y.hat <- apply(dataRW, 1, mean)
# y.obs.v <- y.obs(y.hat = y.hat, tau = tau, intervention_Time = intervention_Time ,alpha = alpha)


#==================================================================================================#

# Compared with y.obs()
# when we regard the observed value is the random varialbe 
# y.obs.random <- y.obs + epsilo
# E (y.obs.random) = y.obs
# Var(y.obs.random) = var(sum(w^2 * sigma(z)  +  rou^2 * sigma(epsilo)))
# sigma(z) = sqrt(1/(1-tau^2))




y.obs.random <- function(rou, tau,
                         controlNumber, timeLength, y.obs.v, seed.obs.random) {
            
            set.seed(seed.obs.random)
            
            if (tau != 1) {
                        var_y.obs <- (1 / controlNumber) * (1 / (1 - tau ^ 2))  +  (rou ^ 2)
                        
                        y.obs.random <- c()
                        for (i in 1:timeLength) {
                                    y.obs.random[i] <- rnorm(1, mean = y.obs.v[i], sd = sqrt(var_y.obs))
                        }
            }else {     
                        y.obs.random <- c()
                        
                        for (k in 1:timeLength) {  
                                    var_y.obs.k <- (1 / controlNumber) * (k)  +  (rou ^ 2)
                                    y.obs.random[k] <- rnorm(1, mean = y.obs.v[k], sd = sqrt(var_y.obs.k))
                        }
            }
            
            return(y.obs.random)
}

# #generate data

seed.obs.random=23


controlNumber = 10; timeLength = 20; intervention_Time = 15;
tau = 0.1; omiga = 0.5; seed.control = 1; alpha = 1
dataRW <- dataRW.sim( controlNumber = controlNumber, timeLength = timeLength,
                      tau = tau, omiga = omiga, seed.control = seed.control )
y.hat <- apply(dataRW, 1, mean)
y.obs.v <- y.obs(y.hat = y.hat, tau = tau, intervention_Time = intervention_Time ,alpha = alpha)
#curren test
rou = 0; tau = 0.6; 

y.obs.ran <- y.obs.random(rou = rou, timeLength = timeLength , tau = tau,
                          seed.obs.random = seed.obs.random,
                          controlNumber = controlNumber , y.obs.v = y.obs.v)

plot(y.obs.ran, type = "l", col = "green")
lines(y.obs.v, col = "black")

# lines(y.hat, col = "dodgerblue4")
#------------------------------------------------------------------------------------------------------#





#==================================================================================================#


#--- Funcion: data.scm() -----------------------------------------------------------------------------#
# data.scm() (1)generate the dataframe which is simulated random walk
#            (2)calculate the difference: alpha = y.Hat - y.Obs
#               (2.1) y.Hat ==> Synthetic value (weighted average) and
#               (2.2) y.Obs ==> True value (follow its own time series value)
#
# return: (1) data.scm as a dataframe


data.scm <- function(timeLength , intervention_Time , controlNumber,
                     tau, omiga , alpha, seed.control) {
            
            # 1. simulate the data
            data <- dataRW.sim(controlNumber = controlNumber, timeLength = timeLength , 
                               tau = tau, omiga = omiga, seed.control = seed.control)
            
            y.Hat <- apply(data, 1, mean)
            
            y.Obs <- y.obs(y.hat = y.Hat, intervention_Time = intervention_Time, 
                           alpha = alpha, tau = tau )
            
            # 2.Merge as the dataframe
            data.scm.frame <- data.frame("Time" = 1:timeLength, data, "Treat.Hat" = y.Hat, "Treat.Obs" = y.Obs) 
            
            return(data.scm.frame)
            
}


##test
# 
# controlNumber = 100; timeLength = 20; intervention_Time = 15;
# tau = 0.1; omiga = 0.5; seed.control = 1; alpha = 1
# 
# 
# data_scm <- data.scm(controlNumber = controlNumber, timeLength = timeLength , intervention_Time = intervention_Time,
#          tau = tau, omiga = omiga, seed.control = seed.control, alpha = alpha
#          )

#--------------------------------------------------------------------------------------------------------#







# script SCM_PlotTrendFunction.R
#
# the purpose of this script is to  plot  how the control units and treated unit varied
#
# first version: 160404
# this version:  160404
# last change by: Yujiao Li

#==================================================================================================#

# raw.data is generated from data.scm()

plot.scm <- function(raw.data, intervention_Time){   
            data.plot <- gather(raw.data, Unit, Value,-Time)
            
            plot <- ggplot(data.plot, aes(x = Time, y = Value, colour = Unit,group = Unit)) +
                        geom_line(lwd = 1, alpha = 0.3,linetype = 1) +
                        
                        
                        
                        geom_line(
                                    data = data.plot[data.plot$Unit == "Treat.Obs",],
                                    aes(x = Time, y = Value), colour = "black",
                                    lwd = 1.5, linetype = "solid" ) +
                        
                        geom_line(
                                    data = data.plot[data.plot$Unit == "Treat.Hat",],
                                    aes(x = Time, y = Value), colour = "royalblue3",
                                    lwd = 1.5,linetype = "longdash" ) +
                        
                        geom_line(
                                    data = data.plot[data.plot$Unit == "Treat.Obs.random",],
                                    aes(x = Time, y = Value), colour = "sienna1",
                                    lwd = 2.2,linetype = "dashed" ) +
                        
                        
                        
                        
                        geom_point(
                                    data = data.plot[data.plot$Unit == "Treat.Hat",],
                                    aes(x = Time, y = Value), colour = "blue",size = 2 ) +
                        geom_point(
                                    data = data.plot[data.plot$Unit == "Treat.Obs",],
                                    aes(x = Time, y = Value), colour = "black",size = 2 ) +
                        geom_point(
                                    data = data.plot[data.plot$Unit == "Treat.Obs.random",],
                                    aes(x = Time, y = Value), colour = "darkorange",size = 2 ) +
                        
                        
                        geom_vline(xintercept = intervention_Time,linetype = "longdash") +
                        
                        ggtitle("Synthetic Control Method\nSimulation") +
                        
                        annotate(
                                    "text", label = paste0("Intervention \n  =  ",intervention_Time),
                                    x = intervention_Time, y = (0.05 + min(data.plot$Value)),
                                    size = 4, colour = "red") +
                        theme(legend.position = "none")
            
            
            if (  "Treat.Hat_SCM" %in%  names(raw.data)  ) {  
                        
                        plot <- plot + geom_line(
                                    data = data.plot[data.plot$Unit == "Treat.Hat_SCM",],
                                    aes(x = Time, y = Value), colour = "darkred",
                                    lwd = 1.8, linetype = "dashed" ) +
                                    ggtitle("Sample in Synthetic Control Method\nSimulation")
            }
            
            
            
            if (  "Treat.Hat_SCM.rand" %in%  names(raw.data)  ) {  
                        
                        plot <- plot + geom_line(
                                    data = data.plot[data.plot$Unit == "Treat.Hat_SCM.rand",],
                                    aes(x = Time, y = Value), colour = "deeppink",
                                    lwd = 1.5, linetype = "solid" ) +
                                    ggtitle("Sample in Synthetic Control Method\nSimulation")
            }
            
            return(plot)
            
}


plot.scm(data_scm, intervention_Time = intervention_Time)


# Function of Sample the data
sampleData <- function(data.scm , seed.sample, sample.percent) {
            set.seed(seed.sample)
            
            control.col <- grep(pattern = "Control", x = colnames(data.scm) , value = FALSE)
            
            sample.size <- (length(control.col) )* sample.percent
            
            sample.r0 <- sample(control.col, size = sample.size)
            
            
            basic.col <- which( !(1:ncol(data.scm))  %in% control.col  )
            
            
            
            data <- data.scm[ , sort(c(basic.col, sample.r0)) ]
            
            
            return(data)
}

#Test
# sample.data <- sampleData(data.scm = data_scm , seed.sample = 10, sample.percent = 0.5)


#==================================================================================================#

## Function of weight paramter estimation by quadratic programming solver

weight.hat <- function(sample.data, intervention_Time, Y.name = NULL,
                       margin.ipop = 0.0005,
                       sigf.ipop = 5,
                       bound.ipop = 10,
                       quadopt = "ipop"
) {
            
            # convert to matrix
            
            if (is.null(Y.name)) {Y.name = "Treat.Obs"}
            
            
            control.col <- grep(pattern = "Control", colnames(sample.data) , value = FALSE)
            Z <- sample.data[,control.col] %>% as.matrix #control
            
            treat.col <- which(colnames(sample.data) == Y.name)
            Y.obs <- sample.data[,treat.col]  %>% as.matrix    #treat
            
            
            
            
            #Z_matrix and Y_obs are used preInterventio value used for estimation 
            Z_matrix <- Z[1:intervention_Time,]
            Y_obs <- Y.obs[1:intervention_Time,] 
            
            
            # set up QP problem
            H <- t(Z_matrix) %*%  (Z_matrix) 
            a <- Y_obs
            c <- -1*c(t(a) %*%  (Z_matrix) )
            A <- t(rep(1, length(c)))
            b <- 1
            l <- rep(0, length(c))
            u <- rep(1, length(c))
            r <- 0
            
            # run QP and obtain w weights
            # ipop
            if (quadopt == "ipop") {
                        
                        res <- ipop(c = c, H = H, A = A, b = b, l = l, u = u, r = r, bound = bound.ipop,
                                    margin = margin.ipop, maxiter = 1000, sigf = sigf.ipop)
                        
                        solution.w <- as.matrix(primal(res))
                        
                        
            } else {
                        
                        # LowRankQP
                        if (quadopt == "LowRankQP") {
                                    
                                    res <- LowRankQP(Vmat = H,dvec = c,Amat = A,
                                                     bvec = 1,uvec = rep(1,length(c)),method = "LU")
                                    solution.w <- as.matrix(res$alpha)
                        } 
            }
            
            
            
            # results    
            
            MSE_pre.int <- as.numeric(t(Y_obs - Z_matrix %*% solution.w) %*% (Y_obs - Z_matrix %*% solution.w )) /nrow(Z_matrix)
            Y.predict <- Z %*% solution.w
            
            
            t.1 <- intervention_Time + 1
            Y_predict.1 <- as.vector(Z %*% solution.w)[t.1] 
            Y_observe.1 <- as.vector(sample.data$Treat.Obs)[t.1] 
            alphaHat.1 <- (Y_observe.1 - Y_predict.1)
            
            
            return(list("w_scm.hat" = solution.w, "MSE_scm_before" = MSE_pre.int, 
                        "Y_scm.hat" = Y.predict, "alpha1_scm.hat" = alphaHat.1 
            )
            )
}


# 
# weight.hat(sample.data, intervention_Time, Y.name = "Treat.Hat",
#                        margin.ipop = 0.0005,
#                        sigf.ipop = 5,
#                        bound.ipop = 10,
#                        quadopt = "ipop"
# )

#==================================================================================================#

## Function of integrate all the results from scm estimation

SCM.sample <- function(data_scm, seed.sample, sample.percent,  intervention_Time){
            
            sample.data <- sampleData(data_scm , seed.sample, sample.percent) 
            
            scm.result <- weight.hat(sample.data, intervention_Time,
                                     margin.ipop = 0.0005,
                                     sigf.ipop = 5,
                                     bound.ipop = 10,
                                     quadopt = "ipop") 
            
            w.Hat_SCM <- scm.result$w_scm.hat
            MSE.Pretime_SCM <- scm.result$MSE_scm_before
            Treat.Hat_SCM <- scm.result$Y_scm.hat 
            alpha1.Hat_SCM <- scm.result$alpha1_scm.hat 
            
            
            sample.data_SCM <- data.frame(sample.data,Treat.Hat_SCM)
            
            return(list(
                        "w.Hat_SCM" = w.Hat_SCM,
                        "MSE.Pretime_SCM" = MSE.Pretime_SCM,
                        "Treat.Hat_SCM" = Treat.Hat_SCM,
                        "alpha1.Hat_SCM" = as.numeric(alpha1.Hat_SCM),
                        "sample.data_SCM" = sample.data_SCM
            ))
            
            
            
}


#==================================================================================================#


# script SCM_ReplicationFunction.R
#
# the purpose of this script is to replicate the estimation of alpha for inference, 
# but every replication is in the same population pool
#
# first version: 160404
# this version:  160404
# last change by: Yujiao Li


rep.scm  <- function(rep.times, seed.replicate,
                     data_scm, 
                     sample.percent,  intervention_Time
) {
            set.seed(seed.replicate)
            seed.sample <- sample(1:1000, rep.times, replace = TRUE)
            
            rep.alpha <- c()
            for (i in 1:rep.times) {
                        seed.sample.i <- seed.sample[i]
                        result.SCM <-
                                    SCM.sample(data_scm, seed.sample.i, sample.percent,  intervention_Time)
                        rep.alpha[i] <- result.SCM$alpha1.Hat_SCM
            }
            
            return(rep.alpha)
            
}





#==================================================================================================#
# alpha_raw is the function considered all of the parameters
# alpha_raw is the integrate function to allow you set any parameter
alpha_raw <- function(seed.control, seed.subsample, 
                      timeLength ,intervention_Time,controlNumber,
                      tau,omiga,alpha,
                      seed.sample,sample.percent){  
            
            
            
            data_scm <- data.scm(
                        timeLength = timeLength,
                        intervention_Time = intervention_Time,
                        controlNumber = controlNumber,
                        
                        tau    = tau,
                        omiga  = omiga,
                        alpha  = alpha,
                        seed.control = seed.control )
            
            seed.sample = seed.sample
            sample.percent = sample.percent
            
            result.SCM <- SCM.sample(data_scm, seed.sample, sample.percent,  intervention_Time)
            alpha_Hat.current <- as.numeric(result.SCM$alpha1.Hat_SCM)
            
            return(alpha_Hat.current)
}


# 
# alpha_raw(
#             seed.control = 3,
#             seed.sample = 10,
#             sample.percent = 0.8,
#             
#             timeLength        = 30,
#             intervention_Time = 20,
#             controlNumber     = 20,
#             
#             
#             tau    = 0.3,
#             omiga  = 0,
#             alpha  = 0.5
# )





#==================================================================================================#

# the purpose of this script is to replicate the estimation of alpha for inference, 
# but every replication is in the different population pool

#repli_alpha_raw  is the function to replicate the simulation of alpha many times
# where the alpha can allow any paramter changed 

repli_alpha_raw <- function(time,  seed.control.rep = 100, 
                            timeLength ,intervention_Time,controlNumber,
                            tau,omiga,alpha,
                            seed.sample,sample.percent){
            
            
            seed.control.rep = seed.control.rep
            set.seed(seed.control.rep) 
            # surely we can change the seed for repli, but for convenience in the app, we set the seed is fixed
            
            seed.control.raw <- sample(1:10000,time)
            replicateRawAlpha <- sapply(seed.control.raw ,   FUN = alpha_raw,     
                                        seed.sample = seed.sample,
                                        sample.percent = sample.percent ,
                                        
                                        timeLength        = timeLength ,
                                        intervention_Time = intervention_Time,
                                        controlNumber     = controlNumber,
                                        
                                        
                                        tau    = tau ,
                                        omiga  = omiga,
                                        alpha  = alpha)
            
            return(as.vector(replicateRawAlpha))
            
}


## Test
# repli_alpha_raw(time = 100,
#             #seed.control.rep=10, # this seed is default as 10, we can change surely
#             seed.sample = 10,
#             sample.percent = 0.8,
#             
#             timeLength        = 30,
#             intervention_Time = 20,
#             controlNumber     = 20,
#             
#             
#             tau    = 0.3,
#             omiga  = 0,
#             alpha  = 0.5
# )

#################################
################################
# random module

###
SCM.sample.random <- function(data.raw, rou, timeLength , seed.obs.random){ 
            
            controlNumber <- length(grep(pattern = "Control", x = colnames(data.raw) , value = FALSE))
            y.obs.v <- data.raw$Treat.Obs
            Treat.Obs.random <- y.obs.random(rou = rou, timeLength = timeLength , tau = tau, 
                                             seed.obs.random = seed.obs.random,
                                             controlNumber = controlNumber , y.obs.v = y.obs.v)
            data.random <- data.frame(data.raw, Treat.Obs.random)
            #plot.scm(data.random, intervention_Time = intervention_Time)
            
            result <- weight.hat(data.random, intervention_Time, Y.name = "Treat.Obs.random",
                                 margin.ipop = 0.0005,
                                 sigf.ipop = 5,
                                 bound.ipop = 10,
                                 quadopt = "ipop")
            
            return(list( "data.random" = data.random, "result" = result))
}


## test
# resu <- SCM.sample.random(data.raw = sample.data, rou = 0.1, timeLength = timeLength, controlNumber = controlNumber)
# plot.scm(resu$data.random, intervention_Time = intervention_Time)

final.plotdata <- function(
            
            raw.data,
            seed.sample,sample.percent,intervention_Time,
            seed.obs.random = seed.obs.random,
            rou,timeLength
            
){   
            samdata <- SCM.sample(raw.data, seed.sample = seed.sample,
                                  sample.percent = sample.percent,  
                                  intervention_Time = intervention_Time)$sample.data_SCM
            result <- SCM.sample.random(
                        data.raw = raw.data , rou = rou, 
                        seed.obs.random = seed.obs.random,
                        timeLength = timeLength)
            
            alpha.hat.random <- result$result$alpha1_scm.hat
            
            data <- data.frame(samdata, "Treat.Obs.random" = result$data.random$Treat.Obs.random)
            final <- data.frame(data, 
                                "Treat.Hat_SCM.rand" = SCM.sample.random(data.raw = data , rou = rou,  
                                                                         seed.obs.random = seed.obs.random,
                                                                         timeLength = timeLength)$result$Y_scm.hat )
            return(list("final.data" = final , "alpha.hat.random" = alpha.hat.random))
            
}



# intervention_Time = 15;timeLength = 20; rou =0.2;
# tau = 0.1; row = 0.4; seed.sample = 10; sample.percent = 0.5
# 
# 
# controlNumber = 100; 
# omiga = 0.5; seed.control = 1; alpha = 1
# 
# seed.obs.random = 21
# 
# data_scm <- data.scm(controlNumber = controlNumber, timeLength = timeLength , intervention_Time = intervention_Time,
#                      tau = tau, omiga = omiga, seed.control = seed.control, alpha = alpha
# )
# 
# result <- final.plotdata(
#             raw.data = data_scm,
#             seed.sample = seed.sample,
#             sample.percent = sample.percent,
#             seed.obs.random = seed.obs.random,
#             intervention_Time = intervention_Time,
#             rou = rou,
#             timeLength = timeLength)
# 
# plotfinal <- result$final.data
# result$alpha.hat.random

# plot.scm(plotfinal, intervention_Time = intervention_Time)



