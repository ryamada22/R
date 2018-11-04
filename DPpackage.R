# Data
      data(rolling)
      y <- cbind(rolling$y1,rolling$y2)


    # Prior information

      prior<-list(a0=1,
                  b0=1,
                  a1=1,
                  b1=1)

    # Initial state
      state <- NULL

    # MCMC parameters

      mcmc <- list(nburn=500,
                   nsave=1000,
                   nskip=3,
                   ndisplay=100)

    # Fitting the model

      fit <- DPbetabinom(y=y,ngrid=100, 
                         prior=prior, 
                         mcmc=mcmc, 
                         state=state, 
                         status=TRUE)

      fit
      summary(fit)

    # density estimate
      plot(fit,output="density")

    # parameters
      plot(fit,output="param")
