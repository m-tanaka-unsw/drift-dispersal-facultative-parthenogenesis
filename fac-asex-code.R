##==========================================================
##
## R script for figures in the paper: 
## 
##  Drift and dispersal hinder the evolution of facultative
##   asexual reproduction
## 
##  by Mark Tanaka and Russell Bonduriansky 
##  in the journal Evolution 
## 
## Required packages: sfsmisc, parallel
## 
##==========================================================

## Set default parameters of model 
setp <- function(N=1000,     # Pop size
                 D=3,        # Dispersal parameter 
                 Tend=3000,  # when to end simulation
                 tol=10^-9,  # tolerance for stopping deterministic process
                 init.type=1,   # which initial condition
                 sims=12000  # num simulations to estimate invasion prob
                 ){
   list(N=N, D=D, init.type=init.type, Tend=Tend, tol=tol, sims=sims)
}


## Set initial conditions based on parameters 
##  The order of the elements: x, z, u, y, v 
setic <- function(p){
   with(p, {
      if (init.type==1) # choice 1: start pop near equal sex ratio
         ## oblig females and males near/at 1/2, heterozygotes at a small fraction
         ic <- list(x=1/2-1/N, z=1/N, u=0, y=1/2, v=0)
      else ## start with many copies 
         ic <- list(x=1/2-10/N, z=10/N, u=0, y=1/2, v=0)
      ic})
}

catn <- function(...) cat(...,"\n") 


##---------------------------------------------------------- 
## Model: deterministic X-linked FP allele
##        (male is hemi-/heterogametic)
## 
## STOPPING METHODS: 
##  1.  stop when UU freq gets to 1/N
##  2.  stop when UU > 1/N AND total amount of change is under tolerance limit
##  3+. don't stop; go to the end -- fixed number of generations, Tend
##---------------------------------------------------------- 
model.det <- function(p, recessive=TRUE, stop.method=1){
   ic <- setic(p)
   x <- ic$x
   z <- ic$z
   u <- ic$u
   y <- ic$y
   v <- ic$v
   traj <- as.matrix(t(unlist(ic)))
   for (i in 1:(p$Tend-1)) { # for-loop to cap the iterations
      ## Prob mating with y, v
      Py <- (1 - exp(-p$D*(y+v)))*y/(v+y)
      Pv <- (1 - exp(-p$D*(y+v)))*v/(v+y)
      
      xp <- Py*(x+z/2)/2
      if (recessive) 
         zp <- Py*(u+z/2)/2 + Pv*(x+z/2)/2
      else # DOMINANT
         zp <- Py*(u+z/2)/2 + Pv*(x+z/2)/2 + exp(-p$D*(y+v))*z
      up <- Pv*(u+z/2)/2 + exp(-p$D*(y+v))*u
      yp <- (Pv+Py)*(x+z/2)/2
      vp <- (Py+Pv)*(u+z/2)/2
      
      wbar <- xp+zp+up+yp+vp
      x <- xp / wbar    # normalise frequency
      z <- zp / wbar
      u <- up / wbar
      y <- yp / wbar
      v <- vp / wbar
      
      traj <- rbind(traj,c(x,z,u,y,v))
      if (i>1) {  # consider stopping after first gen
         if (stop.method==1) { # has UU exceeded some frequency? e.g. 1/N
            if (traj[i,3] > 1/p$N) {
               break
            }
         } else if (stop.method==2) { # look at change rate
            totaldiff <- sum(abs(traj[i-1,] - c(x,z,u,y,v)))
            if (totaldiff < p$tol && traj[i,3] > 1/p$N) {
               break
            }
         }
         else {
            ## Don't stop -- go to the end
         }
      }
   }
   traj
}


##---------------------------------------------------------- 
## Wright-Fisher version of model 
##---------------------------------------------------------- 
model.stoch.recessive <- function(p, earlystop=TRUE){
   N <- p$N 
   ic <- setic(p)
   x <- ic$x
   z <- ic$z
   u <- ic$u
   y <- ic$y
   v <- ic$v
   traj <- as.matrix(t(c(x,z,u,y,v)))
   for (i in 1:(p$Tend-1)) {
      if (v+y==0) { # special boundary: all males extinct
         Py <- 0
         Pv <- 0
      } else {
         Py <- (1 - exp(-p$D*(y+v)))*y/(v+y)
         Pv <- (1 - exp(-p$D*(y+v)))*v/(v+y)
      }
      xp <- Py*(x+z/2)/2
      zp <- Py*(u+z/2)/2 + Pv*(x+z/2)/2
      up <- Pv*(u+z/2)/2 + exp(-p$D*(y+v))*u
      yp <- (Pv+Py)*(x+z/2)/2
      vp <- (Py+Pv)*(u+z/2)/2

      if (xp+zp+up+yp+vp==0) {
         ## This shouldn't happen, but just in case it does
         stop("All frequencies are zero")
      }
      ns <- t(rmultinom(1,size=N, prob=c(xp,zp,up,yp,vp)))
      
      x <- ns[1]/N
      z <- ns[2]/N
      u <- ns[3]/N
      y <- ns[4]/N
      v <- ns[5]/N
      
      traj <- rbind(traj,ns/N)
      if (earlystop) { # stop when things go extinct
         if (x==1 || y==1 || # OS females and males can go to 100% in small pops
             z+u+v==0 || u==1 || v==1)
            break
      }
   }
   traj   
}

##----------------------------------------------------------
## Deterministic only; show effect of increasing D 
##----------------------------------------------------------
plot_dynamics_altD <- function(){
   figpanel <- function(p,dr,fem=TRUE,
                        panel="A.",maintext="D=1"){
      ts <- c(1:p$Tend)  # generations
      if (fem) { # plot female freqs 
         matplot(ts,dr[,1:3], t="l", log="", ylim=c(10^-3,1),
                 lty=1,lwd=1, col=c(1:3), 
                 main=maintext, cex.main=1
               , ylab="Frequency",xlab="" 
                 )
         mtext(panel, 3, adj=-0.25, line=0, cex=0.9) # label panel with letter
         legend("topleft",c("XX","XU","UU"),lty=1, lwd=1, col=c(1:3), cex=0.9)
      } else {  # plot male freqs
         matplot(ts,dr[,4:5], t="l", log="", ylim=c(10^-3,1),
                 lty=1,lwd=1, col=c(4,6), 
                 ylab="Frequency",xlab="", 
                 main=maintext, cex.main=1
                 )
         mtext(panel, 3,  adj=-0.25, line=0, cex=0.9)
         mtext("Time (generations)", 1, line=2.5, cex=0.7)
         legend("topleft",c("XY","UY"),lty=1,lwd=1, col=c(4,6), cex=0.9)
      }
   }
   
   pdf("fig-dyn-alt-D.pdf",width=5,height=4)
   par(mfcol=c(2,2), cex=0.7, 
       mar=c(2.1,   4,   1.4, 0.5),
       oma=c(2.1, 0.1,   0.1, 0.1))
   
   p <- setp(D=1, N=1000, Tend=30, init.type=1)   
   dr <- model.det(p, stop.method=3)
   figpanel(p,dr)
   figpanel(p,dr, fem=FALSE, panel="C.")
   
   p <- setp(D=3, N=1000, Tend=2000, init.type=1)
   dr <- model.det(p, stop.method=3)
   figpanel(p,dr, panel="B.", maintext="D=3")
   figpanel(p,dr, fem=FALSE, panel="D.", maintext="D=3")
   dev.off()   
}


##----------------------------------------------------------
## Plot dynamics for recessive X-linked FP model
##  BOTH deterministic and stochastic -- recessive model 
##----------------------------------------------------------
plot_dynamics_det_stoch <- function(s1=4,s2=2){
   require("sfsmisc")
   p <- setp(D=2, N=1000, Tend=60, init.type=2)
   
   dr <- model.det(p, stop.method=3)  # recessive model 
   set.seed(s1)
   sr <- model.stoch.recessive(p, earlystop=F)
   ts <- c(1:p$Tend)  # generations
   catn("Final UU: ",dr[p$Tend,3])

   ##  plot with success/failure  
   pdf("fig-dynamics-recessive.pdf", width=6, height=4.2)
   par(mfcol=c(2,2), mar=c(2,4, 1, 0.5), oma=c(2.1,0.1, 0, 0),
       cex=0.7)
   mtc <- 0.75
   ##--- First sim: stoch tracks deterministic ---
   ## A: XX, XU, UU
   matplot(ts,dr[,1:3], t="l", log="y", ylim=c(10^-3,1),
           lty=2,lwd=2, col=c(1:3), yaxt="n",
           main="Simulation 1 - females", cex.main=0.9
           , ylab="Frequency",xlab="" #"Time (generations)"
           )
   sfsmisc::eaxis(2, n.axp=1)   
   matplot(ts,sr[,1:3], t="l", add=T,
           lty=1,lwd=1,col=c(1:3))
   legend("bottomright",c("XX","XU","UU"),lty=1,lwd=2, col=c(1:3), cex=0.9)
   mtext("A.", 3,  adj=-0.25, line=0, cex=mtc)
   ## C: XY, UY 
   matplot(ts,dr[,4:5], t="l", log="y", ylim=c(10^-3,1),
           lty=2,lwd=2, col=c(4,6), yaxt="n",
           ylab="Frequency",xlab="", #"Time (generations)"
           main="Simulation 1 - males", cex.main=0.9
           )
   sfsmisc::eaxis(2, n.axp=1)
   matplot(ts,sr[,4:5], t="l", add=T,
           lty=1,lwd=1,col=c(4,6))
   legend("bottomright",c("XY","UY"),lty=1,lwd=2, col=c(4,6), cex=0.9)
   mtext("C.", 3,  adj=-0.25, line=0, cex=mtc)
   mtext("Time (generations)", 1, line=3.0, cex=mtc)
   ## --- Second sim: invasion fails 
   ## Do sim again
   set.seed(s2)
   sr <- model.stoch.recessive(p, earlystop=F)
   ## B: XX, XU, UU
   matplot(ts,dr[,1:3], t="l", log="y", ylim=c(10^-3,1),
           lty=2,lwd=2, col=c(1:3), yaxt="n",
           ylab="" 
         , xlab="" 
         , main="Simulation 2 - females", cex.main=0.9
           )
   sfsmisc::eaxis(2, n.axp=1)   
   matplot(ts,sr[,1:3], t="l", add=T,
           lty=1,lwd=1,col=c(1:3))
   legend("bottomright",c("XX","XU","UU"),lty=1,lwd=2, col=c(1:3), cex=0.9)
   mtext("B.", 3,  adj=-0.25, line=0, cex=mtc)
   ## D: XY, UY 
   matplot(ts,dr[,4:5], t="l", log="y", ylim=c(10^-3,1),
           lty=2,lwd=2, col=c(4,6), yaxt="n",
           ylab="" 
         , xlab="" 
         , main="Simulation 2 - males", cex.main=0.9
           )
   sfsmisc::eaxis(2, n.axp=1)
   matplot(ts,sr[,4:5], t="l", add=T,
           lty=1,lwd=1,col=c(4,6))
   legend("bottomright",c("XY","UY"),lty=1,lwd=2, col=c(4,6), cex=0.9)
   mtext("D.", 3,  adj=-0.25, line=0, cex=mtc)
   mtext("Time (generations)", 1, line=3.0, cex=mtc)
   dev.off()
}


##----------------------------------------------------------
##  Illustrate analytical results 
##----------------------------------------------------------

## Numerical solution to v to obtain v-star
## (the larger of the two solutions and only positive ones) 
vstar <- function(D){
   tiny <- 10^-5
   if (D>2) {
      fn <- function(x,D)
         1-exp(-D*x) - 2*x
      u <- uniroot(fn, c(tiny,1), D)
      return(u$root)
   } else {
      return(0)
   }
}

## Show how D affects: 
##   A. long term UU and
##   B. time to get to 1/N deterministically
plot_D_effects <- function(){
   require("parallel")
   final.state <- function(D){
      dr <- model.det(setp(D=D, Tend=20000), stop.method=2)
      last <- length(dr[,1])
      dr[last,3]  # end state of UU freq
   }
   final.state.stoch <- function(D){
      sr <- model.stoch.recessive(setp(D=D, N=1000, Tend=10000))
      last <- length(sr[,1])
      finalU <- sr[last,3]  # end state of UU freq
      if (finalU > 0)       # Condition on UU invading 
         return(finalU)
      else
         final.state.stoch(D)
   }
   stop.time <- function(D, rec=T){
      ## Under stop.method=1, process stops when 1/N reached
      dr <- model.det(setp(D=D, Tend=20000), recessive=rec, stop.method=1)
      length(dr[,1])  # get the end time 
   }
   nc <- detectCores()
   Ds <- seq(0.2, 5, length.out=50)
   vs <- unlist(mclapply(Ds, vstar, mc.cores=nc))
   Us <- unlist(mclapply(Ds, final.state, mc.cores=nc))
   STR <- unlist(mclapply(Ds, stop.time, mc.cores=nc))  # recessive
   STD <- unlist(mclapply(Ds, stop.time, FALSE, mc.cores=nc)) # dominant
   Ss <- unlist(mclapply(Ds, final.state.stoch, mc.cores=nc))
   catn("Final U")
   print(Us)
   print(Ss)
   catn("End time")
   print(STR)
   print(STD)
   
   pdf("fig-D-effects.pdf", height=4.5, width=3)
   par(mfrow=c(2,1), mar=c(3,4, 2, 0.5), cex=0.8)
   plot(Ds,Us,  t="l"
      , xlab="" 
      , ylab="Long term freq of UU"
      , ylim=c(0,1), lty=1)
   mtext("A.", side=3, line=0.5, adj=0)
   mtext("Dispersal parameter D", side=1, line=2, cex=0.8)
   rd <- adjustcolor("red", alpha.f=0.5)
   bl <- adjustcolor("blue", alpha.f=0.25)
   lines(Ds,1-vs, col=rd, lty=2, lwd=3)
   lines(Ds, Ss, col=bl, lwd=3, lty=1)
   legend("bottomleft", c("Theoretical","Deterministic","Stochastic")
          , lty=c(2,1,1), lwd=c(2,1,2), col=c("red",1,bl), cex=0.8)
   plot(Ds,STR, log="y", t="l"
      , xlab="" 
      , ylab="Time for UU freq to reach 0.001")
   mtext("Dispersal parameter D", side=1, line=2, cex=0.8)
   mtext("B.", side=3, line=0.5, adj=0)
   lines(Ds,STD, lty=2, lwd=2, col="darkorange1")
   legend("topleft",c("Recessive","Dominant"), lty=c(1,2), lwd=c(1,2),
          col=c(1,"darkorange1")
          , cex=0.8)
   dev.off()
}


##----------------------------------------------------------
##  Multiple simulations to get invasion probabilities
##----------------------------------------------------------

## Parallelise for a single set of parameters
parallel.sims.single.set <- function(p){
   simulate <- function(x,p){ # set up function to simulate in parallel 
      ssto <- model.stoch.recessive(p)  # run stochastic model
      last.t <- dim(ssto)[1] # actual end can be <Tend if stopped early
      ssto[last.t,3]*p$N          # use absolute count of u 
   }   
   ncore <- parallel::detectCores()
   x <- c(1:p$sims)
   pendu <- parallel::mclapply(x,simulate,p, mc.cores=ncore)
   sum(pendu>0)/p$sims # proportion with last UU greater than zero
}

## Either run sims or plot sim results 
## Do the sims over N and D, then put data into a file
pinv_sims_over_N_D <- function(p
                        , DD=c(0.5, 2, 5, 10)
                        , NN=c(100, 500, 1000, 5000,10000)
                        , runsims=FALSE ) {
   if (runsims) {  # run the simulations 
      p2 <- p
      PI <- NULL
      print(NN)
      for (k in (1:length(DD))) {
         set.seed(k) 
         p2$D <- DD[k]
         catn(" D is now" , p2$D)   
         pinv <- c()
         for (i in 1:length(NN)) {
            p2$N <- NN[i] 
            catn("N is now", p2$N) 
            pinv <- c(pinv, parallel.sims.single.set(p2))
         }
         PI <- rbind(PI, pinv)
      }
      rownames(PI) <- NULL  # remove rowname "pinv"
      sink("pinv-dat")
      catn("# ", DD)
      catn("# ", NN)
      print(PI)
      sink()
   } else {  # read in sim data and plot results 
      d <- read.table("pinv-dat")
      print(d)
      pdf("fig-prob-invasion-recessive.pdf", width=3,height=3)
      par(mfrow=c(1,1), mar=c(4,4,1,1), cex=0.8)
      cola <- c(2,3,4,5)    # give each D a colour 
      plot(NN, d[1,],  t="l", col=cola[1], log="xy"
           , ylim=c(5e-4, 0.5), yaxt="n"
           , xlab="Population size N"
           , ylab="Probability of invasion")
      points(NN, d[1,], pch=16, col=cola[1])
      lines(NN, d[2,], col=cola[2])
      points(NN, d[2,], pch=16, col=cola[2])
      lines(NN, d[3,], col=cola[3])
      points(NN, d[3,], pch=16, col=cola[3])
      lines(NN, d[4,], col=cola[4])
      points(NN, d[4,], pch=16, col=cola[4])

      sfsmisc::eaxis(2,  n.axp=1)
      legend("topright",
             c(paste("D =",DD[1]), paste("D =",DD[2]), paste("D =",DD[3]),
               paste("D =",DD[4])), 
             lty=1, col=cola, pch=16, cex=0.9)
      dev.off()
   }
}


##----------------------------------------------------------
## Run functions to generate figures 
##----------------------------------------------------------
generate_figures <- function(){
   plot_dynamics_altD()  # figure 1
   plot_D_effects()      # figure 2 
   plot_dynamics_det_stoch()   # figure 3 
   pinv_sims_over_N_D(setp(), runsims=T)  # run sims for fig 4
   pinv_sims_over_N_D(setp())  # plot results (figure 4) 
}
