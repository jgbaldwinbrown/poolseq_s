#!/usr/bin/env Rscript

suppressMessages(library(poolSeq))

print(1)
mySync <- read.sync(file="ex7.sync", gen=c(0, 10, 20, 0, 10, 20), repl=c(1, 1, 1, 2, 2, 2))
print(2)
print(mySync)
mytraj = af.traj(mySync, c("2R", "3R"), c(2304, 2305), repl=c(1,2))
mytraj1 = af.traj(mySync, c("2R", "3R"), c(2304, 2305), repl=c(1))
print(3)
print(mytraj)
print(4)
mycov = coverage(mySync, c("2R", "3R"), c(2304, 2305), repl=c(1,2), gen=c(0,10,20))
mycov1 = coverage(mySync, c("2R", "3R"), c(2304, 2305), repl=c(1), gen=c(0,10,20))
print(5)
print(mycov)
print(6)

estimateNe(p0=mytraj1[,"F0"], pt=mytraj1[,"F10"], cov0=mycov1[,"F0.R1.cov"], covt=mycov1[,"F10.R1.cov"], t=10, Ncensus=1000, poolSize=c(300, 300))
# estimateNe(p0=mytraj[,"F0"], pt=mytraj[,"F20"], cov0=myCov[,"F0"], covt=myCov[,"F20"], t=20, Ncensus=1000, poolSize=c(300, 300))
print(7)

s = estimateSH(mytraj, Ne=1000, t=c(0,10, 20), h=.5, simulate.p.value=TRUE)
print(8)
print(s)
print(9)
confint(s)
print(10)

    
#     simTraj <- wf.traj(p0=0.05, Ne=1000, t=seq(0, 60, by=10), s=0.1, h=0.5)
#     est <- estimateSH(simTraj, Ne=1000, t=seq(0, 60, by=10), h=0.5)
#     estp <- estimateSH(simTraj, Ne=1000, t=seq(0, 60, by=10), h=0.5, simulate.p.value=TRUE)
#     print("simTraj")
#     print(simTraj)
#     print(str(simTraj))
#     print("est")
#     print(est)
#     print("estp")
#     print(estp)
#     print("confint(est)")
#     print(confint(est))
# }

# main()
