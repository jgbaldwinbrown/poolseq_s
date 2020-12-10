#!/usr/bin/env Rscript

suppressMessages(library(poolSeq))

main = function() {
    
    simTraj <- wf.traj(p0=0.05, Ne=1000, t=seq(0, 60, by=10), s=0.1, h=0.5)
    est <- estimateSH(simTraj, Ne=1000, t=seq(0, 60, by=10), h=0.5)
    estp <- estimateSH(simTraj, Ne=1000, t=seq(0, 60, by=10), h=0.5, simulate.p.value=TRUE)
    print("simTraj")
    print(simTraj)
    print(str(simTraj))
    print("est")
    print(est)
    print("estp")
    print(estp)
    print("confint(est)")
    print(confint(est))
}

main()
