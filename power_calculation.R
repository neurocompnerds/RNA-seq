library(PROPER)

# Calculating power for an RNA-seq experiment using LCLs
# The determination of sample sizes is important, however challenges include: 
# Multiple testing problem: Many genes are tested for DE simultaneously
# Coverage depth 
# DE detection procedure --> Different shrinkage procedures are used to stabilize the estimation of within group variation 


### setting up a simulation scenario 
## The Cheung data uses expression of 41 LCLs from unrelated individuals, so the expressions show a large biological variation 

sim.opts.Cheung = RNAseq.SimOptions.2grp(ngenes = 16000, p.DE=0.02,lOD="cheung", lBaselineExpr="cheung")

## The dispersions in this data set wil be very large
## Nreps is number of reps in group 1 --> in this case original control samples, orginal control samples with gtex 
## Nreps2 is number of reps in group 2: in this case, the condition group

simres = runSims(Nreps = c(73,73,7,7),Nreps2 = c(1,3,3,1), sim.opts=sim.opts.Cheung,DEmethod="edgeR", nsims=30)

# Comparing power 
# Method to control for type I error (FDR or raw pvalue)
# Method to stratify genes 
# Method to identify interesting
powers = comparePower(simres, alpha.type="fdr", alpha.nominal=0.1,stratify.by="expr", delta=0.5)

summaryPower(powers)
# Summary table 
# Marginal power 
# TD --> True discovery 
# FD --> False discovery 
# FDC --> False discovery cost defined as the number of FD / TD 

# plot stratified powers 
plotPower(powers)


# plot number of true discoveries 
plotPowerTD(powers)

# plot stratified false discovery cost 
plotFDcost(powers)

# Plot everything in a single finger 
plotAll(powers)

# Power and seq depth 
power.seqDepth(simres,powers)
