##############################################################
# Unsupervised clustering and groupwise characterization of #
# galaxies from the Galaxy And Mass Assembly (GAMA) survey   #
##############################################################

###################
# Galaxy features # 
###################
## stellar mass: mass_stellar_best_fit (in log10 M⊙ - log mass in solar units)
## specific star formation rate: sSFR_0_1Gyr_best_fit (in log10 yr−1 - log SSFR per year)
## u − r colour: uminusr (in mags) 
## half-light radius: R_S_RE (in log10 kpc - log RE converted to kpc) 
## Sérsic index: R_S_NSER (in log10 n) 

#############################
# redshift: 0.05 < z < 0.08 #
#############################

##########################
# Environment parameters #
##########################
## SurfaceDensity (SD): Surface density 
## CountInCyl (CC): Cylindrical count  
## AGEDenPar (AGE): Adaptive Gaussian environment parameter

##################
# Galaxy dataset #
##################
mstar = read.table("gama-galaxies-7187.txt",
                   header = FALSE,
                   sep = "",
                   stringsAsFactors = FALSE)
X = as.matrix(mstar[,-1])

#####################################################################
# Figure: Distribution plots of galaxy features for oringal dataset #
#####################################################################
library(GGally)
type.names = c("Stellar mass", "Star formation rate", "u-r colour", 
               "Half-light radius", "Sersic index") 
dev.off()
pp = ggpairs(data.frame(X),
             upper = list(continuous = wrap("cor", size = 20, color = "brown")),
             diag = list(continuous = wrap("densityDiag", fill="#2171B5")),
             lower = list(continuous = wrap(ggally_points, size = 0.1, color = "#E7298A")),
             columnLabels = type.names)+ theme_bw() +
  scale_y_continuous(n.breaks = 3)+
  scale_x_continuous(n.breaks = 3)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 50),
        #axis.text.x = element_text(angle = 30,hjust = 1),
        strip.background = element_rect(fill="#DEEBF7"),
        strip.text.x = element_text(size = 40),
        strip.text.y = element_text(size = 40))
pp
ggsave(pp, file="~/Desktop/gama-feature-orgdata.png", 
       width=27, height=25)

#########################################
# Environment dataset:                  #
# - values not logged                   #
# - value of -999.9 means not measured  #
#########################################
env.mstar = read.table("gama-galaxies-enviro-7187.txt",
                       header = T,
                       sep = "",
                       stringsAsFactors = FALSE)

# Obtain galaxy samples with complete environment measures
index.flag0 = which((env.mstar$SurfaceDensity > -999.9) & 
                      (env.mstar$CountInCyl > -999.9) &
                      (env.mstar$AGEDenPar > -999.9)
)
mstar.complete = mstar[index.flag0,-1]
dim(mstar.complete)
summary(mstar.complete)

#####################################################################
# Determine the combined environment measure: Optimal density       #
# Reference: https://ui.adsabs.harvard.edu/abs/2023MNRAS.522.4116B/ #
#####################################################################
env.mstar.complete = env.mstar[index.flag0, 3:5]
env.mstar.complete = log10(env.mstar.complete)
log.cc = env.mstar.complete$CountInCyl
log.cc[log.cc<0] = -1
sum(log.cc<0) # 1637 objects with CC = 0 (lowest bin: log.cc <  0) 

# Compute red fraction at the lowest bin
u.r.vals = mstar.complete[which(log.cc<0),3]
red.frac.low = sum((1.5 < u.r.vals) & (u.r.vals< 4))/sum(log.cc<0)
breaks = c(seq(0,1,0.2),max(log.cc))
h <- hist(log.cc[log.cc>=0], breaks = breaks, plot = FALSE)

# Compute red fraction at the highest bin
u.r.vals = mstar.complete[which(log.cc>1),3]
red.frac.high = sum((1.5 < u.r.vals) & (u.r.vals< 4))/sum(log.cc>1)

# Compute the red fraction range for CC 
red.frac.high - red.frac.low

# Compute the red fraction range for SD
log.sigma = env.mstar.complete$SurfaceDensity
sort(log.sigma, decreasing = F)[1637]
breaks = c(seq(-0.6,0.8,0.2),max(log.sigma))
h <- hist(log.sigma[log.sigma>=(-0.6)], breaks = breaks, plot = FALSE)
u.r.vals = mstar.complete[which(log.sigma<(-0.6)),3]
red.frac.low = sum((1.5 < u.r.vals) & (u.r.vals< 4))/sum(log.sigma<(-0.6))
u.r.vals = mstar.complete[which(log.sigma>0.8),3]
red.frac.high = sum((1.5 < u.r.vals) & (u.r.vals< 4))/sum(log.sigma>0.8)
red.frac.high - red.frac.low

# Compute the red fraction range for AGE
log.age = env.mstar.complete$AGEDenPar
sort(log.age, decreasing = F)[1637]
breaks = c(seq(-0.5,0.5,0.2),max(log.age))
h <- hist(log.age[log.age>=(-0.5)], breaks = breaks, plot = FALSE)
u.r.vals = mstar.complete[which(log.age<(-0.5)),3]
red.frac.low = sum((1.5 < u.r.vals) & (u.r.vals< 4))/sum(log.age<(-0.5))
u.r.vals = mstar.complete[which(log.age>0.5),3]
red.frac.high = sum((1.5 < u.r.vals) & (u.r.vals< 4))/sum(log.age>0.5)
red.frac.high - red.frac.low

# Compute the optimal density - log.lambda: 
# - Formula: log.lambda = log.sigma + a*log.cc + b*log.age
# - Find optimal a and b that maximize the red fraction range for log.lambda
step = 0.01
a.seq = seq(0.01,10,step)
b.seq = seq(0.01,10,step)
red.frac.r = matrix(NA,length(a.seq),length(b.seq))
rownames(red.frac.r) = a.seq
colnames(red.frac.r) = b.seq

for (i in 1:length(a.seq)) {
  for (j in 1:length(b.seq)) {
    
    log.lambda = log.sigma + a.seq[i]*log.cc + b.seq[j]*log.age
    
    # Determine the lower bin
    bin.lower = round(sort(log.lambda, decreasing = F)[1640],digits = 1)
    # Compute red fraction at the lowest bin
    u.r.vals = mstar.complete[which(log.lambda < bin.lower),3]
    red.frac.low = sum((1.5 < u.r.vals) & (u.r.vals< 4))/sum(log.lambda < bin.lower)
    
    # Determine the upper bin
    max.bin.upper = bin.lower + 0.2*floor((max(log.lambda)-bin.lower)/0.2)
    breaks =c(seq(bin.lower,max.bin.upper,0.2),max(log.lambda))
    bin.ct = hist(log.lambda[log.lambda >= bin.lower], 
                  breaks = breaks, plot = FALSE)$counts
    # Find the index of the minimum difference
    bin.sum = sum(bin.ct)
    idx = which.min(abs(bin.sum - cumsum(bin.ct) - 436))
    bin.upper = bin.lower+0.2*idx
    # Compute red fraction at the highest bin
    u.r.vals = mstar.complete[which(log.lambda > bin.upper),3]
    red.frac.high = sum((1.5 < u.r.vals) & (u.r.vals< 3))/sum(log.lambda > bin.upper)
    
    # Compute red fraction range 
    red.frac.r[i,j] = abs(red.frac.high - red.frac.low)
  }
}
ind = which(red.frac.r==max(red.frac.r),arr.ind = T)
a.opt = rownames(red.frac.r)[ind[1]]  
b.opt = colnames(red.frac.r)[ind[2]] 
max(red.frac.r)

## Compute log.lambda based on optimal a and b
a.opt = as.numeric(a.opt)
b.opt = as.numeric(b.opt)
log.lambda = log.sigma + a.opt*log.cc + b.opt*log.age

# Create data matrix with the astrophysical features and the optimal density
X.lbd = cbind(mstar.complete,log.lambda)

#######################################################
# Model-based clustering with generalized t-mixtures  #
# Reference: https://arxiv.org/abs/2504.21120.        #
#######################################################
source("mtfad.R")
seed = 75728298
set.seed(seed) 
q1 <- 1:2 # possible factors for each cluster
param_grid1 <- as.matrix(expand.grid(q1 = q1, K = 1))
param_grid2 <- as.matrix(expand.grid(q1 = q1, q2 = q1, K = 2))
param_grid3 <- as.matrix(expand.grid(q1 = q1, q2 = q1, q3 = q1, K = 3))
param_grid4 <- as.matrix(expand.grid(q1 = q1, q2 = q1, q3 = q1, q4 = q1, K = 4))
param_grid5 <- as.matrix(expand.grid(q1 = q1, q2 = q1, q3 = q1, q4 = q1, q5 = q1,
                                     K = 6))
param_grid6 <- as.matrix(expand.grid(q1 = q1, q2 = q1, q3 = q1, q4 = q1, q5 = q1,
                                     q6 = q1, K = 6))
param_grid7 <- as.matrix(expand.grid(q1 = q1, q2 = q1, q3 = q1, q4 = q1, q5 = q1,
                                     q6 = q1, q7 = q1, K = 7))
param_grid8 <- as.matrix(expand.grid(q1 = q1, q2 = q1, q3 = q1, q4 = q1, q5 = q1,
                                     q6 = q1, q7 = q1, q8 = q1, K = 8))
param_grid9 <- as.matrix(expand.grid(q1 = q1, q2 = q1, q3 = q1, q4 = q1, q5 = q1,
                                     q6 = q1, q7 = q1, q8 = q1, q9 = q1, K = 9))
param_grid10 <- as.matrix(expand.grid(q1 = q1, q2 = q1, q3 = q1, q4 = q1, q5 = q1,
                                      q6 = q1, q7 = q1, q8 = q1, q9 = q1, q10 = q1, 
                                      K = 10))
param_grid11 <- as.matrix(expand.grid(q1 = q1, q2 = q1, q3 = q1, q4 = q1, q5 = q1,
                                      q6 = q1, q7 = q1, q8 = q1, q9 = q1, q10 = q1, 
                                      q11 = q1, K = 11))
param_grid12 <- as.matrix(expand.grid(q1 = q1, q2 = q1, q3 = q1, q4 = q1, q5 = q1,
                                      q6 = q1, q7 = q1, q8 = q1, q9 = q1, q10 = q1, 
                                      q11 = q1, q12 = q1, K = 12))
param_grid13 <- as.matrix(expand.grid(q1 = q1, q2 = q1, q3 = q1, q4 = q1, q5 = q1,
                                      q6 = q1, q7 = q1, q8 = q1, q9 = q1, q10 = q1, 
                                      q11 = q1, q12 = q1, q13 = q1, K = 13))
param_grid14 <- as.matrix(expand.grid(q1 = q1, q2 = q1, q3 = q1, q4 = q1, q5 = q1,
                                      q6 = q1, q7 = q1, q8 = q1, q9 = q1, q10 = q1, 
                                      q11 = q1, q12 = q1, q13 = q1, q14 = q1, 
                                      K = 14))
param_grid15 <- as.matrix(expand.grid(q1 = q1, q2 = q1, q3 = q1, q4 = q1, q5 = q1,
                                      q6 = q1, q7 = q1, q8 = q1, q9 = q1, q10 = q1, 
                                      q11 = q1, q12 = q1, q13 = q1, q14 = q1, 
                                      q15 = q1, K = 15))

#### The example below corresponds to K = 8 ####
#### For other values of K, e.g., K = 9, use param_grid9, etc. ####
#### Warning: Executing the code for K = 1 to 15 requires more time. ####
param_grid = param_grid8
BICs <- runtimes <- aris <- loglik  <- niter <- numeric(nrow(param_grid) )
tol = 1e-10
tols = rep(tol, nrow(param_grid) )
for (j in 1:nrow(param_grid) ){
  t1 <- proc.time()
  param_vec <- as.vector(param_grid[j,])
  KK <- param_vec[length(param_vec)]
  qq <- param_vec[-length(param_vec)]
  mtfads <- tryCatch({
    mtfad.q(X.lbd, KK, qq, tol =tol, nstart=20) 
  }, error= function(e) NA)
  runtimes[j]<- (proc.time()-t1)[3]
  est_cl  <- tryCatch({mtfads$clusters }, error = function(e) NA)
  BICs[j] <- tryCatch({mtfads$BIC  }, error= function(e) NA)
  niter[j]  <- tryCatch({mtfads$niter }, error= function(e) NA)
  loglik[j] <- tryCatch({mtfads$loglik },  error= function(e) NA)
  print(Sys.time())
  result_gama_q <- data.frame(
    BIC = BICs[j], runtimes = runtimes[j],
    niter = niter[j], loglik = loglik[j], tolerance= tols[j], 
    K=KK, q1 = qq[1], q2 = qq[2], q3 = qq[3], q4 = qq[4],
    q5 = qq[5], q6 = qq[6], q7 = qq[7], q8 = qq[8])
  print(result_gama_q)
  write.table( result_gama_q, file = paste0("gama_envopt_data_clustering_variedq_results_mtfad_K_",KK,"_seed_",seed,".csv") , 
               append = TRUE, sep = ",",
               col.names = !file.exists(paste0("gama_envopt_data_clustering_variedq_results_mtfad_K_",KK,"_seed_",seed,".csv")
               ), row.names = FALSE )
}
result_gama_q <- data.frame(
  BIC = BICs, runtimes = runtimes,  
  niter = niter, loglik = loglik, tolerance= tols)
result_gama_q <- cbind(result_gama_q, param_grid)
print(result_gama_q)
best_model <- result_gama_q[which.min(result_gama_q$BIC),]
print(best_model)
best_model.8 = best_model

# store the optimal model output
opt_qq <- c(1,  2,  2,  2,  2,  2,  2,  1) # From best_model.8
opt_K <- 8 
tol <- 1e-10
mtfads1 <- mtfad.q(X.lbd, opt_K, opt_qq, tol =tol, nstart=20)
save(mtfads1, X.lbd, opt_qq, file = "gama-data-envopt-mtfad-q-8g.rda")

######################################
# Analysis of the estimated clusters # 
######################################
load("gama-data-envopt-mtfad-q-8g.rda")

######################
# Table: Group sizes # 
######################
cl = factor(mtfads1$clusters)
table(cl) 

##################################################
# Table: Groupwise means and standard deviations # 
##################################################
df.g = data.frame(cl,X.lbd)
K = nlevels(cl)
for (i in 1:K) {
  x.g = df.g[df.g$cl==i,-1]
  print(round(as.numeric(colMeans(x.g)),digits = 2))
  sd = as.numeric(round(apply(x.g,2,sd),digits = 2))
  print(sd)
}

# Rank clusters by averaged u-r color values
u.r.m = NULL
K = nlevels(cl)
for (i in 1:K) {
  x.g = X.lbd[cl==i,]
  u.r.m = c(u.r.m, mean(x.g[,3]))
}
order(u.r.m, decreasing = T)

###################################################
# MOBSynC - Investigate the propensity of merging #
###################################################
library(mvtnorm)
size = 1e6
K = nlevels(cl)
W = diag(K)
p = 6
q.vec = opt_qq

M = t(mtfads1$means) # Means
L = mtfads1$lambda # Loading matrices
D = t(mtfads1$psi) # Uniquenesses
V = mtfads1$v # \nu parameters
prob = mtfads1$weights # Grouping probabilities

for (i in 1:K) {
  for (j in 1:K) {
    
    # density for group i
    ppi = 1
    m = as.matrix(M[,i])
    q = q.vec[i]
    l = array(NA,c(p,q,1))
    l[,,1] = L[[i]]
    d = as.matrix(D[,i])
    v = V[i]
    
    g = dim(d)[2]
    s = array(NA,c(p,p,g))
    for (k in 1:g) {
      s[,,k] = tcrossprod(l[,,k])+diag(d[,k])
    }
    s = s[,,1]
    
    xi = rep(m, each=size) + rmvt(size, sigma=s, df=v)
    
    f <- exp(dmvt(xi, delta = m, sigma = s, df = v) + 
               rep(log(prob[i]), each = nrow(xi)))
    f = as.matrix(f)
    lf.i <- log(rowSums(f))
    
    # density for group j
    ppi = 1
    m = as.matrix(M[,j])
    q = q.vec[j]
    l = array(NA,c(p,q,1))
    l[,,1] = L[[j]]
    d = as.matrix(D[,j])
    v = V[j]
    
    g = dim(d)[2]
    s = array(NA,c(p,p,g))
    for (k in 1:g) {
      s[,,k] = tcrossprod(l[,,k])+diag(d[,k])
    }
    s = s[,,1]
    
    f = exp(dmvt(xi, delta = m, sigma = s, df = v) + 
              rep(log(prob[j]), each = nrow(xi)))
    f = as.matrix(f)
    lf.j = log(rowSums(f))
    
    W[j,i] = sum(lf.j > lf.i)/size
  }
}
O = W + t(W) # Pairwise overlap matrix
diag(O) = rep(1,K)
gomega = (max(eigen(O)$values)-1)/(K-1) # Generalized overlap rate

################################
# Figure: Pairwise overlap map #
################################
library(MixSim)
source("utilis.R")
library(RColorBrewer)
colpal = brewer.pal(10, "RdBu")[c(1:4,7:10)]
colpal.8 = colpal[c(5,7,2,1,8,3,4,6)] 
dev.off()
overlap.map(O, lab.col = colpal.8, 
            map.cex = 1.25, lab.cex = 2,
            legend.cex = 1.5, font = 2, scale.pos = 0.35, legend.width = 0.1)

##################################
# Merge clusters                 #
# - Blue sequence: (1,2,5,6,7,8) #
# - Red sequence: (3,4)          #
##################################
cl.m2 = cl
ind = which(cl==1)
cl.m2[ind] = 1
ind = which(cl==2)
cl.m2[ind] = 1
ind = which(cl==3)
cl.m2[ind] = 2
ind = which(cl==4)
cl.m2[ind] = 2
ind = which(cl==5)
cl.m2[ind] = 1
ind = which(cl==6)
cl.m2[ind] = 1
ind = which(cl==7)
cl.m2[ind] = 1
ind = which(cl==8)
cl.m2[ind] = 1
cl.m2 = factor(cl.m2)
table(cl.m2)

##########################################
# Merge the colors for compound clusters #
##########################################
library(colorspace)
col.rgb = hex2RGB(colpal.8)
col.luv = as(col.rgb,"LUV")
g.size = as.numeric(table(cl))

# Blue sequence
g.ind = c(1,2,5:8) 
wg = g.size[g.ind]/sum(g.size[g.ind])
col.g = LUV(0,0,0,"LUV")
col.g = col.g@coords
count = 1
for (i in g.ind) {
  col.g = col.g + wg[count]*col.luv@coords[i,]
  count = count+1
}
col.g= LUV(col.g[1],col.g[2],col.g[3],"LUV")
col.no34 = hex(col.g)

# Red sequence
g.ind = c(3,4) 
wg = g.size[g.ind]/sum(g.size[g.ind])
col.g = LUV(0,0,0,"LUV")
col.g = col.g@coords
count = 1
for (i in g.ind) {
  col.g = col.g + wg[count]*col.luv@coords[i,]
  count = count+1
}
col.g= LUV(col.g[1],col.g[2],col.g[3],"LUV")
col.34 = hex(col.g)

colpal.m = c("#5F5C84" ,"#8A0B24") 


############################################################
# Figure: Distribution plots for estimated simple clusters #
############################################################
df.g = data.frame(class=factor(cl),X.lbd)
dev.off()
pp = ggpairs(df.g, 
             mapping = aes(color = class),
             upper = list(continuous = wrap("cor", size = 8.5)),
             diag = list(continuous = wrap("densityDiag", alpha = 0.7)),
             lower = list(continuous = wrap("points", alpha = 1, size=0.35)),
             columns = colnames(df.g),
             columnLabels = c("Cluster", "Stellar mass", "Star formation rate",
                              "u-r colour", "Half-light radius", "Sersic index",
                              "Optimal density")
) +
  scale_fill_manual(values = colpal.8) + 
  scale_colour_manual(values = colpal.8) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 60,hjust = 1),
        strip.background = element_rect(fill="#DEEBF7"),
        strip.text.x = element_text(size = 27),
        strip.text.y = element_text(size = 27))
for (i in 2:7) { 
  pp[i,1] = pp[i,1] + scale_y_continuous(n.breaks = 3)
}
for (i in 2:7) {
  pp[7,i] = pp[7,i] + scale_x_continuous(n.breaks = 3)
}
pp
ggsave(pp, file="gama-envopt-mtfad-q-initial-dist.png", 
       width=25, height=25)

####################################################################
# Figure: 3D star coordinates visualization of the simple clusters #
####################################################################
colnames(X.lbd) = c("Stellar mass", "Star formation rate","u-r colour", 
                    "Half-light radius", "Sersic index",
                    "Optimal density") 
source("starcoord3d.R")
plot.starcoords3D(data = X.lbd, pradius = 1, class = cl, 
                  cex = 0.5, axes.cex = 1.35, lwd = 1.5,
                  colors = colpal.8, pch = as.character(1:8))

####################################################
# Figure: Distribution plots for compound clusters #
####################################################
df.g = data.frame(class=factor(cl.m2),X.lbd)
pp = ggpairs(df.g, 
             mapping = aes(color = class),
             upper = list(continuous = wrap("cor", size = 8.5)),
             diag = list(continuous = wrap("densityDiag", alpha = 0.7)),
             lower = list(continuous = wrap("points", alpha = 1, size=0.35)),
             columns = colnames(df.g),
             columnLabels = c("Cluster", "Stellar mass", "Star formation rate",
                              "u-r colour", "Half-light radius", "Sersic index",
                              "Optimal density") 
) +
  scale_fill_manual(values = colpal.m) + 
  scale_colour_manual(values = colpal.m) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 60,hjust = 1),
        strip.background = element_rect(fill="#DEEBF7"),
        strip.text.x = element_text(size = 27),
        strip.text.y = element_text(size = 27))
for (i in 2:7) { 
  pp[i,1] = pp[i,1] + scale_y_continuous(n.breaks = 3)
}
for (i in 2:7) {
  pp[7,i] = pp[7,i] + scale_x_continuous(n.breaks = 3)
}
pp
ggsave(pp, file="gama-envopt-mtfad-q-merge1-dist.png", 
       width=25, height=25)

#################################################
# Estimated factor loadings for simple clusters #
#################################################
K = nlevels(cl)
p = 6
q.vec = opt_qq

ll = mtfads1$lambda # Groupwise factor loadings
dd = t(mtfads1$psi) # Uniquenesses

library(GPArotation)
library(ggplot2)
L = list()
for (k in 1:K) {
  sd = sqrt(rowSums(ll[[k]]^2) + dd[,k])
  L[[k]] = 1/sd*ll[[k]]
  if(dim(L[[k]])[2] > 1){
    L[[k]] = oblimin(L[[k]], maxit = 5000)$loadings
  }
}
for (k in 1:K) {
  mL = L[[k]]
  class(mL) = "loadings"
  print(mL)
}

##################################################################
# Figure: Heatmaps of the estimated loadings for simple clusters # 
##################################################################
lvalues = fa = feature = group = NULL

for (k in 1:K) {
  group = c(group, rep(paste0("Cluster ",k-1),p*q.vec[k]))
  for (q in 1:q.vec[k]) {
    lvalues = c(lvalues, L[[k]][,q])
    fa = c(fa, rep(q,p))
    feature = c(feature, factor(1:p))
  }
}
df = data.frame(feature,lvalues,fa,group)
df$feature = factor(df$feature)
df$fa = factor(df$fa)
df$group = factor(df$group)

lab_col  <- "#2717B5"
colors <- colorRampPalette(RColorBrewer::brewer.pal(name = "RdYlGn", n = 11)[1:11])(40)
label <- round(x = df$lvalues, digits = 3)
theme_set(theme_bw(base_size = 10))
feature.names = c("Stellar mass", "Star formation rate","u-r colour", 
                  "Half-light radius", "Sersic index",
                  "Optimal density") 
group.names <- list(
  'Cluster 0'="Cluster 1",
  'Cluster 1'="Cluster 2",
  'Cluster 2'="Cluster 3",
  'Cluster 3'="Cluster 4",
  'Cluster 4'="Cluster 5",
  'Cluster 5'="Cluster 6",
  'Cluster 6'="Cluster 7",
  'Cluster 7'="Cluster 8" 
)
group.labeller <- function(variable,value){
  return(group.names[value])
}

dev.off()
p = ggplot(data = df, aes(feature, fa, fill = lvalues)) + 
  ggplot2::geom_point(
    shape = 21,
    colour = "black",      
    stroke = 0.5,          
    ggplot2::aes_string(size=40),
  ) +
  ggplot2::scale_size(range = c(8.5, 9.5)) +
  ggplot2::guides(scale = 'none') +
  labs(fill='') +
  guides(size = "none",fill=guide_colorbar(ticks.colour = NA)) + 
  facet_wrap(.~group, ncol = 1, scales = "free_y",
             labeller = group.labeller) + 
  ylab('') + xlab('') + 
  scale_x_discrete(labels=feature.names) + 
  scale_y_discrete(limits=rev) +
  scale_fill_gradientn(colours = colors, values = c(0, 1),
                       breaks=c(min(df$lvalues),0.5,0,-0.5,max(df$lvalues)),
                       labels=c(-1,0.5,0,-0.5,1)) +
  theme(axis.text.x = element_text(angle = 90,
                                   vjust = 1, 
                                   size = 14.5, hjust = 1),
        legend.position = 'right',
        strip.text.x = element_text(size = 15),
        strip.background = element_rect(fill = "#deebf7",colour = "#000000"), 
        legend.key.width = unit(1.75, "mm"),
        legend.key.height = unit(43.5, "mm"),
        legend.text=element_text(size=14.5),
        legend.justification="right", 
        legend.box.spacing = unit(0, "pt"),
        legend.margin=margin(-30,0,0,5),
        axis.text.y = element_text(size = 14.5),
        panel.margin = unit(0.5, "pt"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  guides(
    size = "none",
    fill = guide_colorbar(
      ticks.colour = NA,
      frame.colour = "black",
      frame.linewidth = 0.5
    )
  )
p
ggsave(p, file="gama-envopt-mtfad-q-initial-loadings.png", 
       width=3.25, height=10.75)

###########################
# Estimated factor scores #
###########################
Z = gr <- NULL
K = nlevels(cl)
p = 6

ll = mtfads1$lambda 
dd = t(mtfads1$psi)
mm = t(mtfads1$means)

# For 1-factor model
for (k in c(1,8)) { # 1,3,4,6,7,8
  X.k <- X.lbd[cl==k,]
  dd.inv <- as.numeric(1/dd[,k])
  svd.H <- svd(crossprod(ll[[k]],as.numeric(1/dd[,k])*ll[[k]]))
  Z <- rbind(Z, t(t(X.k) - mm[,k]) %*% (as.numeric(1/dd[,k])*ll[[k]]) %*% (svd.H$u %*% (1/svd.H$d * t(svd.H$v))))
  gr <- c(gr,rep(k,sum(cl==k)))
}
grb.score = data.frame(Z, "group" = gr)
# Averaged groupwise factor scores
for (k in c(1,8)) {
  xx = grb.score[grb.score$group==k,-3]
  print(round(colMeans(xx),digits = 3))
}

# For 2-factor model
for (k in c(1,3,4,6,7,8)) {
  X.k <- X.lbd[cl==k,]
  dd.inv <- as.numeric(1/dd[,k])
  Ga <- solve(crossprod(ll[[k]],dd.inv*ll[[k]]))
  for (i in 1:nrow(X.k)) {
    Z <- rbind(Z, t(Ga%*%t(ll[[k]])%*%diag(dd.inv)%*%t(X.k[i,]-mm[,k])))
   }
  gr <- c(gr,rep(k,sum(cl==k)))
}
grb.score = data.frame(Z, "group" = gr)
# Averaged groupwise factor scores
for (k in c(1,8)) {
  xx = grb.score[grb.score$group==k,-3]
  print(round(colMeans(xx),digits = 3))
}

#########################################
# Characterization of compound clusters #
#########################################
source("Gtrans-sync.R")
index1 = c(1,2,5:8) 
index2 = c(3,4) 
K = nlevels(cl)
p = 6
mm = t(mtfads1$means) 
ll = mtfads1$lambda 
d = t(mtfads1$psi) 
V = mtfads1$v 
prob = mtfads1$weights 

# For blue sequence - compound cluster (1,2,5,6,7,8)
mu = mm[,index1]
L = ll[index1] 
D = dd[,index1]
prob = prob[index1]
kk = length(prob)
p = nrow(mu)
sigma = array(NA,c(p,p,kk))
for (k in 1:kk) {
  sigma[,,k] = tcrossprod(L[[k]]) + diag(D[,k]) 
}
df.blue = X.lbd[(cl!=3) & (cl!=5), ]
# Gaussianization
df.blue.gdt = Gtrans.sync(data = df.blue, mu, sigma, prob) 
# Factor analysis
library(fad)
bic.blue = rep(NA,2)
for (q in 1:2) {
  g.blue.fad = fad(as.matrix(df.blue.gdt),q)
  bic.blue[q] = g.blue.fad$BIC
}
q.blue = which.min(bic.blue) 
g.blue.fad = fad(as.matrix(df.blue.gdt),q.blue,scores = "Bartlett")
# Rotate factor loadings
mL.blue = oblimin(g.blue.fad$loadings, maxit = 5000)$loadings
class(mL.blue) = "loadings"
print(mL.blue)

# For red sequence - compound cluster (3,4)
mu = mm[,index2]
L = ll[index2] 
D = dd[,index2]
prob = prob[index2]
kk = length(prob)
p = nrow(mu)
sigma = array(NA,c(p,p,kk))
for (k in 1:kk) {
  sigma[,,k] = tcrossprod(L[[k]]) + diag(D[,k]) 
}
df.red = X.lbd[(cl==3) | (cl==4), ] # 3, 5
df.red.gdt = Gtrans.sync(data = df.red, mu, sigma, prob)
bic.red = rep(NA,2)
for (q in 1:2) {
  g.red.fad = fad(as.matrix(df.red.gdt),q)
  bic.red[q] = g.red.fad$BIC
}
q.red = which.min(bic.red) 
g.red.fad = fad(as.matrix(df.red.gdt),q.red,scores = "Bartlett")
mL.red = oblimin(g.red.fad$loadings, maxit = 5000)$loadings
class(mL.red) = "loadings"
print(mL.red)

###########################################################################
# Figure: Heatmaps of the estimated factor loadings for compound clusters #
###########################################################################
p = 6
lvalues = fa = feature = group = NULL
qq = q.blue
group = c(group, rep("(1,2,5,6,7,8)",p*qq)) 
for (q in 1:qq) {
  lvalues = c(lvalues, mL.blue[,q])
  fa = c(fa, rep(q,p))
  feature = c(feature, factor(1:p))
}
qq = q.red
group = c(group, rep("(3,4)",p*qq)) 
for (q in 1:qq) {
  lvalues = c(lvalues, mL.red[,q])
  fa = c(fa, rep(q,p))
  feature = c(feature, factor(1:p))
}
df = data.frame(feature,lvalues,fa,group)
df$feature = factor(df$feature)
df$fa = factor(df$fa)
df$group = factor(df$group)

lab_col  <- "#2717B5"
colors <- colorRampPalette(RColorBrewer::brewer.pal(name = "RdYlGn", n = 11)[1:11])(40)
label <- round(x = df$lvalues, digits = 3)
theme_set(theme_bw(base_size = 10))
fa.names = c("1"="Factor 1","2"="Factor 2")
feature.names = c("Stellar mass", "Star formation rate", "u-r colour", 
                  "Half-light radius", "Sersic index",
                  "Optimal density") 
dev.off()
p = ggplot(data = df, aes(feature, fa, fill = lvalues)) + 
  ggplot2::geom_point(
    shape = 21,
    colour = "black",      
    stroke = 0.5,          
    ggplot2::aes_string(size=40),
  ) +
  ggplot2::scale_size(range = c(10.5, 11.5)) +
  ggplot2::guides(scale = 'none') +
  labs(fill='') +
  guides(size = "none",fill=guide_colorbar(ticks.colour = NA)) + 
  facet_wrap(.~group,ncol = 1, scales = "free_y") + 
  ylab('') + xlab('') + 
  scale_x_discrete(labels=feature.names) + 
  scale_y_discrete(limits=rev) +
  scale_fill_gradientn(colours = colors, values = c(0, 1),
                       breaks=c(min(df$lvalues),0.5,0,-0.5,max(df$lvalues)),
                       labels=c(-1,0.5,0,-0.5,1)) +
  theme(axis.text.x = element_text(angle = 90,
                                   vjust = 1, 
                                   size = 18.5, hjust = 1),
        legend.position = 'right',
        strip.text.x = element_text(size = 18),
        strip.background = element_rect(fill = "#deebf7",colour = "#000000"), 
        legend.key.width = unit(1.75, "mm"),
        legend.key.height = unit(11.75, "mm"),
        legend.text=element_text(size=18.5),
        legend.justification="right", 
        legend.box.spacing = unit(0, "pt"),
        legend.margin=margin(-30,0,0,5),
        axis.text.y = element_text(size = 18.5),
        panel.margin = unit(0.5, "pt"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  guides(
    size = "none",
    fill = guide_colorbar(
      ticks.colour = NA,
      frame.colour = "black",
      frame.linewidth = 0.5
    )
  )
p
ggsave(p, file="gama-envopt-mtfad-q-merge1-loadings.png", 
       width=3.75, height=5) 
