# For Bayes Factor post
# Zijun Zhang, 12.31.2016

library(ggplot2)
library(gridExtra)

myTheme= theme(legend.text = element_text(size = 10),
						plot.title = element_text(size=11, face="bold"), axis.title.y = element_text(size=10), axis.title.x = element_text(size=10),
						axis.text.y = element_text(size=10, angle = 90), axis.text.x = element_text(size=10, angle=0),
						panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1), legend.background = element_rect(fill = "transparent", colour = "transparent"))

# predicted probability
pred_prob = function(y, N, theta)  choose(N,y) * theta^y * (1-theta)^(N-y)

# truncated normal density
trunc_dnorm = function(x, mu, sd) dnorm(x, mean=mu,sd=sd)/(pnorm(1,mu,sd)-pnorm(0,mu,sd))

plot_pred_prob = function(N, y.obs, p1=0.6, p2=0.7)
# plot predicted probability under two models
# N: total counts; y.obs: observed count of heads; p1, p2: two specified parameter (theta) value
{
	all_outcome = 1:N
	pred1 = sapply(all_outcome, pred_prob, N=N, theta=p1)
	pred2 = sapply(all_outcome, pred_prob, N=N, theta=p2)
	df = data.frame(x=rep(all_outcome, 2), y=c(pred1, pred2), H=c(rep('1',N), rep('2',N)))
	p = ggplot() + geom_bar(data=df, aes(x=x, y=y, fill=H),
		stat='identity', position='dodge') + myTheme +
		geom_rect(mapping=aes(xmin=y.obs-0.5, xmax=y.obs+0.5, ymin=0, ymax=max(df$y), fill="3"), alpha=0.5) + labs(x='No. of heads', y='Predicted prob.') + xlim(40,90) + 
		scale_fill_manual(values=c('red',"steelblue", 'grey'), labels=c('M1','M2', 'observed'))
	return(p)
}

plot_prior_pred_dist = function(p1, p2, tau1, tau2, y.obs, N, dx=0.001)
# use a small dx to approximate the density
# p1, p2: mean for two model; tau1, tau2: standard deviation
# y.obs: observed y/heads; N: total counts
# dx: precision of approximation
{
	dx.step = seq(0, 1, by=dx)
	p1_pmf = sapply(dx.step, trunc_dnorm, mu=p1, sd=tau1) * dx
	p2_pmf = sapply(dx.step, trunc_dnorm, mu=p2, sd=tau2) * dx
	df1 = data.frame(x=rep(1:length(dx.step),2), y=c(p1_pmf, p2_pmf), H=c(rep('1',length(dx.step)),rep('2',length(dx.step))))
	pg1 = ggplot() + geom_bar(data=df1, aes(x=x, y=y, fill=H),
		stat='identity', position='dodge') + myTheme + 
		scale_fill_manual(values=c('red',"steelblue"), labels=c('M1','M2')) + labs(title='Prior distribution for theta',x=paste0('theta x ', 1/dx), y='prob.')
	
	all_outcome = 1:N
	compute_post_pmf = function(pmf)
	{
		post_pmf = rep(0,N)
		for(i in 1:length(all_outcome)) 
		{
			y = all_outcome[i]
			p.y = 0
			for(j in 1:length(dx.step))
			{
				theta = dx.step[j]
				weight = pmf[j]
				p.y = p.y + pred_prob(y, N, theta) * weight
			}
			post_pmf[i] = p.y
		}
		return(post_pmf)
	}
	p1_post_pmf = compute_post_pmf(p1_pmf)
	p2_post_pmf = compute_post_pmf(p2_pmf)
	df2 = data.frame(x=rep(all_outcome, 2), y=c(p1_post_pmf, p2_post_pmf), H=c(rep('1',N), rep('2',N)))
	pg2 = ggplot() + geom_bar(data=df2, aes(x=x, y=y, fill=H),
		stat='identity', position='dodge') + myTheme +
		geom_rect(mapping=aes(xmin=y.obs-0.5, xmax=y.obs+0.5, ymin=0, ymax=max(df2$y), fill="3"), alpha=0.5) + 
		labs(title='Predicted prob. for y',x='No. of heads', y='Predicted prob.') + xlim(0,100) + 
		scale_fill_manual(values=c('red',"steelblue", 'grey'), labels=c('M1','M2', 'obs.'))
	return(list(p1=pg1, p2=pg2))
}


binom_norm_pdf = function(theta, stuff)
# theta: parameter to be integrate
# stuff: a list consists of N, y, p, tau
{
	N = stuff$N
	y = stuff$y
	p = stuff$p
	tau = stuff$tau
	dbinom(y, size=N, prob=theta, log=F) * 
		dnorm(theta, mean=p, sd=tau)/(pnorm(1, mean=p, sd=tau)-pnorm(0, mean=p, sd=tau))
}

bayes.factor =function(p1, p2, tau1, tau2, N, y)
# a wrapper for computing bayes factor of binomial with truncated normal prior
# p: prior mean; tau: prior standard deviation
# N: total trials; y: No. of successes
{
	m1=integrate(binom_norm_pdf, 0, 1, stuff=list(N=N, y=y, p=p1, tau=tau1))
	m2=integrate(binom_norm_pdf, 0, 1, stuff=list(N=N, y=y, p=p2, tau=tau2))
	m1$value / m2$value
}

## command-line usage
# generating figure 1. 
f1 = plot_pred_prob(N=100, y.obs=62, p1=0.6, p2=0.7)
png('f1.png', 500,250)
f1
dev.off()

# generating figure 2.
f2.list =  plot_prior_pred_dist(p1=0.6, p2=0.7, tau1=0.2, tau2=0.1, y.obs=62, N=100)
png('f2.png', 500, 500)
grid.arrange(f2.list$p1, f2.list$p2, nrow=2, ncol=1)
dev.off()

# testing our hand-made bayes factor calculator
# reproduce the example, second scenario in the first section
bayes.factor(0.6, 0.7, 0.2, 0.1, 100, 62)
### 0.722
# shrink the variance to approximate the first scenario 
bayes.factor(0.6, 0.7, 0.01, 0.01, 100, 62)
### 3.726
