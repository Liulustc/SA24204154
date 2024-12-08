## -----------------------------------------------------------------------------
# 绘制正弦曲线
x <- seq(0, 2*pi, length.out = 100)
y <- sin(x)

plot(x, y, type = "l", col = "blue", lwd = 2, 
     main = "正弦函数 y = sin(x)", xlab = "x", ylab = "y")

## -----------------------------------------------------------------------------
set.seed(123)
n <- 6
students <- c("A", "B", "C", "D", "E", "F")
Chinese_scores <- rnorm(n, mean = 80, sd = 10)  
Math_scores <- rnorm(n, mean = 80, sd = 10)     
English_scores <- rnorm(n, mean = 80, sd = 10) 

print(round(Chinese_scores, 1))
print(round(Math_scores, 1))
print(round(English_scores, 1))

## -----------------------------------------------------------------------------
# 定义x的范围
x <- seq(0.01, 5, length.out = 100)  # 避免x为0
# n=2时
plot(x, x^2, type="l", col="skyblue", lwd=2, ylim=c(-1, 10), 
     xlab="x", ylab="y", main="y = x^n, n = -1, 1/2, 2")
# 添加n=1/2和-1
lines(x, x^(-1), col="pink", lwd=2)  # n = -1
lines(x, x^(1/2), col="lightgreen", lwd=2)  # n = 1/2
# 添加图例
legend("topright", legend=c("n = -1", "n = 1/2", "n = 2"), 
       col=c("pink", "lightgreen", "skyblue"), lwd=2)

## -----------------------------------------------------------------------------
generate_Rayleigh <- function(n, sigma){
  U <- runif(n)
  X <- sigma*sqrt(-2*log(U))
  return(X)
}

set.seed(789)
sigma <- 1
n <- 10000
rayleigh <- generate_Rayleigh(n, sigma)

hist(rayleigh, prob = "T", main = paste("Rayleigh 分布 (σ =", sigma, ")"),breaks = 20)

x <- seq(0, max(rayleigh), length = 10000)
f_x <- (x / sigma^2) * exp(-x^2 / (2 * sigma^2))
#理论函数图像
lines(x, f_x, col = "red")  
abline(v = sigma, col = "blue", lty = 2) 

## -----------------------------------------------------------------------------
generate_Mixture <- function(n, p1){
  # 生成0 1,0代表N(0,1), 1代表N(3,1)
  mix <- rbinom(n, 1, p1)
  samples <- ifelse(mix == 1, 
                    rnorm(n, mean = 0, sd = 1), 
                    rnorm(n, mean = 3, sd = 1))  
  return(samples)
}

set.seed(789)
p1_values <- c(0.25, 0.35, 0.5, 0.65, 0.75)  
n <- 1000  

for (p1 in p1_values) {
  samples <- generate_Mixture(n, p1)
  hist(samples, breaks = 30, freq = FALSE, col = "lightblue", border = "black", 
         main = paste("正态混合分布 (p1 =", p1, ")"), xlab = "x", ylab = "密度")
  x <- seq(min(samples), max(samples), length.out = 1000)
  f_x <- p1 * dnorm(x, mean = 0, sd = 1) + (1 - p1) * dnorm(x, mean = 3, sd = 1)
  lines(x, f_x, col = "red")
}

## -----------------------------------------------------------------------------
set.seed(789)
generate_Poisson_Gamma <- function(lambda, alpha, beta, t){
  N_t <- rpois(1, lambda * t)
  if (N_t > 0) {
    Y <- rgamma(N_t, alpha, beta) 
    X_t <- sum(Y)
  } else {
    X_t <- 0  
  }
  return(X_t)
}

lambda <- 2  
alpha <- 2 
beta <- 1   
t <- 10    
n <- 10000 
X_t_values <- 1:n

for (i in 1:n) {
    X_t_values[i] <- generate_Poisson_Gamma(lambda, alpha, beta, t)
}

# 理论结果
mean1 <- lambda * t * (alpha / beta)
var1 <- lambda * t * (alpha * (alpha + 1) / beta^2)
# 实验结果
mean2 <- mean(X_t_values)
var2 <- var(X_t_values)

cat("模拟均值: ", mean2, "\n")
cat("理论均值: ", mean1, "\n")
cat("模拟方差: ", var2, "\n")
cat("理论方差: ", var1, "\n")
hist(X_t_values, breaks = 30, col = "lightblue", border = "black",
     main = paste("复合泊松-伽玛过程"),
     xlab = "X(t)")

## -----------------------------------------------------------------------------
beta <- function(x, n = 10000){
  beta_samples <- rbeta(n, shape1 = 3, shape2 = 3)
  estimate <- mean(beta_samples <= x)
  return(estimate)
}
x <- seq(0.1, 0.9, 0.1)
# monte carlo
monte_carlo_estimates <- sapply(x, beta)
# pbeta
pbeta_values <- pbeta(x, shape1 = 3, shape2 = 3)

results <- data.frame(x, monte_carlo_estimates, pbeta_values)
print(results)

## -----------------------------------------------------------------------------
# 对偶变量法
rayleigh <- function(sigma, n) {
  U <- runif(n)
  X <- sigma * sqrt(-2 * log(U))
  X_prime <- sigma * sqrt(-2 * log(1 - U))
  return((X + X_prime) / 2)  
}

sigma <- 1  
n <- 10000
# 对偶
samples_1 <- rayleigh(sigma, n)
# 独立
samples_2_1 <- sigma * sqrt(-2 * log(runif(n)))
samples_2_2 <- sigma * sqrt(-2 * log(runif(n)))
# 方差
var_1 <- var(samples_1)
print(var_1)
var_2 <- (var(samples_2_1) + var(samples_2_2)) / 2
print(var_2)
# 方差缩减百分比
var_reduce <- (var_2 - var_1) / var_2 * 100
print(var_reduce)

## -----------------------------------------------------------------------------
# g(x)
g <- function(x) {
  (x^2 / sqrt(2 * pi)) * exp(-x^2 / 2)
}
# 求C1和C2的值
C1 <- 1 / integrate(function(x) exp(-x), 1, Inf)$value
C2 <- 1 / integrate(function(x) x * exp(-x), 1, Inf)$value
# 重要性函数 f1 和 f2
f1 <- function(x) {
  C1 * exp(-x)
}
f2 <- function(x) {
  C2 * x * exp(-x)
}
# 计算方差
var_estimate <- function(f, g, n) {
  samples <- rexp(n, rate = 1) + 1  
  weights <- g(samples) / f(samples)
  return(var(weights))
}

n <- 10000
var_f1 <- var_estimate(f1, g, n)
var_f2 <- var_estimate(f2, g, n)
print(var_f1)
print(var_f2)

## -----------------------------------------------------------------------------
# 快速排序
fast_sort <- function(x) {
  if (length(x) <= 1) {
    return(x)
  }
  pivot_index <- sample(length(x), 1)
  pivot <- x[pivot_index]
  left <- x[x < pivot]
  right <- x[x > pivot]
  return(c(fast_sort(left), pivot, fast_sort(right)))
}
# 计算平均时间
compute_average_time <- function(n, simulations = 100) {
  times <- numeric(simulations)
  for (i in 1:simulations) {
    numbers <- sample(1:n)
    # 开始时间
    start_time <- Sys.time()
    sorted_numbers <- fast_sort(numbers)
    # 结束时间
    end_time <- Sys.time()
    # 时间差
    times[i] <- as.numeric(end_time - start_time, units = "secs")
  }
  return(mean(times))
}

n_values <- c(10^4, 2 * 10^4, 4 * 10^4, 6 * 10^4, 8 * 10^4)
average_times <- sapply(n_values, compute_average_time)
#  t_n = n log(n)
t_n <- n_values * log(n_values)
# 回归分析
model <- lm(average_times ~ t_n)
summary(model)
# 绘制散点图和回归线
plot(t_n, average_times, 
     main = "平均计算时间与 n log(n) 的关系",
     xlab = "t_n = n log(n)", 
     ylab = "平均计算时间 (秒)",
     pch = 19, col = "skyblue")
abline(model, col = "pink")

## -----------------------------------------------------------------------------
set.seed(456)
n <- 1000    # 样本大小
N <- 10000   # 模拟次数
b1 <- numeric(N) # 偏度
for (i in 1:N) {
  sample_data <- rnorm(n)
  sample_sd <- sd(sample_data)
  if (sample_sd > 0) {
    b1[i] <- sum((sample_data - mean(sample_data))^3) / (n * sample_sd^3)
  } else {
    b1[i] <- NA  
  }
}
# 删除NA值
b1 <- na.omit(b1)
quantiles <- quantile(b1, probs = c(0.025, 0.05, 0.95, 0.975))

# 输出蒙特卡洛结果
print(quantiles)

# 输出理论结果
se_sqrt_b1 <- sqrt(6 / n)
# 输出标准误差
print(se_sqrt_b1)
theoretical_quantiles <- qnorm(c(0.025, 0.05, 0.95, 0.975), mean = 0, sd = sqrt(6 / n))
print(theoretical_quantiles)


## -----------------------------------------------------------------------------
library(MASS)
n <- 1000      
N <- 1000      
alpha <- 0.05  
# 拒绝原假设的次数
pearson_reject <- 0
spearman_reject <- 0
kendall_reject <- 0

# 模拟二元正态分布
for (i in 1:N) {
  mu <- c(0, 0)
  sigma <- matrix(c(1, 0.03, 0.03, 6), 2, 2)  
  data <- mvrnorm(n, mu, sigma)
  # 进行相关性检验
  pearson_test <- cor.test(data[,1], data[,2], method = "pearson")
  spearman_test <- cor.test(data[,1], data[,2], method = "spearman")
  kendall_test <- cor.test(data[,1], data[,2], method = "kendall")
  # 检查是否拒绝原假设
  if (pearson_test$p.value < alpha) pearson_reject <- pearson_reject + 1
  if (spearman_test$p.value < alpha) spearman_reject <- spearman_reject + 1
  if (kendall_test$p.value < alpha) kendall_reject <- kendall_reject + 1
}
cat("二元正态分布下，拒绝原假设的次数：\n")
cat("Pearson:", pearson_reject, "\n")
cat("Spearman:", spearman_reject, "\n")
cat("Kendall:", kendall_reject, "\n")


## ----warning=FALSE------------------------------------------------------------

set.seed(789)  

N <- 1000
null_hypotheses <- 950
alt_hypotheses <- 50
alpha <- 0.1
m <- 10000

# 存储不同方法下的FWER、FDR 和 TPR
results <- matrix(0, nrow = 3, ncol = 2)
rownames(results) <- c("FWER", "FDR", "TPR")
colnames(results) <- c("Bonferroni correction", "B-H correction")
# 10000次模拟
for (sim in 1:m) {
  # 生成 p 值
  p_values <- c(runif(null_hypotheses), rbeta(alt_hypotheses, 0.1, 1))
  
  # Bonferroni 
  bonferroni_p <- p.adjust(p_values, method = "bonferroni")
  bonferroni_significant <- bonferroni_p < alpha
  
  # B-H 
  bh_p <- p.adjust(p_values, method = "BH")
  bh_significant <- bh_p < alpha
  
  # 计算 FWER,FDR,TPR
  FWER_bonferroni <- mean(bonferroni_significant[1:null_hypotheses])
  FDR_bonferroni <- mean(bonferroni_significant) / max(1, sum(bonferroni_significant))
  TPR_bonferroni <- sum(bonferroni_significant[(null_hypotheses + 1):N]) / alt_hypotheses
  
  FWER_BH <- mean(bh_significant[1:null_hypotheses])
  FDR_BH <- mean(bh_significant) / max(1, sum(bh_significant))
  TPR_BH <- sum(bh_significant[(null_hypotheses + 1):N]) / alt_hypotheses
  
  # 累加
  results[1, 1] <- results[1, 1] + FWER_bonferroni
  results[2, 1] <- results[2, 1] + FDR_bonferroni
  results[3, 1] <- results[3, 1] + TPR_bonferroni
  
  results[1, 2] <- results[1, 2] + FWER_BH
  results[2, 2] <- results[2, 2] + FDR_BH
  results[3, 2] <- results[3, 2] + TPR_BH
}

results <- results / m
results


## -----------------------------------------------------------------------------

library(boot)

data <- c(3, 5, 7, 18, 43, 85, 91, 98, 100, 130, 230, 487)
# 计算指数分布的MLE函数
lambda_mle <- function(data) {
  1 / mean(data) 
}

lambda_hat <- lambda_mle(data)

# bootstrap估计偏差和标准误
bootstrap_function <- function(data, indices) {
  sample_data <- data[indices]  
  return(lambda_mle(sample_data))  
}

set.seed(789)  
boot_results <- boot(data, bootstrap_function, R = 10000)

# 计算偏差和标准误
bias <- mean(boot_results$t) - lambda_hat
std_error <- sd(boot_results$t)
# 输出
cat("MLE of λ:", lambda_hat, "\n")
cat("Bias:", bias, "\n")
cat("Standard Error:", std_error, "\n")

## -----------------------------------------------------------------------------
# 自助法函数，用于估计平均故障时间
mean_time_between_failures <- function(data, indices) {
  sample_data <- data[indices]  # 抽样
  lambda_sample <- lambda_mle(sample_data)  # 计算抽样的 λ
  return(1 / lambda_sample)  # 返回平均故障时间
}

# 执行自助法
set.seed(123)  # 设置随机种子以保证可重复性
boot_results <- boot(data, mean_time_between_failures, R = 10000)

# 计算 95% 自助法置信区间

# 1. 标准正态法
mean_estimate <- mean(boot_results$t)  # 自助样本均值
std_error <- sd(boot_results$t)  # 自助样本标准差
z_value <- qnorm(0.975)  # 正态分布的临界值
ci_normal <- mean_estimate + c(-z_value, z_value) * std_error  # 置信区间

# 2. 基本法
ci_basic <- boot.ci(boot_results, type = "basic")$basic[4:5]

# 3. 百分位法
ci_percentile <- quantile(boot_results$t, c(0.025, 0.975))

# 4. BCa 法
ci_bca <- boot.ci(boot_results, type = "bca")$bca[4:5]


# 输出置信区间
cat("95% Confidence Interval (Normal):", ci_normal, "\n")
cat("95% Confidence Interval (Basic):", ci_basic, "\n")
cat("95% Confidence Interval (Percentile):", ci_percentile, "\n")
cat("95% Confidence Interval (BCa):", ci_bca, "\n")

## ----warning=FALSE------------------------------------------------------------
set.seed(456)
d1 <- c(-2.961, 0.478, -0.391, -0.869, -0.460, -0.937, 0.779, -1.409, 0.027, -1.569)
d2 <- c(1.608, 1.009, 0.878, 1.600, -0.263, 0.680, 2.280, 2.390, 1.793, 8.091, 1.468)

mean <- mean(d2) - mean(d1)

# bootstrap方法
bootstrap_diff <- function(data1, data2, R=10000) {
  n1 <- length(data1)
  n2 <- length(data2)
  boot_diff <- numeric(R)
  
  for (i in 1:R) {
    resample1 <- sample(data1, n1, replace = TRUE)
    resample2 <- sample(data2, n2, replace = TRUE)
    boot_diff[i] <- mean(resample2) - mean(resample1)
  }
  
  return(boot_diff)
}

# bootstrap抽样
boot_diff <- bootstrap_diff(d1, d2)

# bootstrap bias
bootstrap_bias <- mean(boot_diff) - mean

# bootstrap标准差
bootstrap_se <- sd(boot_diff)
# 样本标准差
obs_se <- sqrt(var(d1)/length(d1) + var(d2)/length(d2))

cat("样本均值:", mean, "\n")
cat("bootstrap偏差:", bootstrap_bias, "\n")
cat("原始样本标准差:", obs_se, "\n")
cat("bootstrap标准差:", bootstrap_se, "\n")

## -----------------------------------------------------------------------------
library(bootstrap)
hat_Sigma <- cov(scor)
n <- nrow(scor)

eigen_values <- eigen(hat_Sigma)$values
hat_theta <- eigen_values[1] / sum(eigen_values)

jackknife_theta <- function(data) {
  theta_jack <- numeric(n)
  for (i in 1:n) {
    # 删除第 i 个观测
    temp_data <- data[-i, ]
    temp_Sigma <- cov(temp_data)
    temp_eigen_values <- eigen(temp_Sigma)$values
    theta_jack[i] <- temp_eigen_values[1] / sum(temp_eigen_values)
  }
  return(theta_jack)
}

jackknife_estimates <- jackknife_theta(scor)

bias <- (n-1)*(mean(jackknife_estimates) - hat_theta)
se <- sqrt((n - 1) * mean((jackknife_estimates - mean(jackknife_estimates))^2))

bias
se

## -----------------------------------------------------------------------------
library(DAAG)
attach(ironslag)
n <- length(magnetic) #in DAAG ironslag
e1 <- e2 <- e3 <- e4 <- numeric(n)

for (k in 1:n){
    y <- magnetic[-k]
    x <- chemical[-k]
    J1 <- lm(y ~ x)
    yhat1 <- J1$coef[1] + J1$coef[2] * chemical[k]
    e1[k] <- magnetic[k] - yhat1
    
    J2 <- lm(y ~ x + I(x^2))
    yhat2 <- J2$coef[1] + J2$coef[2] * chemical[k] +
        J2$coef[3] * chemical[k]^2
    e2[k] <- magnetic[k] - yhat2
    
    J3 <- lm(log(y) ~ x)
    logyhat3 <- J3$coef[1] + J3$coef[2] * chemical[k]
    yhat3 <- exp(logyhat3)
    e3[k] <- magnetic[k] - yhat3
    
    # 将双对数模型换为三次多项式模型
    J4 <- lm(y~x + I(x^2) + I(x^3))
    yhat4 <- J4$coef[1] + J4$coef[2]*chemical[k] + J4$coef[3]*
      chemical[k]^2 + J4$coef[4]*chemical[k]^3
    e4[k] <- magnetic[k]- yhat4
}

c(mean(e1^2), mean(e2^2), mean(e3^2), mean(e4^2))

## -----------------------------------------------------------------------------
y <- magnetic
x <- chemical
lm(y ~ x + I(x^2))

## -----------------------------------------------------------------------------
library(DAAG)
 n <- length(magnetic) 
 r1 <- r2 <- r3 <- r4 <- numeric(n)
 for (k in 1:n) {
   y <- ironslag$magnetic[-k]
   x <- ironslag$chemical[-k]
   
   J1 <- lm(y ~ x)
   r1[k] <- summary(J1)$adj.r.squared
   
   J2 <- lm(y ~ x + I(x^2))
   r2[k] <- summary(J2)$adj.r.squared
   
   J3 <- lm(log(y) ~ x)
   r3[k] <- summary(J3)$adj.r.squared
   
   # 将双对数模型换为三次多项式模型
   J4<-lm(y~x+I(x^2)+I(x^3))
   r4[k] <- summary(J4)$adj.r.squared
 
 }
 c(mean(r1^2), mean(r2^2), mean(r3^2), mean(r4^2))

## -----------------------------------------------------------------------------
y <- magnetic
x <- chemical
lm(y ~ x + I(x^2))

## -----------------------------------------------------------------------------
attach(chickwts)
x <- sort(as.vector(chickwts$weight[chickwts$feed == "soybean"]))
y <- sort(as.vector(chickwts$weight[chickwts$feed == "linseed"]))
detach(chickwts)
# Cramér-von Mises 统计量
cramer_von_mises_stat <- function(x, y) {
  n1 <- length(x)
  n2 <- length(y)
  z <- c(x, y)
  n <- n1 + n2
  
  ecdf_x <- ecdf(x)(sort(z))
  ecdf_y <- ecdf(y)(sort(z))
  
  W <- sum((ecdf_x - ecdf_y)^2) * (n1 * n2) / n
  return(W)
}

R <- 999 
z <- c(x, y) 
n1 <- length(x)
n2 <- length(y)
reps <- numeric(R) 

# 计算原始统计量
W0 <- cramer_von_mises_stat(x, y)

# 置换检验
for (i in 1:R) {
  permuted <- sample(z)
  x1 <- permuted[1:n1]
  y1 <- permuted[(n1 + 1):(n1 + n2)]
  
  reps[i] <- cramer_von_mises_stat(x1, y1)
}

p_value <- mean(c(W0, reps) >= W0)
p_value


## -----------------------------------------------------------------------------
data("chickwts")
x <- as.vector(chickwts$weight[chickwts$feed == "sunflower"])
y <- as.vector(chickwts$weight[chickwts$feed == "linseed"])

spearman_stat <- cor(x, y, method = "spearman")

R <- 999  
z <- c(x, y) 
n1 <- length(x)  
n2 <- length(y) 
reps <- numeric(R)  

set.seed(123) 
for (i in 1:R) {
  permuted_z <- sample(z)
  x1 <- permuted_z[1:n1]
  y1 <- permuted_z[(n1 + 1):(n1 + n2)]
  reps[i] <- cor(x1, y1, method = "spearman")
}

p_value_permutation <- mean(c(spearman_stat, reps) >= spearman_stat)
p_value_permutation


# cor.test p 值
cor_test_result <- cor.test(x, y, method = "spearman")
p_value_cor_test <- cor_test_result$p.value
p_value_cor_test

## -----------------------------------------------------------------------------
set.seed(456)

# 柯西分布的密度函数
dcauchy_target <- function(x) {
  return(1 / (pi * (1 + x^2)))
}
q_proposal <- function(x_current, x_proposed, sigma = 1) {
  # 建议分布选择正态分布
  return(dnorm(x_proposed, mean = x_current, sd = sigma))
}
# Metropolis-Hastings算法生成样本
mh_cauchy <- function(n, burn_in = 1000, sigma = 1) {
  x <- numeric(n + burn_in)
  x[1] <- 0  
  
  for (i in 2:(n + burn_in)) {
    proposal <- rnorm(1, mean = x[i - 1], sd = sigma)
    target_ratio <- dcauchy_target(proposal) / dcauchy_target(x[i - 1])
    proposal_ratio <- q_proposal(proposal, x[i - 1], sigma) / q_proposal(x[i - 1], proposal, sigma)
    alpha <- min(1, target_ratio * proposal_ratio)
    if (runif(1) < alpha) {
      x[i] <- proposal
    } else {
      x[i] <- x[i - 1]
    }
  }
  return(x[(burn_in + 1):(n + burn_in)])
}

n <- 10000
samples <- mh_cauchy(n)

# 生成样本的分位数
sample_deciles <- quantile(samples, probs = seq(0.1, 0.9, by = 0.1))

# 标准柯西分布的理论分位数
theoretical_deciles <- qcauchy(seq(0.1, 0.9, by = 0.1))

# 分位数
comparison <- data.frame(
  'Decile' = seq(0.1, 0.9, by = 0.1),
  'Sample_Deciles' = sample_deciles,
  'Theoretical_Deciles' = theoretical_deciles
)

print(comparison)

# 密度函数
plot(seq(0.1, 0.9, by = 0.1), sample_deciles, type = 'b', col = 'blue',
     xlab = 'Deciles', ylab = 'Values', main = 'Sample vs Theoretical Deciles')
points(seq(0.1, 0.9, by = 0.1), theoretical_deciles, type = 'b', col = 'red')
legend('topleft', legend = c('Sample Deciles', 'Theoretical Deciles'),
       col = c('blue', 'red'), lty = 1, bty = 'n')


## -----------------------------------------------------------------------------
set.seed(12345)
library(coda)

run_multiple_chains <- function(num_chains, iter, burn_in, sigma) {
  chains <- list()
  
  for (i in 1:num_chains) {
    chains[[i]] <- mh_cauchy(iter, burn_in, sigma)
  }
  
  mcmc_chains <- mcmc.list(lapply(chains, mcmc))
  return(mcmc_chains)
}

num_chains <- 4       
burn_in <- 1000       
sigma <- 1           
converged_n <- NA     
max_iterations <- 1000
r_values <- numeric()  # 每次迭代的 R2 
n_values <- numeric()  

# 寻找最先收敛的n
for (n in seq(100, max_iterations, by = 1)) {
  chains <- run_multiple_chains(num_chains, n, burn_in, sigma)
  gelman_results <- gelman.diag(chains)

  if (all(gelman_results$psrf[, "Point est."] < 1.2)) {
    cat("Chains have converged at n =", n, "(R < 1.2)!\n")
    converged_n <- n
    break  
  } 
}


if (!is.na(converged_n)) {
  cat("The chains first converged at n =", converged_n, "\n")
} else {
  cat("The chains did not converge within the maximum iterations.\n")
}


## -----------------------------------------------------------------------------
set.seed(456)
gibbs_sampler <- function(a, b, n, iter, burn_in) {
 
  x_chain <- numeric(iter)
  y_chain <- numeric(iter)
  
  x_chain[1] <- rbinom(1, n, 0.5)  # 随机选择一个初值
  y_chain[1] <- rbeta(1, x_chain[1] + a, n - x_chain[1] + b)
  
  for (i in 2:iter) {
    x_chain[i] <- rbinom(1, n, y_chain[i - 1])
    y_chain[i] <- rbeta(1, x_chain[i] + a, n - x_chain[i] + b)
  }
  return(list(x = x_chain[(burn_in + 1):iter], y = y_chain[(burn_in + 1):iter]))
}

a <- 2
b <- 2
n <- 10
iterations <- 10000
burn_in <- 1000

samples <-gibbs_sampler (a, b, n, iterations, burn_in)

head(samples$x)
head(samples$y)

plot(samples$x, samples$y, 
     xlab = "x", ylab = "y", pch = 20, col = rgb(0.2, 0.4, 0.6, 0.5))


## -----------------------------------------------------------------------------
library(ggplot2)
set.seed(2456)
# Gibbs 采样函数
gibbs <- function(n_samples, a, b, n) {
  x <- numeric(n_samples)
  y <- numeric(n_samples)
  
  y[1] <- 0.5
  x[1] <- 0.5
  
  for (i in 2:n_samples) {
    x1 <- x[i - 1]
    y1 <- y[i - 1]
    x1 <- rbinom(1, n, y1)
    x[i] <- x1
    y[i] <- rbeta(1, x1 + a, n - x1 + b)
  }
  
  return(data.frame(x = x, y = y))
}

# Gelman-Rubin 诊断函数
gelman_rubin <- function(chains) {
  m <- ncol(chains)  # 链数
  n <- nrow(chains)  # 样本数
  means <- colMeans(chains)
  overall_mean <- mean(means)
  
  B <- n * var(means)  # 序列间方差
  W <- mean(apply(chains, 2, var))  # 序列内方差
  var_hat <- (n - 1) / n * W + B / n
  R_hat <- sqrt(var_hat / W)  # 使用 sqrt 提供更标准的 Rhat
  return(R_hat)
}

# 生成多个链
generate_chains <- function(n_samples, n_chains, a, b, n) {
  chains_x <- matrix(NA, n_samples, n_chains)
  chains_y <- matrix(NA, n_samples, n_chains)
  
  for (i in 1:n_chains) {
    samples <- gibbs(n_samples, a, b, n)
    chains_x[, i] <- samples$x
    chains_y[, i] <- samples$y
  }
  
  return(list(chains_x = chains_x, chains_y = chains_y))
}

# 计算 R_hat 并绘图
calculate_rhat_and_plot <- function(n_samples, n_chains, a, b, n) {
  chains <- generate_chains(n_samples, n_chains, a, b, n)
  chains_x <- chains$chains_x
  chains_y <- chains$chains_y
  
  R_hat_x <- numeric(n_samples)
  R_hat_y <- numeric(n_samples)
  
  for (k in 2:n_samples) {
    R_hat_x[k] <- gelman_rubin(chains_x[1:k, ])
    R_hat_y[k] <- gelman_rubin(chains_y[1:k, ])
  }
  
  # 创建数据框以用于 ggplot
  df <- data.frame(
    Iteration = 1:n_samples,
    Rhat_x = R_hat_x,
    Rhat_y = R_hat_y
  )
  
  # 绘制 R_hat 曲线
  ggplot(df, aes(x = Iteration)) +
    geom_line(aes(y = Rhat_x, color = "Rhat_x")) +
    geom_line(aes(y = Rhat_y, color = "Rhat_y")) +
    labs(
      title = "Gelman-Rubin Diagnostic (Rhat) over Iterations",
      y = "Rhat",
      color = "Chains"
    ) +
    theme_minimal() +
    geom_hline(yintercept = 1.2, linetype = "dashed", color = "red")  # 加入收敛标准线
}

# 参数设置
a <- 2     
b <- 3 
n_samples <- 10000
n_chains <- 4
n <- 10

# 计算并绘制 Rhat
calculate_rhat_and_plot(n_samples, n_chains, a, b, n)


## -----------------------------------------------------------------------------
library(pracma)

compute_k <- function(k, a, d) {
  
  numerator_1 <- (-1)^k / (factorial(k) * 2^k)  
  numerator_2 <- norm(a,"2")^(2*k+2) / ((2 * k + 1) * (2 * k + 2))  
  gamma_part <- exp(lgamma((d+1)/2)+lgamma(k+3/2)-lgamma(k+d/2+1))
  
  kth <- numerator_1 * numerator_2 * gamma_part
  
  return(kth)
}

## -----------------------------------------------------------------------------
compute_sum <- function(max_k, a, d) {
  sum <- 0  
  for (k in 0:max_k) {
    term_k <- compute_k(k, a, d)
    sum <- sum + term_k
  }
  return(sum)
}

## -----------------------------------------------------------------------------
a <- c(1,2)
d <- 1

compute_sum(100, a, d)

## -----------------------------------------------------------------------------
left_integral <- function(ck1, k) {
  integrand <- function(u) (1 + u^2 / (k - 1))^(-k / 2)
  integrate(integrand, 0, ck1)$value
}

right_integral <- function(ck, k) {
  integrand <- function(u) (1 + u^2 / k)^(-(k + 1) / 2)
  integrate(integrand, 0, ck)$value
}

equation <- function(a, k) {
  ck1 <- sqrt(a^2 * (k - 1) / (k - a^2))
  ck  <- sqrt(a^2 * k / (k + 1 - a^2))
  left <- (2 * exp(lgamma(k / 2) - lgamma((k - 1) / 2)) / sqrt(pi * (k - 1))) * left_integral(ck1, k)
  right <- (2 * exp(lgamma((k + 1) / 2) - lgamma(k / 2)) / sqrt(pi * k)) * right_integral(ck, k)
  return(left - right)
}

solve_for_a <- function(k) {
  solution <- uniroot(equation, interval = c(0.1, sqrt(k) - 1), k = k)
  return(solution$root)
}

# 运行并输出结果
k <- 10
a_solution <- solve_for_a(k)
cat("当 k =", k, "时，a 的解为：", a_solution, "\n")


## -----------------------------------------------------------------------------
y <- c(0.54, 0.48, 0.33, 0.43, 1.00, 1.00, 0.91, 1.00, 0.21, 0.85)
tau <- 1  
uncensored <- y < tau
censored <- y == tau
n_u <- sum(uncensored)  # 未截尾数据个数
n_c <- sum(censored)     # 截尾数据个数
# 初始lambda值，使用未截尾数据的均值作为初始值
lambda <- mean(y[uncensored])
# 迭代条件
tol <- 1e-6
max_iter <- 1000
iter <- 0
lambda_diff <- Inf
# E-M 算法
while (iter < max_iter && lambda_diff > tol) {
  iter <- iter + 1
  # E步
  E_T_censored <- tau + lambda 
  # M步
  lambda_new <- (sum(y[uncensored]) + n_c * E_T_censored) / (n_u + n_c)
  # 计算收敛条件
  lambda_diff <- abs(lambda_new - lambda)
  lambda <- lambda_new
}
cat("通过E-M算法估计的lambda值为:", lambda, "\n")
cat("迭代次数:", iter, "\n")

## -----------------------------------------------------------------------------
# MLE
# 计算未截尾数据的和
sum_uncensored <- sum(y[uncensored])

lambda_mle <- (sum_uncensored + n_c * tau) / (n_u )
cat("lambda值的MLE估计为:", lambda_mle, "\n")

## -----------------------------------------------------------------------------
library(lpSolve)

ob <- c(4, 2, 9)
constraints <- matrix(c(
  2, 1, 1,  
  1, -1, 3 
), nrow = 2, byrow = TRUE)

rhs <- c(2, 3)

directions <- c("<=", "<=")

result <- lp("min", ob, constraints, directions, rhs, compute.sens = TRUE)

# 输出结果
cat("最优解:\n")
cat("x =", result$solution[1],";")
cat("y =", result$solution[2],";")
cat("z =", result$solution[3],";")
cat("最小值为：", result$objval, "\n")


## -----------------------------------------------------------------------------
data(mtcars)

formulas <- list(
    mpg ~ disp,
    mpg ~ I(1 / disp),
    mpg ~ disp + wt,
    mpg ~ I(1 / disp) + wt
)
# for loop
loop_for <- list()  
for (i in seq_along(formulas)) {
    loop_for[[i]] <- lm(formulas[[i]], data = mtcars)
}
loop_for
# lapply() 
loop_lapply <- lapply(formulas, function(f) lm(f, data = mtcars))
loop_lapply

## -----------------------------------------------------------------------------
set.seed(456)  
bootstraps <- lapply(1:10, function(i) {
    rows <- sample(1:nrow(mtcars), replace = TRUE)
    mtcars[rows, ]
})
# for loop
loop_for4 <- list()  
for (i in seq_along(bootstraps)) {
    loop_for4[[i]] <- lm(mpg ~ disp, data = bootstraps[[i]])
}
loop_for4

fit_lm <- function(data) {
    lm(mpg ~ disp, data = data)
}
# lapply() 
lapply4 <- lapply(bootstraps, fit_lm)
lapply4

## -----------------------------------------------------------------------------
rsq <- function(mod) summary(mod)$r.squared

rsq_for3 <- numeric(length(loop_for))  
for (i in seq_along(loop_for)) {
    rsq_for3[i] <- rsq(loop_for[[i]])
}
rsq_for3

rsq_lapply3 <- sapply(loop_lapply, rsq)
rsq_lapply3

rsq_for4 <- numeric(length(loop_for4))  
for (i in seq_along(loop_for4)) {
    rsq_for4[i] <- rsq(loop_for4[[i]])
}
rsq_for4

rsq_lapply4 <- sapply(lapply4, rsq)
rsq_lapply4

## -----------------------------------------------------------------------------
set.seed(45678)
trials <- replicate(
    100,
    t.test(rpois(10, 10), rpois(7, 10)),
    simplify = FALSE
)

p_values <- sapply(trials, function(trial) trial$p.value)
head(p_values)

p_values_direct <- sapply(trials, `[[`, "p.value")
p_values_direct

# 比较
all.equal(p_values, p_values_direct) 

## -----------------------------------------------------------------------------
map_vapply <- function(FUN, ..., FUN.VALUE) {
    inputs <- list(...)  
    mapped <- Map(FUN, ...)  
    vapply(mapped, identity, FUN.VALUE)  
}

x <- 1:8
y <- c(1, 2, 5, 8, 7, 9, 5, 7)

add <- function(a, b) a + b

result <- map_vapply(add, x, y, FUN.VALUE = numeric(1))
print(result)  

add_diff <- function(a, b) c(sum = a + b, diff = a - b)

result_matrix <- map_vapply(add_diff, x, y, FUN.VALUE = numeric(2))
print(result_matrix)

## ----warning=FALSE------------------------------------------------------------
fast_chisq_test1 <- function(x, y) {
    tbl <- table(x, y)
    row_sums <- rowSums(tbl)
    col_sums <- colSums(tbl)
    total <- sum(tbl)
    expected <- outer(row_sums, col_sums) / total
    chisq_fast1 <- sum((tbl - expected)^2 / expected)
    return(chisq_fast1)
}

fast_chisq_test2 <- function(x, y) {
    tbl <- table(x, y)
    row_sums <- rowSums(tbl)
    col_sums <- colSums(tbl)
    total <- sum(tbl)
    expected <- outer(row_sums, col_sums) / total
    observed <- as.numeric(tbl)
    expected <- as.numeric(expected)
    chisq_fast2 <- sum((observed  - expected)^2 / expected)
    return(chisq_fast2)
}

x <- c(3, 1, 2, 2, 3, 2, 3)
y <- c(3, 2, 3, 3, 1, 3, 1)

chisq_fast1 <- fast_chisq_test1(x, y)
print(chisq_fast1) 

chisq_fast2 <- fast_chisq_test2(x, y)
print(chisq_fast2) 

chisq_test <- chisq.test(table(x, y))
print(chisq_test$statistic)  


## ----warning=FALSE------------------------------------------------------------
fast_table <- function(x, y) {
    stopifnot(is.integer(x), is.integer(y))
    unique_x <- unique(x)
    unique_y <- unique(y)
    tbl <- matrix(0, nrow = length(unique_x), ncol = length(unique_y))
    for (i in 1:length(x)) {
        x_idx <- which(unique_x == x[i])
        y_idx <- which(unique_y == y[i])
        tbl[x_idx, y_idx] <- tbl[x_idx, y_idx] + 1
    }
    dimnames(tbl) <- list(unique_x, unique_y)
    return(tbl)
}
fast_chisq_test <- function(x, y) {
    tbl <- fast_table(x, y)
    row_sums <- rowSums(tbl)
    col_sums <- colSums(tbl)
    total <- sum(tbl)
    expected <- outer(row_sums, col_sums) / total
    chisq_stat <- sum((tbl - expected)^2 / expected)
    return(chisq_stat)
}

x <- as.integer(c(3, 1, 2, 2, 3, 2, 3))
y <- as.integer(c(3, 2, 3, 3, 1, 3, 1))

tbl_fast <- fast_table(x, y)
print(tbl_fast)

chisq_stat <- fast_chisq_test(x, y)
print(chisq_stat)

chisq_test <- chisq.test(table(x, y))
print(chisq_test$statistic)  

## ----warning=FALSE------------------------------------------------------------
 library(Rcpp) 
cppFunction('
NumericMatrix gibbs_sampler(int n, double a, double b, int num_samples, int burn_in) {
    NumericMatrix samples(num_samples, 2); 
    int x = 0; 
    double y = 0.5; 
    // Gibbs 采样
    for (int i = 0; i < num_samples + burn_in; i++) {
        x = R::rbinom(n, y);
        y = R::rbeta(x + a, n - x + b);
        if (i >= burn_in) {
            samples(i - burn_in, 0) = x;
            samples(i - burn_in, 1) = y;
        }
    }
    return samples;
}
')
set.seed(456)
n <- 10
a <- 2
b <- 3
num_samples <- 10000
burn_in <- 1000

cppsamples <- gibbs_sampler(n, a, b, num_samples, burn_in)
par(mfrow = c(1, 2))
plot(cppsamples[, 2], type = 'l', ylab = 'y')
hist(cppsamples[, 2], breaks = 30, main = 'y', xlab = 'y')

## -----------------------------------------------------------------------------
rgibbs_sampler <- function(a, b, n, iter, burn_in) {
 
  x_chain <- numeric(iter)
  y_chain <- numeric(iter)
  
  x_chain[1] <- rbinom(1, n, 0.5)  
  y_chain[1] <- rbeta(1, x_chain[1] + a, n - x_chain[1] + b)
  
  for (i in 2:iter) {
    x_chain[i] <- rbinom(1, n, y_chain[i - 1])
    y_chain[i] <- rbeta(1, x_chain[i] + a, n - x_chain[i] + b)
  }
  return(list(x = x_chain[(burn_in + 1):iter], y = y_chain[(burn_in + 1):iter]))
}
set.seed(456)
a <- 2
b <- 3
n <- 10
iterations <- 10000
burn_in <- 1000

rsamples <-rgibbs_sampler (a, b, n, iterations, burn_in)
qqplot(cppsamples[, 1], rsamples[[1]])
qqplot(cppsamples[, 2], rsamples[[2]])

## -----------------------------------------------------------------------------
library('microbenchmark')
# c++ gibbs_sampler(n, a, b, num_samples, burn_in)
# r rgibbs_sampler (a, b, n, iterations, burn_in)
time <- microbenchmark(gibbs_sampler(n, a, b, num_samples, burn_in), rgibbs_sampler (a, b, n, iterations, burn_in))
summary(time)

