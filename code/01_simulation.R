##
## 01_simulation.R
## 一元稳健回归课程报告 —— 数值模拟（蒙特卡洛实验）
## 模型: y = beta0 + beta1 * x + eps,  beta0 = 1, beta1 = 2, n = 100
## 情景: 正态 / t(3) 厚尾 / Cauchy 厚尾 / 污染正态 / 纵向离群点 / 坏杠杆点
## 方法: OLS, LAD, Huber-M, MM, LTS
## 输出: figures/ 下的图形 + output/sim_summary.csv
## 

pkgs <- c("MASS", "robustbase", "quantreg")
invisible(lapply(pkgs, library, character.only = TRUE))

base_dir <- "f:/27 研一上/2.回归分析/稳健报告"
fig_dir  <- file.path(base_dir, "figures")
out_dir  <- file.path(base_dir, "output")

## 中文字体（Windows / R>=4.2 原生 UTF-8）
tryCatch({
  windowsFonts(yahei = windowsFont("Microsoft YaHei"))
  par(family = "yahei")
}, error = function(e) NULL)

## 可复现性：固定随机数生成器
## "Rejection" 与 R >= 3.6 的默认一致；显式指定可避免不同 R 版本下
## sample() 行为变化导致的差异。
RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion",
        sample.kind = "Rejection")

SEED_DESIGN <- 2026      ## 设计点 x 的种子
SEED_BASE   <- 2026000   ## 每个重复的种子基数: SEED_BASE + 1000*s + r
SEED_FIG    <- 2027      ## 绘图示例数据的种子（与模拟互不干扰）


#### 1. 数据生成机制 (DGP) #####

beta0 <- 1
beta1 <- 2
n <- 100

## 固定设计点（便于不同重复之间比较方差）
set.seed(SEED_DESIGN)
x_design <- runif(n, 0, 10)

## 各情景的误差 / 污染机制
gen_data <- function(scenario, n, x) {
  x_new <- x
  if (scenario == "leverage") {
    ## 10% 坏杠杆点: x 推到 [15,20], y 偏离真实直线 -40
    idx <- sample(n, ceiling(0.10 * n))
    x_new[idx] <- runif(length(idx), 15, 20)
    y <- beta0 + beta1 * x_new + rnorm(n)
    y[idx] <- beta0 + beta1 * x_new[idx] - 40
    return(list(x = x_new, y = y))
  }
  eps <- switch(scenario,
    normal   = rnorm(n),
    t3       = rt(n, df = 3),
    cauchy   = rcauchy(n),
    contam   = {                              ## 90% N(0,1) + 10% N(0,10^2)
      m <- rbinom(n, 1, 0.10)
      rnorm(n, mean = 0, sd = ifelse(m == 1, 10, 1))
    },
    vertical = rnorm(n)
  )
  y <- beta0 + beta1 * x + eps
  if (scenario == "vertical") {               ## 10% 纵向离群点: y 上移 20
    idx <- sample(n, ceiling(0.10 * n))
    y[idx] <- y[idx] + 20
  }
  list(x = x_new, y = y)
}

scen_ids   <- c("normal", "t3", "cauchy", "contam", "vertical", "leverage")
scen_names <- c("S1 正态误差", "S2 t(3) 厚尾", "S3 Cauchy 厚尾",
                "S4 污染正态", "S5 纵向离群点", "S6 坏杠杆点")


################ 2. 估计方法  ################

methods <- c("OLS", "LAD", "Huber", "MM", "LTS")
method_cn <- c(OLS = "OLS 最小二乘", LAD = "LAD 最小绝对偏差",
               Huber = "Huber M", MM = "MM 估计", LTS = "LTS 最小截尾平方")

fit_all <- function(x, y, seed = NULL) {
  res <- matrix(NA_real_, length(methods), 2,
                dimnames = list(methods, c("b0", "b1")))
  try(res["OLS", ]      <- coef(lm(y ~ x)), silent = TRUE)
  try(res["LAD", ]      <- coef(quantreg::rq(y ~ x)), silent = TRUE)
  try(res["Huber", ]    <- coef(MASS::rlm(y ~ x, psi = psi.huber,
                                          k = 1.345, maxit = 100)), silent = TRUE)
  ## MM 与 LTS 内部含随机子采样：各自调用前重置随机流，
  ## 使其不受前序拟合消耗的随机数影响（ltsReg 的 seed 参数需 .Random.seed，故用 set.seed）
  if (!is.null(seed)) set.seed(seed)
  try(res["MM", ] <- coef(robustbase::lmrob(y ~ x)), silent = TRUE)
  if (!is.null(seed)) set.seed(seed + 1L)
  try(res["LTS", ] <- coef(robustbase::ltsReg(y ~ x, alpha = 0.5, nsamp = 500)),
      silent = TRUE)
  res
}


######################## 3. 蒙特卡洛主循环  ######################## 

R <- 500   ## 重复次数

slope_sims <- array(NA_real_, dim = c(R, length(methods), length(scen_ids)),
                    dimnames = list(NULL, methods, scen_ids))
intc_sims  <- slope_sims

t0 <- Sys.time()
for (s in seq_along(scen_ids)) {
  sc <- scen_ids[s]
  for (r in 1:R) {
    ## 关键：每个 (情景, 重复) 都重新设种子。
    ## 这样 lmrob()/ltsReg() 内部随机子采样所消耗的随机数，不会影响下一次重复。
    set.seed(SEED_BASE + 1000L * s + r)
    dat <- gen_data(sc, n, x_design)
    est <- fit_all(dat$x, dat$y, seed = SEED_BASE + 1000L * s + r)
    slope_sims[r, , s] <- est[, "b1"]
    intc_sims[r, , s]  <- est[, "b0"]
  }
  cat(sprintf("[%s] 完成情景 %s, 累计用时 %.1f s\n",
              format(Sys.time(), "%H:%M:%S"), sc,
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}


###############  4. 汇总指标 ################ 

mk_summary <- function(sims, true_val, param_label) {
  out <- do.call(rbind, lapply(scen_ids, function(sc) {
    mat <- sims[, , sc, drop = FALSE][, , 1]
    mse_ols <- mean((mat[, "OLS"] - true_val)^2, na.rm = TRUE)
    do.call(rbind, lapply(methods, function(m) {
      v  <- mat[, m]
      ok <- is.finite(v)
      data.frame(
        scenario = sc, parameter = param_label, method = m,
        bias    = mean(v[ok]) - true_val,
        var     = var(v[ok]),
        mse     = mean((v[ok] - true_val)^2),
        rel_eff = mse_ols / mean((v[ok] - true_val)^2),  ## 相对 OLS 的 MSE 效率
        n_fail  = sum(!ok),
        row.names = NULL
      )
    }))
  }))
  out
}

sum_b1 <- mk_summary(slope_sims, beta1, "slope")
sum_b0 <- mk_summary(intc_sims,  beta0, "intercept")
sim_summary <- rbind(sum_b1, sum_b0)
write.csv(sim_summary, file.path(out_dir, "sim_summary.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

## 原始抽样分布也存一份(供报告画密度图等)
saveRDS(list(slope = slope_sims, intercept = intc_sims,
             methods = methods, scenarios = scen_ids),
        file.path(out_dir, "sim_raw.rds"))

## ---- 复现清单：记录随机种子、RNG 算法与软件版本 ----
manifest <- c(
  paste0("生成时间: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("R 版本: ", R.version.string),
  paste0("RNGkind: ", paste(RNGkind(), collapse = " / ")),
  paste0("SEED_DESIGN = ", SEED_DESIGN),
  paste0("SEED_BASE   = ", SEED_BASE, "（第 s 情景第 r 次重复的种子 = SEED_BASE + 1000*s + r）"),
  paste0("SEED_FIG    = ", SEED_FIG),
  paste0("重复次数 R = ", R),
  paste0("样本量 n   = ", n),
  paste0("robustbase = ", as.character(packageVersion("robustbase"))),
  paste0("MASS       = ", as.character(packageVersion("MASS"))),
  paste0("quantreg   = ", as.character(packageVersion("quantreg"))),
  paste0("sim_summary.csv 的 MD5: ",
         unname(tools::md5sum(file.path(out_dir, "sim_summary.csv"))))
)
writeLines(manifest, file.path(out_dir, "sim_manifest.txt"), useBytes = TRUE)
cat(paste(manifest, collapse = "\n"), "\n")


########################## 5. 图形 ######################## 


## 图1: 每个情景一个示例数据集 + 真值线 / OLS / MM / LTS
png(file.path(fig_dir, "fig1_scenarios.png"), width = 1800, height = 1200,
    res = 200, family = "yahei")
par(mfrow = c(2, 3), mar = c(4, 4, 2.5, 1), family = "yahei")
set.seed(SEED_FIG)
for (s in seq_along(scen_ids)) {
  dat <- gen_data(scen_ids[s], n, x_design)
  plot(dat$x, dat$y, pch = 16, cex = 0.7, col = grey(0.35),
       main = scen_names[s], xlab = "x", ylab = "y")
  xx <- seq(min(dat$x), max(dat$x), length.out = 200)
  lines(xx, beta0 + beta1 * xx, lwd = 2, lty = 2)
  lines(xx, coef(lm(dat$y ~ dat$x))[1] + coef(lm(dat$y ~ dat$x))[2] * xx,
        col = "red", lwd = 2)
  set.seed(SEED_FIG + s)
  cf <- coef(robustbase::lmrob(dat$y ~ dat$x))
  lines(xx, cf[1] + cf[2] * xx, col = "blue", lwd = 2)
  set.seed(SEED_FIG + s + 1L)
  cf <- coef(robustbase::ltsReg(dat$y ~ dat$x, alpha = 0.5, nsamp = 500))
  lines(xx, cf[1] + cf[2] * xx, col = "darkgreen", lwd = 2)
  if (s == 1)
    legend("topleft", bty = "n", cex = 0.85, lwd = 2,
           lty = c(2, 1, 1, 1), col = c("black", "red", "blue", "darkgreen"),
           legend = c("真实直线", "OLS", "MM", "LTS"))
}
dev.off()

## 图2: 各情景下斜率估计抽样分布箱线图
png(file.path(fig_dir, "fig2_boxplot_slope.png"), width = 1800, height = 1200,
    res = 200, family = "yahei")
par(mfrow = c(2, 3), mar = c(4.5, 4, 2.5, 1), family = "yahei")
box_col <- c("red", "orange", "skyblue3", "blue", "darkgreen")
for (s in seq_along(scen_ids)) {
  mat <- slope_sims[, , s]
  lim <- range(mat[is.finite(mat)])
  lim <- c(max(lim[1], beta1 - 3), min(lim[2], beta1 + 3))
  boxplot(mat, col = box_col, names = methods, las = 2, cex.axis = 0.8,
          main = scen_names[s], ylab = expression(hat(beta)[1]),
          ylim = lim, outline = FALSE)
  abline(h = beta1, lty = 2, lwd = 2)
}
dev.off()

## 图3: 效率-稳健性权衡散点
## x: 正态情景下相对 OLS 的效率(渐近相对效率的蒙特卡洛估计, 用方差比)
## y: 最坏情景(坏杠杆点)下的 MSE (log 轴)
are_norm <- sapply(methods, function(m)
  var(slope_sims[, "OLS", "normal"]) / var(slope_sims[, m, "normal"]))
mse_lev  <- sapply(methods, function(m)
  mean((slope_sims[, m, "leverage"] - beta1)^2, na.rm = TRUE))

png(file.path(fig_dir, "fig3_tradeoff.png"), width = 1400, height = 1100,
    res = 200, family = "yahei")
par(mar = c(4.5, 4.5, 2.5, 1), family = "yahei")
plot(are_norm, mse_lev, log = "y", pch = 19, cex = 1.6, col = box_col,
     xlab = "正态误差下的相对效率(方差比, 相对 OLS)",
     ylab = "坏杠杆点情景下的 MSE(对数刻度)",
     main = "效率—稳健性权衡(鱼与熊掌)", xlim = c(0, 1.15))
lab_pos <- c(OLS = 2, LAD = 4, Huber = 1, MM = 4, LTS = 3)  ## 手动调整标签方位避免重叠
for (m in methods)
  text(are_norm[m], mse_lev[m], labels = m, pos = lab_pos[m],
       cex = 0.9, offset = 0.8)
abline(v = 1, lty = 3, col = grey(0.6))
dev.off()

## 图4: 正态情景下斜率估计的抽样分布密度(效率可视化)
png(file.path(fig_dir, "fig4_density_normal.png"), width = 1500, height = 1000,
    res = 200, family = "yahei")
par(mar = c(4.5, 4.5, 2.5, 1), family = "yahei")
mat <- slope_sims[, , "normal"]
d <- lapply(methods, function(m) density(mat[, m], na.rm = TRUE))
xlim <- range(unlist(lapply(d, function(z) z$x)))
ylim <- range(unlist(lapply(d, function(z) z$y)))
plot(d[[1]], col = box_col[1], lwd = 2, xlim = xlim, ylim = ylim,
     main = "正态误差下斜率估计的抽样分布", xlab = expression(hat(beta)[1]))
for (i in 2:length(methods)) lines(d[[i]], col = box_col[i], lwd = 2)
abline(v = beta1, lty = 2)
legend("topright", bty = "n", col = box_col, lwd = 2, cex = 0.85,
       legend = paste0(methods, " (SD=",
                       formatC(apply(mat, 2, sd, na.rm = TRUE), digits = 3,
                               format = "f"), ")"))
dev.off()

cat("数值模拟完成, 结果已写入:", out_dir, "\n")
