## 02_real_data.R
## 一元稳健回归课程报告 —— 实例分析
## 数据: MASS::phones (比利时电话呼叫, 经典稳健回归教材数据, 含离群点)
##       datasets::cars (干净对照数据)
##       datasets::anscombe (安斯科姆四联图, 诊断教学)
## 输出: figures/ 下的图形 + output/real_phones.csv / real_cars.csv
##
## 【可复现性说明】自助法同样采用"每个重复独立设种子"：
## 第 b 次自助抽样的种子为 SEED_BOOT + b，且 MM/LTS 的内部随机子采样也被固定，
## 因此结果不依赖运行顺序、R 版本或包版本。


pkgs <- c("MASS", "robustbase", "quantreg")

invisible(lapply(pkgs, library, character.only = TRUE))

base_dir <- "f:/27 研一上/2.回归分析/稳健报告"
fig_dir  <- file.path(base_dir, "figures")
out_dir  <- file.path(base_dir, "output")

tryCatch({
  windowsFonts(yahei = windowsFont("Microsoft YaHei"))
  par(family = "yahei")
}, error = function(e) NULL)

## ---------- 可复现性：固定随机数生成器 ----------
RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion",
        sample.kind = "Rejection")
SEED_BOOT <- 2026

methods <- c("OLS", "LAD", "Huber", "MM", "LTS")
box_col <- c(OLS = "red", LAD = "orange", Huber = "skyblue3",
             MM = "blue", LTS = "darkgreen")

fit_all <- function(x, y, seed = NULL) {
  res <- matrix(NA_real_, length(methods), 2,
                dimnames = list(methods, c("b0", "b1")))
  try(res["OLS", ]      <- coef(lm(y ~ x)), silent = TRUE)
  try(res["LAD", ]      <- coef(quantreg::rq(y ~ x)), silent = TRUE)
  try(res["Huber", ]    <- coef(MASS::rlm(y ~ x, psi = psi.huber,
                                          k = 1.345, maxit = 100)), silent = TRUE)
  if (!is.null(seed)) set.seed(seed)
  try(res["MM", ] <- coef(robustbase::lmrob(y ~ x)), silent = TRUE)
  if (!is.null(seed)) set.seed(seed + 1L)
  try(res["LTS", ] <- coef(robustbase::ltsReg(y ~ x, alpha = 0.5, nsamp = 500)),
      silent = TRUE)
  res
}

#### 4.1 MASS::phones —— 比利时电话呼叫数据 (1950-1973) ####

data(phones, package = "MASS")
x <- phones$year
y <- phones$calls   ## 单位: 百万次

set.seed(SEED_BOOT)
est <- fit_all(x, y, seed = SEED_BOOT)

## ---- 自助法 (pairs bootstrap) 估计斜率标准误 ----
B <- 1000
boot_b1 <- matrix(NA_real_, B, length(methods),
                  dimnames = list(NULL, methods))
for (b in 1:B) {
  ## 每次自助抽样独立设种子，避免 MM/LTS 内部随机消耗影响后续批次
  set.seed(SEED_BOOT + b)
  idx <- sample(length(x), replace = TRUE)
  boot_b1[b, ] <- fit_all(x[idx], y[idx], seed = SEED_BOOT + b)[, "b1"]
}
se_b1 <- apply(boot_b1, 2, sd, na.rm = TRUE)

tab_phones <- data.frame(
  方法 = rownames(est),
  截距 = round(est[, "b0"], 3),
  斜率 = round(est[, "b1"], 4),
  斜率自助标准误 = round(se_b1, 4),
  斜率t统计量 = round(est[, "b1"] / se_b1, 2)
)
write.csv(tab_phones, file.path(out_dir, "real_phones.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
print(tab_phones)

## ---- 图5: phones 散点 + 各方法拟合直线 ----
png(file.path(fig_dir, "fig5_phones_fit.png"), width = 1500, height = 1100,
    res = 200, family = "yahei")
par(mar = c(4.5, 4.5, 2.5, 1), family = "yahei")
plot(x, y, pch = 16, col = grey(0.3),
     xlab = "年份", ylab = "国际电话呼叫量(百万次)",
     main = "比利时电话呼叫数据: 各方法拟合直线")
xx <- seq(min(x), max(x), length.out = 200)
for (m in methods)
  lines(xx, est[m, "b0"] + est[m, "b1"] * xx, col = box_col[m], lwd = 2)
legend("topleft", bty = "n", lwd = 2, col = box_col[methods], cex = 0.85,
       legend = c("OLS 最小二乘", "LAD 最小绝对偏差", "Huber M",
                  "MM 估计", "LTS 最小截尾平方"))
dev.off()

## ---- 图6: phones 诊断图 ----
fit_ols <- lm(y ~ x)
set.seed(SEED_BOOT)
fit_mm  <- robustbase::lmrob(y ~ x)
png(file.path(fig_dir, "fig6_phones_diag.png"), width = 1800, height = 1200,
    res = 200, family = "yahei")
par(mfrow = c(2, 2), mar = c(4.5, 4.5, 2.5, 1), family = "yahei")

## (a) OLS 残差-拟合值
plot(fitted(fit_ols), resid(fit_ols), pch = 16, col = grey(0.3),
     xlab = "拟合值", ylab = "OLS 残差", main = "(a) OLS: 残差-拟合值")
abline(h = 0, lty = 2, col = "red")

## (b) OLS 残差正态 QQ
qqnorm(resid(fit_ols), pch = 16, col = grey(0.3),
       main = "(b) OLS: 残差正态 Q-Q 图")
qqline(resid(fit_ols), col = "red", lwd = 2)

## (c) MM 残差-拟合值
plot(fitted(fit_mm), resid(fit_mm), pch = 16, col = grey(0.3),
     xlab = "拟合值", ylab = "MM 残差", main = "(c) MM: 残差-拟合值")
abline(h = 0, lty = 2, col = "blue")

## (d) MM 稳健权重(离群点自动降权/剔除)
plot(x, fit_mm$rweights, pch = 16, col = grey(0.3), ylim = c(0, 1.05),
     xlab = "年份", ylab = "稳健权重", main = "(d) MM: 稳健权重")
abline(h = 1, lty = 3, col = grey(0.6))
dev.off()


###### 4.2 datasets::cars —— 干净数据对照 ####

data(cars)
xc <- cars$speed
yc <- cars$dist

SEED_CARS <- 3026
set.seed(SEED_CARS)
est_c <- fit_all(xc, yc, seed = SEED_CARS)
B2 <- 1000
boot_c <- matrix(NA_real_, B2, length(methods),
                 dimnames = list(NULL, methods))
for (b in 1:B2) {
  set.seed(SEED_CARS + b)
  idx <- sample(length(xc), replace = TRUE)
  boot_c[b, ] <- fit_all(xc[idx], yc[idx], seed = SEED_CARS + b)[, "b1"]
}
tab_cars <- data.frame(
  方法 = rownames(est_c),
  截距 = round(est_c[, "b0"], 3),
  斜率 = round(est_c[, "b1"], 4),
  斜率自助标准误 = round(apply(boot_c, 2, sd, na.rm = TRUE), 4)
)
write.csv(tab_cars, file.path(out_dir, "real_cars.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
print(tab_cars)

## ---- 图7: cars 散点 + 拟合直线 ----
png(file.path(fig_dir, "fig7_cars_fit.png"), width = 1500, height = 1100,
    res = 200, family = "yahei")
par(mar = c(4.5, 4.5, 2.5, 1), family = "yahei")
plot(xc, yc, pch = 16, col = grey(0.3),
     xlab = "车速 (mph)", ylab = "刹车距离 (ft)",
     main = "cars 数据(近似干净): 各方法拟合直线")
xx <- seq(min(xc), max(xc), length.out = 200)
for (m in methods)
  lines(xx, est_c[m, "b0"] + est_c[m, "b1"] * xx, col = box_col[m], lwd = 2)
legend("topleft", bty = "n", lwd = 2, col = box_col[methods], cex = 0.85,
       legend = c("OLS 最小二乘", "LAD 最小绝对偏差", "Huber M",
                  "MM 估计", "LTS 最小截尾平方"))
dev.off()


###### 4.3 anscombe 四联图 —— 仅凭 OLS 摘要统计无法发现结构差异 ####

data(anscombe)
png(file.path(fig_dir, "fig8_anscombe.png"), width = 1600, height = 1200,
    res = 200, family = "yahei")
par(mfrow = c(2, 2), mar = c(4.5, 4.5, 2.5, 1), family = "yahei")
for (i in 1:4) {
  xi <- anscombe[[paste0("x", i)]]
  yi <- anscombe[[paste0("y", i)]]
  plot(xi, yi, pch = 16, col = grey(0.3), xlim = c(4, 20), ylim = c(2, 14),
       xlab = "x", ylab = "y",
       main = paste0("Anscombe 数据集 ", i))
  abline(lm(yi ~ xi), col = "red", lwd = 2)
  set.seed(SEED_BOOT + i)
  cf <- coef(robustbase::lmrob(yi ~ xi))
  abline(cf, col = "blue", lwd = 2, lty = 2)
  if (i == 1)
    legend("topleft", bty = "n", lwd = 2, lty = c(1, 2),
           col = c("red", "blue"), cex = 0.9, legend = c("OLS", "MM"))
}
dev.off()

## ---- 复现清单 ----
manifest <- c(
  paste0("生成时间: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("R 版本: ", R.version.string),
  paste0("RNGkind: ", paste(RNGkind(), collapse = " / ")),
  paste0("SEED_BOOT = ", SEED_BOOT, "（第 b 次自助抽样种子 = SEED_BOOT + b）"),
  paste0("SEED_CARS = ", SEED_CARS),
  paste0("自助重复数 B = ", B, " / ", B2),
  paste0("robustbase = ", as.character(packageVersion("robustbase"))),
  paste0("real_phones.csv 的 MD5: ",
         unname(tools::md5sum(file.path(out_dir, "real_phones.csv")))),
  paste0("real_cars.csv   的 MD5: ",
         unname(tools::md5sum(file.path(out_dir, "real_cars.csv"))))
)
writeLines(manifest, file.path(out_dir, "real_manifest.txt"), useBytes = TRUE)

cat("实例分析完成, 结果已写入:", out_dir, "\n")
