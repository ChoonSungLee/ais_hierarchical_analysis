# ============================================================
# STEP 1. 모델 검증 및 진단 (Model Validation & Diagnostics)
# ============================================================
# 목적: 사후분포를 신뢰하기 전에 MCMC가 제대로 수렴했는지,
#       모델이 실제 데이터를 잘 재현하는지 확인한다.
# ============================================================

library(rstan)
library(bayesplot)
library(ggplot2)
library(dplyr)
library(loo)        # LOO-CV 계산용
library(patchwork)  # 그래프 배치용

# fit_pedicle 객체가 이미 있다고 가정합니다.
# (없다면 main_hierarchical.R을 먼저 실행하세요)

# ------------------------------------------------------------
# 1-1. MCMC 수렴 진단 (Convergence Diagnostics)
# ------------------------------------------------------------

## (A) Rhat 확인: 모든 모수가 1.01 미만이어야 안전
rhat_vals <- summary(fit_pedicle)$summary[, "Rhat"]
cat("=== Rhat 요약 ===\n")
cat("최댓값:", max(rhat_vals, na.rm = TRUE), "\n")
cat("1.01 초과 모수 수:", sum(rhat_vals > 1.01, na.rm = TRUE), "\n\n")

# Rhat 히스토그램
rhat_df <- data.frame(rhat = rhat_vals[!is.na(rhat_vals)])
p_rhat <- ggplot(rhat_df, aes(x = rhat)) +
  geom_histogram(binwidth = 0.001, fill = "#2C7BB6", color = "white") +
  geom_vline(xintercept = 1.01, color = "red", linetype = "dashed", linewidth = 1) +
  labs(
    title = "1-1(A). Rhat Distribution",
    subtitle = "모든 값이 점선(1.01) 왼쪽에 있어야 수렴 성공",
    x = "Rhat", y = "Count"
  ) +
  theme_minimal(base_size = 12)

## (B) n_eff (유효 샘플 수) 확인: 총 샘플(4000)의 10% 이상 권장
neff_vals <- summary(fit_pedicle)$summary[, "n_eff"]
cat("=== n_eff 요약 ===\n")
cat("최솟값:", min(neff_vals, na.rm = TRUE), "\n")
cat("400 미만 모수 수:", sum(neff_vals < 400, na.rm = TRUE), "\n\n")

## (C) Trace Plot: 핵심 모수들의 체인 혼합 확인
# mu (T1~T6, Lt/Rt) + sigma 모수들만 추출
posterior_array <- as.array(fit_pedicle)

p_trace_mu <- mcmc_trace(
  posterior_array,
  pars = c("mu[1,1]", "mu[1,2]", "mu[3,1]", "mu[3,2]", "mu[6,1]", "mu[6,2]"),
  facet_args = list(ncol = 2)
) +
  labs(title = "1-1(C). Trace Plot: mu (T1, T3, T6 × Lt/Rt)",
       subtitle = "4개 체인이 서로 잘 겹쳐야 수렴 성공") +
  theme_minimal(base_size = 11)

p_trace_sigma <- mcmc_trace(
  posterior_array,
  pars = c("sigma_obs", "sigma_patient")
) +
  labs(title = "1-1(C). Trace Plot: sigma 모수들") +
  theme_minimal(base_size = 11)

## (D) Pairs Plot: 주요 모수 간 상관관계 (funnel 확인)
p_pairs <- mcmc_pairs(
  posterior_array,
  pars = c("sigma_patient", "mu[1,1]", "mu[3,1]", "mu[6,1]"),
  off_diag_args = list(size = 0.3, alpha = 0.3)
)

# ------------------------------------------------------------
# 1-2. 사후예측검증 (Posterior Predictive Check, PPC)
# ------------------------------------------------------------
# 모델이 만들어낸 데이터(width_rep)가 실제 데이터(width)와
# 얼마나 비슷한지 시각적으로 비교

# width_rep 추출 (행: MCMC 샘플, 열: 관측치)
width_rep_matrix <- as.matrix(fit_pedicle, pars = "width_rep")
y_obs <- stan_data$width

## (A) 전체 밀도 비교 (100개 샘플 오버레이)
p_ppc_dens <- ppc_dens_overlay(
  y     = y_obs,
  yrep  = width_rep_matrix[1:100, ]
) +
  labs(
    title    = "1-2(A). PPC: Density Overlay",
    subtitle = "파란 선(예측)이 검은 선(실제)을 잘 감싸면 모델 적합",
    x        = "Pedicle Width (mm)"
  ) +
  theme_minimal(base_size = 12)

## (B) 통계량 비교: mean, sd, min, max
p_ppc_stat <- ppc_stat_2d(
  y    = y_obs,
  yrep = width_rep_matrix,
  stat = c("mean", "sd")
) +
  labs(
    title    = "1-2(B). PPC: Mean vs SD",
    subtitle = "점(실제값)이 분포 중앙에 있으면 모델이 평균과 분산을 잘 포착"
  ) +
  theme_minimal(base_size = 12)

## (C) 극단값 확인: 실제 최솟값이 예측 범위 안에 있는지
p_ppc_min <- ppc_stat(y_obs, width_rep_matrix, stat = "min") +
  labs(title = "1-2(C). PPC: Minimum Value",
       subtitle = "극단적으로 좁은 척추경(0에 가까운 값)을 모델이 재현하는지 확인") +
  theme_minimal(base_size = 12)

## (D) 분위수 비교
p_ppc_ecdf <- ppc_ecdf_overlay(
  y    = y_obs,
  yrep = width_rep_matrix[1:50, ]
) +
  labs(
    title    = "1-2(D). PPC: ECDF Overlay",
    subtitle = "누적분포 전체에서 예측(파랑)과 실제(검정)의 일치 확인",
    x        = "Pedicle Width (mm)"
  ) +
  theme_minimal(base_size = 12)

# ------------------------------------------------------------
# 1-3. LOO-CV (Leave-One-Out Cross Validation)
# ------------------------------------------------------------
# 모델의 예측력을 수치로 평가; 나중에 모델 비교 시 기준값이 됨

log_lik_matrix <- extract_log_lik(fit_pedicle, parameter_name = "log_lik",
                                   merge_chains = FALSE)
loo_result <- loo(log_lik_matrix, relative_eff = relative_eff(log_lik_matrix))

cat("=== LOO-CV 결과 ===\n")
print(loo_result)

# Pareto-k 진단: k > 0.7이면 해당 관측치가 영향력이 너무 큰 것
p_pareto <- plot(loo_result, diagnostic = "k", label_points = FALSE) +
  labs(
    title    = "1-3. LOO-CV Pareto-k Diagnostics",
    subtitle = "k > 0.7(빨간 선) 관측치가 없어야 LOO 추정이 신뢰 가능",
    y        = "Pareto Shape k"
  ) +
  theme_minimal(base_size = 12)

# ------------------------------------------------------------
# 결과 저장
# ------------------------------------------------------------
ggsave("diag_rhat.png",         p_rhat,       width = 7,  height = 4, dpi = 300)
ggsave("diag_trace_mu.png",     p_trace_mu,   width = 10, height = 7, dpi = 300)
ggsave("diag_trace_sigma.png",  p_trace_sigma,width = 8,  height = 4, dpi = 300)
ggsave("diag_ppc_density.png",  p_ppc_dens,   width = 8,  height = 5, dpi = 300)
ggsave("diag_ppc_stat2d.png",   p_ppc_stat,   width = 7,  height = 6, dpi = 300)
ggsave("diag_ppc_min.png",      p_ppc_min,    width = 7,  height = 4, dpi = 300)
ggsave("diag_ppc_ecdf.png",     p_ppc_ecdf,   width = 8,  height = 5, dpi = 300)
ggsave("diag_pareto_k.png",     p_pareto,     width = 8,  height = 5, dpi = 300)

cat("\n[STEP 1 완료] 모든 진단 그래프가 저장되었습니다.\n")
cat("다음 단계: step2_bayesian_reanalysis.R\n")
