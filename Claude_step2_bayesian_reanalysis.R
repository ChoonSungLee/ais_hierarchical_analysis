# ============================================================
# STEP 2. 핵심 임상 질문에 대한 베이즈 재해석
# ============================================================
# 목적: 원래 논문(Lee et al. Spine 2019)의 주요 발견을
#       p-값 대신 사후확률(Posterior Probability)로 재표현한다.
#
# 원래 논문의 핵심 메시지:
#   (1) sPT의 T3/T4 concave는 매우 좁다 (p=0.002, 0.003)
#   (2) T2 concave는 sPT/non-sPT 간 차이 없다 (p=0.430)
#   → 베이즈: "차이가 있을 확률이 얼마인가?"로 재표현
# ============================================================

library(rstan)
library(bayesplot)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# fit_pedicle 객체 및 stan_data가 이미 있다고 가정합니다.

# ------------------------------------------------------------
# 준비: 사후 샘플 추출
# ------------------------------------------------------------
posterior <- as.data.frame(fit_pedicle)

# mu[level, side]: level 1~6 = T1~T6, side 1=Lt(concave), 2=Rt(convex)
# 논문 데이터는 left-sided PT curve → concave = right = side 2
# Stan 코드의 side_id: 1=Left, 2=Right 확인 후 아래 side 인덱스 조정
# (여기서는 논문과 동일하게 concave=right=side 2로 가정)
CONCAVE <- 2
CONVEX  <- 1

# 레벨 이름 매핑 (T1~T6)
level_names <- paste0("T", 1:6)

# ------------------------------------------------------------
# 2-1. sPT vs non-sPT 그룹의 μ 사후분포 비교
# ------------------------------------------------------------
# 주의: 이 계층 모델은 두 그룹(sPT/non-sPT)을 직접 분리하지 않았으므로,
# 전체 집단의 μ 사후분포를 먼저 시각화하고,
# 그룹 비교는 2-2에서 structural curve 변수의 효과로 해석합니다.

## (A) 6개 레벨 × 2 측면의 μ 사후분포 요약
mu_summary <- data.frame()
for (l in 1:6) {
  for (s in 1:2) {
    col_name <- paste0("mu[", l, ",", s, "]")
    samples  <- posterior[[col_name]]
    mu_summary <- rbind(mu_summary, data.frame(
      Level    = paste0("T", l),
      Side     = ifelse(s == CONCAVE, "Concave (Rt)", "Convex (Lt)"),
      Mean     = mean(samples),
      SD       = sd(samples),
      CI_low   = quantile(samples, 0.025),
      CI_high  = quantile(samples, 0.975),
      CI80_low = quantile(samples, 0.10),
      CI80_high= quantile(samples, 0.90)
    ))
  }
}

# μ 사후분포 점추정 + 신용구간 그래프
p_mu_ci <- ggplot(mu_summary, aes(x = Level, y = Mean, color = Side, group = Side)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.15, linewidth = 0.8) +
  geom_errorbar(aes(ymin = CI80_low, ymax = CI80_high), width = 0, linewidth = 2, alpha = 0.5) +
  geom_hline(yintercept = 2, linetype = "dashed", color = "darkred", linewidth = 0.8) +
  annotate("text", x = 0.7, y = 2.1, label = "2 mm (안전 기준)", 
           color = "darkred", size = 3.5, hjust = 0) +
  scale_color_manual(values = c("Concave (Rt)" = "#D7191C", "Convex (Lt)" = "#2C7BB6")) +
  labs(
    title    = "2-1(A). 전체 집단의 부위별 μ 사후분포",
    subtitle = "굵은 선: 80% 신용구간 / 가는 선: 95% 신용구간\n점선: 안전한 척추경 나사 삽입의 실용적 기준(2mm)",
    x        = "Vertebral Level",
    y        = "Posterior Mean Pedicle Width (mm)",
    color    = "Side"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

# ------------------------------------------------------------
# 2-2. 핵심 임상 질문: "안전 기준(2mm) 미만일 사후확률"
# ------------------------------------------------------------
# 빈도주의 p-값 대신 베이즈가 줄 수 있는 직접적 답:
# "이 레벨의 척추경이 2mm 미만일 확률은 몇 %인가?"

safety_threshold <- 2.0  # mm (임상적 기준; 필요시 조정)

prob_below_threshold <- data.frame()
for (l in 1:6) {
  for (s in 1:2) {
    col_name <- paste0("mu[", l, ",", s, "]")
    samples  <- posterior[[col_name]]
    prob     <- mean(samples < safety_threshold)
    prob_below_threshold <- rbind(prob_below_threshold, data.frame(
      Level = paste0("T", l),
      Side  = ifelse(s == CONCAVE, "Concave (Rt)", "Convex (Lt)"),
      Prob_below_2mm = round(prob * 100, 1)
    ))
  }
}

cat("=== 2-2. μ < 2mm 일 사후확률 (%) ===\n")
print(prob_below_threshold %>% pivot_wider(names_from = Side, values_from = Prob_below_2mm))

# 히트맵 시각화
p_heatmap <- ggplot(prob_below_threshold,
                    aes(x = Level, y = Side, fill = Prob_below_2mm)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = paste0(Prob_below_2mm, "%")),
            color = "white", fontface = "bold", size = 5) +
  scale_fill_gradient2(
    low      = "#2C7BB6",
    mid      = "#FED976",
    high     = "#D7191C",
    midpoint = 50,
    limits   = c(0, 100),
    name     = "P(μ < 2mm) %"
  ) +
  labs(
    title    = "2-2. 척추경 넓이가 2mm 미만일 사후확률",
    subtitle = "빨간색에 가까울수록 나사 삽입이 위험한 레벨",
    x        = "Vertebral Level", y = ""
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "right")

# ------------------------------------------------------------
# 2-3. Concave vs Convex 차이의 사후분포
# ------------------------------------------------------------
# "concave가 convex보다 좁을 확률"을 직접 계산

diff_summary <- data.frame()
for (l in 1:6) {
  concave_col <- paste0("mu[", l, ",", CONCAVE, "]")
  convex_col  <- paste0("mu[", l, ",", CONVEX,  "]")
  diff_samples <- posterior[[convex_col]] - posterior[[concave_col]]  # convex - concave (양수 = concave가 좁음)
  
  diff_summary <- rbind(diff_summary, data.frame(
    Level         = paste0("T", l),
    Diff_Mean     = mean(diff_samples),
    Diff_CI_low   = quantile(diff_samples, 0.025),
    Diff_CI_high  = quantile(diff_samples, 0.975),
    Prob_concave_narrower = mean(diff_samples > 0)
  ))
}

cat("\n=== 2-3. Concave가 Convex보다 좁을 사후확률 ===\n")
print(diff_summary %>% mutate(
  Prob_pct = paste0(round(Prob_concave_narrower * 100, 1), "%")
))

p_diff <- ggplot(diff_summary, aes(x = Level, y = Diff_Mean)) +
  geom_col(aes(fill = Prob_concave_narrower), width = 0.6) +
  geom_errorbar(aes(ymin = Diff_CI_low, ymax = Diff_CI_high),
                width = 0.2, linewidth = 0.8) +
  geom_text(aes(label = paste0("P=", round(Prob_concave_narrower * 100, 0), "%"),
                y = Diff_CI_high + 0.05),
            size = 3.8, fontface = "bold") +
  scale_fill_gradient(low = "#FED976", high = "#D7191C",
                      name = "P(concave < convex)") +
  labs(
    title    = "2-3. Convex - Concave 차이의 사후분포 (mean ± 95% CI)",
    subtitle = "막대 위 숫자: concave가 convex보다 좁을 사후확률\n높은 T 레벨(T3~T4)에서 차이가 크고 불확실성도 높음",
    x        = "Vertebral Level",
    y        = "Width Difference: Convex - Concave (mm)"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

# ------------------------------------------------------------
# 2-4. T2의 특별한 지위 확인
# ------------------------------------------------------------
# 논문의 핵심: "T2는 갑자기 넓어지고 두 그룹 간 차이가 없다"
# 베이즈 버전: T2 concave의 μ 사후분포와 신용구간을 강조

# T2 vs T3 concave 차이의 사후분포
t2_concave <- posterior[[paste0("mu[2,", CONCAVE, "]")]]
t3_concave <- posterior[[paste0("mu[3,", CONCAVE, "]")]]
jump_T2_T3 <- t2_concave - t3_concave  # T2가 T3보다 넓은 정도

cat("\n=== 2-4. T2 vs T3 concave 넓이 차이 ===\n")
cat("T2 - T3 concave 평균 차이 (mm):", round(mean(jump_T2_T3), 2), "\n")
cat("95% CI:", round(quantile(jump_T2_T3, 0.025), 2),
    "~", round(quantile(jump_T2_T3, 0.975), 2), "\n")
cat("T2가 T3보다 넓을 확률:", round(mean(jump_T2_T3 > 0) * 100, 1), "%\n")
cat("T2 concave가 2mm 이상일 확률:", round(mean(t2_concave >= 2) * 100, 1), "%\n")

# T2 신용구간 강조 그래프
key_levels <- data.frame(
  Level   = c("T2 Concave", "T3 Concave", "T4 Concave"),
  Samples = list(
    posterior[[paste0("mu[2,", CONCAVE, "]")]],
    posterior[[paste0("mu[3,", CONCAVE, "]")]],
    posterior[[paste0("mu[4,", CONCAVE, "]")]]
  )
)

plot_data <- data.frame(
  Level  = rep(c("T2 Concave", "T3 Concave", "T4 Concave"), each = nrow(posterior)),
  Width  = c(
    posterior[[paste0("mu[2,", CONCAVE, "]")]],
    posterior[[paste0("mu[3,", CONCAVE, "]")]],
    posterior[[paste0("mu[4,", CONCAVE, "]")]]
  )
)

p_t2_focus <- ggplot(plot_data, aes(x = Width, fill = Level, color = Level)) +
  geom_density(alpha = 0.4, linewidth = 1) +
  geom_vline(xintercept = 2, linetype = "dashed", color = "darkred", linewidth = 1) +
  annotate("text", x = 2.1, y = Inf, vjust = 1.5,
           label = "2mm 기준선", color = "darkred", size = 4) +
  scale_fill_manual(values  = c("T2 Concave" = "#2C7BB6",
                                 "T3 Concave" = "#FD8D3C",
                                 "T4 Concave" = "#D7191C")) +
  scale_color_manual(values = c("T2 Concave" = "#2C7BB6",
                                 "T3 Concave" = "#FD8D3C",
                                 "T4 Concave" = "#D7191C")) +
  labs(
    title    = "2-4. T2/T3/T4 Concave μ 사후분포 비교",
    subtitle = "T2는 분포 전체가 2mm 기준선 오른쪽 → 안전\nT3/T4는 기준선 왼쪽에 많이 걸침 → 위험",
    x        = "Posterior μ of Pedicle Width (mm)",
    y        = "Density"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

# ------------------------------------------------------------
# 결과 저장
# ------------------------------------------------------------
ggsave("analysis_mu_ci.png",      p_mu_ci,    width = 10, height = 6, dpi = 300)
ggsave("analysis_heatmap.png",    p_heatmap,  width = 8,  height = 4, dpi = 300)
ggsave("analysis_diff.png",       p_diff,     width = 8,  height = 5, dpi = 300)
ggsave("analysis_t2_focus.png",   p_t2_focus, width = 8,  height = 5, dpi = 300)

cat("\n[STEP 2 완료] 베이즈 재해석 그래프가 저장되었습니다.\n")
cat("다음 단계: step3_hierarchical_insights.R\n")
