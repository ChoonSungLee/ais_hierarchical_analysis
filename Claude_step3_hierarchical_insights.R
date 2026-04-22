# ============================================================
# STEP 3. 계층 모델만이 줄 수 있는 새로운 임상 인사이트
# ============================================================
# 목적: 빈도주의 분석에서는 불가능했던 세 가지를 보여준다.
#
#   3-1. 환자별 random effect (η_p) 시각화
#        → "척추경이 전반적으로 좁은 환자 유형"을 데이터로 포착
#   3-2. Shrinkage 효과 시각화
#        → 데이터가 적은 환자의 추정치가 전체 평균으로 당겨지는 현상
#   3-3. 새 환자 예측: 임상적 핵심 결과
#        → "처음 보는 sPT 환자의 T2/T3/T4 척추경이
#           안전 범위일 사후 예측 확률"
# ============================================================

library(rstan)
library(bayesplot)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# fit_pedicle, stan_data, df_long 객체가 있다고 가정합니다.

CONCAVE <- 2
CONVEX  <- 1

posterior <- as.data.frame(fit_pedicle)

# ------------------------------------------------------------
# 3-1. 환자별 Random Effect (η_p) 시각화
# ------------------------------------------------------------
# patient_effect[p]: p번 환자의 전반적 척추경 수준이
# 집단 평균보다 얼마나 다른지를 나타내는 값

# 사후 중앙값 + 95% CI 추출
n_patient <- stan_data$N_patient
patient_effects <- data.frame(
  Patient  = 1:n_patient,
  Median   = sapply(1:n_patient, function(p)
    median(posterior[[paste0("patient_effect[", p, "]")]])),
  CI_low   = sapply(1:n_patient, function(p)
    quantile(posterior[[paste0("patient_effect[", p, "]")]], 0.025)),
  CI_high  = sapply(1:n_patient, function(p)
    quantile(posterior[[paste0("patient_effect[", p, "]")]], 0.975))
)

# 중앙값 기준으로 정렬
patient_effects <- patient_effects %>%
  arrange(Median) %>%
  mutate(Rank = row_number(),
         Risk_group = case_when(
           Median < -0.5 ~ "좁은 척추경 환자 (η < -0.5mm)",
           Median >  0.5 ~ "넓은 척추경 환자 (η > +0.5mm)",
           TRUE           ~ "평균 범위 환자"
         ))

p_patient_effect <- ggplot(patient_effects, aes(x = Rank, y = Median, color = Risk_group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +
  geom_hline(yintercept = c(-0.5, 0.5), linetype = "dotted",
             color = "gray60", linewidth = 0.6) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high),
                width = 0, alpha = 0.3, linewidth = 0.5) +
  geom_point(size = 1.5) +
  scale_color_manual(values = c(
    "좁은 척추경 환자 (η < -0.5mm)" = "#D7191C",
    "평균 범위 환자"                  = "#4575B4",
    "넓은 척추경 환자 (η > +0.5mm)"  = "#1A9641"
  )) +
  labs(
    title    = "3-1. 환자별 Random Effect (η_p) 사후분포",
    subtitle = paste0(
      "전체 ", n_patient, "명 환자를 η_p 크기 순으로 정렬\n",
      "빨간 환자: 집단 평균보다 전반적으로 척추경이 좁음 → 수술 시 주의 필요"
    ),
    x     = "환자 순위 (η_p 오름차순)",
    y     = "η_p: 집단 평균 대비 개인 편차 (mm)",
    color = "환자 분류"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 10))

# η_p 분포의 요약 통계 출력
cat("=== 3-1. 환자별 Random Effect 요약 ===\n")
cat("η_p < -0.5 mm (척추경 좁음 위험군):", sum(patient_effects$Median < -0.5), "명\n")
cat("η_p > +0.5 mm (척추경 넓음):",        sum(patient_effects$Median >  0.5), "명\n")
cat("σ_patient 사후 중앙값:", round(median(posterior$sigma_patient), 3), "mm\n\n")

# sigma_patient 사후분포
p_sigma_patient <- ggplot(data.frame(sigma = posterior$sigma_patient), aes(x = sigma)) +
  geom_density(fill = "#4575B4", alpha = 0.5, color = "#2C7BB6", linewidth = 1) +
  geom_vline(xintercept = median(posterior$sigma_patient),
             color = "darkblue", linetype = "dashed") +
  labs(
    title    = "3-1(B). σ_patient 사후분포 (환자 간 변동성)",
    subtitle = "이 값이 크면 환자마다 척추경 수준의 차이가 크다는 의미",
    x        = "σ_patient (mm)", y = "Density"
  ) +
  theme_minimal(base_size = 12)

# ------------------------------------------------------------
# 3-2. Shrinkage 효과 시각화
# ------------------------------------------------------------
# 각 환자의 "날 것 데이터 평균(no-pooling)"과
# "계층 모델 추정치(partial pooling)"를 비교

# 환자별 실제 측정 평균 계산
raw_patient_mean <- df_long %>%
  mutate(patient_idx = as.integer(factor(ID))) %>%
  group_by(patient_idx) %>%
  summarise(
    raw_mean = mean(width, na.rm = TRUE),
    n_obs    = n(),
    .groups  = "drop"
  )

# 계층 모델의 환자별 mu_patient 사후 중앙값
# mu_patient[n]은 각 관측치별이므로, 환자별로 평균
# patient_id와 mu_patient를 매핑
n_obs <- stan_data$N
patient_id_vec <- stan_data$patient_id

mu_patient_posterior_median <- sapply(1:n_obs, function(n)
  median(posterior[[paste0("mu_patient[", n, "]")]]))

shrinkage_df <- data.frame(
  patient_idx   = patient_id_vec,
  mu_patient_est = mu_patient_posterior_median
) %>%
  group_by(patient_idx) %>%
  summarise(hierarchical_mean = mean(mu_patient_est), .groups = "drop") %>%
  left_join(raw_patient_mean, by = "patient_idx") %>%
  mutate(
    global_mean = mean(posterior$`mu[1,1]`) # 전체 평균의 대리값 (T1 Lt 기준)
  )

# Shrinkage plot
p_shrinkage <- ggplot(shrinkage_df, aes(x = raw_mean, y = hierarchical_mean)) +
  geom_abline(slope = 1, intercept = 0, color = "gray70",
              linetype = "dashed", linewidth = 1) +
  geom_hline(yintercept = mean(shrinkage_df$hierarchical_mean),
             color = "red", linetype = "dotted", linewidth = 0.8) +
  geom_point(aes(size = n_obs, alpha = n_obs), color = "#2C7BB6") +
  geom_segment(aes(xend = raw_mean, yend = raw_mean),
               arrow = arrow(length = unit(0.08, "cm")),
               color = "gray40", alpha = 0.4, linewidth = 0.3) +
  scale_size_continuous(name = "관측치 수", range = c(1, 5)) +
  scale_alpha_continuous(name = "관측치 수", range = c(0.3, 0.9)) +
  labs(
    title    = "3-2. Shrinkage 효과",
    subtitle = "점선(y=x): shrinkage 없을 때\n빨간선: 전체 평균\n점이 y=x 선에서 빨간선 쪽으로 이동한 거리 = shrinkage 양",
    x        = "No-pooling 추정 (환자별 데이터 평균, mm)",
    y        = "Partial-pooling 추정 (계층 모델, mm)"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right")

# ------------------------------------------------------------
# 3-3. 새 환자 예측: 핵심 임상 결과
# ------------------------------------------------------------
# width_new_patient[l, s]: 처음 보는 환자의 레벨별 척추경 예측 분포

new_patient_samples <- data.frame()
for (l in 1:6) {
  for (s in 1:2) {
    col_name <- paste0("width_new_patient[", l, ",", s, "]")
    samples  <- posterior[[col_name]]
    new_patient_samples <- rbind(new_patient_samples, data.frame(
      Level         = paste0("T", l),
      Side          = ifelse(s == CONCAVE, "Concave", "Convex"),
      Width         = samples
    ))
  }
}

## (A) 새 환자 예측 분포 (violin + 95% CI)
pred_summary <- new_patient_samples %>%
  group_by(Level, Side) %>%
  summarise(
    Median   = median(Width),
    CI_low   = quantile(Width, 0.025),
    CI_high  = quantile(Width, 0.975),
    P_below2 = mean(Width < 2) * 100,
    P_below0 = mean(Width < 0) * 100,   # 음수 예측(비현실적) 비율 점검
    .groups  = "drop"
  )

cat("=== 3-3. 새 환자 예측 요약 ===\n")
print(pred_summary %>% select(Level, Side, Median, CI_low, CI_high, P_below2))

p_new_patient <- ggplot(
  new_patient_samples %>% filter(Width > 0),  # 음수 제거
  aes(x = Level, y = Width, fill = Side)
) +
  geom_violin(trim = TRUE, alpha = 0.6, scale = "width", position = position_dodge(0.8)) +
  geom_boxplot(width = 0.15, position = position_dodge(0.8),
               outlier.shape = NA, alpha = 0.8) +
  geom_hline(yintercept = 2, linetype = "dashed", color = "darkred", linewidth = 1) +
  annotate("text", x = 0.6, y = 2.15,
           label = "2 mm 안전 기준", color = "darkred", size = 3.5) +
  scale_fill_manual(values = c("Concave" = "#D7191C", "Convex" = "#2C7BB6")) +
  labs(
    title    = "3-3(A). 새 환자(Future Patient) 예측 사후분포",
    subtitle = "바이올린 넓이: 예측 불확실성의 크기\n점선 아래 영역: 나사 삽입 위험 구간",
    x        = "Vertebral Level",
    y        = "Predicted Pedicle Width (mm)",
    fill     = "Side"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

## (B) 핵심 임상 메시지: "새 환자에서 2mm 미만일 예측 확률"
p_risk_table <- ggplot(
  pred_summary,
  aes(x = Level, y = Side, fill = P_below2)
) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = paste0(round(P_below2, 0), "%")),
            color = "white", fontface = "bold", size = 6) +
  scale_fill_gradient2(
    low      = "#2166AC",
    mid      = "#F7F7F7",
    high     = "#D7191C",
    midpoint = 30,
    limits   = c(0, 100),
    name     = "P(width < 2mm) %"
  ) +
  labs(
    title    = "3-3(B). 새 환자에서 척추경이 2mm 미만일 예측 확률",
    subtitle = "이 수치가 높을수록 해당 레벨에서 나사 삽입을 피하거나\n특수 기법(extrapedicular 등)을 고려해야 할 가능성이 높음",
    x        = "Vertebral Level", y = "Side"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "right",
        axis.text = element_text(size = 12))

## (C) mcmc_areas를 이용한 레벨별 예측 분포 (concave 집중)
concave_samples_wide <- new_patient_samples %>%
  filter(Side == "Concave") %>%
  select(Level, Width) %>%
  group_by(Level) %>%
  mutate(iter = row_number()) %>%
  pivot_wider(names_from = Level, values_from = Width) %>%
  select(-iter)

pred_matrix_concave <- as.matrix(concave_samples_wide)

p_areas <- mcmc_areas(
  pred_matrix_concave,
  prob       = 0.80,
  prob_outer = 0.95,
  point_est  = "median"
) +
  geom_vline(xintercept = 2, color = "darkred",
             linetype = "dashed", linewidth = 1) +
  labs(
    title    = "3-3(C). Concave 척추경 예측 분포 (레벨별)",
    subtitle = "진한 영역: 80% 예측구간 / 선: 95% 예측구간\n점선: 2mm 안전 기준",
    x        = "Predicted Pedicle Width for a New Patient (mm)"
  ) +
  theme_minimal(base_size = 12)

# ------------------------------------------------------------
# 3-4. 종합 요약: 베이즈 vs 빈도주의 비교표 출력
# ------------------------------------------------------------
cat("\n")
cat("===========================================================\n")
cat("   STEP 3 종합: 계층 베이즈 모델의 추가 인사이트 요약\n")
cat("===========================================================\n")
cat("\n[1] 환자 개인 수준 정보\n")
cat("  - σ_patient:", round(median(posterior$sigma_patient), 2),
    "mm → 환자 간 척추경 수준의 개인차 크기\n")
cat("  - '위험 환자' (η < -0.5mm):",
    sum(patient_effects$Median < -0.5), "명 /", n_patient, "명 중\n")

cat("\n[2] T2 Safety 확인 (새 환자 예측 기준)\n")
t2_concave_new <- new_patient_samples %>%
  filter(Level == "T2", Side == "Concave")
cat("  - T2 Concave 예측 중앙값:",
    round(median(t2_concave_new$Width), 2), "mm\n")
cat("  - T2 Concave P(width < 2mm):",
    round(mean(t2_concave_new$Width < 2) * 100, 1), "%\n")

cat("\n[3] T3/T4 Risk 확인 (새 환자 예측 기준)\n")
for (lv in c("T3", "T4")) {
  d <- new_patient_samples %>% filter(Level == lv, Side == "Concave")
  cat(sprintf("  - %s Concave P(width < 2mm): %.1f%%\n",
              lv, mean(d$Width < 2) * 100))
}

# ------------------------------------------------------------
# 결과 저장
# ------------------------------------------------------------
ggsave("insight_patient_effect.png",  p_patient_effect, width = 10, height = 6, dpi = 300)
ggsave("insight_sigma_patient.png",   p_sigma_patient,  width = 7,  height = 4, dpi = 300)
ggsave("insight_shrinkage.png",       p_shrinkage,      width = 8,  height = 6, dpi = 300)
ggsave("insight_new_patient_violin.png", p_new_patient, width = 10, height = 6, dpi = 300)
ggsave("insight_risk_heatmap.png",    p_risk_table,     width = 8,  height = 4, dpi = 300)
ggsave("insight_pred_areas.png",      p_areas,          width = 8,  height = 5, dpi = 300)

cat("\n[STEP 3 완료] 모든 인사이트 그래프가 저장되었습니다.\n")
cat("=== 전체 분석 파이프라인 완료 ===\n")
