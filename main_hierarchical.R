library(rstan)
library(bayesplot)  # 시각화 패키지 추가
library(ggplot2)    # 베이즈플롯은 ggplot2 기반이므로 함께 로드하면 좋습니다.

# 병렬 처리 설정 [cite: 783]
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# 모델 적합 [cite: 592, 1009]
fit_pedicle <- stan(
  file = "hierarchical_pedicle_width.stan", 
  data = stan_data,
  chains = 4, 
  iter = 2000, 
  warmup = 1000, 
  seed = 42
)

# 결과 확인 [cite: 607, 1012]
print(fit_pedicle, pars = c("mu", "sigma_obs", "sigma_patient"))

# 1. GQB에서 만든 width_rep(복제 데이터) 추출
posterior_samples <- as.matrix(fit_pedicle)
width_rep_samples <- posterior_samples[, grepl("width_rep", colnames(posterior_samples))]

# 2. 사후예측검증(PPC) 시각화
# 실제 데이터(stan_data$width)와 모델이 예측한 분포를 비교합니다.
ppc_dens_overlay(y = stan_data$width, yrep = width_rep_samples[1:100, ])

# 2. GQB 결과물 추출
# width_new_patient: 새로운 환자의 [레벨, 측면]별 예측 사후분포
# Stan 결과에서 해당 모수만 추출합니다.
post_pred_matrix <- as.matrix(fit_pedicle, pars = "width_new_patient")

# 3. 데이터 구조 정리 (T1~T9 레벨별 시각화를 위해)
# 예: T1-Lt, T1-Rt 등을 하나의 레벨(T1)로 묶거나 각각 표시할 수 있습니다.
# 여기서는 '레벨'의 차이를 보여주는 데 집중합니다.

# bayesplot의 mcmc_areas 함수 사용: 사후분포의 밀도와 신용구간을 동시에 보여줌
p <- mcmc_areas(
  post_pred_matrix,
  prob = 0.8,       # 안쪽 짙은 영역 (80% 신용구간)
  prob_outer = 0.95, # 바깥쪽 선 (95% 신용구간)
  point_est = "median"
) +
  # 논문용 스타일 입히기
  labs(
    title = "Posterior Predictive Distributions of Pedicle Width",
    subtitle = "Predicted widths for a new patient across vertebral levels (T1-T9)",
    x = "Pedicle Width (mm)",
    y = "Vertebral Level & Side"
  ) +
  theme_minimal(base_family = "sans") + # 학술지에 적합한 폰트 설정
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(face = "bold", size = 14),
    panel.grid.minor = element_blank()
  )

# 4. 그래프 출력
print(p)

# 5. 논문용 고해상도 이미지 저장 (PDF나 TIFF 추천)
# ggsave("figure_pedicle_width_prediction.pdf", plot = p, width = 8, height = 10, units = "in", dpi = 300)