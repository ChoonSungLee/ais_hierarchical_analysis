
source("data_hierarchical.R")

library(rstan)
library(bayesplot)  # 시각화 패키지 추가
library(ggplot2)    # 베이즈플롯은 ggplot2 기반

# 병렬 처리 설정 
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# 모델 적합 
fit_pedicle <- stan(
  file = "hierarchical_pedicle_width.stan", 
  data = stan_data,
  chains = 4, 
  iter = 2000, 
  warmup = 1000, 
  seed = 42
)

# 결과 확인 
print(fit_pedicle, pars = c("mu", "sigma_obs", "sigma_patient"))

# 1. GQB에서 만든 width_rep(복제 데이터) 추출
posterior_samples <- as.matrix(fit_pedicle)
width_rep_samples <- posterior_samples[, grepl("width_rep", colnames(posterior_samples))]

# 2. 사후예측검증(PPC) 시각화
# 실제 데이터(stan_data$width)와 모델이 예측한 분포를 비교
ppc_dens_overlay(y = stan_data$width, yrep = width_rep_samples[1:100, ])

# 2. GQB 결과물 추출
# width_new_patient: 새로운 환자의 [레벨, 측면]별 예측 사후분포
# Stan 결과에서 해당 모수만 추출
post_pred_matrix <- as.matrix(fit_pedicle, pars = "width_new_patient")

# 3. 데이터 구조 정리 (T1~T9 레벨별 시각화를 위해)
# 예: T1-Lt, T1-Rt 등을 하나의 레벨(T1)로 묶거나 각각 표시할 수 있음.
# 여기서는 '레벨'의 차이를 보여주는 데 집중함.

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


# 민감도 분석용 R 코드
# 1. 두 가지 시나리오 설정
# Case 1: 기존 설정 (평균 5, 표준편차 3)
# Case 2: 느슨한 설정 (평균 5, 표준편차 10) - 데이터의 영향력을 더 크게 만듬

# Case 1 실행
fit_orig <- stan(
  file = "hierarchical_pedicle_width.stan",
  data = stan_data, # 기존 레시피 코드의 데이터
  iter = 2000, chains = 4, seed = 42
)

# Case 2를 위해 데이터 리스트의 사전분포 파라미터만 수정 (Stan 파일 내에서 변수화했을 경우)
# 만약 Stan 파일에 숫자를 직접 적으셨다면, Stan 파일을 복사해 숫자를 수정 후 실행해야 함.
# 여기서는 시각적 비교를 위해 fit_orig와 fit_weak(가칭)가 있다고 가정함.

# 2. 결과 추출 및 비교 데이터프레임 생성
library(dplyr)
library(tidyr)

# 주요 모수(mu[1,1] ~ mu[9,2])의 평균값 추출 예시
summ_orig <- as.data.frame(summary(fit_orig, pars = "mu")$summary) %>% 
  mutate(Parameter = rownames(.), Model = "Informative (SD=3)")

# (동일한 방식으로 fit_weak에서도 추출)
# summ_weak <- as.data.frame(summary(fit_weak, pars = "mu")$summary) %>% 
#   mutate(Parameter = rownames(.), Model = "Weakly Informative (SD=10)")

# 3. 시각적 비교 (Posterior Overlay)
library(bayesplot)

# 특정 레벨(예: T1)의 사후분포가 사전분포 변경에 따라 얼마나 변하는지 확인
posterior_orig <- as.matrix(fit_orig, pars = "mu[1,1]")
# posterior_weak <- as.matrix(fit_weak, pars = "mu[1,1]")

# 두 분포를 겹쳐서 그리기
# m_list <- list(Informative = posterior_orig, Weakly_Informative = posterior_weak)
# ppc_dens_overlay(y = as.vector(posterior_orig), yrep = posterior_weak[1:100, ]) +
#   labs(title = "Sensitivity Analysis: Posterior of mu[1,1]")


