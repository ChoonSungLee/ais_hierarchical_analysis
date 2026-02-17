data {
  int<lower=1> N;                // 전체 관측치 개수
  int<lower=1> N_patient;        // 전체 환자 수
  int<lower=1> N_level;          // 척추 레벨 수 (T1~T9 등)
  int<lower=1> N_side;           // 좌우 구분 (2)
  
  vector[N] width;               // 종속변수: 척추경 넓이
  int<lower=1> patient_id[N];    // 각 관측치의 환자 번호
  int<lower=1> level_id[N];      // 각 관측치의 레벨 번호
  int<lower=1> side_id[N];       // 각 관측치의 좌우 번호
}

parameters {
  matrix[N_level, N_side] mu;    // 레벨/측면별 평균 넓이
  real<lower=0> sigma_obs;       // 관측 오차 (Residual)
  
  // 계층적 구조를 위한 파라미터
  vector[N_patient] patient_effect_raw; // 비중심화 파라미터화용
  real<lower=0> sigma_patient;          // 환자 간 편차 (초모수)
}

transformed parameters {
  // 환자별 랜덤 효과 계산 (평균 0, 표준편차 sigma_patient인 분포를 따름)
  vector[N_patient] patient_effect;
  patient_effect = patient_effect_raw * sigma_patient;
}

model {
  // 1. 사전분포 및 초사전분포 (Priors & Hyperpriors)
  to_vector(mu) ~ normal(5, 3);
  sigma_obs ~ exponential(1);
  sigma_patient ~ cauchy(0, 2.5);       // 환자 간 변동성에 대한 사전분포
  patient_effect_raw ~ normal(0, 1);    // 효율적인 샘플링을 위한 표준정규분포
  
  // 2. 가능도 (Likelihood)
  for (n in 1:N) {
    // 평균 = (레벨/측면별 평균) + (해당 환자의 고유 편차)
    width[n] ~ normal(mu[level_id[n], side_id[n]] + patient_effect[patient_id[n]], sigma_obs);
  }
}

// hierarchical_pedicle_width.stan 파일의 맨 아래에 추가

generated quantities {
  vector[N] log_lik;             // 모델 비교를 위한 로그 가능도
  vector[N] width_rep;           // 사후예측검증(PPC)을 위한 복제 데이터
  matrix[N_level, N_side] width_new_patient; // 새로운 환자의 레벨별 예측 넓이

  // 1. log_lik 및 width_rep 계산
  for (n in 1:N) {
    real current_mu = mu[level_id[n], side_id[n]] + patient_effect[patient_id[n]];
    log_lik[n] = normal_lpdf(width[n] | current_mu, sigma_obs);
    width_rep[n] = normal_rng(current_mu, sigma_obs);
  }

  // 2. 새로운 환자(Future Patient)에 대한 예측
  // 선생님의 원고 777페이지(academy_model.stan)와 같은 논리입니다.
  // 새로운 환자의 고유 편차를 sigma_patient 분포에서 무작위로 추출합니다.
  for (l in 1:N_level) {
    for (s in 1:N_side) {
      real random_patient_effect = normal_rng(0, sigma_patient);
      width_new_patient[l, s] = normal_rng(mu[l, s] + random_patient_effect, sigma_obs);
    }
  }
}

