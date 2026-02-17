# 0. 필요한 패키지 로드
library(tidyverse)
library(here)

# 1. 데이터 불러오기 (이미 불러온 df가 있다고 가정)
# 예: df <- read.csv("pedicle_data.csv")

# 2. 전처리: Wide to Long (지난번에 논의한 코드 활용)
df_long <- df %>%
  # ID 제외하고 나머지 열에서 '...숫자' 제거 [cite: 1, 125, 439]
  rename_with(~str_remove(., "\\.\\.\\..*"), -ID) %>%
  # 세로로 긴 형태로 변환 [cite: 1, 107, 44]
  pivot_longer(cols = -ID, names_to = "level_side", values_to = "width") %>%
  # level(T1~T9)과 side(Lt/Rt) 분리
  mutate(
    level = str_extract(level_side, "T[1-9]"),
    side = str_extract(level_side, "Lt|Rt"),
    # Stan 입력을 위한 숫자 인덱스 생성
    patient_idx = as.integer(as.factor(ID)),
    level_idx = as.integer(as.factor(level)),
    side_idx = as.integer(as.factor(side))
  ) %>%
  filter(!is.na(width)) # 결측치 제거

# 3. Stan 전달용 데이터 리스트 생성 [cite: 593, 1002]
stan_data <- list(
  N = nrow(df_long),
  N_patient = max(df_long$patient_idx),
  N_level = max(df_long$level_idx),
  N_side = max(df_long$side_idx),
  width = df_long$width,
  patient_id = df_long$patient_idx,
  level_id = df_long$level_idx,
  side_id = df_long$side_idx
)