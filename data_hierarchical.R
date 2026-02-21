# data_hierarchical.R
library(readxl)
library(tidyverse)
library(here)

# 1. 데이터 불러오기: 이전 연구와 동일한 엑셀 파일 로드
raw_df <- read_excel(here("kbpark_modify.xlsm"), skip = 1)

# 2. 전처리: Wide to Long (이전의 clean_df 로직을 df_long으로 계승)
df_long <- raw_df %>%
  # 33번부터 66번 컬럼(T1~T12 척추경 넓이 등) 선택
  select(ID, 33:66) %>% 
  # 컬럼 이름에서 불필요한 마침표(...) 등 제거
  rename_with(~str_remove(., "\\.\\.\\..*"), -ID) %>%
  # 긴 형식으로 변환
  pivot_longer(cols = -ID, names_to = "key", values_to = "width") %>%
  # 레벨(T1, T2 등)과 좌우(Left, Right) 정보 추출
  mutate(
    Level = str_extract(key, "[TL]\\d+"), 
    Side_Label = ifelse(str_detect(key, "Lt"), "Left", "Right")
  ) %>%
  # 결측치 및 0 이하의 값 제외
  filter(!is.na(width), width > 0) %>%
  # 현재 분석 대상인 T1~T6 레벨만 필터링 (필요시 T1~T9 등으로 수정 가능)
  filter(Level %in% c("T1", "T2", "T3", "T4", "T5", "T6")) %>%
  # Stan 모델 입력을 위한 숫자 인덱스 생성
  mutate(level_idx = as.integer(factor(Level, levels = c("T1", "T2", "T3", "T4", "T5", "T6"))))

# 3. Stan 데이터 리스트 생성
# main_hierarchical.R에서 사용할 stan_data 객체를 만듭니다.
# data_hierarchical.R의 마지막 부분
stan_data <- list(
  # 1. 차원 정의 (Scalars)
  N = nrow(df_long),                          # 전체 데이터 개수
  N_patient = length(unique(df_long$ID)),       # 고유한 환자 수
  N_level = 6,                                  # PT curve 범위 (T1~T6)
  N_side = 2,                                   # 좌(1), 우(2)
  
  # 2. 관측값 및 인덱스 (Vectors) - Stan의 이름과 100% 일치
  width = as.numeric(df_long$width),            # 종속변수
  patient_id = as.integer(factor(df_long$ID)),  # 환자 식별 번호 (1 ~ N_patient)
  level_id = df_long$level_idx,                 # 척추 레벨 번호 (1 ~ 6)
  side_id = ifelse(df_long$Side_Label == "Left", 1, 2) # 좌우 번호 (1 or 2)
)