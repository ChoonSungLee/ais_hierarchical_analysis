# AIS Scoliosis Analysis (Complete Pooling Model)

이 프로젝트는 2,700례의 청소년기 특발성 척추 측만증(AIS) 데이터를 활용하여 척추경 넓이를 분석하는 **완전 풀링 모델** 연구입니다.

## 📌 현재 진행 상태 (2026-02-17)
- [x] 데이터 전처리 로직 완성 (data_recipe.R)
- [x] 완전 풀링 모델 Stan 코드 작성
- [ ] MCMC 샘플링 및 수렴 진단 (Rhat 확인 필요)
- [ ] 사전분포 민감도 분석 수행

## 💡 주요 메모
- 현재는 모든 환자를 하나의 그룹으로 보는 **Pooling 모델**에 집중하고 있음.
- 이 단계가 완벽히 정리된 후, 별도의 프로젝트(`ais_hierarchical_analysis`)에서 계층 모델 작업을 재개할 예정임.
- 사전분포 `mu ~ normal(5, 3)` 설정의 근거를 문헌과 대조해 볼 것.

## 📂 파일 구성
- `main_analysis.R`: 전체 분석 실행 스크립트
- `pedicle_pooling.stan`: 풀링 모델용 Stan 코드
- `data_recipe.R`: 2,700례 데이터 로딩 및 Stan 입력용 변환