# AIS Scoliosis Analysis (True Hierarchical Model)

이 프로젝트는 2,700례의 AIS 데이터를 활용하여 환자별 변동성을 고려하는 **진정한 계층적 모델(True Hierarchical Model)** 연구입니다.

## ⚠️ 잠정 중단 안내 (2026-02-17)
- 현재 **완전 풀링 모델(Complete Pooling Model)**의 완벽한 숙달을 위해 이 프로젝트 작업을 잠시 멈춥니다.
- 풀링 모델에서 데이터의 특성을 충분히 파악한 후 다시 이 폴더로 돌아올 예정입니다.

## 📌 현재까지 완료된 작업
- [x] 계층적 구조를 반영한 Stan 코드 작성 (GQB 포함)
- [x] `bayesplot`을 활용한 레벨별 예측 분포 시각화 코드 준비
- [x] GitHub 원격 저장소(`origin`) 연결 및 초기 푸시 완료

## 💡 복귀 시 검토 사항
- 풀링 모델의 결과와 계층 모델의 **수축 효과(Shrinkage effect)** 비교 분석.
- 사전분포(Hyperprior) 설정이 사후분포에 미치는 민감도 분석 실행.
- GQB 블록에서 생성된 `width_new_patient` 결과 검증.

## 📂 주요 파일
- `hierarchical_pedicle_width.stan`: 환자 효과가 포함된 계층 모델
- `main_hierarchical.R`: 분석 실행 및 시각화 스크립트