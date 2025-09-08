:: ============================================
:: Plink genotype data QC (run in RStudio Terminal on Windows CMD)
:: ============================================

:: 0) 출력 폴더 생성
mkdir 00_import 01_autosomesX 02_reports 03_qc 04_split

:: 1) PED/MAP → PGEN 변환
::    --cow        : 소(cattle) 빌트/염색체 규칙 사용
::    --pedmap     : 입력 prefix (NIAS_ibv3_296ea.ped/.map)
::    --sort-vars  : 변이를 유전체 좌표 순으로 정렬
::    --make-pgen  : 출력 형식은 PLINK2(.pgen/.pvar/.psam)
::    --out        : 출력 접두어
plink2 --cow --pedmap NIAS_ibv3_296ea --sort-vars --make-pgen --out 00_import\NIAS_ibv3_296ea

:: 2) PVAR에서 REF/ALT가 둘 다 '.'인 문제 변이 목록 추출 및 제거
::    (이 상태는 "Duplicate allele code" 오류를 유발함)

:: 2-1) REF(4열)과 ALT(5열)가 모두 '.'인 전체 행 저장(헤더 제외)
awk -F '\t' "NR>1 && $4==\".\" && $5==\".\"" 00_import\NIAS_ibv3_296ea.pvar > bad_refalt_ids_both_strict.txt

:: 2-2) 같은 조건의 변이 ID(3열)만 추출
awk -F '\t' "NR>1 && $4==\".\" && $5==\".\" {print $3}" 00_import\NIAS_ibv3_296ea.pvar > bad_refalt_ids.txt

:: 2-3) 문제 변이 제외 후 깨끗한 pfile로 저장
plink2 --pfile 00_import/NIAS_ibv3_296ea --exclude bad_refalt_ids.txt --make-pgen --out 00_import/NIAS_ibv3_296ea.clean

:: 3) 클린 필터(필수): 오토솜+X염색체, biallelic, 단일형 제거, 중복 ID 제거 → 이후 단계에서 사용하는 autosomesX 생성
::    --chr 1-29,X         : 염색체 1–29 + X만 사용
::    --max-alleles 2      : 대립유전자 ≤ 2 (biallelic만 유지)
::    --mac 1              : 최소 소수 대립유전자 수 ≥ 1 (단일형 변이 제거)
::    --rm-dup force-first : 중복 변이 발견 시 첫 번째 것만 유지 (좌표/allele 기준)
::    --sort-vars          : 좌표 순으로 정렬
::    --make-pgen          : pgen/pvar/psam 출력
::    --out                : autosomesX 세트 출력 위치

plink2 --pfile 00_import/NIAS_ibv3_296ea.clean --chr 1-29,X --max-alleles 2 --mac 1 --rm-dup force-first --sort-vars --make-pgen --out 01_autosomesX/NIAS_ibv3_296ea.autosomesX

:: 3-추가) 품종별 데이터셋 분리 준비: Holstein/Jersey IID 목록을 FID IID 두 열로 변환
::        (holstein.txt / jersey.txt에는 IID만 있다고 가정, .psam에서 FID를 매칭해온다)

:: Holstein: 첫 번째 파일(holstein.txt)의 IID를 해시에 저장하고, .psam에서 일치하는 IID 행의 FID IID를 출력
awk "NR==FNR {iid[$1]=1; next} FNR>1 && ($2 in iid) {print $1, $2}" holstein.txt 01_autosomesX\NIAS_ibv3_296ea.autosomesX.psam > holstein_ids.txt

:: Jersey: 동일 처리
awk "NR==FNR {iid[$1]=1; next} FNR>1 && ($2 in iid) {print $1, $2}" jersey.txt 01_autosomesX\NIAS_ibv3_296ea.autosomesX.psam > jersey_ids.txt

:: 4) 품종별 서브셋 pfile 만들기

:: Holstein 서브셋 pfile 생성 (--keep: FID IID 목록의 샘플만 유지)
plink2 --pfile 01_autosomesX/NIAS_ibv3_296ea.autosomesX --keep holstein_ids.txt --make-pgen --out 04_split/NIAS_ibv3_296ea.holstein

:: Jersey 서브셋 pfile 생성
plink2 --pfile 01_autosomesX/NIAS_ibv3_296ea.autosomesX --keep jersey_ids.txt --make-pgen --out 04_split/NIAS_ibv3_296ea.jersey

:: 5) 품종별 기본 QC 리포트 만들기
::    --missing : 샘플/변이 결측률(missing genotype rate)
::    --hardy   : HWE(하디–바인베르그) 검정
::    --freq    : 대립유전자 빈도

:: 5-1) 품종별 QC 리포트 폴더 생성
mkdir 02_reports\holstein 02_reports\jersey

:: Holstein QC 리포트
plink2 --pfile 04_split\NIAS_ibv3_296ea.holstein --missing --hardy --freq --out 02_reports\holstein\NIAS_ibv3_296ea.holstein.autosomesX

:: Jersey QC 리포트
plink2 --pfile 04_split\NIAS_ibv3_296ea.jersey --missing --hardy --freq --out 02_reports\jersey\NIAS_ibv3_296ea.jersey.autosomesX

:: 6) RStudio에서 품종별 QC 리포트 확인 및 하드 필터 값 결정
:: check_genotype_data_vF.R script 실행

:: 7) 품종별 하드 필터 적용(값은 연구 목적/샘플 크기에 맞게 조정)
::    --mind 0.05 : 개체 결측률 ≤ 5%만 유지(0.05 초과 개체 제외)
::    --geno 0.05 : 변이 결측률 ≤ 5%만 유지(0.05 초과 변이 제외)
::    --maf  0.05(Holstein) 및 0.01(Jersey) : MAF ≥ 0.05/0.01 변이만 유지
::    --hwe  1e-6 : HWE p < 1e-6 변이 제외
::    --make-bed  : PLINK1 바이너리(.bed/.bim/.fam) 파일로 저장

:: Holstein
plink2 --pfile 04_split\NIAS_ibv3_296ea.holstein --mind 0.05 --geno 0.05 --maf 0.05 --hwe 1e-6 --make-bed --out 03_qc/NIAS_ibv3_296ea.holstein

:: Jersey 
plink2 --pfile 04_split\NIAS_ibv3_296ea.jersey --mind 0.05 --geno 0.05 --maf 0.01 --hwe 1e-6 --make-bed --out 03_qc/NIAS_ibv3_296ea.jersey
