# OpenACC_CG_sample
GPU移植のための実装例（共役勾配法によるポアソン方程式求解）

## 概要
* 共役勾配法のコードを、OpenACCでGPU対応する実例です。
* 同一のアルゴリズムをC言語(`src_c/`)およびFORTRAN(`src_f/`)で実装しています。
* 直方体格子上での熱伝導方程式を共役勾配法で解いており、主に`src_c/solver_CG.c`および`src_f/solver_CG.f`で求解しています。
* OpenACCでは、元のCPU専用コードに指示文を挿入することでGPU対応するため、下記の場合に特に有用です。
  * 大量の既存コードを最小限の手間でGPU化したい場合
  * 同一コードでCPU機・NVIDIA社製GPU搭載機の両方をサポートしたい場合
* 本レポジトリのプログラムのベースは、[東京大学情報基盤センターお試しアカウント付き並列プログラミング講習会「並列有限要素法で学ぶ並列プログラミング徹底入門」](http://nkl.cc.u-tokyo.ac.jp/FEMall/)で解説されている、CPU用のプログラムです。共役勾配法のアルゴリズムや、MPI・OpenMPによる並列化の方法に関しては、上記URLをご覧ください。
* 本レポジトリではサンプルプログラムを掲載していますが、現時点では解説は付けていません。東京大学情報基盤センターでは[お試しアカウント付き並列プログラミング講習会](https://www.cc.u-tokyo.ac.jp/events/lectures/)を多数開催しており、過去のGPU化に関する講習会の資料も公開されていますので、そちらもご参照ください。
  * [第208回講習会「OpenMPで並列化されたC++プログラムのGPU移植手法」](https://www.cc.u-tokyo.ac.jp/events/lectures/208/)
  * [第209回講習会「OpenACCとMPIによるマルチGPUプログラミング入門」](https://www.cc.u-tokyo.ac.jp/events/lectures/209/)
  * [第211回講習会「MPI+OpenMPで並列化されたFortranプログラムのGPUへの移行手法」](https://www.cc.u-tokyo.ac.jp/events/lectures/211/)

## 実行までの流れ

本サンプルコードは直方体格子上でポアソン方程式を解くプログラムです。直方体格子の定義ファイルを事前に生成しておく必要があるため、実行は下記の流れで行います。
1. 本レポジトリをローカルにコピーする
1. 主プログラムと格子生成プログラムをコンパイルする
1. 格子生成プログラムを実行する
1. 主プログラムを実行し、ポアソン方程式を解く

ここでは設定ファイルを変更せず、$`128^3`$個の格子を1GPUで処理する手順を示しています。[格子サイズや使用GPU数を変える](#格子サイズや使用GPU数を変える)こともできます。

### 本レポジトリをローカルにコピーする
任意の作業ディレクトリで下記のコマンドを実行します。
```
git clone https://github.com/kazuya-yamazaki/OpenACC_CG_sample.git
```
ディレクトリOpenACC_CG_sampleのフルパスを以下ではSAMPLEROOTと表記します。
```
SAMPLEROOT=`readlink -f OpenACC_CG_sample`
```

### 主プログラムと格子生成プログラムのコンパイル
Wisteria Aqariusで実行する場合は、下記コマンドでNVIDIA HPC SDKを利用可能な状態にします。
```
module load nvidia nvmpi cuda
```
まず、ディレクトリ`${SAMPLEROOT}/pmesh/`に移動し、`mpif90 -o pmesh pmesh.f`を実行します。

次に、主プログラムのディレクトリに移動します。C版を使う場合は`${SAMPLEROOT}/src_c/`、FORTRAN版を使う場合は`${SAMPLEROOT}/src_f/`に移動し、`make`を実行します。
ファイル`${SAMPLEROOT}/run/prog_gpu`が生成されていれば成功です。

### 格子生成プログラムを実行する
```
cd ${SAMPLEROOT}/pmesh/

# Wisteria Aquariusで実行する場合
pjsub pjsub_pmesh.sh # スクリプト内のyour_group_nameをグループ名に置き換えてください
# ワークステーションなどで実行する場合
mpirun -n 1 ./a.out
```
### 主プログラムを実行する
```
cd ${SAMPLEROOT}/run/

# Wisteria Aquariusで実行する場合
pjsub pjsub_run.sh # スクリプト内のyour_group_nameをグループ名に置き換えてください
# ワークステーションなどで実行する場合
mpirun -n 1 ./prog_gpu
```
上記の流れを一度行った後に、コードや設定を変えずに再実行する場合は、コンパイルと格子生成を省略して主プログラムをすぐに実行できます。

格子の隣接状況(connectivity)の算出に要した時間は`*** matrix conn. ? sec`、ポアソン方程式求解のための行列要素設定に要した時間は`*** matrix ass. ? sec`、CG法によるポアソン方程式求解に要した時間は`*** real  COMP. ? sec`として出力されます。また、熱伝導方程式を解いた結果得られた直方体領域の端での無次元温度が`Value on edge=?`という形で出力されます。

## 格子サイズや使用GPU数を変える
主プログラムを実行する前に、下記の流れで設定ファイルを編集し格子を再生成する必要があります。
1. 格子サイズと使用GPU数を決める
1. 格子生成プログラムと主プログラムの設定ファイルを編集する
1. 格子生成プログラムを再実行する
1. 主プログラムを実行する

### 直方体格子のサイズと使用GPU数を決める
格子数と使用GPU数(並列数)は下記のフォーマットで指定します。
```
(x方向格子数) (y方向格子数) (z方向格子数)
(x方向並列数) (y方向並列数) (z方向並列数)
```
合計で、(x方向並列数)×(y方向並列数)×(z方向並列数)個のGPUを使用します。また、各軸の格子数はその軸の並列数で割り切れる必要があります。

例えば、
```
256 256 256
2 2 2
```
とすると、Aquariusで1ノードに搭載された8個のGPUを全て使って$`256^3`$格子を扱うことができます。

また、6GPUのマシンを使う場合は、
```
240 256 256
3 2 1
```
のように並列数を割り振ることができます。ここでは、x軸を3並列に分割しているため、x方向の格子数を3の倍数としています。

問題サイズを固定し、並列数を変えて速度を比較するstrong scalingの性能を測定する場合は、並列数に様々な値を設定することになるので、そのいずれでも割り切れるように格子数を指定しておく必要があります。

### 格子生成プログラムと主プログラムの設定ファイルを編集する
格子生成プログラムの設定ファイル`${SAMPLEROOT}/pmesh/mesh.inp`の1行目と2行目を、決めた格子数・並列数に従って書き換えます：
```
(x方向格子数) (y方向格子数) (z方向格子数)
(x方向並列数) (y方向並列数) (z方向並列数)
```
3行目で格子ファイル名を指定します。pcubeのままにしておくと古い設定での格子ファイルを上書きします。設定毎に別々のファイル名で格子を生成しておくと、以後使いまわすことができます。

主プログラムの設定ファイルは`${SAMPLEROOT}/run/INPUT.DAT`です。1行目が格子ファイルへの相対パスなので、`${SAMPLEROOT}/pmesh/mesh.inp`の3行目の格子ファイル名を変更した場合はそれに合わせて変更します。2行目以降はポアソン方程式の収束条件に関わるパラメータで、結果の精度や実行時間に影響しますが、変えなくても実行可能です。

### 格子生成プログラムを再実行する
```
cd ${SAMPLEROOT}/pmesh/

# Wisteria Aquariusで8並列で実行する場合
pjsub pjsub_pmesh8.sh # スクリプト内のyour_group_nameをグループ名に置き換えてください
# ワークステーションなどで実行する場合
mpirun -n 並列数 ./a.out
```
Aquariusで実行する場合、ジョブスクリプトのノード数(1箇所)・プロセス数(2箇所)を並列数に合わせて書き換えてください。8GPUで実行するジョブスクリプトを`pmesh/pjsub_pmesh8.sh`に置いています。

### 主プログラムを実行する
```
cd ${SAMPLEROOT}/run/

# Wisteria Aquariusで8並列で実行する場合
pjsub pjsub_run8.sh # スクリプト内のyour_group_nameをグループ名に置き換えてください
# ワークステーションなどで実行する場合
mpirun -n 並列数 ./prog_gpu
```
Aquariusで実行する場合、ジョブスクリプトのノード数(1箇所)・プロセス数(2箇所)を並列数に合わせて書き換えてください。8GPUで実行するジョブスクリプトを`run/pjsub_run8.sh`に置いています。主プログラムはGPUを使用するため、Aquariusの1ノード(8GPU搭載)で実行できるのは8並列までです。2ノード使うと16並列まで実行することができます。

## 行列生成部分をGPU化する
サブルーチンMAT_CON0の行列生成処理は、変数INLUの同一インデックスにループから複数回アクセスする可能性があるため、そのまま並列化するとアクセス競合が起きて誤った結果になるおそれがあります。そのため、標準コードではGPU化せず、CPU上で逐次処理しているため時間がかかります。

`!$acc atomic`による排他アクセスを利用すると、こうした競合を避けてGPUで並列化することができます。そのコードはmat_con0_atomic.f (C版ではmat_con0_atomic.c)に記載しています。Makefileの中で`mat_con0.o`とある部分を`mat_con0_atomic.o`に置き換え、`make clean`したうえで`make`すると行列生成GPU版を利用することができます。

C言語では、外側（左側）の次元が非常に長い配列に限って、OpenACCでのCPU・GPU間の転送が遅いという制約があります。行列生成部分をGPU化する際に、この種の配列の転送が増えたため、所要時間はGPU化でかえって増えています。ただ、行列生成処理のループそのものは`#pragma acc atomic`によって高速化しています。また、Fortranでは上記のような性能低下は特に起きないため、行列生成処理のGPU化によって所要時間が減少します。
