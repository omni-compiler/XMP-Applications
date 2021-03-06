MODYLAS-MINI
============

* version: 1.0.0 (based on MODYLAS 1.0.1)
* date: 2014/10/02
* contact: miniapp@riken.jp


MODYLASについて
---------------

本ミニアプリは，名古屋大学および分子科学研究所で開発が進められている
汎用古典分子動力学シミュレーションプログラムMODYLASをベースにしています．
MODYLASは，大規模並列計算に対応するために，長距離クーロン相互作用計算に
FMM(Fast Multipole Method)を使用しています．
MODYLASの詳細については，[MODYLAS Web site](http://www.modylas.org/)を参照ください．

オリジナルMODYLAS連絡先: 安藤嘉倫<yoshimichi.andoh@apchem.nagoya-u.ac.jp>


インストール
------------

外部ライブラリとして，MPIが必要．
ビルドツールとして，GNU Make, OpenMPに対応したFortran90コンパイラ，Cコンパイラが必要．

 1. 本パッケージの入手とファイルの展開

 2. srcディレクトリに移動し、"make_setting"ファイルの内容を実行環境に合わせて編集する．
    同ディレクトリには，以下の設定例が含まれる:
    * make_setting.intel : Intelコンパイラ
    * make_setting.gcc   : GCCコンパイラ
    * make_setting.fx10  : 富士通コンパイ(京/FX10)

 3. srcディレクトリでmakeを実行．
    正しくmakeが完了すると、srcディレクトリ内に実行プログラム`modylas_mini`が作成される．


テスト
------

インタラクティブにMPIジョブが実行できる環境用に，簡単なテストスクリプトをtestディレクトリに用意してある．
テストを実行するには， testディレクトリでシェルスクリプトgo.shを実行するか，srcディレクトリで「make test」を実行する．


プログラム実行方法
-------------------

### 入力ファイル

modylas_miniの実行には，以下の3つの入力ファイルが必要．
(以下で，「セッション名」は，実行条件を識別するための任意の文字列)

* 設定ファイル

    プログラムの各種実行パラメータを指定するテキストファイル．
    ファイル名は「セッション名.mdconf」

* 力場パラメータファイル

    分子動力学計算で参照される力場パラメータなどの情報を納めたバイナリファイル．
    ファイル名は「セッション名.mdff.bin」

* 座標ファイル

    分子動力学計算の初期条件となる系全体の原子座標情報を納めたバイナリファイル．
    ファイル名は「セッション名.mdxyz.bin」

### 実行

実行時に，実行プログラムの引数としてセッション名を指定する．
以下の例は，セッションwat111を，8スレッド，64ノードで実行する．

    $ export OMP_NUM_THREADS=8
    $ mpiexec -n 64 path_to_src_directory/modylas_mini wat111   


### 出力ファイル

* 物理量モニタファイル

    物理量の時間発展を出力したテキストファイル．
    ファイル名は「セッション名.mdmntr」

* 軌道ファイル
    原子の座標および速度の時間発展を出力したバイナリファイル．
    ファイル名は「セッション名.mdtrj.bin」

* リスタートファイル
    リスタート用バイナリファイル．
    ファイル名は「セッション名.mdrestart.bin」

物理量モニタファイル，軌道ファイル，リスタートファイルについては，
出力開始ステップ(start)と出力ステップ間隔(interval)を設定ファイル中で指定する．

### 注意 ###

#### 設定ファイルのフルアプリとの非互換性

本ミニアプリの設定ファイルは，フルアプリ(オリジナルMODYLAS)の設定ファイル(拡張子mddef)とは書式および指定内容が異なっている．そのため，ミニアプリの設定ファイルをフルアプリ実行時に，フルアプリ設定ファイルをミニアプリ実行時に使用することはできない．

####  並列プロセス数(ノード数)に関する制限

 - セル数(設定ファイル内で指定したncellの3乗)以下であること
 - 8以上であること
 - 2のべき乗であること

#### Big endian環境での実行

付属するサンプル入力データでは，力場パラメータファイルと座標ファイルはlittle endian環境で作成してある．
そのため，京などのSPARCベースのbig endian環境で実行するには，バイナリファイル読み込み時のendian変換機能を使用する必要がある(京の場合は環境変数「FORT90L='-Wl,-T'」を指定)．


設定ファイル
------------

プログラムの実行パラメータを，ファイル「セッション名.mdconf」の中で指定する．
一行に一つのパラメータを，等号(=)の左辺に「key(パラメータ名)を右辺に「 value(パラメータ値)」を記述することにより指定する．
key, valueの両端には，任意個のスペースを入れることができる．
また，空行および各行の「#」以降は無視される．

以下に，指定可能なkeyとvalueの型を示す．
bool型のパラメータには，yes/no, on/off, true/false, およびこれらの全体または一部を大文字にした文字列が指定できる．

### 出力ファイルパラメータ

    mntr_start(int):       物理量モニタファイル出力開始ステップ(ディフォルト 0)
    mntr_interval(int):    物理量モニタファイル出力ステップ間隔(ディフォルト 1)
    trj_start(int):        軌道ファイル出力開始ステップ(ディフォルト 0)
    trj_interval(int):     軌道ファイル出力ステップ間隔(ディフォルト 1)
    restart_start(int):    リスタートファイル出力開始ステップ(ディフォルト 0)
    restart_interval(int): リスタートファイル出力ステップ間隔(ディフォルト 1)

### 初期値設定パラメータ

    maxwell_velocities(bool): 原子速度をMaxwell分布に再設定(ディフォルト no)
    temperature(real):        温度(K)(ディフォルト 300.0)
    randomseed(int):          擬似乱数種(ディフォルト 1235)

### 時間ステップパラメータ

    dt(real):                時間ステップ間隔(sec)(ディフォルトなし)
    step(int):               時間ステップ数(ディフォルトなし)
    nstep_skip_middle(int):  詳細情報「マルチ時間ステップ」を参照(ディフォルト 1)
    nstep_skip_long(int):    詳細情報「マルチ時間ステップ」を参照(ディフォルト 1)

### 領域分割パラメータ

    manual_division(bool):  マニュアル分割フラグ(ディフォルト no)
    nxdiv(int):             x方向ノード数(マニュアル分割時には必須)
    nydiv(int):             y方向ノード数(マニュアル分割時には必須)
    nzdiv(int):             z方向ノード数(マニュアル分割時には必須)

### FMMパラメータ

    ncell(int):             セル数(各辺のセル分割数)(ディフォルトなし)
    nmax(int):              多重極子展開/局所展開での展開次数(ディフォルト 4)
    ULswitch(int):          詳細情報「MPIプロセス間通信」を参照(ディフォルト 1)

### SHAKE/RATTLEパラメータ

    shake_max_iteration(int): 最大反復回数(ディフォルト 100)
    shake_tolerance(real):    収束判定値(ディフォルト 1.0e-8)

### その他のMDパラメータ

    cutoff(real):             カットオフ半径(A)(ディフォルトなし)
    ewald_surface_term(bool): Ewald表面項の有無(ディフォルト off)


入力データサンプル
------------------

以下の3ケースの水分子系のデータをdataディレクトリに用意した．
各ディレクトリには，入力データファイルの他に，京でのジョブスクリプト
サンプルgo.shと，実行結果例refvaluesを納めてある．

### wat111

 - 原子数: 19,530
 - セル分割: 8x8x8
 - FMMツリーレベル: 3

### wat222

 - 原子数: 156,240
 - セル分割: 16x16x16
 - FMMツリーレベル: 4

### wat444

 - 原子数: 1,249,920
 - セル分割: 32x32x32
 - FMMツリーレベル: 5


性能計測時の計算条件設定については，"performance.md"を参照すること．

他サイズの入力データを希望する場合は，オリジナルMODYLAS開発者に連絡すること．


オリジナルMODYLASからの変更点
-----------------------------

* 機能を限定することにより，コードを削減
 
    - NVEアンサンブル(原子数，体積，エネルギー一定)の分子動力学計算に限定
   
    - 対象を水分子系に限定(bond, angle, ub, dihedral, CMAP関連のコードを削除)
   
* 未使用なファイル，関数，配列，コードの削除

* その他のコード整備


エクサスケールでの想定計算規模
------------------------------

10億原子系に対して，10^9ステップの計算を150時間で完了したい．



詳細情報
--------

### マルチ時間ステップ ####

入力パラメータ:

    config. parameter    valiable name
    --------------------+----------------------------
    steps                md_condition__howmany_steps
    nstep_skip_middle    maxMTm
    nstep_skip_long      maxMTl
    dt                   dt


力の計算を以下の手順で行っている:

  - 毎ステップごとに

     短距離成分力(bond, angle, ...)を計算
     (水原子系に限定しているので，実際には何も計算しない)

  - maxMTmステップごとに

     中距離成分(FMMの直接二体計算部分，ファンデルワールス力)を計算

  - maxMTm * maxMTlステップごとに

     長距離成分(FMMの多重極展開〜局所展開による計算部分)を計算


実際の実行ステップ数は md_condition__howmany_steps * maxMTm * maxMTl回，
ステップの時間刻み幅は dt / maxMTm / maxMTl.


### MPIプロセス間通信 

主な通信は，以下の4部分からなる．

- COMM_DIRECT:

    力計算の前に，袖領域に属する原子の座標データを隣接ノードからコピー．
    毎時間ステップ実行，サブルーチンcomm_direct_3.

- MIGRATION:

    原子のノード担当領域間移動にともなう，原子データの隣接ノード間でのコピー．
    maxMTm*maxMTlステップごとに実行．サブルーチンcomm_bound.

- COMM_FMM:

    M2M(多重極展開の展開中心シフト)およびM2L(多重極展開から局所展開への変換)
    計算時に参照される多重極展開係数データのノード間コピー．
    maxMTm*maxMTlステップごとに，ツリーレベルを上昇しながら，各レベルで実行．
    サブルーチンcomm_fmm_local_multi(低ツリーレベル用),
    comm_fmm_local_top(高ツリーレベル用).

- ENE_REDUCTON:

    全系のエネルギー，熱力学量の集約．
    maxMTm*maxMTlステップごとに実行．サブルーチンcalc_hamiltonian_nve.

現在の実装では，通信担当ルーチンは3-Dトーラス構造を持ったネットワークに最適化されている．
それらのシステムにおいてベストな性能を得るためには，"performance.md"を参照のこと．

なお，MODYLAS-MINIでは二体力の直接計算に作用・反作用則を利用していないため，
力計算後に袖領域に属する原子の力データを隣接ノードへ通信する必要はない．


References
----------

*  Andoh, Y. et al., "MODYLAS: A Highly Parallelized General-Purpose
   Molecular Dynamics Simulation Program for Large-Scale Systems with
   Long-Range Forces Calculated by Fast Multipole Method (FMM) and
   Highly Scalable Fine-Grained New Parallel Processing Algorithms",
   J. Chem. Theory Comput., 2013, 9 (7), pp 3201-3209,
   DOI: 10.1021/ct400203a.
