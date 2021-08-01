/***********************************/
/* 遺伝的アルゴリズムによるTSPの解 */
/*      coded by Y.Suganuma        */
/***********************************/

/***********************/
/* クラスSpeciesの定義 */
/***********************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "MT.h"

class Species {

	protected:

		double max;   // 最大適応度
		double mean;   // 平均適応度
		double *pi;   // 適応度
		double *ro;   // ルーレット板
		int allele_u;   // 対立遺伝子上限
		int allele_l;   // 対立遺伝子下限
		int size;   // 個体総数
		int max_ch;   // 子供の数の最大値
		int max_len;   // 最大遺伝子長
		int min_len;   // 最小遺伝子長（負の時は，最大遺伝子長で固定）
		int max_n;   // 最大適応度の個体番号
		int dup_a;   // 遺伝子の重複
                     //   =0 : 重複を許さない
                     //   =1 : 重複を許す
		int dup_s;   // 個体の重複（同じ染色体の個体）
                     //   =0 : 重複を許さない
                     //   =1 : 重複を許す
		int **ind;   // 集団（個体の集まり）
		int *len;   // 各個体の遺伝子長
		int *kou1;   // 交叉・突然変異用作業場所１
		int *kou2;   // 交叉・突然変異用作業場所２
		int *s_w;   // 淘汰用指標（選択された回数）
		int **edge;   // エッジ組み替え交叉用ワークエリア
		char *pi_w;   // 適応度計算指標
                     //   =0 : 未使用
                     //   =1 : 適応度計算前（突然変異はこの個体だけに適用）
                     //   =2 : 適応度計算済み（交叉時に親とみなす）

		int Position(int);
		int Select(int method=-1, double bias=0.0, double step=1.0);

	public:
					// コンストラクタ
		Species(char *, long);
					// デストラクタ
		~Species();
					// 標準的な初期設定
		void Init_std();
					// 標準的な出力
		void Out_std(int, int, int, char *);
					// 交叉
						// 親のコピー
		void C_copy(int method=2, int pair=0, int k_method=-1,
                    double k_bias=0.0, double k_step=1.0);
						// 多点交叉
		void C_point(double, int k_point=1, int k_vr=0, int k_method=-1,
                     double k_bias=0.0, double k_step=1.0);
						// 一様交叉
		void C_uniform(double, int k_method=-1, double k_bias=0.0,
                       double k_step=1.0);
						// 平均化交叉
		void C_mean(double, int k_method=-1, double k_bias=0.0,
                    double k_step=1.0);
						// 循環交叉
		void C_cycle(double, int k_method=-1, double k_bias=0.0,
                     double k_step=1.0);
						// 部分的交叉
		void C_part(double, int k_method=-1, double k_bias=0.0,
                    double k_step=1.0);
						// 順序交叉
		void C_seq(double, int k_method=-1, double k_bias=0.0,
                   double k_step=1.0);
						// 一様順序交叉
		void C_useq(double, int k_method=-1, double k_bias=0.0,
                    double k_step=1.0);
						// 一様位置交叉
		void C_upos(double, int k_method=-1, double k_bias=0.0,
                    double k_step=1.0);
						// エッジ組み替え交叉
		void C_edge(double, int k_method=-1, double k_bias=0.0,
                    double k_step=1.0);
						// サブツアー交叉
		void C_sub(double, int count=10);
					// 突然変異
						// 対立遺伝子への置換
		void M_alle(double);
						// 移動
		void M_move(double);
						// 逆位
		void M_inv(double, int wd=0);
						// スクランブル
		void M_scram(double pr, int wd=0);
						// 転座
		void M_chg(double pr, int wd=0);
						// 重複
		void M_dup(double pr, int wd=0);
						// 摂動
		void M_per(double pr, int method=0, double m=0.0, double s=1.0);
						// 挿入
		void M_ins(double pr, int wd=0);
						// 削除
		void M_del(double pr, int wd=0);
                    // エリート・ルーレット選択
		void S_roul(int elite=0, int s_method=1, double s_bias=0.0,
                    double s_step=1.0);
};

double norm_d(double, double);

/****************************/
/* コンストラクタ           */
/*      name : ファイル名   */
/*      seed : 乱数の初期値 */
/****************************/
Species::Species(char *name, long seed)
{
	int i1, kind, num;
	FILE *in;
/*
     データの入力
*/
	in = fopen(name, "r");

	fscanf(in,"%*s %d %*s %d", &allele_u, &allele_l);
	fscanf(in,"%*s %d %*s %d", &max_len, &min_len);
	fscanf(in,"%*s %d %*s %d", &dup_a, &dup_s);
	fscanf(in,"%*s %d %*s %d", &size, &max_ch);
/*
     データのチェック
*/
	if (size <= 0) {
		printf("***error  個体総数≦０ (Constructor)\n");
		exit(1);
	}

	if (max_ch < 0) {
		printf("***error  子供の数＜０ (Constructor)\n");
		exit(1);
	}

	if (max_len <= 0 || min_len == 0) {
		printf("***error  遺伝子長≦０ (Constructor)\n");
		exit(1);
	}

	if (max_len < min_len) {
		printf("***error  最大遺伝子長＜最小遺伝子長 (Constructor)\n");
		exit(1);
	}

	if (allele_u <= allele_l) {
		printf("***error  対立遺伝子上限≦対立遺伝子下限 (Constructor)\n");
		exit(1);
	}

	kind = allele_u - allele_l + 1;
	if (dup_a == 0 && max_len > kind) {
		printf("***error  遺伝子の重複を防ぐことはできない (Constructor)\n");
		exit(1);
	}
/*
     領域の確保
*/
	num = size + max_ch;

	ind = new int * [num];
	for (i1 = 0; i1 < num; i1++)
		ind[i1] = new int [max_len];

	edge = new int * [max_len];
	for (i1 = 0; i1 < max_len; i1++)
		edge[i1] = new int [5];

	pi   = new double [num];
	ro   = new double [num];
	len  = new int [num];
	kou1 = new int [max_len];
	kou2 = new int [max_len];
	s_w  = new int [num];
	pi_w = new char [num];
/*
     乱数の初期設定
*/
	init_genrand(seed);
}

/****************/
/* デストラクタ */
/****************/
Species::~Species()
{
	int i1;

	for (i1 = 0; i1 < size+max_ch; i1++)
		delete [] ind[i1];
	delete [] ind;

	for (i1 = 0; i1 < max_len; i1++)
		delete [] edge[i1];
	delete [] edge;

	delete [] pi;
	delete [] len;
	delete [] kou1;
	delete [] kou2;
	delete [] pi_w;
	delete [] s_w;
	delete [] ro;
}

/**************************************************/
/* 場所を探す                                     */
/*      n : >=0 : n番目の親を捜す                 */
/*          -1 : 空いている場所を探す             */
/*      return : 親の場所，または，空いている場所 */
/*               （存在しないときは負の値）       */
/**************************************************/
int Species::Position(int n)
{
	int i1, k = -1, sw = 0;
/*
     空いている場所を探す
*/
	if (n < 0) {
		for (i1 = 0; i1 < size+max_ch && k < 0; i1++) {
			if (pi_w[i1] == 0)
				k = i1;
		}
		if (k < 0) {
			printf("***error  空いている場所がない --Position--\n");
			exit(1);
		}
	}
/*
     ｎ番目の親（pi_w[i]=2）を捜す
*/
	else {
		for (i1 = 0; i1 < size+max_ch && sw == 0; i1++) {
			if (pi_w[i1] == 2) {
				k++;
				if (k == n) {
					sw = 1;
					k  = i1;
				}
			}
		}
	}

	return k;
}

/*******************************************************************/
/* 個体の選択                                                      */
/*      method : 選択方法                                          */
/*                 =-1 : ランダム(default)                         */
/*                 =0 : 適応度をそのまま使用                       */
/*                 =1 : 最小値からの差（ただし，α以下の場合はα） */
/*                 =2 : 評価値に順位をつけ，減少率βで線形化       */
/*      bias : α，または，method=2の場合は初期値(default=0)       */
/*      step : β(default=1)                                       */
/*      return : 個体番号                                          */
/*******************************************************************/
int Species::Select(int method, double bias, double step)
{
	double sum = 0.0, x;
	int i1, k, min, n, sw;
					// ルーレット板の用意
	switch (method) {
						// ランダム
		case -1:
			n = 0;
			for (i1 = 0; i1 < size+max_ch; i1++) {
				if (pi_w[i1] > 1)
					n++;
			}
			sum = 1.0 / n;
			for (i1 = 0; i1 < size+max_ch; i1++) {
				if (pi_w[i1] > 1)
					ro[i1] = sum;
			}
			break;
						// 評価値をそのまま利用
		case 0:
			n = 0;
			for (i1 = 0; i1 < size+max_ch; i1++) {
				if (pi_w[i1] > 1) {
					sum += pi[i1];
					n++;
				}
			}
			if (fabs(sum) > 1.0e-10) {
				sum = 1.0 / fabs(sum);
				for (i1 = 0; i1 < size+max_ch; i1++) {
					if (pi_w[i1] > 1)
						ro[i1] = pi[i1] * sum;
				}
			}
			else {
				sum = 1.0 / n;
				for (i1 = 0; i1 < size+max_ch; i1++) {
					if (pi_w[i1] > 1)
						ro[i1] = sum;
				}
			}
			break;
						// 最小値からの差
		case 1:
			min = -1;
			n   = 0;
			for (i1 = 0; i1 < size+max_ch; i1++) {
				if (pi_w[i1] > 1) {
					n++;
					if (min < 0 || pi[i1] < pi[min])
						min = i1;
				}
			}
			for (i1 = 0; i1 < size+max_ch; i1++) {
				if (pi_w[i1] > 1) {
					ro[i1] = pi[i1] - pi[min];
					if (ro[i1] < bias)
						ro[i1] = bias;
					sum += ro[i1];
				}
			}
			if (sum > 1.0e-10) {
				sum = 1.0 / sum;
				for (i1 = 0; i1 < size+max_ch; i1++) {
					if (pi_w[i1] > 1)
						ro[i1] *= sum;
				}
			}
			else {
				sum = 1.0 / n;
				for (i1 = 0; i1 < size+max_ch; i1++) {
					if (pi_w[i1] > 1)
						ro[i1] = sum;
				}
			}
			break;
						// 線形化
		case 2:
			n = 0;
			for (i1 = 0; i1 < size+max_ch; i1++) {
				if (pi_w[i1] > 1) {
					ro[i1] = -1.0;
					n++;
				}
				else
					ro[i1] = 1.0;
			}
			sw  = 0;
			sum = bias;
			while (sw == 0) {
				min = -1;
				for (i1 = 0; i1 < size+max_ch; i1++) {
					if (ro[i1] < 0.0 && (min < 0 || pi[i1] < pi[min]))
						min = i1;
				}
				if (min < 0)
					sw = 1;
				else {
					ro[min]  = sum;
					sum     += step;
				}
			}
			sum = 1.0 / (0.5 * (2.0 * bias + step * (n - 1)) * n);
			for (i1 = 0; i1 < size+max_ch; i1++) {
				if (pi_w[i1] > 1)
					ro[i1] *= sum;
			}
			break;
	}

	sum = 0.0;
	for (i1 = 0; i1 < size+max_ch; i1++) {
		if (pi_w[i1] > 1) {
			sum    += ro[i1];
			ro[i1]  = sum;
		}
	}
					// 選択
	x  = genrand_real3();
	sw = 0;
	k  = 0;
	for (i1 = 0; i1 < size+max_ch && sw == 0; i1++) {
		if (pi_w[i1] > 1) {
			if (x <= ro[i1]) {
				sw = 1;
				k  = i1;
			}
		}
	}

	return k;
}

/********************/
/* 標準的な初期設定 */
/********************/
void Species::Init_std()
{
	int i1, i2, i3, length, lid, sw1, sw2;
/*
     初期設定
*/
	for (i1 = 0; i1 < size+max_ch; i1++) {
		if (i1 < size)
			pi_w[i1] = 1;   // 適応度の計算前
		else
			pi_w[i1] = 0;   // 未使用
	}
/*
     遺伝子の決定
*/
	for (i1 = 0; i1 < size; i1++) {

		sw1 = 0;

		while (sw1 == 0) {
					// 遺伝子長の決定
			if (min_len < 0)
				length = max_len;
			else {
				length = (int)(genrand_real3() * (max_len - min_len + 1) + min_len);
				if (length > max_len)
					length = max_len;
			}
			len[i1] = length;
					// 遺伝子の決定
			for (i2 = 0; i2 < length; i2++) {
				sw2 = 0;
				while (sw2 == 0) {
					lid = (int)(genrand_real3() * (allele_u - allele_l + 1) + allele_l);
					if (lid > allele_u)
						lid = allele_u;
					ind[i1][i2] = lid;
						// 重複遺伝子のチェック
					sw2 = 1;
					if (dup_a == 0) {
						for (i3 = 0; i3 < i2 && sw2 > 0; i3++) {
							if (lid == ind[i1][i3])
								sw2 = 0;
						}
					}
				}
			}
					// 重複個体のチェック
			sw1 = 1;
			if (dup_s == 0) {
				for (i2 = 0; i2 < i1 && sw1 > 0; i2++) {
					if (len[i1] == len[i2]) {
						sw2 = 0;
						for (i3 = 0; i3 < len[i1] && sw2 == 0; i3++) {
							if (ind[i1][i3] != ind[i2][i3])
								sw2 = 1;
						}
						if (sw2 == 0)
							sw1 = 0;
					}
				}
			}
		}
	}
}

/****************************************************/
/* 標準的な出力                                     */
/*      sw : 出力レベル                             */
/*             =0 : 最終出力だけ                    */
/*             n>0 : ｎ世代毎に出力（負はファイル） */
/*      out_m : 出力方法                            */
/*                =0 : すべての個体を出力           */
/*                =1 : 最大適応度の個体だけを出力   */
/*      gen : 現在の世代番号                        */
/*      name : 出力ファイル名                       */
/****************************************************/
void Species::Out_std(int sw, int out_m, int gen, char *name)
{
	int i1, i2, k = 0, pr;
	char *now;
	time_t aclock;
	FILE *out;

	if (sw >= 0) {
		printf("   出力先は（0:出力なし，n:画面にｎ個づつ，-1:ファイル）？ ");
		scanf("%d", &pr);
	}
	else
		pr = -1;

	if (pr != 0) {
					// 出力先の決定と評価値の出力
		if (pr > 0) {
			out = stdout;
			getchar();
		}
		else {
			time(&aclock);
			now = ctime(&aclock);
			out = fopen(name, "a");
			fprintf(out, "***世代 %d 適応度 max %f (%d) mean %f 時間 %s\n",
                    gen, max, max_n, mean, now);
		}
					// 詳細出力
		for (i1 = 0; i1 < size+max_ch; i1++) {
			if ((pi_w[i1] > 1) && (out_m ==0 || out_m == 1 && i1 == max_n)) {
				fprintf(out, "%d allele", i1);
				for (i2 = 0; i2 < len[i1]; i2++)
					fprintf(out, " %d", ind[i1][i2]);
				fprintf(out, " value %f\n", pi[i1]);
				if (pr > 0) {
					k++;
					if (k == pr) {
						getchar();
						k = 0;
					}
				}
			}
		}

		if (pr < 0)
			fclose(out);
	}
}

/*******************************************************************/
/* 交叉（親のコピー）                                              */
/*      method : =2 : 有性（２つの親から２つの子供）(default)      */
/*               =1 : １つの親から１つの子供                       */
/*      pair : method=2 の時は親のペア数(default=max_ch/2)         */
/*             method=1 の時は親の数（＝子供の数）                 */
/*      k_method : 選択方法                                        */
/*                 =-1 : ランダム(default)                         */
/*                 =0 : 適応度をそのまま使用                       */
/*                 =1 : 最小値からの差（ただし，α以下の場合はα） */
/*                 =2 : 評価値に順位をつけ，減少率βで線形化       */
/*      k_bias : α，または，method=2の場合は初期値(default=0)     */
/*      k_step : β(default=1)                                     */
/*******************************************************************/
void Species::C_copy(int method, int pair, int k_method,
                     double k_bias, double k_step)
{
   int i1, i2, i3, k, p, p1, p2, sw;
/*
     初期設定とデータチェック
*/
	if (method != 1)
		method = 2;

	if (pair <= 0)
		pair = (method==2) ? max_ch/2 : max_ch;
	else {
		if (method == 2 && 2*pair > max_ch || method == 1 && pair > max_ch) {
			printf("***error  子供が多すぎる (C_copy)\n");
			exit(1);
		}
	}
/*
     実行
*/
	for (i1 = 0; i1 < pair; i1++) {
					// 親の選択
		p1 = Select(k_method, k_bias, k_step);
		sw = 0;

		while (sw == 0) {
			p2 = Select(k_method, k_bias, k_step);
			if (p1 != p2)
				sw = 1;
		}
					// コピー
		for (i2 = 0; i2 < method; i2++) {
			p       = (i2 == 0) ? p1 : p2;
			k       = Position(-1);
			len[k]  = len[p];
			pi_w[k] = 1;
			for (i3 = 0; i3 < len[k]; i3++)
				ind[k][i3] = ind[p][i3];
		}
	}
}

/*******************************************************************/
/* 交叉（多点交叉）                                                */
/*      kosa : 交叉確率                                            */
/*      k_point : 交叉点の数(default=1)                            */
/*                （負の時は，1から-k_point間のランダム）          */
/*      k_vr : =0 : 両親とも同じ位置で交叉(default)                */
/*             =1 : 両親が異なる位置で交叉（遺伝子長は可変）       */
/*      k_method : 選択方法                                        */
/*                 =-1 : ランダム(default)                         */
/*                 =0 : 適応度をそのまま使用                       */
/*                 =1 : 最小値からの差（ただし，α以下の場合はα） */
/*                 =2 : 評価値に順位をつけ，減少率βで線形化       */
/*      k_bias : α，または，method=2の場合は初期値(default=0)     */
/*      k_step : β(default=1)                                     */
/*******************************************************************/
void Species::C_point(double kosa, int k_point, int k_vr, int k_method,
                      double k_bias, double k_step)
{
	int abs_p, c1, c2, i1, i2, i3, k1, k2, mn = 0, num, p1, p2, pair,
        sw, t11, t12, t21, t22;
/*
     初期設定とデータのチェック
*/
	pair = max_ch / 2;

	if (dup_a == 0) {
		printf("***error  交叉方法が不適当 (C_point)\n");
		exit(1);
	}

	abs_p = abs(k_point);
	if (abs_p == 0 || abs_p > max_len-1 || min_len > 0 && abs_p > min_len-1) {
		printf("***error  交叉点の数が不適当 (C_point)\n");
		exit(1);
	}

	if (k_vr > 0 && min_len < 0) {
		printf("***error  遺伝子長は可変でなければならない (C_point)\n");
		exit(1);
	}
/*
     交叉
*/
	num = k_point;

	for (i1 = 0; i1 < pair; i1++) {
					// 交叉しない場合
		if (genrand_real3() > kosa)
			C_copy(2, 1);
					// 交叉する場合
		else {
						// 親の選択
			p1 = Select(k_method, k_bias, k_step);
			sw = 0;
			while (sw == 0) {
				p2 = Select(k_method, k_bias, k_step);
				if (p1 != p2)
					sw = 1;
			}
						// 交叉位置の数の決定
			if (k_point < 0) {
				num = (int)(genrand_real3() * abs_p + 1);
				if (num > abs_p)
					num = abs_p;
			}
						// 交叉位置の決定（点の後ろで交叉）
			for (i2 = 0; i2 < num; i2++) {
							// 親１の交叉位置
				sw = 0;
				while (sw == 0) {
					sw       = 1;
					kou1[i2] = (int)(genrand_real3() * (len[p1] - 1));
					if (kou1[i2] > len[p1]-2)
						kou1[i2] = len[p1] - 2;
					if (k_vr == 0 && kou1[i2] > len[p2]-2)
						kou1[i2] = len[p2] - 2;
					for (i3 = 0; i3 < i2 && sw > 0; i3++) {
						if (kou1[i3] == kou1[i2])
							sw = 0;
					}
				}
							// 親２の交叉位置
				if (k_vr > 0) {
					sw = 0;
					while (sw == 0) {
						sw       = 1;
						kou2[i2] = (int)(genrand_real3() * (len[p2] - 1));
						if (kou2[i2] > len[p2]-2)
							kou2[i2] = len[p2] - 2;
						for (i3 = 0; i3 < i2 && sw > 0; i3++) {
							if (kou2[i3] == kou2[i2])
								sw = 0;
						}
					}
				}
			}
						// 交叉の実行
						//   親１のt11からt12を子１のc1へコピー
						//   親２のt21からt22を子２のc2へコピー
						//     次は，
						//   親１のt11からt12を子２のc2へコピー
						//   親２のt21からt22を子１のc1へコピー
						//     ・・・・・
			c1  = 0;
			c2  = 0;
			t11 = 0;
			t21 = 0;
							// 遺伝子長
			k1       = Position(-1);
			pi_w[k1] = 1;
			len[k1]  = len[p1];
			k2       = Position(-1);
			pi_w[k2] = 1;
			len[k2]  = len[p2];

			for (i2 = 0; i2 < num+1; i2++ ) {
							// 次の交叉位置を求める
				if (i2 == num) {            // 最後
					t12 = len[p1];
					t22 = len[p2];
				}
				else {
								// 親１
					t12 = max_len;
					for (i3 = 0; i3 < num; i3++) {
						if (kou1[i3] >= 0 && kou1[i3] <= t12) {
							t12 = kou1[i3];
							mn  = i3;
						}
					}
					kou1[mn] = -1;
					t12++;
								// 親２
					if (k_vr == 0)
						t22 = t12;
					else {
						t22 = max_len;
						for (i3 = 0; i3 < num; i3++) {
							if (kou2[i3] >= 0 && kou2[i3] <= t22) {
								t22 = kou2[i3];
								mn  = i3;
							}
						}
						kou2[mn] = -1;
						t22++;
					}
				}
							// 指定箇所のコピー
				for (i3 = t11; i3 < t12; i3++) {
					if (i2%2 == 0) {
						if (c1 < max_len) {
							ind[k1][c1] = ind[p1][i3];
							c1++;
						}
					}
					else {
						if (c2 < max_len) {
							ind[k2][c2] = ind[p1][i3];
							c2++;
						}
					}
				}

				for (i3 = t21; i3 < t22; i3++) {
					if (i2%2 == 0) {
						if (c2 < max_len) {
							ind[k2][c2] = ind[p2][i3];
							c2++;
						}
					}
					else {
						if (c1 < max_len) {
							ind[k1][c1] = ind[p2][i3];
							c1++;
						}
					}
				}
							// 交叉位置の移動
				t11 = t12;
				t21 = t22;
			}
		}
	}
}

/*******************************************************************/
/* 交叉（一様交叉．[0,1]を等確率で発生させ，1であれば，            */
/*       親１，0であれば親２の遺伝子を子１が受け継ぐ）             */
/*      kosa : 交叉確率                                            */
/*      k_method : 選択方法                                        */
/*                 =-1 : ランダム(default)                         */
/*                 =0 : 適応度をそのまま使用                       */
/*                 =1 : 最小値からの差（ただし，α以下の場合はα） */
/*                 =2 : 評価値に順位をつけ，減少率βで線形化       */
/*      k_bias : α，または，method=2の場合は初期値(default=0)     */
/*      k_step : β(default=1)                                     */
/*******************************************************************/
void Species::C_uniform(double kosa, int k_method, double k_bias,
                        double k_step)
{
	int i1, i2, k1, k2, p1, p2, pair, sw;
/*
     初期設定とデータのチェック
*/
	pair = max_ch / 2;

	if (dup_a == 0) {
		printf("***error  交叉方法が不適当 (C_uniform)\n");
		exit(1);
	}

	if (min_len > 0) {
		printf("***error  遺伝子長は固定長でなければならない (C_uniform)\n");
		exit(1);
	}
/*
     交叉
*/
	for (i1 = 0; i1 < pair; i1++) {
					// 交叉しない場合
		if (genrand_real3() > kosa)
			C_copy(2, 1);
					// 交叉する場合
		else {
						// 親の選択
			p1 = Select(k_method, k_bias, k_step);
			sw = 0;
			while (sw == 0) {
				p2 = Select(k_method, k_bias, k_step);
				if (p1 != p2)
					sw = 1;
			}
						// 遺伝子長
			k1       = Position(-1);
			pi_w[k1] = 1;
			len[k1]  = len[p1];
			k2       = Position(-1);
			pi_w[k2] = 1;
			len[k2]  = len[p2];
						// 交叉
			for (i2 = 0; i2 < len[p1]; i2++) {
				if (genrand_real3() > 0.5) {
					ind[k1][i2] = ind[p1][i2];
					ind[k2][i2] = ind[p2][i2];
				}
				else {
					ind[k1][i2] = ind[p2][i2];
					ind[k2][i2] = ind[p1][i2];
				}
			}
		}
	}
}

/*******************************************************************/
/* 交叉（平均化交叉．２つの親の平均値を受け継ぐ）                  */
/*      kosa : 交叉確率                                            */
/*      k_method : 選択方法                                        */
/*                 =-1 : ランダム(default)                         */
/*                 =0 : 適応度をそのまま使用                       */
/*                 =1 : 最小値からの差（ただし，α以下の場合はα） */
/*                 =2 : 評価値に順位をつけ，減少率βで線形化       */
/*      k_bias : α，または，method=2の場合は初期値(default=0)     */
/*      k_step : β(default=1)                                     */
/*******************************************************************/
void Species::C_mean(double kosa, int k_method, double k_bias,
                     double k_step)
{
	int i1, i2, k, p1, p2, sw;
/*
     初期設定とデータのチェック
*/
	if (min_len > 0) {
		printf("***error  遺伝子長は固定長でなければならない (C_mean)\n");
		exit(1);
	}
/*
     交叉
*/
	for (i1 = 0; i1 < max_ch; i1++) {
					// 交叉しない場合
		if (genrand_real3() > kosa)
			C_copy(1, 1);
					// 交叉する場合
		else {
						// 親の選択
			p1 = Select(k_method, k_bias, k_step);
			sw = 0;
			while (sw == 0) {
				p2 = Select(k_method, k_bias, k_step);
				if (p1 != p2)
					sw = 1;
			}
						// 遺伝子長
			k       = Position(-1);
			len[k]  = len[p1];
			pi_w[k] = 1;
						// 交叉
			for (i2 = 0; i2 < len[k]; i2++)
				ind[k][i2] = (ind[p1][i2] + ind[p2][i2]) / 2;
		}
	}
}

/*******************************************************************/
/* 交叉（循環交叉．ランダムに１点を選択し，その位置にある遺伝子を  */
/*       そのまま各子供が選択する．その位置にある親２（１）の遺伝  */
/*       子を，その遺伝子の親１（２）の場所に，子１（２）が受け継  */
/*       ぐ（ただし，doubleの場合は，この手続きをのぞく）．この手  */
/*       続きを，すでに受け継いだ遺伝子の位置が選択されるまで繰り  */
/*       返し，残りの遺伝子については，子１（２）は，親２（１）の  */
/*       遺伝子をその順番通りに受け継ぐ）                          */
/*         2 4 1 3 6 5    + + 1 + + 5    3 2 1 4 6 5               */
/*             *       →             →                           */
/*         3 2 5 4 1 6    + + 5 + 1 +    2 4 5 3 1 6               */
/*      kosa : 交叉確率                                            */
/*      k_method : 選択方法                                        */
/*                 =-1 : ランダム(default)                         */
/*                 =0 : 適応度をそのまま使用                       */
/*                 =1 : 最小値からの差（ただし，α以下の場合はα） */
/*                 =2 : 評価値に順位をつけ，減少率βで線形化       */
/*      k_bias : α，または，method=2の場合は初期値(default=0)     */
/*      k_step : β(default=1)                                     */
/*******************************************************************/
void Species::C_cycle(double kosa, int k_method, double k_bias,
                      double k_step)
{
   int i1, i2, i3, k1, k2, p, pair, p1, p2, sw;
/*
     初期設定とデータのチェック
*/
	pair = max_ch / 2;

	if (dup_a != 0) {
		printf("***error  交叉方法が不適当 (C_cycle)\n");
		exit(1);
	}

	if (min_len > 0) {
		printf("***error  遺伝子長は固定長でなければならない (C_cycle)\n");
		exit(1);
	}
/*
     交叉
*/
	for (i1 = 0; i1 < pair; i1++) {
					// 交叉しない場合
		if (genrand_real3() > kosa)
			C_copy(2, 1);
					// 交叉する場合
		else {
						// 親の選択
			p1 = Select(k_method, k_bias, k_step);
			sw = 0;
			while (sw == 0) {
				p2 = Select(k_method, k_bias, k_step);
				if (p1 != p2)
					sw = 1;
			}
						// 初期設定
			for (i2 = 0; i2 < len[p1]; i2++) {
				kou1[i2] = 0;
				kou2[i2] = 0;
			}
						// 遺伝子長
			k1       = Position(-1);
			pi_w[k1] = 1;
			len[k1]  = len[p1];
			k2       = Position(-1);
			pi_w[k2] = 1;
			len[k2]  = len[p2];
						// 交叉
			sw = 0;

			while (sw == 0) {
				sw = 1;
				p  = (int)(genrand_real3() * len[p1]);
				if (p >= len[p1])
					p = len[p1] - 1;
				if (kou1[p] == 0 && kou2[p] == 0) {
					kou1[p]    = 1;
					kou2[p]    = 1;
					ind[k1][p] = ind[p1][p];
					ind[k2][p] = ind[p2][p];
					for (i2 = 0; i2 < len[p1] && sw > 0; i2++) {
						if (ind[p2][p] == ind[p1][i2]) {
							ind[k1][i2] = ind[p1][i2];
							kou1[i2]    = 1;
							sw          = 0;
						}
					}
					sw = 1;
					for (i2 = 0; i2 < len[p2] && sw > 0; i2++) {
						if (ind[p1][p] == ind[p2][i2]) {
							ind[k2][i2] = ind[p2][i2];
							kou2[i2]    = 1;
							sw          = 0;
						}
					}
				}
			}

			sw = 0;
			i2 = 0;
			i3 = 0;
			while (sw == 0) {
				while (sw == 0 && i2 < len[p1]) {
					if (kou1[i2] == 0)
						sw = 1;
					else
						i2++;
				}
				sw = 0;
				while (sw == 0 && i3 < len[p2]) {
					if (kou2[i3] == 0)
						sw = 1;
					else
						i3++;
				}
				if (i2 < len[p1] && i3 < len[p2]) {
					ind[k1][i2] = ind[p2][i3];
					ind[k2][i3] = ind[p1][i2];
					sw          = 0;
					i2++;
					i3++;
				}
				else
					sw = 1;
			}
		}
	}
}

/*******************************************************************/
/* 交叉（部分的交叉．ランダムに１点を選択し，その位置にある親１と  */
/*       親２の遺伝子を取り出す．次に，親１と親２の染色体上で，こ  */
/*       の２つの遺伝子の位置を交換する．この操作を，選択した点よ  */
/*       り右にあるすべての遺伝子に対して実施する                  */
/*         2 4 1 3 6 5    2 4 5 3 6 1                              */
/*             *       →             → ･････                     */
/*         3 2 5 4 1 6    3 2 1 4 5 6                              */
/*      kosa : 交叉確率                                            */
/*      k_method : 選択方法                                        */
/*                 =-1 : ランダム(default)                         */
/*                 =0 : 適応度をそのまま使用                       */
/*                 =1 : 最小値からの差（ただし，α以下の場合はα） */
/*                 =2 : 評価値に順位をつけ，減少率βで線形化       */
/*      k_bias : α，または，method=2の場合は初期値(default=0)     */
/*      k_step : β(default=1)                                     */
/*******************************************************************/
void Species::C_part(double kosa, int k_method, double k_bias,
                     double k_step)
{
	int i1, i2, i3, k1, k2, lv, p, pair, p1, p2, sw;
/*
     初期設定とデータのチェック
*/
	pair = max_ch / 2;

	if (dup_a != 0) {
		printf("***error  交叉方法が不適当 (C_part)\n");
		exit(1);
	}

	if (min_len > 0) {
		printf("***error  遺伝子長は固定長でなければならない (C_part)\n");
		exit(1);
	}
/*
     交叉
*/
	for (i1 = 0; i1 < pair; i1++) {
					// 交叉しない場合
		if (genrand_real3() > kosa)
			C_copy(2, 1);
					// 交叉する場合
		else {
						// 親の選択
			p1 = Select(k_method, k_bias, k_step);
			sw = 0;
			while (sw == 0) {
				p2 = Select(k_method, k_bias, k_step);
				if (p1 != p2)
					sw = 1;
			}
						// 遺伝子長
			k1       = Position(-1);
			pi_w[k1] = 1;
			len[k1]  = len[p1];
			k2       = Position(-1);
			pi_w[k2] = 1;
			len[k2]  = len[p2];
						// 交叉
			p = (int)(genrand_real3() * len[p1]);
			if (p >= len[p1])
				p = len[p1] - 1;

			for (i2 = 0; i2 < len[p1]; i2++) {
				ind[k1][i2] = ind[p1][i2];
				ind[k2][i2] = ind[p2][i2];
			}

			for (i2 = p; i2 < len[p1]; i2++) {
				sw = 0;
				lv = ind[k1][i2];
				for (i3 = 0; i3 < len[p1] && sw == 0; i3++) {
					if (ind[k2][i2] == ind[k1][i3]) {
						ind[k1][i2] = ind[k1][i3];
						ind[k1][i3] = lv;
						sw          = 1;
					}
				}
				sw = 0;
				for (i3 = 0; i3 < len[p1] && sw == 0; i3++) {
					if (lv == ind[k2][i3]) {
						ind[k2][i3] = ind[k2][i2];
						ind[k2][i2] = lv;
						sw          = 1;
					}
				}
			}
		}
	}
}

/*******************************************************************/
/* 交叉（順序交叉．ランダムに切れ目を決定し，子１に対し，切れ目の  */
/*       左側では，親１の遺伝子をそのまま受け継ぎ，右側では，親１  */
/*       の遺伝子を親２の遺伝子の出現順序に並べ替える．            */
/*         2 4 1 3 6 5    2 4 1 3 5 6                              */
/*             *       →                                          */
/*         3 2 5 4 1 6    3 2 5 4 1 6                              */
/*      kosa : 交叉確率                                            */
/*      k_method : 選択方法                                        */
/*                 =-1 : ランダム(default)                         */
/*                 =0 : 適応度をそのまま使用                       */
/*                 =1 : 最小値からの差（ただし，α以下の場合はα） */
/*                 =2 : 評価値に順位をつけ，減少率βで線形化       */
/*      k_bias : α，または，method=2の場合は初期値(default=0)     */
/*      k_step : β(default=1)                                     */
/*******************************************************************/
void Species::C_seq(double kosa, int k_method, double k_bias,
                    double k_step)
{
	int i1, i2, i3, i4, k1, k2, p, pair, pp, p1, p2, sw;
/*
     初期設定とデータのチェック
*/
	pair = max_ch / 2;

	if (dup_a != 0) {
		printf("***error  交叉方法が不適当 (C_seq)\n");
		exit(1);
	}

	if (min_len > 0) {
		printf("***error  遺伝子長は固定長でなければならない (C_seq)\n");
		exit(1);
	}
/*
     交叉
*/
	for (i1 = 0; i1 < pair; i1++) {
					// 交叉しない場合
		if (genrand_real3() > kosa)
			C_copy(2, 1);
					// 交叉する場合
		else {
						// 親の選択
			p1 = Select(k_method, k_bias, k_step);
			sw = 0;
			while (sw == 0) {
				p2 = Select(k_method, k_bias, k_step);
				if (p1 != p2)
					sw = 1;
			}
						// 遺伝子長
			k1       = Position(-1);
			pi_w[k1] = 1;
			len[k1]  = len[p1];
			k2       = Position(-1);
			pi_w[k2] = 1;
			len[k2]  = len[p2];
						// 交叉
			p = (int)(genrand_real3() * (len[p1] - 1));
			if (p >= len[p1]-1)
				p = len[p1] - 2;

			for (i2 = 0; i2 <= p; i2++) {
				ind[k1][i2] = ind[p1][i2];
				ind[k2][i2] = ind[p2][i2];
			}

			pp = 0;
			for (i2 = p+1; i2 < len[p1]; i2++) {
				sw = 0;
				for (i3 = pp; i3 < len[p2] && sw == 0; i3++) {
					for (i4 = p+1; i4 < len[p1] && sw == 0; i4++) {
						if (ind[p2][i3] == ind[p1][i4]) {
							sw          = 1;
							pp          = i3 + 1;
							ind[k1][i2] = ind[p1][i4];
						}
					}
				}
			}
			pp = 0;
			for (i2 = p+1; i2 < len[p2]; i2++) {
				sw = 0;
				for (i3 = pp; i3 < len[p1] && sw == 0; i3++) {
					for (i4 = p+1; i4 < len[p2] && sw == 0; i4++) {
						if (ind[p1][i3] == ind[p2][i4]) {
							sw          = 1;
							pp          = i3 + 1;
							ind[k2][i2] = ind[p2][i4];
						}
					}
				}
			}
		}
	}
}

/*******************************************************************/
/* 交叉（一様順序交叉．位置の集合をランダムに選択し，一方の親の選  */
/*       択された位置における遺伝子の順序に従って，他の親の遺伝子  */
/*       を並べ替える                                              */
/*         2 4 1 3 6 5    2 4 1 3 6 5                              */
/*           *   *     →                                          */
/*         3 2 5 4 1 6    4 2 5 3 1 6                              */
/*      kosa : 交叉確率                                            */
/*      k_method : 選択方法                                        */
/*                 =-1 : ランダム(default)                         */
/*                 =0 : 適応度をそのまま使用                       */
/*                 =1 : 最小値からの差（ただし，α以下の場合はα） */
/*                 =2 : 評価値に順位をつけ，減少率βで線形化       */
/*      k_bias : α，または，method=2の場合は初期値(default=0)     */
/*      k_step : β(default=1)                                     */
/*******************************************************************/
void Species::C_useq(double kosa, int k_method, double k_bias,
                     double k_step)
{
   int i1, i2, i3, i4, k1, k2, p, pair, p1, p2, sw;
/*
     初期設定とデータのチェック
*/
	pair = max_ch / 2;

	if (dup_a != 0) {
		printf("***error  交叉方法が不適当 (C_useq)\n");
		exit(1);
	}

	if (min_len > 0) {
		printf("***error  遺伝子長は固定長でなければならない (C_useq)\n");
		exit(1);
	}
/*
     交叉
*/
	for (i1 = 0; i1 < pair; i1++) {
					// 交叉しない場合
		if (genrand_real3() > kosa)
			C_copy(2, 1);
					// 交叉する場合
		else {
						// 親の選択
			p1 = Select(k_method, k_bias, k_step);
			sw = 0;
			while (sw == 0) {
				p2 = Select(k_method, k_bias, k_step);
				if (p1 != p2)
					sw = 1;
			}
						// 遺伝子長
			k1       = Position(-1);
			pi_w[k1] = 1;
			len[k1]  = len[p1];
			k2       = Position(-1);
			pi_w[k2] = 1;
			len[k2]  = len[p2];
						// 交叉
			for (i2 = 0; i2 < len[p1]; i2++) {
				ind[k1][i2] = ind[p1][i2];
				ind[k2][i2] = ind[p2][i2];
				kou1[i2]    = (genrand_real3() < 0.5) ? 0 : 1;
			}

			p = 0;
			for (i2 = 0; i2 < len[p1]; i2++) {
				if (kou1[i2] > 0) {
					sw = 0;
					for (i3 = p; i3 < len[p2] && sw == 0; i3++) {
						for (i4 = 0; i4 < len[p1] && sw == 0; i4++) {
							if (ind[p2][i3] == ind[p1][i4] && kou1[i4] > 0) {
								sw          = 1;
								p           = i3 + 1;
								ind[k1][i2] = ind[p1][i4];
							}
						}
					}
				}
			}
			p = 0;
			for (i2 = 0; i2 < len[p2]; i2++) {
				if (kou1[i2] > 0) {
					sw = 0;
					for (i3 = p; i3 < len[p1] && sw == 0; i3++) {
						for (i4 = 0; i4 < len[p2] && sw == 0; i4++) {
							if (ind[p1][i3] == ind[p2][i4] && kou1[i4] > 0) {
								sw          = 1;
								p           = i3 + 1;
								ind[k2][i2] = ind[p2][i4];
							}
						}
					}
				}
			}
		}
	}
}

/*******************************************************************/
/* 交叉（一様位置交叉．位置の集合をランダムに選択し，一方の親の選  */
/*       択された位置における遺伝子の位置に，他の親の同じ遺伝子を  */
/*       配置する．残りの遺伝子は，親と同じ順序に配置する．        */
/*         2 4 1 3 6 5    + + 5 + 1 +    2 4 5 3 1 6               */
/*             *   *   →             →                           */
/*         3 2 5 4 1 6    + + 1 + 6 +    3 2 1 5 6 4               */
/*      kosa : 交叉確率                                            */
/*      k_method : 選択方法                                        */
/*                 =-1 : ランダム(default)                         */
/*                 =0 : 適応度をそのまま使用                       */
/*                 =1 : 最小値からの差（ただし，α以下の場合はα） */
/*                 =2 : 評価値に順位をつけ，減少率βで線形化       */
/*      k_bias : α，または，method=2の場合は初期値(default=0)     */
/*      k_step : β(default=1)                                     */
/*******************************************************************/
void Species::C_upos(double kosa, int k_method, double k_bias,
                     double k_step)
{
   int i1, i2, i3, k1, k2, p, pair, p1, p2, sw;
/*
     初期設定とデータのチェック
*/
	pair = max_ch / 2;

	if (dup_a != 0) {
		printf("***error  交叉方法が不適当 (C_upos)\n");
		exit(1);
	}

	if (min_len > 0) {
		printf("***error  遺伝子長は固定長でなければならない (C_upos)\n");
		exit(1);
	}
/*
     交叉
*/
	for (i1 = 0; i1 < pair; i1++) {
					// 交叉しない場合
		if (genrand_real3() > kosa)
			C_copy(2, 1);
					// 交叉する場合
		else {
						// 親の選択
			p1 = Select(k_method, k_bias, k_step);
			sw = 0;
			while (sw == 0) {
				p2 = Select(k_method, k_bias, k_step);
				if (p1 != p2)
					sw = 1;
			}
						// 遺伝子長
			k1       = Position(-1);
			pi_w[k1] = 1;
			len[k1]  = len[p1];
			k2       = Position(-1);
			pi_w[k2] = 1;
			len[k2]  = len[p2];
						// 交叉
			for (i2 = 0; i2 < len[p1]; i2++) {
				kou1[i2] = (genrand_real3() < 0.5) ? 0 : 1;
				if (kou1[i2] > 0) {
					ind[k1][i2] = ind[p2][i2];
					ind[k2][i2] = ind[p1][i2];
				}
			}

			p = 0;
			for (i2 = 0; i2 < len[p1]; i2++) {
				sw = 0;
				for (i3 = 0; i3 < len[p1] && sw == 0; i3++) {
					if (kou1[i3] > 0 && ind[p1][i2] == ind[k1][i3])
						sw = 1;
				}
				if (sw == 0) {
					for (i3 = p; i3 < len[p1] && sw == 0; i3++) {
						if (kou1[i3] == 0) {
							ind[k1][i3] = ind[p1][i2];
							p           = i3 + 1;
							sw          = 1;
						}
					}
				}
			}
			p = 0;
			for (i2 = 0; i2 < len[p2]; i2++) {
				sw = 0;
				for (i3 = 0; i3 < len[p2] && sw == 0; i3++) {
					if (kou1[i3] > 0 && ind[p2][i2] == ind[k2][i3])
						sw = 1;
				}
				if (sw == 0) {
					for (i3 = p; i3 < len[p2] && sw == 0; i3++) {
						if (kou1[i3] == 0) {
							ind[k2][i3] = ind[p2][i2];
							p           = i3 + 1;
							sw          = 1;
						}
					}
				}
			}
		}
	}
}

/*******************************************************************/
/* 交叉（エッジ組み替え交叉．以下の手順に従って行う．対立遺伝子は  */
/*       0～(max_len-1)である必要がある）                          */
/*         (0) エッジマップを作成する．エッジマップとは，２つの親  */
/*             を見て，ノードがどこに接続されているのかを表すもの  */
/*             であり，例えば，２つの親が，                        */
/*                 [A B C D E F]                                   */
/*                 [B D C A E F]                                   */
/*             である場合は，                                      */
/*                 A : B F C E                                     */
/*                 B : A C D F                                     */
/*                 C : B D A                                       */
/*                 D : C E B                                       */
/*                 E : D F A                                       */
/*                 F : A E B                                       */
/*             となる．                                            */
/*         (1) 両親の２つの出発点の内１つで初期化する．ランダムま  */
/*             たはステップ(4)の基準に従って選ぶ（現在のノード）   */
/*         (2) エッジマップから，現在のノードを除く                */
/*         (3) 現在のノードが接続先のノードを持っていたら，(4)に   */
/*             進む．さもなければ，(5)に進む                       */
/*         (4) 現在のノードが持っている接続先ノードの内，最も少な  */
/*             い接続先ノードを持ったノードを選択し（同じ条件の場  */
/*             合は，ランダム），それを現在のノードとし，(2)に進む */
/*         (5) 未接続のノードが残っていればランダムに選択し，(2)に */
/*             戻る．さもなければ，終了する                        */
/*      kosa : 交叉確率                                            */
/*      k_method : 選択方法                                        */
/*                 =-1 : ランダム(default)                         */
/*                 =0 : 適応度をそのまま使用                       */
/*                 =1 : 最小値からの差（ただし，α以下の場合はα） */
/*                 =2 : 評価値に順位をつけ，減少率βで線形化       */
/*      k_bias : α，または，method=2の場合は初期値(default=0)     */
/*      k_step : β(default=1)                                     */
/*******************************************************************/
void Species::C_edge(double kosa, int k_method, double k_bias,
                     double k_step)
{
	int e[2], i1, i2, i3, i4, i5, k, kk, k0 = 0, k1, k2, min, num,
        p, pair, p1, p2, sw;
/*
     初期設定とデータのチェック
*/
	pair = max_ch;

	if (dup_a != 0) {
		printf("***error  交叉方法が不適当 (C_edge)\n");
		exit(1);
	}

	if (min_len > 0) {
		printf("***error  遺伝子長は固定長でなければならない (C_edge)\n");
		exit(1);
	}
/*
     交叉
*/
	for (i1 = 0; i1 < pair; i1++) {
					// 交叉しない場合
		if (genrand_real3() > kosa)
			C_copy(1, 1);
					// 交叉する場合
		else {
						// 親の選択
			p1 = Select(k_method, k_bias, k_step);
			sw = 0;
			while (sw == 0) {
				p2 = Select(k_method, k_bias, k_step);
				if (p1 != p2)
					sw = 1;
			}
						// 遺伝子長
			k       = Position(-1);
			pi_w[k] = 1;
			len[k]  = len[p1];
						// エッジマップの初期化
			for (i2 = 0; i2 < len[k]; i2++) {
				edge[i2][0] = 0;
				for (i3 = 1; i3 <= 4; i3++)
					edge[i2][i3] = -1;
			}
						// 交叉
							// エッジマップの作成
			for (i2 = 0; i2 < len[k]; i2++) {

				sw = 0;
				for (i3 = 0; i3 < len[k] && sw == 0; i3++) {
					if (i2 == ind[p1][i3]) {
						sw = 1;
						if (i3 == 0) {
							e[0] = ind[p1][len[k]-1];
							e[1] = ind[p1][1];
						}
						else {
							if (i3 == len[k]-1) {
								e[0] = ind[p1][i3-1];
								e[1] = ind[p1][0];
							}
							else {
								e[0] = ind[p1][i3-1];
								e[1] = ind[p1][i3+1];
							}
						}
						for (i4 = 0; i4 < 2; i4++) {
							edge[i2][0]++;
							edge[i2][edge[i2][0]] = e[i4];
						}
					}
				}

				sw = 0;
				for (i3 = 0; i3 < len[k] && sw == 0; i3++) {
					if (i2 == ind[p2][i3]) {
						sw = 1;
						if (i3 == 0) {
							e[0] = ind[p2][len[k]-1];
							e[1] = ind[p2][1];
						}
						else {
							if (i3 == len[k]-1) {
								e[0] = ind[p2][i3-1];
								e[1] = ind[p2][0];
							}
							else {
								e[0] = ind[p2][i3-1];
								e[1] = ind[p2][i3+1];
							}
						}
						for (i4 = 0; i4 < 2; i4++) {
							sw = 1;
							for (i5 = 1; i5 <= edge[i2][0] && sw == 1; i5++) {
								if (edge[i2][i5] == e[i4])
									sw = 2;
							}
							if (sw == 1) {
								edge[i2][0]++;
								edge[i2][edge[i2][0]] = e[i4];
							}
						}
					}
				}
			}
							// 交叉の実行
								// 出発点の決定
			k1 = ind[p1][0];
			k2 = ind[p2][0];
			if (edge[k1][0] == edge[k2][0])
				kk = (genrand_real3() > 0.5) ? k2 : k1;
			else
				kk = (edge[k1][0] < edge[k2][0]) ? k1 : k2;
			ind[k][0] = kk;
			p         = 1;

			while (p < len[k]) {
								// ノードの除去
				for (i2 = 0; i2 < len[k]; i2++) {
					sw = 0;
					if (edge[i2][0] > 0) {
						for (i3 = 1; i3 <= 4 && sw == 0; i3++) {
							if (edge[i2][i3] == kk) {
								sw           = 1;
								edge[i2][i3] = -1;
								edge[i2][0]--;
							}
						}
					}
				}
								// 次の現在ノードの選択
				min = 10;
				num = 0;
				for (i2 = 1; i2 <= 4; i2++) {
					if (edge[kk][i2] >= 0) {
						k1 = edge[kk][i2];
						if (edge[k1][0] >= 0 && edge[k1][0] < min) {
							num = 1;
							min = edge[k1][0];
							k0  = k1;
						}
						else {
							if (edge[k1][0] == min)
								num++;
						}
					}
				}
				if (num > 1) {
					k1 = (int)(genrand_real3() * num) + 1;
					if (k1 > num)
						k1 = num;
					k2 = 0;
					k0 = -1;
					for (i2 = 1; i2 <= 4 && k0 < 0; i2++) {
						if (edge[kk][i2] >= 0) {
							if (edge[edge[kk][i2]][0] == min) {
								k2++;
								if (k1 == k2)
									k0 = edge[kk][i2];
							}
						}
					}
				}
				else {
					if (num <= 0) {
						num = 0;
						for (i2 = 0; i2 < len[k]; i2++) {
							if (i2 != kk && edge[i2][0] >= 0)
								num++;
						}
						if (num <= 0) {
							printf("***error  invalid data (C_edge)\n");
							exit(1);
						}
						else {
							k1 = (int)(genrand_real3() * num) + 1;
							if (k1 > num)
								k1 = num;
							k2 = 0;
							k0 = -1;
							for (i2 = 0; i2 < len[k] && k0 < 0; i2++) {
								if (i2 != kk && edge[i2][0] >= 0) {
									k2++;
									if (k1 == k2)
										k0 = i2;
								}
							}
						}
					}
				}
				edge[kk][0] = -1;
				ind[k][p]   = k0;
				kk          = k0;
				p++;
			}
		}
	}
}

/*************************************************************/
/* 交叉（サブツアー交叉．２点交叉の拡張である．ただし，相手に*/
/*       同じ遺伝子のグループがない限り実行されない．たとえば*/
/*         ***abcd**                                         */
/*         *cdab****                                         */
/*       のような両親の時実行され，以下の４つの子供が生成され*/
/*       る）                                                */
/*         ***cdab**                                         */
/*         *abcd****                                         */
/*         ***badc**                                         */
/*         *dcba****                                         */
/*       最大，４＊交叉回数＊個体総数＊(個体総数－１) 個の子 */
/*       供が生成される可能性があるので，子供の数としてこの値*/
/*       以上のデータを入力しておく必要がある．              */
/*      kosa : 交叉確率                                      */
/*      count : １つのペアーに対する交差回数(default=10)     */
/*************************************************************/
void Species::C_sub(double kosa, int count)
{
	int i1, i2, i3, i4, i5, k1, k2, k3, k4, p1, p2,
        t11, t12, t21, t22 = 0, sw;
/*
     初期設定とデータのチェック
*/
	if ((4*count*size*(size-1)) > max_ch) {
		printf("***error  子供が多すぎる (C_sub)\n");
		exit(1);
	}
/*
     交叉
*/
	for (i1 = 0; i1 < size-1; i1++) {
					// 親１
		p1 = Position(i1);

		if (p1 >= 0) {

			for (i2 = i1; i2 < size; i2++) {
					// 親２
				p2 = Position(i2);

				if (p2 >= 0) {
					// 交叉しない場合
					if (genrand_real3() > kosa)
						C_copy(2, 1);
					// 交叉する場合
					else {
						// 交叉回数の制御
						for (i3 = 0; i3 < count; i3++) {
							// 交叉位置の決定（点の後ろで交叉）
								// 親１の交叉位置
							t11 = (int)(genrand_real3() * len[p1]);
							if (t11 > (len[p1]-1))
								t11 = len[p1] - 1;
							sw = 0;
							while (sw == 0) {
								t12 = (int)(genrand_real3() * len[p1]);
								if (t12 > (len[p1]-1))
									t12 = len[p1] - 1;
								if (t12 != t11)
									sw = 1;
							}
							if (t11 > t12) {
								k1  = t11;
								t11 = t12;
								t12 = k1;
							}
								// 親２の交叉位置
							sw  = 0;
							t21 = -1;
							for (i4 = 0; i4 < len[p2] && t21 < 0; i4++) {
								for (i5 = t11; i5 <= t12 && t21 < 0; i5++) {
									if (ind[p2][i4] == ind[p1][i5])
										t21 = i4;
								}
							}
							if (t21 >= 0) {
								t22 = t21 + t12 - t11;
								if (t22 < len[p2]) {
									sw = 1;
									for (i4 = t21+1; i4 <= t22 && sw > 0; i4++) {
										sw = 0;
										for (i5 = t11; i5 <= t12 && sw == 0; i5++) {
											if (ind[p2][i4] == ind[p1][i5])
												sw = 1;
										}
									}
								}
							}
								// 交叉の実行
							if (sw > 0) {

								k1       = Position(-1);
								pi_w[k1] = 1;
								len[k1]  = len[p1];
								k2       = Position(-1);
								pi_w[k2] = 1;
								len[k2]  = len[p1];
								k3       = Position(-1);
								pi_w[k3] = 1;
								len[k3]  = len[p2];
								k4       = Position(-1);
								pi_w[k4] = 1;
								len[k4]  = len[p2];

								for (i4 = 0; i4 < t11; i4++) {
									ind[k1][i4] = ind[p1][i4];
									ind[k2][i4] = ind[p1][i4];
								}
								for (i4 = t11; i4 <= t12; i4++) {
									ind[k1][i4] = ind[p2][t21+i4-t11];
									ind[k2][i4] = ind[p2][t22-i4+t11];
								}
								for (i4 = t12+1; i4 < len[p1]; i4++) {
									ind[k1][i4] = ind[p1][i4];
									ind[k2][i4] = ind[p1][i4];
								}
								for (i4 = 0; i4 < t21; i4++) {
									ind[k3][i4] = ind[p2][i4];
									ind[k4][i4] = ind[p2][i4];
								}
								for (i4 = t21; i4 <= t22; i4++) {
									ind[k3][i4] = ind[p1][t11+i4-t21];
									ind[k4][i4] = ind[p1][t12-i4+t21];
								}
								for (i4 = t22+1; i4 < len[p2]; i4++) {
									ind[k3][i4] = ind[p2][i4];
									ind[k4][i4] = ind[p2][i4];
								}
							}
						}
					}
				}
			}
		}
	}
}

/**************************************/
/* 突然変異（対立遺伝子との置き換え） */
/*      pr : 突然変異率               */
/**************************************/
void Species::M_alle(double pr)
{
	int i1, i2, lid;
/*
     データのチェックと初期設定
*/
	if (dup_a == 0) {
		printf("***error  突然変異方法が不適当 (M_alle)\n");
		exit(1);
	}
/*
     実行
*/
	for (i1 = 0; i1 < size+max_ch; i1++) {
		if (pi_w[i1] == 1) {
			for (i2 = 0; i2 < len[i1]; i2++) {
				if (genrand_real3() <= pr) {
					lid = (int)(genrand_real3() * (allele_u - allele_l + 1) + allele_l);
					if (lid > allele_u)
						lid = allele_u;
					if (lid != ind[i1][i2])
						ind[i1][i2] = lid;
				}
			}
		}
	}
}

/**********************************************************************/
/* 突然変異（移動．２点を選択し，２番目の遺伝子を１番目の遺伝子の前に */
/*           移動する）                                               */
/*      pr : 突然変異率                                               */
/**********************************************************************/
void Species::M_move(double pr)
{
	int i1, i2, ld, p1, p2, sw;

	for (i1 = 0; i1 < size+max_ch; i1++) {

		if (pi_w[i1] == 1 && genrand_real3() <= pr) {
/*
     位置の決定
*/
					// p1
			p1 = (int)(genrand_real3() * len[i1]);
			if (p1 >= len[i1])
				p1 = len[i1] - 1;
					// p2
			sw = 0;
			while (sw == 0) {
				p2 = (int)(genrand_real3() * len[i1]);
				if (p2 >= len[i1])
					p2 = len[i1] - 1;
				if (p2 != p1)
					sw = 1;
			}
/*
     実行
*/
			if (p2 > p1) {
				ld = ind[i1][p2];
				for (i2 = p2; i2 > p1; i2--)
					ind[i1][i2] = ind[i1][i2-1];
				ind[i1][p1] = ld;
			}
			else {
				ld = ind[i1][p2];
				for (i2 = p2; i2 < p1-1; i2++)
					ind[i1][i2] = ind[i1][i2+1];
				ind[i1][p1-1] = ld;
			}
		}
	}
}

/********************************************************/
/* 突然変異（逆位．２点間の遺伝子順序を逆に並べ替える） */
/*      pr : 突然変異率                                 */
/*      wd : >0 : 幅を固定                              */
/*           =0 : 幅をランダム(default)                 */
/********************************************************/
void Species::M_inv(double pr, int wd)
{
	int i1, lid, p, p1, p2, sw;

	for (i1 = 0; i1 < size+max_ch; i1++) {

		if (pi_w[i1] == 1 && genrand_real3() <= pr) {
/*
     区間の決定
*/
			if (wd == 0) {
				p1 = (int)(genrand_real3() * len[i1]);
				if (p1 >= len[i1])
					p1 = len[i1] - 1;
				sw = 0;
				while (sw == 0) {
					p2 = (int)(genrand_real3() * len[i1]);
					if (p2 >= len[i1])
						p2 = len[i1] - 1;
					if (p2 != p1)
						sw = 1;
				}
				if (p1 > p2) {
					p  = p1;
					p1 = p2;
					p2 = p;
				}
			}

			else {
				p1 = len[i1];
				while (p1 > len[i1]-2)
					p1 = (int)(genrand_real3() * len[i1]);
				p2 = p1 + wd - 1;
				if (p2 >= len[i1])
					p2 = len[i1] - 1;
			}
/*
     実行
*/
			sw = 0;
			while (sw == 0) {
				lid         = ind[i1][p1];
				ind[i1][p1] = ind[i1][p2];
				ind[i1][p2] = lid;
				p1++;
				p2--;
				if (p1 >= p2)
					sw = 1;
			}
		}
	}
}

/**********************************************************************/
/* 突然変異（スクランブル．２点間の遺伝子順序をランダムに並べ替える） */
/*      pr : 突然変異率                                               */
/*      wd : >0 : 幅を固定                                            */
/*           =0 : 幅をランダム(default)                               */
/**********************************************************************/
void Species::M_scram(double pr, int wd)
{
	int i1, i2, ld, p, p1, p2, sw;

	for (i1 = 0; i1 < size+max_ch; i1++) {

		if (pi_w[i1] == 1 && genrand_real3() <= pr) {
/*
     区間の決定
*/
			if (wd == 0) {
				p1 = (int)(genrand_real3() * len[i1]);
				if (p1 >= len[i1])
					p1 = len[i1] - 1;
				sw = 0;
				while (sw == 0) {
					p2 = (int)(genrand_real3() * len[i1]);
					if (p2 >= len[i1])
						p2 = len[i1] - 1;
					if (p2 != p1)
						sw = 1;
				}
				if (p1 > p2) {
					p  = p1;
					p1 = p2;
					p2 = p;
				}
			}

			else {
				p1 = len[i1];
				while (p1 > len[i1]-2)
					p1 = (int)(genrand_real3() * len[i1]);
				p2 = p1 + wd - 1;
				if (p2 >= len[i1])
					p2 = len[i1] - 1;
			}
/*
     実行
*/
			for (i2 = p1; i2 <= p2; i2++) {
				p = (int)(genrand_real3() * (p2 - p1 + 1) + p1);
				if (p > p2)
					p = p2;
				ld          = ind[i1][i2];
				ind[i1][i2] = ind[i1][p];
				ind[i1][p]  = ld;
			}
		}
	}
}

/**********************************************************************/
/* 突然変異（転座．２点間の遺伝子を他の位置のものと置き換える．ただし */
/*           重複部分はそのままとする）                               */
/*      pr : 突然変異率                                               */
/*      wd : >0 : 幅を固定                                            */
/*           =0 : 幅をランダム(default)                               */
/**********************************************************************/
void Species::M_chg(double pr, int wd)
{
	int i1, i2, ld, p, p1, p2, p3, p4, sw;

	for (i1 = 0; i1 < size+max_ch; i1++) {

		if (pi_w[i1] == 1 && genrand_real3() <= pr) {
/*
     区間等の決定（[p1,p2]と[p3,p4]の入れ替え）
*/
					// p1
			p1 = (int)(genrand_real3() * len[i1]);
			if (p1 >= len[i1])
				p1 = len[i1] - 1;
					// p3
			sw = 0;
			while (sw == 0) {
				p3 = (int)(genrand_real3() * len[i1]);
				if (p3 >= len[i1])
					p3 = len[i1] - 1;
				if (p3 != p1)
					sw = 1;
			}
					// 小さい方をp1,p2にする
			if (p1 > p3) {
				p  = p1;
				p1 = p3;
				p3 = p;
			}
					// p4, p2
			p4 = (wd == 0) ? (int)(genrand_real3() * (len[i1] - p3)) + p3 :
                             p1 + wd - 1;
			if (p4 >= len[i1])
				p4 = len[i1] - 1;
			p2 = p1 + (p4 - p3);
					// 重複部分のチェック
			if (p2 >= p3) {
				p  = p3 - 1;
				p3 = p2 + 1;
				p2 = p;
				p4 = p3 + (p2 - p1);
			}
/*
     実行
*/
			p = p3;
			for (i2 = p1; i2 <= p2; i2++) {
				ld          = ind[i1][i2];
				ind[i1][i2] = ind[i1][p];
				ind[i1][p]  = ld;
				p++;
			}
		}
	}
}

/**********************************************************************/
/* 突然変異（重複．２点間の遺伝子を他の位置にコピーする               */
/*      pr : 突然変異率                                               */
/*      wd : >0 : 幅を固定                                            */
/*           =0 : 幅をランダム(deafult)                               */
/**********************************************************************/
void Species::M_dup(double pr, int wd)
{
	int i1, i2, p, p1, p2, p3, p4, sw;
/*
     データのチェック
*/
	if (dup_a == 0) {
		printf("***error  突然変異方法が不適当 (M_dup)\n");
		exit(1);
	}
/*
     実行
*/
	for (i1 = 0; i1 < size+max_ch; i1++) {

		if (pi_w[i1] == 1 && genrand_real3() <= pr) {
					// 区間の決定（[p1,p2]を[p3,p4]にコピー）
						// p1
			p1 = (int)(genrand_real3() * len[i1]);
			if (p1 >= len[i1])
				p1 = len[i1] - 1;
						// p3
			sw = 0;
			while (sw == 0) {
				p3 = (int)(genrand_real3() * len[i1]);
				if (p3 >= len[i1])
					p3 = len[i1] - 1;
				if (p3 != p1)
					sw = 1;
			}
						// 区間を決める
			if (p3 > p1) {
				p4 = (wd == 0) ? (int)(genrand_real3() * (len[i1] - p3)) + p3 :
                                 p3 + wd - 1;
				if (p4 >= len[i1])
					p4 = len[i1] - 1;
				p2 = p1 + (p4 - p3);
			}
			else {
				p2 = (wd == 0) ? (int)(genrand_real3() * (len[i1] - p1)) + p1 :
                                 p1 + wd - 1;
				if (p2 >= len[i1])
					p2 = len[i1] - 1;
				p4 = p3 + (p2 - p1);
			}
					// 実行
			p = p4;
			for (i2 = p2; i2 >= p1; i2--) {
				ind[i1][p] = ind[i1][i2];
				p--;
			}
		}
	}
}

/******************************************************/
/* 突然変異（摂動．値をある量だけ変化させる）         */
/*      pr : 突然変異率                               */
/*      method : =0 : 正規分布(default)               */
/*               =1 : 一様分布                        */
/*      m : 平均または一様分布の下限(default=0.0)     */
/*      s : 標準偏差または一様分布の上限(default=1.0) */
/******************************************************/
void Species::M_per(double pr, int method, double m, double s)
{
	double w, wd = 0.0, x1;
	int i1, i2;
/*
     データのチェックと初期設定
*/
	if (dup_a == 0) {
		printf("***error  突然変異方法が不適当 (M_per)\n");
		exit(1);
	}

	if (method > 0)
		wd = s - m;
/*
     実行
*/
	for (i1 = 0; i1 < size+max_ch; i1++) {
		if (pi_w[i1] == 1) {
			for (i2 = 0; i2 < len[i1]; i2++) {
				if (genrand_real3() <= pr) {
					if (method == 0)
						w = norm_d(m, s);
					else {
						w = genrand_real3() * wd;
						if (genrand_real3() < 0.5)
							w = -w;
					}
					x1 = (double)ind[i1][i2] + w;
					if (x1 > allele_u)
						x1 = allele_u;
					else {
						if (x1 < allele_l)
							x1 = allele_l;
					}
					ind[i1][i2] = (int)x1;
				}
			}
		}
	}
}

/**********************************************************************/
/* 突然変異（挿入．ある長さの遺伝子を挿入する）                       */
/*      pr : 突然変異率                                               */
/*      wd : >0 : 幅を固定                                            */
/*           =0 : 幅をランダム(default)                               */
/**********************************************************************/
void Species::M_ins(double pr, int wd)
{
	int i1, i2, l, ld, p;
/*
     データのチェック
*/
	if (dup_a == 0 || min_len < 0) {
		printf("***error  突然変異方法が不適当 (M_ins)\n");
		exit(1);
	}
/*
     実行
*/
	for (i1 = 0; i1 < size+max_ch; i1++) {

		if (pi_w[i1] == 1 && genrand_real3() <= pr) {
					// 挿入位置の決定
			p = (int)(genrand_real3() * (len[i1]+1));
			if (p > len[i1])
				p = len[i1];
					// 挿入する遺伝子長の決定
			l = (wd == 0) ? (int)(genrand_real3() * (max_len - len[i1] + 1)) : wd;
			if (l > max_len-len[i1])
				l = max_len - len[i1];
			else {
				if (l <= 0)
					l = 1;
			}
					// 実行
						// 挿入場所の確保
			if (p < len[i1]) {
				for (i2 = len[i1]+l-1; i2 >= p; i2--)
					ind[i1][i2] = ind[i1][i2-l];
			}
						// 挿入場所の遺伝子の決定
			for (i2 = p; i2 < p+l; i2++) {
				ld = (int)(genrand_real3() * (allele_u - allele_l + 1) + allele_l);
				if (ld > allele_u)
					ld = allele_u;
				ind[i1][i2] = ld;
			}

			len[i1]  += l;
		}
	}
}

/**********************************************************************/
/* 突然変異（削除．ある長さの遺伝子を削除する）                       */
/*      pr : 突然変異率                                               */
/*      wd : >0 : 幅を固定                                            */
/*           =0 : 幅をランダム(default)                               */
/**********************************************************************/
void Species::M_del(double pr, int wd)
{
	int i1, i2, l, max, p;
/*
     データのチェック
*/
	if (dup_a == 0 || min_len < 0) {
		printf("***error  突然変異方法が不適当 (M_del)\n");
		exit(1);
	}
/*
     実行
*/
	for (i1 = 0; i1 < size+max_ch; i1++) {

		if (pi_w[i1] == 1 && genrand_real3() <= pr) {
					// 削除位置の決定
			p = (int)(genrand_real3() * len[i1]);
			if (p >= len[i1])
				p = len[i1] - 1;
					// 削除する遺伝子長の決定
			max = (len[i1]-min_len < len[i1]-p) ? len[i1] - min_len : len[i1] - p;
			l   = (wd == 0) ? (int)(genrand_real3() * max + 1) : wd;
			if (l > max)
				l = max;
					// 実行
			for (i2 = 0; i2 < len[i1]-p-l; i2++)
				ind[i1][p+i2] = ind[i1][p+i2+l];

			len[i1]  -= l;
		}
	}
}

/*********************************************************************/
/* 淘汰（エリート・ルーレット選択）                                  */
/*      elite : エリートで残す個体数(default=0)                      */
/*      s_method : ルーレット板の作成方法(default=1)                 */
/*                   =0 : 適応度をそのまま使用                       */
/*                   =1 : 最小値からの差（ただし，α以下の場合はα） */
/*                   =2 : 評価値に順位をつけ，減少率βで線形化       */
/*      s_bias : α，または，method=2の場合は初期値(default=0)       */
/*      s_step : β(default=1)                                       */
/*********************************************************************/
void Species::S_roul(int elite, int s_method, double s_bias, double s_step)
{
	int count = 0, i1, i2, i3, k = 0, max, n = 0, p, sw;
/*
     値のチェックと初期設定
*/
	if (s_method != 0 && s_method != 2)
		s_method = 1;

	if (elite > size) {
		printf("***error  エリートで残す数が多すぎる (S_roul)\n");
		exit(1);
	}

	if (s_method == 2 && s_step <= 0.0)
		s_step = 1.0;

	for (i1 = 0; i1 < size+max_ch; i1++)
		s_w[i1] = 0;
/*
     重複個体を削除
*/
	if (dup_s == 0) {
		for (i1 = 0; i1 < size+max_ch; i1++) {
			if (pi_w[i1] > 0) {
				for (i2 = i1+1; i2 < size+max_ch; i2++) {
					if (pi_w[i2] > 0 && len[i1] == len[i2]) {
						sw = 0;
						for (i3 = 0; i3 < len[i1] && sw == 0; i3++) {
							if (ind[i1][i3] != ind[i2][i3])
								sw = 1;
						}
						if (sw == 0)
							pi_w[i2] = 0;
					}
				}
			}
		}
	}

	for (i1 = 0; i1 < size+max_ch; i1++) {
		if (pi_w[i1] > 1)
			n++;
	}

	if (n < 0 || dup_s == 0 && n < size) {
		printf("***error  残す個体がない (S_roul)\n");
		exit(1);
	}
/*
     淘汰して残す個体を選ぶ
*/
					// エリートの選択
	sw = 0;

	while (k < elite && k < n && sw == 0) {
		max = -1;
		for (i1 = 0; i1 < size+max_ch; i1++) {
			if (pi_w[i1] > 1 && s_w[i1] == 0) {
				if (max < 0 || pi[i1] > pi[max])
					max = i1;
			}
		}
		if (max < 0)
			sw = 1;
		else {
			s_w[max] = 1;
			k++;
		}
	}
					// ルーレット選択
	while (count < size+max_ch && k < size) {
		p = Select(s_method, s_bias, s_step);
		if (dup_s == 0 && s_w[p] > 0)
			count++;
		else {
			count = 0;
			s_w[p]++;
			k++;
		}
	}
						// 選択に失敗した場合の処理
	if (dup_s == 0 && k < size) {
		for (i1 = 0; i1 < size+max_ch && k < size; i1++) {
			if (pi_w[i1] > 1 && s_w[i1] == 0) {
				s_w[i1] = 1;
				k++;
			}
		}
	}
						// 複数回選択されたものの処理
	for (i1 = 0; i1 < size+max_ch; i1++) {
		if (s_w[i1] == 0)
			pi_w[i1] = 0;
	}

	for (i1 = 0; i1 < size+max_ch; i1++) {
		if (s_w[i1] > 0) {
			if (s_w[i1] > 1) {
				for (i2 = 2; i2 <= s_w[i1]; i2++) {
					k       = Position(-1);
					len[k]  = len[i1];
					pi_w[k] = 2;
					pi[k]   = pi[i1];
					for (i3 = 0; i3 < len[i1]; i3++)
						ind[k][i3] = ind[i1][i3];
				}
			}
		}
	}
}

/***********************************/
/* 正規分布変量の発生              */
/*      m : 平均                   */
/*      s : 標準偏差               */
/*           return : 正規分布変量 */
/***********************************/
double norm_d(double m, double s)
{
	double x = 0.0;
	int i1;

	for (i1 = 0; i1 < 12; i1++)
		x += genrand_real3();

	x = s * (x - 6.0) + m;

	return x;
}

/*******************/
/* クラスTSPの定義 */
/*******************/

class TSP : public Species {

		int max_gen;   // 最大世代交代数
		int kosa_m;   // 交叉方法
                      //   =-1 : 交叉を使用しない
                      //   =0 : 親のコピー
                      //   =1 : 循環交叉
                      //   =2 : 部分的交叉
                      //   =3 : 順序交叉
                      //   =4 : 一様順序交叉
                      //   =5 : 一様位置交叉
                      //   =6 : エッジ組み替え交叉
                      //   =7 : サブツアー交叉
		double kosa;   // 交叉確率
		int k_point;   // 交差点の数（負の時は，1から-k_point間のランダム）
		int k_vr;   // =0 : 両親とも同じ位置で交叉
                    // =1 : 両親が異なる位置で交叉（遺伝子長は可変）
		int k_method;   // 交叉の時の親の選択方法
                        //   =-1 : ランダム
                        //   =0 : 適応度をそのまま使用
                        //   =1 : 最小値からの差（ただし，α以下の場合はα）
                        //   =2 : 評価値に順位をつけ，減少率βで線形化
		double k_bias;   // α，または，method=2の場合は初期値
		double k_step;   // β
		int mute_m;   // 突然変異方法
                      //   =-1 : 突然変異を使用しない
                      //   =0 : 移動
                      //   =1 : 逆位
                      //   =2 : スクランブル
                      //   =3 : 転座
		double mute;   // 突然変異率
		int wd;   // 突然変異に使用する部分遺伝子長
		double m_mean;   // 摂動の平均値
		double m_std;   // 摂動の標準偏差
		int elite;   // エリート選択で残す数
		int s_method;   // ルーレット板の作成方法
                        //   =0 : 適応度をそのまま使用
                        //   =1 : 最小値からの差（ただし，α以下の場合はα）
                        //   =2 : 評価値に順位をつけ，減少率βで線形化
		double s_bias;   // α，または，s_method=2の場合は初期値
		double s_step;   // β
		int out_d;   // 表示間隔
		int out_lvl;   // 出力レベル
                       //   =0 : 最終出力だけ
                       //   n>0 : ｎ世代毎に出力（負の時はファイル）
		int out_m;   // 出力方法
                     //   =0 : すべてを出力
                     //   =1 : 最大適応度の個体だけを出力
		char o_file[100];   // 出力ファイル名
		int **city;   //都市の位置データ
		int n_city;   // 都市の数
		int **rg;   // 都市間の距離
		int kinbo;   // 近傍探索（０：行わない，１：行う）
		int neib;   // 近傍（2 or 3）
		int sel;   // エッジの選択方法
                   //   =0 : 最良のものを選択
                   //   =1 : 最初のものを選択

	public:
					// コンストラクタ
		TSP (char *, char *, long);
					// デストラクタ
		~TSP ();
					// 全体の実行制御
		void Control();
					// 距離の計算
		int Kyori(int, int*);
					// 適応度の計算
		void Adap();
					// 枝の入れ替え
		int Change(int, int *, int *);
					// 近傍の探索
		void Kinbo();
					// 出力
		void Output(int);
};

/***************************************/
/* コンストラクタ                      */
/*      name1 : Species定義ファイル名  */
/*      name2 : TSP定義ファイル名      */
/*      seed : 乱数の初期値            */
/***************************************/
TSP::TSP (char *name1, char *name2, long seed) : Species (name1, seed)
{
	double x, y;
	int i1, i2;
	FILE *in;
					// 基本データの入力
	in = fopen(name2, "r");

	fscanf(in, "%*s %d %*s %d", &out_lvl, &out_m);
	fscanf(in, "%*s %s %*s %d", o_file, &out_d);
	fscanf(in, "%*s %d %*s %lf %*s %d %*s %d %*s %d %*s %lf %*s %lf",
           &kosa_m, &kosa, &k_point, &k_vr, &k_method, &k_bias, &k_step);
	fscanf(in, "%*s %d %*s %lf %*s %d %*s %lf %*s %lf",
           &mute_m, &mute, &wd, &m_mean, &m_std);
	fscanf(in, "%*s %d %*s %d %*s %lf %*s %lf",
           &elite, &s_method, &s_bias, &s_step);
	fscanf(in, "%*s %d %*s %d", &n_city, &max_gen);
	fscanf(in, "%*s %d %*s %d", &kinbo, &neib);
	fscanf(in, "%*s %d", &sel);

	if (kinbo > 0 && neib != 2 && neib != 3) {
		printf("***error  近傍の値が不適当 \n");
		exit(1);
	}

	if (n_city != max_len) {
		printf("***error  都市数が不適当 \n");
		exit(1);
	}
					// 都市の位置データ
	city = new int * [n_city];
	for (i1 = 0; i1 < n_city; i1++) {
		city[i1] = new int [2];
		fscanf(in, "%d %d", &city[i1][0], &city[i1][1]);
	}
					// 距離テーブル
	rg = new int * [n_city];

	for (i1 = 0; i1 < n_city; i1++) {
		rg[i1] = new int [n_city];
		for (i2 = i1+1; i2 < n_city; i2++) {
			x          = city[i2][0] - city[i1][0];
			y          = city[i2][1] - city[i1][1];
			rg[i1][i2] = (int)(sqrt(x * x + y * y) + 0.5);
		}
	}

	for (i1 = 1; i1 < n_city; i1++) {
		for (i2 = 0; i2 < i1; i2++)
			rg[i1][i2] = rg[i2][i1];
	}

	fclose(in);
}

/****************/
/* デストラクタ */
/****************/
TSP::~TSP()
{
	int i1;

	for (i1 = 0; i1 < size+max_ch; i1++)
		delete [] ind[i1];
	delete [] ind;

	for (i1 = 0; i1 < max_len; i1++)
		delete [] edge[i1];
	delete [] edge;

	delete [] pi;
	delete [] len;
	delete [] kou1;
	delete [] kou2;
	delete [] pi_w;
	delete [] s_w;
	delete [] ro;

	for (i1 = 0; i1 < n_city; i1++) {
		delete [] city[i1];
		delete [] rg[i1];
	}
	delete [] city;
	delete [] rg;
}

/**************/
/* 全体の制御 */
/**************/
void TSP::Control()
{
	int gen = 1, k1;
					// 初期集団の発生
	Init_std();
					// 評価
	if (kinbo > 0)
		Kinbo();
	else
		Adap();
					// 出力
	printf("***世代 %d 適応度 max %f (%d) mean %f\n",
           gen, max, max_n, mean);

	if (abs(out_lvl) > 0)
		Output(gen);
					// 世代交代
	for (gen = 2; gen <= max_gen; gen++) {
						// 交叉
		switch (kosa_m) {
			case -1:
				break;
			case 0:
				C_copy();   // 親のコピー
				break;
			case 1:
				C_cycle(kosa);   // 循環交叉
				break;
			case 2:
				C_part(kosa);   // 部分的交叉
				break;
			case 3:
				C_seq(kosa);   // 順序交叉
				break;
			case 4:
				C_useq(kosa);   // 一様順序交叉
				break;
			case 5:
				C_upos(kosa);   // 一様位置交叉
				break;
			case 6:
				C_edge(kosa);   // エッジ組み替え交叉
				break;
			case 7:
				C_sub(kosa, k_point);   // サブツアー交叉
				break;
			default:
				break;
		}
						// 突然変異
		switch (mute_m) {
			case -1:
				break;
			case 0:
				M_move(mute);   // 移動
				break;
			case 1:
				M_inv(mute);   // 逆位
				break;
			case 2:
				M_scram(mute);   // スクランブル
				break;
			case 3:
				M_chg(mute);   // 転座
				break;
			default:
				break;
		}
						// 適応度
		if (kinbo > 0)
			Kinbo();
		else
			Adap();
						// 淘汰
		S_roul(elite);
						// 出力
		if (gen%out_d == 0)
			printf("***世代 %d 適応度 max %f (%d) mean %f\n",
                   gen, max, max_n, mean);

		if (abs(out_lvl) > 0) {
			if (gen%abs(out_lvl) == 0)
				Output(gen);
		}
	}

	gen--;
	k1    = out_m;
	out_m = 0;
	printf("***世代 %d 適応度 max %f (%d) mean %f\n",
            gen, max, max_n, mean);
	Output(gen);
	out_m = k1;
}

/*********************************/
/* 距離の計算                    */
/*      n_c : 都市の数           */
/*      p : 都市番号             */
/*      return : 距離（負）      */
/*********************************/
int TSP::Kyori(int n_c, int *p)
{
	int i1, n1, n2, range = 0;

	n1 = p[0];

	for (i1 = 1; i1 < n_c; i1++) {
		n2     = p[i1];
		range -= rg[n1][n2];
		n1     = n2;
	}

	n2     = p[0];
	range -= rg[n1][n2];

	return range;
}

/****************/
/* 適応度の計算 */
/****************/
void TSP::Adap()
{
	int i1, k = 0;

	mean  = 0.0;
	max   = 0.0;
	max_n = -1;

	for (i1 = 0; i1 < size+max_ch; i1++) {
		if (pi_w[i1] == 1) {
			pi_w[i1] = 2;
			pi[i1]   = Kyori(len[i1], ind[i1]);
		}
		if (pi_w[i1] > 0) {
			k++;
			mean += pi[i1];
			if (max_n < 0 || pi[i1] > max) {
				max   = pi[i1];
				max_n = i1;
			}
		}
	}

	if (k > 0)
		mean /= k;
}

/**************************************/
/* エッジの入れ替え                   */
/*      n_city : 都市の数             */
/*      seq : 訪問する順番            */
/*      r_m : 距離の負値              */
/*      return : =0 : 改善がなかった  */
/*               =1 : 改善があった    */
/**************************************/
int TSP::Change(int n_city, int *seq, int *r_m)
{
	int ch = 0, i1, i2, i3, i4, k, k1, k2, max, n1, n2, n3, nn, r, sw = 0;

	max = *r_m;

	n3  = (int)(genrand_real3() * (n_city - 2));
	if (n3 > n_city-3)
		n3 = n_city - 3;
                         // ２近傍
	for (i1 = 0; i1 <= n_city-3 && ch == 0; i1++) {

		if (n3 == 0)
			n1 = n_city - 2;
		else
			n1 = n_city - 1;

		for (i2 = n3+2; i2 <= n1 && ch == 0; i2++) {
                              // 枝の場所（(n3,n3+1), (k1,k2)）
			k1 = i2;
			if (i2 == n_city-1)
				k2 = 0;
			else
				k2 = i2 + 1;
                              // 枝の入れ替え
			kou1[0] = seq[n3];
			k       = 1;
			for (i3 = k1; i3 >= n3+1; i3--) {
				kou1[k] = seq[i3];
				k++;
			}

			nn = k2;
			while (nn != n3) {
				kou1[k] = seq[nn];
				k++;
				nn++;
				if (nn > n_city-1)
					nn = 0;
			}
                              // 評価
			r = Kyori(n_city, kou1);

			if (r > max) {
				max = r;
				sw  = 1;
				for (i3 = 0; i3 < n_city; i3++)
					kou2[i3] = kou1[i3];
				if (sel > 0)
					ch = 1;
			}
		}

		n3++;
		if (n3 > n_city-3)
			n3 = 0;
	}
                         // ３近傍
	if (neib == 3 && ch == 0) {

		for (i1 = 0; i1 <= n_city-3 && ch == 0; i1++) {

			n1 = n_city - 2;
			n2 = n_city - 1;

			for (i2 = n3+1; i2 <= n1 && ch == 0; i2++) {

				for (i3 = i2+1; i3 <= n2 && ch == 0; i3++) {
                              // 枝の場所（(n3,n3+1), (i2,i2+1), (k1,k2)）
					k1 = i3;
					if (i3 == n_city-1)
						k2 = 0;
					else
						k2 = i3 + 1;
                              // 枝の入れ替えと評価
                                   // 入れ替え（その１）
					kou1[0] = seq[n3];
					k       = 1;
					for (i4 = i2; i4 >= n3+1; i4--) {
						kou1[k] = seq[i4];
						k++;
					}

					for (i4 = k1; i4 >= i2+1; i4--) {
						kou1[k] = seq[i4];
						k++;
					}

					nn = k2;
					while (nn != n3) {
						kou1[k] = seq[nn];
						k++;
						nn++;
						if (nn > n_city-1)
							nn = 0;
					}
                                   // 評価（その１）
					r = Kyori(n_city, kou1);

					if (r > max) {
						max = r;
						sw  = 1;
						for (i3 = 0; i3 < n_city; i3++)
							kou2[i3] = kou1[i3];
						if (sel > 0)
							ch = 1;
					}
                                   // 入れ替え（その２）
					kou1[0] = seq[n3];
					k       = 1;
					for (i4 = k1; i4 >= i2+1; i4--) {
						kou1[k] = seq[i4];
						k++;
					}

					for (i4 = n3+1; i4 <= i2; i4++) {
						kou1[k] = seq[i4];
						k++;
					}

					nn = k2;
					while (nn != n3) {
						kou1[k] = seq[nn];
						k++;
						nn++;
						if (nn > n_city-1)
							nn = 0;
					}
                                   // 評価（その２）
					r = Kyori(n_city, kou1);

					if (r > max) {
						max = r;
						sw  = 1;
						for (i3 = 0; i3 < n_city; i3++)
							kou2[i3] = kou1[i3];
						if (sel > 0)
							ch = 1;
					}
                                   // 入れ替え（その３）
					kou1[0] = seq[n3];
					k       = 1;
					for (i4 = i2+1; i4 <= k1; i4++) {
						kou1[k] = seq[i4];
						k++;
					}

					for (i4 = i2; i4 >= n3+1; i4--) {
						kou1[k] = seq[i4];
						k++;
					}

					nn = k2;
					while (nn != n3) {
						kou1[k] = seq[nn];
						k++;
						nn++;
						if (nn > n_city-1)
							nn = 0;
					}
                                   // 評価（その３）
					r = Kyori(n_city, kou1);

					if (r > max) {
						max = r;
						sw  = 1;
						for (i3 = 0; i3 < n_city; i3++)
							kou2[i3] = kou1[i3];
						if (sel > 0)
							ch = 1;
					}
                                   // 入れ替え（その４）
					kou1[0] = seq[n3];
					k       = 1;
					for (i4 = i2+1; i4 <= k1; i4++) {
						kou1[k] = seq[i4];
						k++;
					}

					for (i4 = n3+1; i4 <= i2; i4++) {
						kou1[k] = seq[i4];
						k++;
					}

					nn = k2;
					while (nn != n3) {
						kou1[k] = seq[nn];
						k++;
						nn++;
						if (nn > n_city-1)
							nn = 0;
					}
                                   // 評価（その４）
					r = Kyori(n_city, kou1);

					if (r > max) {
						max = r;
						sw  = 1;
						for (i3 = 0; i3 < n_city; i3++)
							kou2[i3] = kou1[i3];
						if (sel > 0)
							ch = 1;
					}
				}
			}

			n3++;
			if (n3 > n_city-3)
				n3 = 0;
		}
	}
                         // 設定
	if (sw > 0) {
		*r_m = max;
		for (i1 = 0; i1 < n_city; i1++)
			seq[i1] = kou2[i1];
	}

	return sw;
}

/**************/
/* 近傍の探索 */
/**************/
void TSP::Kinbo()
{
	int i1, k = 0, r, sw;

	max   = 0.0;
	max_n = -1;
	mean  = 0.0;

	for (i1 = 0; i1 < size+max_ch; i1++) {
		if (pi_w[i1] == 1) {
			pi_w[i1] = 2;
			sw       = 1;
			r        = Kyori(len[i1], ind[i1]);
			while (sw > 0)
				sw = Change(len[i1], ind[i1], &r);
			pi[i1] = r;
		}
		if (pi_w[i1] > 0) {
			k++;
			mean += pi[i1];
			if (max_n < 0 || pi[i1] > max) {
				max   = pi[i1];
				max_n = i1;
			}
		}
	}

	if (k > 0)
		mean /= k;
}

/*****************************/
/* 結果の出力                */
/*      gen : 現在の世代番号 */
/*****************************/
void TSP::Output(int gen)
{
	int i1, k = 0, n, pr;
	char *now;
	time_t aclock;
	FILE *out;

	if (out_lvl >= 0) {
		printf("   出力先は（0:出力なし，n:画面にｎ個づつ，-1:ファイル）？ ");
		scanf("%d", &pr);
	}
	else
		pr = -1;

	if (pr != 0) {
					// 出力先の決定と評価値の出力
		if (pr > 0) {
			out = stdout;
			getchar();
		}
		else {
			time(&aclock);
			now = ctime(&aclock);
			out = fopen(o_file, "a");
			fprintf(out, "***世代 %d 適応度 max %f (%d) mean %f 時間 %s\n",
                    gen, max, max_n, mean, now);
		}
					// 巡回順序の出力
		if (out_m == 0) {
			for (i1 = 0; i1 < len[max_n]; i1++) {
				n = ind[max_n][i1];
				fprintf(out, "%d %d %d\n", n, city[n][0], city[n][1]);
				if (pr > 0) {
					k++;
					if (k == pr) {
						getchar();
						k = 0;
					}
				}
			}
		}

		if (pr < 0)
			fclose(out);
	}
}

/****************/
/* main program */
/****************/
int main(int argc, char *argv[])
{
	long *seed;
	int i1, n;
	char **i_file1, **i_file2;
	FILE *in;
	TSP *tsp;
					// 入力ミス
	if (argc <= 1) {
		printf("***error  ファイル名を入力して下さい\n");
		exit(1);
	}
					// 入力ＯＫ
	else {
						// ファイルのオープン
		in = fopen(argv[1], "r");
		if (in == NULL) {
			printf("***error  ファイル名が不適当です\n");
			exit(1);
		}
						// 入力データファイル名の入力
		fscanf(in, "%d", &n);   // データの数

		seed    = new long [n];
		i_file1 = new char * [n];
		i_file2 = new char * [n];

		for (i1 = 0; i1 < n; i1++) {
			i_file1[i1] = new char [100];
			i_file2[i1] = new char [100];
			seed[i1]    = 1000 * i1 + 1234567;
			fscanf(in, "%s %s", i_file1[i1], i_file2[i1]);
		}

		fclose(in);
						// 実行（乱数の初期値を変える）
		for (i1 = 0; i1 < n; i1++) {

			printf("\n+++++ケース %d+++++\n", i1+1);

			tsp = new TSP (i_file1[i1], i_file2[i1], seed[i1]);

			tsp->Control();
		}
	}

	return 0;
}

