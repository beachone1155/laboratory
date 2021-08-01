# -*- coding: UTF-8 -*-
import numpy as np
import sys
from math import *
from random import *
from datetime import *

#***********************************#
#* 遺伝的アルゴリズムによるTSPの解 *#
#*      coded by Y.Suganuma        *#
#***********************************#

######################
# クラスSpeciesの定義
######################

class Species :
	
	#########################
	# コンストラクタ
	#      name : ファイル名
	#########################

	def __init__(self, name) :
	
				# データの入力
		inn = open(name, "r")
	
		s = inn.readline().split()
		self.allele_u = int(s[1])   # 対立遺伝子上限
		self.allele_l = int(s[3])   # 対立遺伝子下限

		s = inn.readline().split()
		self.max_len = int(s[1])   # 最大遺伝子長
		self.min_len = int(s[3])   # 最小遺伝子長（負の時は，最大遺伝子長で固定）

		s = inn.readline().split()
		self.dup_a = int(s[1])  # 遺伝子の重複
		                        #   =0 : 重複を許さない
		                        #   =1 : 重複を許す
		self.dup_s = int(s[3])   # 個体の重複（同じ染色体の個体）
		                         #   =0 : 重複を許さない
		                         #   =1 : 重複を許す

		s = inn.readline().split()
		self.size   = int(s[1])   # 個体総数
		self.max_ch = int(s[3])   # 子供の数の最大値
				# データのチェック
		if self.size <= 0 :
			print("***error  個体総数≦０ (Constructor)")
	
		if self.max_ch < 0 :
			print("***error  子供の数＜０ (Constructor)")
	
		if self.max_len <= 0 or self.min_len == 0 :
			print("***error  遺伝子長≦０ (Constructor)")
	
		if self.max_len < self.min_len :
			print("***error  最大遺伝子長＜最小遺伝子長 (Constructor)")
	
		if self.allele_u <= self.allele_l :
			print("***error  対立遺伝子上限≦対立遺伝子下限 (Constructor)")
	
		kind = self.allele_u - self.allele_l + 1
		if self.dup_a == 0 and self.max_len > kind :
			print("***error  遺伝子の重複を防ぐことはできない (Constructor)")
				# 領域の確保
		num       = self.size + self.max_ch
		self.ind  = np.empty((num, self.max_len), np.int)   # 集団（個体の集まり）
		self.edge = np.empty((self.max_len, 5), np.int)   # エッジ組み替え交叉用ワークエリア
		self.pi   = np.empty(num, np.float)   # 適応度
		self.ro   = np.empty(num, np.float)   # ルーレット板
		self.len  = np.empty(num, np.int)   # 各個体の遺伝子長
		self.kou1 = np.empty(self.max_len, np.int)   # 交叉・突然変異用作業場所１
		self.kou2 = np.empty(self.max_len, np.int)   # 交叉・突然変異用作業場所２
		self.s_w  = np.empty(num, np.int)   # 淘汰用指標（選択された回数）
		self.pi_w = np.empty(num, np.int)   # 適応度計算指標
		                                    #   =0 : 未使用
		                                    #   =1 : 適応度計算前（突然変異はこの個体だけに適用）
		                                    #   =2 : 適応度計算済み（交叉時に親とみなす）
		self.max   = -inf   # 最大適応度
		self.mean  = 0.0   # 平均適応度
		self.max_n = -1   # 最大適応度の個体番号

	##################################################
	# 場所を探す
	#      n : >=0 : n番目の親を捜す
	#          -1 : 空いている場所を探す
	#      return : 親の場所，または，空いている場所
	#               （存在しないときは負の値）
	##################################################

	def Position(self, n) :
	
		k  = -1
		sw = 0
				# 空いている場所を探す
		if n < 0 :
			for i1 in range(0, self.size+self.max_ch) :
				if self.pi_w[i1] == 0 :
					k = i1
					break
			if k < 0 :
				print("***error  空いている場所がない --Position--")
				# ｎ番目の親（pi_w[i]=2）を捜す
		else :
			for i1 in range(0, self.size+self.max_ch) :
				if self.pi_w[i1] == 2 :
					k += 1
					if k == n :
						sw = 1
						k  = i1
						break
	
		return k
	
	###################################################################
	# 個体の選択
	#      method : 選択方法
	#                 =-1 : ランダム(default)
	#                 =0 : 適応度をそのまま使用
	#                 =1 : 最小値からの差（ただし，α以下の場合はα）
	#                 =2 : 評価値に順位をつけ，減少率βで線形化
	#      bias : α，または，method=2の場合は初期値(default=0)
	#      step : β(default=1)
	#      return : 個体番号
	###################################################################

	def Select(self, method, bias, step) :
	
		sum = 0.0
					# ルーレット板の用意
							# ランダム
		if method == -1 :
			n = 0
			for i1 in range(0, self.size+self.max_ch) :
				if self.pi_w[i1] > 1 :
					n += 1
			sum = 1.0 / n
			for i1 in range(0, self.size+self.max_ch) :
				if self.pi_w[i1] > 1 :
					self.ro[i1] = sum
							# 評価値をそのまま利用
		elif method == 0 :
			n = 0
			for i1 in range(0, self.size+self.max_ch) :
				if self.pi_w[i1] > 1 :
					sum += self.pi[i1]
					n   += 1
			if fabs(sum) > 1.0e-10 :
				sum = 1.0 / fabs(sum)
				for i1 in range(0, self.size+self.max_ch) :
					if self.pi_w[i1] > 1 :
						self.ro[i1] = self.pi[i1] * sum
			else :
				sum = 1.0 / n
				for i1 in range(0, self.size+self.max_ch) :
					if self.pi_w[i1] > 1 :
						self.ro[i1] = sum
							# 最小値からの差
		elif method == 1 :
			min = -1
			n   = 0
			for i1 in range(0, self.size+self.max_ch) :
				if self.pi_w[i1] > 1 :
					n += 1
					if min < 0 or self.pi[i1] < self.pi[min] :
						min = i1
			for i1 in range(0, self.size+self.max_ch) :
				if self.pi_w[i1] > 1 :
					self.ro[i1] = self.pi[i1] - self.pi[min]
					if self.ro[i1] < bias :
						self.ro[i1] = bias
					sum += self.ro[i1]
			if sum > 1.0e-10 :
				sum = 1.0 / sum
				for i1 in range(0, self.size+self.max_ch) :
					if self.pi_w[i1] > 1 :
						self.ro[i1] *= sum
			else :
				sum = 1.0 / n
				for i1 in range(0, self.size+self.max_ch) :
					if self.pi_w[i1] > 1 :
						self.ro[i1] = sum
							# 線形化
		elif method == 2 :
			n = 0
			for i1 in range(0, self.size+self.max_ch) :
				if self.pi_w[i1] > 1 :
					self.ro[i1] = -1.0
					n += 1
				else :
					self.ro[i1] = 1.0
			sw  = 0
			sum = bias
			while sw == 0 :
				min = -1
				for i1 in range(0, self.size+self.max_ch) :
					if self.ro[i1] < 0.0 and (min < 0 or self.pi[i1] < self.pi[min]) :
						min = i1
				if min < 0 :
					sw = 1
				else :
					self.ro[min]  = sum
					sum          += self.step
			sum = 1.0 / (0.5 * (2.0 * bias + step * (n - 1)) * n)
			for i1 in range(0, self.size+self.max_ch) :
				if self.pi_w[i1] > 1 :
					self.ro[i1] *= sum
	
		sum = 0.0
		for i1 in range(0, self.size+self.max_ch) :
			if self.pi_w[i1] > 1 :
				sum         += self.ro[i1]
				self.ro[i1]  = sum
						# 選択
		x  = random()
		sw = 0
		k  = 0
		for i1 in range(0, self.size+self.max_ch) :
			if self.pi_w[i1] > 1 :
				if x <= self.ro[i1] :
					sw = 1
					k  = i1
					break
	
		return k
	
	###################
	# 標準的な初期設定
	###################

	def Init_std(self) :
	
				# 初期設定
		for i1 in range(0, self.size+self.max_ch) :
			if i1 < self.size :
				self.pi_w[i1] = 1   # 適応度の計算前
			else :
				self.pi_w[i1] = 0   # 未使用
				# 遺伝子の決定
		for i1 in range(0, self.size) :
	
			sw1    = 0
			length = 0

			while sw1 == 0 :
						# 遺伝子長の決定
				if self.min_len < 0 :
					length = self.max_len
				else :
					length = int(random() * (self.max_len - self.min_len + 1) + self.min_len)
					if length > self.max_len :
						length = self.max_len
				self.len[i1] = length
						# 遺伝子の決定
				for i2 in range(0, length) :
					sw2 = 0
					while sw2 == 0 :
						lid = int(random() * (self.allele_u - self.allele_l + 1) + self.allele_l)
						if lid > self.allele_u :
							lid = self.allele_u
						self.ind[i1][i2] = lid
							# 重複遺伝子のチェック
						sw2 = 1
						if self.dup_a == 0 :
							for i3 in range(0, i2) :
								if lid == self.ind[i1][i3] :
									sw2 = 0
									break
						# 重複個体のチェック
				sw1 = 1
				if self.dup_s == 0 :
					for i2 in range(0, i1) :
						if self.len[i1] == self.len[i2] :
							sw2 = 0
							for i3 in range(0, self.len[i1]) :
								if self.ind[i1][i3] != self.ind[i2][i3] :
									sw2 = 1
									break
							if sw2 == 0 :
								sw1 = 0
								break
	
	####################################################
	# 標準的な出力
	#      sw : 出力レベル
	#             =0 : 最終出力だけ
	#             n>0 : ｎ世代毎に出力（負はファイル）
	#      out_m : 出力方法
	#                =0 : すべての個体を出力
	#                =1 : 最大適応度の個体だけを出力
	#      gen : 現在の世代番号
	#      name : 出力ファイル名
	####################################################

	def Out_std(self, sw, out_m, gen, name) :
	
		k  = 0
		pr = -1

		if sw >= 0 :
			pr = int(input("   出力先は（0:出力なし，n:画面にｎ個づつ，-1:ファイル）？ "))

		if pr != 0 :
						# 出力先の決定と評価値の出力
			if pr > 0 :
				out = sys.stdout
				input("")
			else :
				now = datetime.today().time().isoformat()
				out = open(name, "a")
				out.write("***世代 " + str(gen) + " 適応度 max " + str(self.max) + " (" + str(self.max_n) + ") mean " + str(self.mean) + " 時間 " + now + "\n")
						# 詳細出力
			for i1 in range(0, self.size+self.max_ch) :
				if (self.pi_w[i1] > 1) and (out_m == 0 or (out_m == 1 and i1 == self.max_n)) :
					out.write(str(i1) + " allele")
					for i2 in range(0, self.len[i1]) :
						out.write(" " + str(self.ind[i1][i2]))
					out.write(" value " + str(self.pi[i1]) + "\n")
					if pr > 0 :
						k += 1
						if k == pr :
							input("")
							k = 0
	
			if pr < 0 :
				out.close()
	
	###################################################################
	# 交叉（親のコピー）
	#      method : =2 : 有性（２つの親から２つの子供）(default)
	#               =1 : １つの親から１つの子供
	#      pair : method=2 の時は親のペア数(default=max_ch/2)
	#             method=1 の時は親の数（＝子供の数）
	#      k_method : 選択方法
	#                 =-1 : ランダム(default)
	#                 =0 : 適応度をそのまま使用
	#                 =1 : 最小値からの差（ただし，α以下の場合はα）
	#                 =2 : 評価値に順位をつけ，減少率βで線形化
	#      k_bias : α，または，method=2の場合は初期値(default=0)
	#      k_step : β(default=1)
	###################################################################

	def C_copy(self, method = 2, pair = 0, k_method = -1, k_bias = 0.0, k_step = 1.0) :
	
				# 初期設定とデータチェック
		if method != 1 :
			method = 2
	
		if pair <= 0 :
			if method == 2 :
				pair = int(self.max_ch / 2)
			else :
				pair = self.max_ch
		else :
			if method == 2 and 2*pair > self.max_ch or method == 1 and pair > self.max_ch :
				print("***error  子供が多すぎる (C_copy)")
				# 実行
		for i1 in range(0, pair) :
						# 親の選択
			p1 = self.Select(k_method, k_bias, k_step)
			p2 = p1
			sw = 0
	
			while sw == 0 :
				p2 = self.Select(k_method, k_bias, k_step)
				if p1 != p2 :
					sw = 1
						# コピー
			for i2 in range(0, method) :
				p = p2
				if i2 == 0 :
					p = p1
				k            = self.Position(-1)
				self.len[k]  = self.len[p]
				self.pi_w[k] = 1
				for i3 in range(0, self.len[k]) :
					self.ind[k][i3] = self.ind[p][i3]
	
	###################################################################
	# 交叉（多点交叉）
	#      kosa : 交叉確率
	#      k_point : 交叉点の数(default=1)
	#                （負の時は，1から-k_point間のランダム）
	#      k_vr : =0 : 両親とも同じ位置で交叉(default) 
	#             =1 : 両親が異なる位置で交叉（遺伝子長は可変）
	#      k_method : 選択方法
	#                 =-1 : ランダム(default) 
	#                 =0 : 適応度をそのまま使用
	#                 =1 : 最小値からの差（ただし，α以下の場合はα）
	#                 =2 : 評価値に順位をつけ，減少率βで線形化
	#      k_bias : α，または，method=2の場合は初期値(default=0)
	#      k_step : β(default=1)
	###################################################################

	def C_point(self, kosa, k_point = 1, k_vr = 0, k_method = -1, k_bias = 0.0, k_step = 1.0) :
	
		mn = 0
				# 初期設定とデータのチェック
		pair = int(self.max_ch / 2)
	
		if self.dup_a == 0 :
			print("***error  交叉方法が不適当 (C_point)")
	
		abs_p = abs(k_point)
		if abs_p == 0 or abs_p > self.max_len-1 or self.min_len > 0 and abs_p > self.min_len-1 :
			print("***error  交叉点の数が不適当 (C_point)")
	
		if k_vr > 0 and self.min_len < 0 :
			print("***error  遺伝子長は可変でなければならない (C_point)")
				# 交叉
		num = k_point
	
		for i1 in range(0, pair) :
						# 交叉しない場合
			if random() > kosa :
				self.C_copy(2, 1)
						# 交叉する場合
			else :
							# 親の選択
				p1 = self.Select(k_method, k_bias, k_step)
				p2 = p1
				sw = 0
				while sw == 0 :
					p2 = self.Select(k_method, k_bias, k_step)
					if p1 != p2 :
						sw = 1
							# 交叉位置の数の決定
				if k_point < 0 :
					num = int(random() * abs_p + 1)
					if num > abs_p :
						num = abs_p
							# 交叉位置の決定（点の後ろで交叉）
				for i2 in range(0, num) :
								# 親１の交叉位置
					sw = 0
					while sw == 0 :
						sw            = 1
						self.kou1[i2] = int(random() * (self.len[p1] - 1))
						if self.kou1[i2] > self.len[p1]-2 :
							self.kou1[i2] = self.len[p1] - 2
						if k_vr == 0 and self.kou1[i2] > self.len[p2]-2 :
							self.kou1[i2] = self.len[p2] - 2
						for i3 in range(0, i2) :
							if self.kou1[i3] == self.kou1[i2] :
								sw = 0
								break
								# 親２の交叉位置
					if k_vr > 0 :
						sw = 0
						while sw == 0 :
							sw            = 1
							self.kou2[i2] = int(random() * (self.len[p2] - 1))
							if self.kou2[i2] > self.len[p2]-2 :
								self.kou2[i2] = self.len[p2] - 2
							for i3 in range(0, i2) :
								if self.kou2[i3] == self.kou2[i2] :
									sw = 0
									break
							# 交叉の実行
							#   親１のt11からt12を子１のc1へコピー
							#   親２のt21からt22を子２のc2へコピー
							#     次は，
							#   親１のt11からt12を子２のc2へコピー
							#   親２のt21からt22を子１のc1へコピー
							#     ・・・・・
				c1  = 0
				c2  = 0
				t11 = 0
				t21 = 0
								# 遺伝子長
				k1            = self.Position(-1)
				self.pi_w[k1] = 1
				self.len[k1]  = self.len[p1]
				k2            = self.Position(-1)
				self.pi_w[k2] = 1
				self.len[k2]  = self.len[p2]
	
				for i2 in range(0, num+1) :
								# 次の交叉位置を求める
					if i2 == num :            # 最後
						t12 = self.len[p1]
						t22 = self.len[p2]
					else :
									# 親１
						t12 = self.max_len
						for i3 in range(0, num) :
							if self.kou1[i3] >= 0 and self.kou1[i3] <= t12 :
								t12 = self.kou1[i3]
								mn  = i3
						self.kou1[mn] = -1
						t12 += 1
									# 親２
						if k_vr == 0 :
							t22 = t12
						else :
							t22 = self.max_len
							for i3 in range(0, num) :
								if self.kou2[i3] >= 0 and self.kou2[i3] <= t22 :
									t22 = self.kou2[i3]
									mn  = i3
							self.kou2[mn] = -1
							t22 += 1
								# 指定箇所のコピー
					for i3 in range(t11, t12) :
						if i2%2 == 0 :
							if c1 < self.max_len :
								self.ind[k1][c1] = self.ind[p1][i3]
								c1 += 1
						else :
							if c2 < self.max_len :
								self.ind[k2][c2] = self.ind[p1][i3]
								c2 += 1
	
					for i3 in range(t21, t22) :
						if i2%2 == 0 :
							if c2 < self.max_len :
								self.ind[k2][c2] = self.ind[p2][i3]
								c2 += 1
						else :
							if c1 < self.max_len :
								self.ind[k1][c1] = self.ind[p2][i3]
								c1 += 1
								# 交叉位置の移動
					t11 = t12
					t21 = t22
	
	###################################################################
	# 交叉（一様交叉．[0,1]を等確率で発生させ，1であれば，
	#       親１，0であれば親２の遺伝子を子１が受け継ぐ）
	#      kosa : 交叉確率
	#      k_method : 選択方法
	#                 =-1 : ランダム(default)
	#                 =0 : 適応度をそのまま使用
	#                 =1 : 最小値からの差（ただし，α以下の場合はα）
	#                 =2 : 評価値に順位をつけ，減少率βで線形化
	#      k_bias : α，または，method=2の場合は初期値(default=0)
	#      k_step : β(default=1)
	###################################################################

	def C_uniform(self, kosa, k_method = -1, k_bias = 0.0, k_step = 1.0) :
	
				# 初期設定とデータのチェック
		pair = int(self.max_ch / 2)
	
		if self.dup_a == 0 :
			print("***error  交叉方法が不適当 (C_uniform)")
	
		if self.min_len > 0 :
			print("***error  遺伝子長は固定長でなければならない (C_uniform)")
				# 交叉
		for i1 in range(0, pair) :
						# 交叉しない場合
			if random() > kosa :
				self.C_copy(2, 1)
						# 交叉する場合
			else :
							# 親の選択
				p1 = self.Select(k_method, k_bias, k_step)
				p2 = p1
				sw = 0
				while sw == 0 :
					p2 = self.Select(k_method, k_bias, k_step)
					if p1 != p2 :
						sw = 1
							# 遺伝子長
				k1            = self.Position(-1)
				self.pi_w[k1] = 1
				self.len[k1]  = self.len[p1]
				k2            = self.Position(-1)
				self.pi_w[k2] = 1
				self.len[k2]  = self.len[p2]
							# 交叉
				for i2 in range(0, self.len[p1]) :
					if random() > 0.5 :
						self.ind[k1][i2] = self.ind[p1][i2]
						self.ind[k2][i2] = self.ind[p2][i2]
					else :
						self.ind[k1][i2] = self.ind[p2][i2]
						self.ind[k2][i2] = self.ind[p1][i2]
	
	###################################################################
	# 交叉（平均化交叉．２つの親の平均値を受け継ぐ）
	#      kosa : 交叉確率
	#      k_method : 選択方法
	#                 =-1 : ランダム(default)
	#                 =0 : 適応度をそのまま使用
	#                 =1 : 最小値からの差（ただし，α以下の場合はα）
	#                 =2 : 評価値に順位をつけ，減少率βで線形化
	#      k_bias : α，または，method=2の場合は初期値(default=0)
	#      k_step : β(default=1)
	###################################################################

	def C_mean(self, kosa, k_method = -1, k_bias = 0.0, k_step = 1.0) :
	
				# 初期設定とデータのチェック
		if self.min_len > 0 :
			print("***error  遺伝子長は固定長でなければならない (C_mean)")
				# 交叉
		for i1 in range(0, self.max_ch) :
						# 交叉しない場合
			if random() > kosa :
				self.C_copy(1, 1)
						# 交叉する場合
			else :
							# 親の選択
				p1 = self.Select(k_method, k_bias, k_step)
				p2 = p1
				sw = 0
				while sw == 0 :
					p2 = self.Select(k_method, k_bias, k_step)
					if p1 != p2 :
						sw = 1
							# 遺伝子長
				k            = self.Position(-1)
				self.len[k]  = self.len[p1]
				self.pi_w[k] = 1
							# 交叉
				for i2 in range(0, self.len[k]) :
					self.ind[k][i2] = int((self.ind[p1][i2] + self.ind[p2][i2]) / 2)
	
	###################################################################
	# 交叉（循環交叉．ランダムに１点を選択し，その位置にある遺伝子を
	#       そのまま各子供が選択する．その位置にある親２（１）の遺伝
	#       子を，その遺伝子の親１（２）の場所に，子１（２）が受け継
	#       ぐ（ただし，doubleの場合は，この手続きをのぞく）．この手
	#       続きを，すでに受け継いだ遺伝子の位置が選択されるまで繰り
	#       返し，残りの遺伝子については，子１（２）は，親２（１）の
	#       遺伝子をその順番通りに受け継ぐ）
	#         2 4 1 3 6 5    + + 1 + + 5    3 2 1 4 6 5
	#             *       →             →
	#         3 2 5 4 1 6    + + 5 + 1 +    2 4 5 3 1 6
	#      kosa : 交叉確率
	#      k_method : 選択方法
	#                 =-1 : ランダム(default)
	#                 =0 : 適応度をそのまま使用
	#                 =1 : 最小値からの差（ただし，α以下の場合はα）
	#                 =2 : 評価値に順位をつけ，減少率βで線形化
	#      k_bias : α，または，method=2の場合は初期値(default=0)
	#      k_step : β(default=1)
	###################################################################
	def C_cycle(self, kosa, k_method = -1, k_bias = 0.0, k_step = 1.0) :
	
				# 初期設定とデータのチェック
		pair = int(self.max_ch / 2)
	
		if self.dup_a != 0 :
			print("***error  交叉方法が不適当 (C_cycle)")
	
		if self.min_len > 0 :
			print("***error  遺伝子長は固定長でなければならない (C_cycle)")
				# 交叉
		for i1 in range(0, pair) :
						# 交叉しない場合
			if random() > kosa :
				self.C_copy(2, 1)
						# 交叉する場合
			else :
							# 親の選択
				p1 = self.Select(k_method, k_bias, k_step)
				p2 = p1
				sw = 0
				while sw == 0 :
					p2 = self.Select(k_method, k_bias, k_step)
					if p1 != p2 :
						sw = 1
							# 初期設定
				for i2 in range(0, self.len[p1]) :
					self.kou1[i2] = 0
					self.kou2[i2] = 0
							# 遺伝子長
				k1            = self.Position(-1)
				self.pi_w[k1] = 1
				self.len[k1]  = self.len[p1]
				k2            = self.Position(-1)
				self.pi_w[k2] = 1
				self.len[k2]  = self.len[p2]
							# 交叉
				sw = 0
	
				while sw == 0 :
					sw = 1
					p  = int(random() * self.len[p1])
					if p >= self.len[p1] :
						p = self.len[p1] - 1
					if self.kou1[p] == 0 and self.kou2[p] == 0 :
						self.kou1[p]    = 1
						self.kou2[p]    = 1
						self.ind[k1][p] = self.ind[p1][p]
						self.ind[k2][p] = self.ind[p2][p]
						for i2 in range(0, self.len[p1]) :
							if self.ind[p2][p] == self.ind[p1][i2] :
								self.ind[k1][i2] = self.ind[p1][i2]
								self.kou1[i2]    = 1
								sw               = 0
								break
						sw = 1
						for i2 in range(0, self.len[p2]) :
							if self.ind[p1][p] == self.ind[p2][i2] :
								self.ind[k2][i2] = self.ind[p2][i2]
								self.kou2[i2]    = 1
								sw               = 0
								break
	
				sw = 0
				i2 = 0
				i3 = 0
				while sw == 0 :
					while sw == 0 and i2 < self.len[p1] :
						if self.kou1[i2] == 0 :
							sw = 1
						else :
							i2 += 1
					sw = 0
					while sw == 0 and i3 < self.len[p2] :
						if self.kou2[i3] == 0 :
							sw = 1
						else :
							i3 += 1
					if i2 < self.len[p1] and i3 < self.len[p2] :
						self.ind[k1][i2] = self.ind[p2][i3]
						self.ind[k2][i3] = self.ind[p1][i2]
						sw          = 0
						i2 += 1
						i3 += 1
					else :
						sw = 1
	
	###################################################################
	# 交叉（部分的交叉．ランダムに１点を選択し，その位置にある親１と
	#       親２の遺伝子を取り出す．次に，親１と親２の染色体上で，こ
	#       の２つの遺伝子の位置を交換する．この操作を，選択した点よ
	#       り右にあるすべての遺伝子に対して実施する
	#         2 4 1 3 6 5    2 4 5 3 6 1
	#             *       →             → ･････
	#         3 2 5 4 1 6    3 2 1 4 5 6
	#      kosa : 交叉確率
	#      k_method : 選択方法
	#                 =-1 : ランダム(default)
	#                 =0 : 適応度をそのまま使用
	#                 =1 : 最小値からの差（ただし，α以下の場合はα）
	#                 =2 : 評価値に順位をつけ，減少率βで線形化
	#      k_bias : α，または，method=2の場合は初期値(default=0)
	#      k_step : β(default=1)
	#******************************************************************/
	def C_part(self, kosa, k_method = -1, k_bias = 0.0, k_step = 1.0) :
	
				# 初期設定とデータのチェック
		pair = int(self.max_ch / 2)
	
		if self.dup_a != 0 :
			print("***error  交叉方法が不適当 (C_part)")
	
		if self.min_len > 0 :
			print("***error  遺伝子長は固定長でなければならない (C_part)")
				# 交叉
		for i1 in range(0, pair) :
						# 交叉しない場合
			if random() > kosa :
				self.C_copy(2, 1)
						# 交叉する場合
			else :
							# 親の選択
				p1 = self.Select(k_method, k_bias, k_step)
				p2 = p1
				sw = 0
				while sw == 0 :
					p2 = self.Select(k_method, k_bias, k_step)
					if p1 != p2 :
						sw = 1
							# 遺伝子長
				k1            = self.Position(-1)
				self.pi_w[k1] = 1
				self.len[k1]  = self.len[p1]
				k2            = self.Position(-1)
				self.pi_w[k2] = 1
				self.len[k2]  = self.len[p2]
							# 交叉
				p = int(random() * self.len[p1])
				if p >= self.len[p1] :
					p = self.len[p1] - 1
	
				for i2 in range(0, self.len[p1]) :
					self.ind[k1][i2] = self.ind[p1][i2]
					self.ind[k2][i2] = self.ind[p2][i2]
	
				for i2 in range(p, self.len[p1]) :
					sw = 0
					lv = self.ind[k1][i2]
					for i3 in range(0, self.len[p1]) :
						if self.ind[k2][i2] == self.ind[k1][i3] :
							self.ind[k1][i2] = self.ind[k1][i3]
							self.ind[k1][i3] = lv
							sw               = 1
							break
					sw = 0
					for i3 in range(0, self.len[p1]) :
						if lv == self.ind[k2][i3] :
							self.ind[k2][i3] = self.ind[k2][i2]
							self.ind[k2][i2] = lv
							sw               = 1
							break
	
	###################################################################
	# 交叉（順序交叉．ランダムに切れ目を決定し，子１に対し，切れ目の
	#       左側では，親１の遺伝子をそのまま受け継ぎ，右側では，親１
	#       の遺伝子を親２の遺伝子の出現順序に並べ替える．
	#         2 4 1 3 6 5    2 4 1 3 5 6
	#             *       →
	#         3 2 5 4 1 6    3 2 5 4 1 6
	#      kosa : 交叉確率
	#      k_method : 選択方法
	#                 =-1 : ランダム(default)
	#                 =0 : 適応度をそのまま使用
	#                 =1 : 最小値からの差（ただし，α以下の場合はα）
	#                 =2 : 評価値に順位をつけ，減少率βで線形化
	#      k_bias : α，または，method=2の場合は初期値(default=0)
	#      k_step : β(default=1)
	###################################################################
	def C_seq(self, kosa, k_method = -1, k_bias = 0.0, k_step = 1.0) :
	
				# 初期設定とデータのチェック
		pair = int(self.max_ch / 2)
	
		if self.dup_a != 0 :
			print("***error  交叉方法が不適当 (C_seq)")
	
		if self.min_len > 0 :
			print("***error  遺伝子長は固定長でなければならない (C_seq)")
				# 交叉
		for i1 in range(0, pair) :
						# 交叉しない場合
			if random() > kosa :
				self.C_copy(2, 1)
						# 交叉する場合
			else :
							# 親の選択
				p1 = self.Select(k_method, k_bias, k_step)
				p2 = p1
				sw = 0
				while sw == 0 :
					p2 = self.Select(k_method, k_bias, k_step)
					if p1 != p2 :
						sw = 1
							# 遺伝子長
				k1            = self.Position(-1)
				self.pi_w[k1] = 1
				self.len[k1]  = self.len[p1]
				k2            = self.Position(-1)
				self.pi_w[k2] = 1
				self.len[k2]  = self.len[p2]
							# 交叉
				p = int(random() * (self.len[p1] - 1))
				if p >= self.len[p1]-1 :
					p = self.len[p1] - 2
	
				for i2 in range(0, p+1) :
					self.ind[k1][i2] = self.ind[p1][i2]
					self.ind[k2][i2] = self.ind[p2][i2]
	
				pp = 0
				for i2 in range(p+1, self.len[p1]) :
					sw = 0
					i3 = pp
					while i3 < self.len[p2] and sw == 0 :
						i4 = p + 1
						while i4 < self.len[p1] and sw == 0 :
							if self.ind[p2][i3] == self.ind[p1][i4] :
								sw               = 1
								pp               = i3 + 1
								self.ind[k1][i2] = self.ind[p1][i4]
							i4 += 1
						i3 += 1
				pp = 0
				for i2 in range(p+1, self.len[p2]) :
					sw = 0
					i3 = pp
					while i3 < self.len[p1] and sw == 0 :
						i4 = p + 1
						while i4 < self.len[p2] and sw == 0 :
							if self.ind[p1][i3] == self.ind[p2][i4] :
								sw               = 1
								pp               = i3 + 1
								self.ind[k2][i2] = self.ind[p2][i4]
							i4 += 1
						i3 += 1
	
	###################################################################
	# 交叉（一様順序交叉．位置の集合をランダムに選択し，一方の親の選
	#       択された位置における遺伝子の順序に従って，他の親の遺伝子
	#       を並べ替える
	#         2 4 1 3 6 5    2 4 1 3 6 5
	#           *   *     →
	#         3 2 5 4 1 6    4 2 5 3 1 6
	#      kosa : 交叉確率
	#      k_method : 選択方法
	#                 =-1 : ランダム(default)
	#                 =0 : 適応度をそのまま使用
	#                 =1 : 最小値からの差（ただし，α以下の場合はα）
	#                 =2 : 評価値に順位をつけ，減少率βで線形化
	#      k_bias : α，または，method=2の場合は初期値(default=0)
	#      k_step : β(default=1)
	###################################################################
	
	def C_useq(self, kosa, k_method = -1, k_bias = 0.0, k_step = 1.0) :
	
				# 初期設定とデータのチェック
		pair =int(self.max_ch / 2)
	
		if self.dup_a != 0 :
			print("***error  交叉方法が不適当 (C_useq)")
	
		if self.min_len > 0 :
			print("***error  遺伝子長は固定長でなければならない (C_useq)")
				# 交叉
		for i1 in range(0, pair) :
						# 交叉しない場合
			if random() > kosa :
				self.C_copy(2, 1)
						# 交叉する場合
			else :
							# 親の選択
				p1 = self.Select(k_method, k_bias, k_step)
				p2 = p1
				sw = 0
				while sw == 0 :
					p2 = self.Select(k_method, k_bias, k_step)
					if p1 != p2 :
						sw = 1
							# 遺伝子長
				k1            = self.Position(-1)
				self.pi_w[k1] = 1
				self.len[k1]  = self.len[p1]
				k2            = self.Position(-1)
				self.pi_w[k2] = 1
				self.len[k2]  = self.len[p2]
							# 交叉
				for i2 in range(0, self.len[p1]) :
					self.ind[k1][i2] = self.ind[p1][i2]
					self.ind[k2][i2] = self.ind[p2][i2]
					if random() < 0.5 :
						self.kou1[i2] = 0
					else :
						self.kou1[i2] = 1
	
				p = 0
				for i2 in range(0, self.len[p1]) :
					if self.kou1[i2] > 0 :
						sw = 0
						i3 = p
						while i3 < self.len[p2] and sw == 0 :
							i4 = 0
							while i4 < self.len[p1] and sw == 0 :
								if self.ind[p2][i3] == self.ind[p1][i4] and self.kou1[i4] > 0 :
									sw               = 1
									p                = i3 + 1
									self.ind[k1][i2] = self.ind[p1][i4]
								i4 += 1
							i3 += 1
				p = 0
				for i2 in range(0, self.len[p3]) :
					if self.kou1[i2] > 0 :
						sw = 0
						i3 = p
						while i3 < self.len[p1] and sw == 0 :
							i4 = 0
							while i4 < self.len[p2] and sw == 0 :
								if self.ind[p1][i3] == self.ind[p2][i4] and self.kou1[i4] > 0 :
									sw               = 1
									p                = i3 + 1
									self.ind[k2][i2] = self.ind[p2][i4]
								i4 += 1
							i3 += 1
	
	###################################################################
	# 交叉（一様位置交叉．位置の集合をランダムに選択し，一方の親の選
	#       択された位置における遺伝子の位置に，他の親の同じ遺伝子を
	#       配置する．残りの遺伝子は，親と同じ順序に配置する．
	#         2 4 1 3 6 5    + + 5 + 1 +    2 4 5 3 1 6
	#             *   *   →             →
	#         3 2 5 4 1 6    + + 1 + 6 +    3 2 1 5 6 4
	#      kosa : 交叉確率
	#      k_method : 選択方法
	#                 =-1 : ランダム(default)
	#                 =0 : 適応度をそのまま使用
	#                 =1 : 最小値からの差（ただし，α以下の場合はα）
	#                 =2 : 評価値に順位をつけ，減少率βで線形化
	#      k_bias : α，または，method=2の場合は初期値(default=0)
	#      k_step : β(default=1)
	###################################################################

	def C_upos(self, kosa, k_method = -1, k_bias = 0.0, k_step = 1.0) :
	
				# 初期設定とデータのチェック
		pair = int(self.max_ch / 2)
	
		if self.dup_a != 0 :
			print("***error  交叉方法が不適当 (C_upos)")
	
		if self.min_len > 0 :
			print("***error  遺伝子長は固定長でなければならない (C_upos)")
				# 交叉
		for i1 in range(0, pair) :
						# 交叉しない場合
			if random() > kosa :
				self.C_copy(2, 1)
						# 交叉する場合
			else :
							# 親の選択
				p1 = self.Select(k_method, k_bias, k_step)
				p2 = p1
				sw = 0
				while sw == 0 :
					p2 = self.Select(k_method, k_bias, k_step)
					if p1 != p2 :
						sw = 1
							# 遺伝子長
				k1            = self.Position(-1)
				self.pi_w[k1] = 1
				self.len[k1]  = self.len[p1]
				k2            = self.Position(-1)
				self.pi_w[k2] = 1
				self.len[k2]  = self.len[p2]
							# 交叉
				for i2 in range(0, self.len[p1]) :
					self.kou1[i2] = 1
					if random() < 0.5 :
						self.kou1[i2] = 1
					if self.kou1[i2] > 0 :
						self.ind[k1][i2] = self.ind[p2][i2]
						self.ind[k2][i2] = self.ind[p1][i2]
	
				p = 0
				for i2 in range(0, self.len[p1]) :
					sw = 0
					for i3 in range(0, self.len[p1]) :
						if self.kou1[i3] > 0 and self.ind[p1][i2] == self.ind[k1][i3] :
							sw = 1
							break
					if sw == 0 :
						for i3 in range(p, self.len[p1]) :
							if self.kou1[i3] == 0 :
								self.ind[k1][i3] = self.ind[p1][i2]
								p                = i3 + 1
								sw               = 1
								break
				p = 0
				for i2 in range(0, self.len[p2]) :
					sw = 0
					for i3 in range(0, self.len[p2]) :
						if self.kou1[i3] > 0 and self.ind[p2][i2] == self.ind[k2][i3] :
							sw = 1
							break
					if sw == 0 :
						for i3 in range(p, self.len[p2]) :
							if self.kou1[i3] == 0 :
								self.ind[k2][i3] = self.ind[p2][i2]
								p                = i3 + 1
								sw               = 1
								break
	
	###################################################################
	# 交叉（エッジ組み替え交叉．以下の手順に従って行う．対立遺伝子は
	#       0～(max_len-1)である必要がある）
	#         (0) エッジマップを作成する．エッジマップとは，２つの親
	#             を見て，ノードがどこに接続されているのかを表すもの
	#             であり，例えば，２つの親が，
	#                 [A B C D E F]
	#                 [B D C A E F]
	#             である場合は，
	#                 A : B F C E
	#                 B : A C D F
	#                 C : B D A
	#                 D : C E B
	#                 E : D F A
	#                 F : A E B
	#             となる． 
	#         (1) 両親の２つの出発点の内１つで初期化する．ランダムま
	#             たはステップ(4)の基準に従って選ぶ（現在のノード）
	#         (2) エッジマップから，現在のノードを除く
	#         (3) 現在のノードが接続先のノードを持っていたら，(4)に
	#             進む．さもなければ，(5)に進む
	#         (4) 現在のノードが持っている接続先ノードの内，最も少な
	#             い接続先ノードを持ったノードを選択し（同じ条件の場
	#             合は，ランダム），それを現在のノードとし，(2)に進む
	#         (5) 未接続のノードが残っていればランダムに選択し，(2)に
	#             戻る．さもなければ，終了する
	#      kosa : 交叉確率
	#      k_method : 選択方法
	#                 =-1 : ランダム(default)
	#                 =0 : 適応度をそのまま使用
	#                 =1 : 最小値からの差（ただし，α以下の場合はα）
	#                 =2 : 評価値に順位をつけ，減少率βで線形化
	#      k_bias : α，または，method=2の場合は初期値(default=0)
	#      k_step : β(default=1)
	###################################################################

	def C_edge(self, kosa, k_method = -1, k_bias = 0.0, k_step = 1.0) :
	
		e  = np.empty(2, np.int)
		k0 = 0
				# 初期設定とデータのチェック
		pair = self.max_ch
	
		if self.dup_a != 0 :
			print("***error  交叉方法が不適当 (C_edge)")
	
		if self.min_len > 0 :
			print("***error  遺伝子長は固定長でなければならない (C_edge)")
				# 交叉
		for i1 in range(0, pair) :
						# 交叉しない場合
			if random() > kosa :
				self.C_copy(1, 1)
						# 交叉する場合
			else :
							# 親の選択
				p1 = self.Select(k_method, k_bias, k_step)
				p2 = p1
				sw = 0
				while sw == 0 :
					p2 = self.Select(k_method, k_bias, k_step)
					if p1 != p2 :
						sw = 1
							# 遺伝子長
				k            = self.Position(-1)
				self.pi_w[k] = 1
				self.len[k]  = self.len[p1]
							# エッジマップの初期化
				for i2 in range(0, self.len[k]) :
					self.edge[i2][0] = 0
					for i3 in range(1, 5) :
						self.edge[i2][i3] = -1
							# 交叉
								# エッジマップの作成
				for i2 in range(0, self.len[k]) :
	
					sw = 0
					for i3 in range(0, self.len[k]) :
						if i2 == self.ind[p1][i3] :
							sw = 1
							if i3 == 0 :
								e[0] = self.ind[p1][self.len[k]-1]
								e[1] = self.ind[p1][1]
							else :
								if i3 == self.len[k]-1 :
									e[0] = self.ind[p1][i3-1]
									e[1] = self.ind[p1][0]
								else :
									e[0] = self.ind[p1][i3-1]
									e[1] = self.ind[p1][i3+1]
							for i4 in range(0, 2) :
								self.edge[i2][0] += 1
								self.edge[i2][self.edge[i2][0]] = e[i4]
							break
	
					sw = 0
					for i3 in range(0, self.len[k]) :
						if i2 == self.ind[p2][i3] :
							sw = 1
							if i3 == 0 :
								e[0] = self.ind[p2][self.len[k]-1]
								e[1] = self.ind[p2][1]
							else :
								if i3 == self.len[k]-1 :
									e[0] = self.ind[p2][i3-1]
									e[1] = self.ind[p2][0]
								else :
									e[0] = self.ind[p2][i3-1]
									e[1] = self.ind[p2][i3+1]
							for i4 in range(0, 2) :
								sw = 1
								for i5 in range(1, self.edge[i2][0]+1) :
									if self.edge[i2][i5] == e[i4] :
										sw = 2
										break
								if sw == 1 :
									self.edge[i2][0] += 1
									self.edge[i2][self.edge[i2][0]] = e[i4]
							break
								# 交叉の実行
									# 出発点の決定
				k1 = self.ind[p1][0]
				k2 = self.ind[p2][0]
				if self.edge[k1][0] == self.edge[k2][0] :
					kk = k1
					if random() > 0.5 :
						kk = k2
				else :
					kk = k2
					if self.edge[k1][0] < self.edge[k2][0] :
						kk = k2
				self.ind[k][0] = kk
				p              = 1
	
				while p < self.len[k] :
									# ノードの除去
					for i2 in range(0, self.len[k]) :
						sw = 0
						if self.edge[i2][0] > 0 :
							for i3 in range(1, 5) :
								if self.edge[i2][i3] == kk :
									sw                 = 1
									self.edge[i2][i3]  = -1
									self.edge[i2][0]  -= 1
									break
									# 次の現在ノードの選択
					min = 10
					num = 0
					for i2 in range(1, 5) :
						if self.edge[kk][i2] >= 0 :
							k1 = self.edge[kk][i2]
							if self.edge[k1][0] >= 0 and self.edge[k1][0] < min :
								num = 1
								min = self.edge[k1][0]
								k0  = k1
							else :
								if self.edge[k1][0] == min :
									num += 1
					if num > 1 :
						k1 = int(random() * num) + 1
						if k1 > num :
							k1 = num
						k2 = 0
						k0 = -1
						i2 = 1
						while i2 <= 4 and k0 < 0 :
							if self.edge[kk][i2] >= 0 :
								if self.edge[self.edge[kk][i2]][0] == min :
									k2 += 1
									if k1 == k2 :
										k0 = self.edge[kk][i2]
							i2 += 1
					else :
						if num <= 0 :
							num = 0
							for i2 in range(0, self.len[k]) :
								if i2 != kk and self.edge[i2][0] >= 0 :
									num += 1
							if num <= 0 :
								print("***error  invalid data (C_edge)")
							else :
								k1 = int(random() * num) + 1
								if k1 > num :
									k1 = num
								k2 = 0
								k0 = -1
								i2 = 0
								while i2 < self.len[k] and k0 < 0 :
									if i2 != kk and self.edge[i2][0] >= 0 :
										k2 += 1
										if k1 == k2 :
											k0 = i2
									i2 += 1
					self.edge[kk][0]  = -1
					self.ind[k][p]    = k0
					kk                = k0
					p                += 1
	
	##############################################################
	# 交叉（サブツアー交叉．２点交叉の拡張である．ただし，相手に
	#       同じ遺伝子のグループがない限り実行されない．たとえば
	#         ***abcd**
	#         *cdab****
	#       のような両親の時実行され，以下の４つの子供が生成され
	#       る）
	#         ***cdab**
	#         *abcd****
	#         ***badc**
	#         *dcba****
	#       最大，４＊交叉回数＊個体総数＊(個体総数－１) 個の子
	#       供が生成される可能性があるので，子供の数としてこの値
	#       以上のデータを入力しておく必要がある．
	#      kosa : 交叉確率
	#      count : １つのペアーに対する交差回数(default=10)
	##############################################################

	def C_sub(self, kosa, count = 10) :
	
		t22 = 0
				# 初期設定とデータのチェック
		if (4*count*self.size*(self.size-1)) > self.max_ch :
			print("***error  子供が多すぎる (C_sub)")
				# 交叉
		for i1 in range(0, self.size-1) :
						# 親１
			p1 = self.Position(i1)
	
			if p1 >= 0 :
	
				for i2 in range(i1, self.size) :
						# 親２
					p2 = self.Position(i2)
	
					if p2 >= 0 :
						# 交叉しない場合
						if random() > kosa :
							self.C_copy(2, 1)
						# 交叉する場合
						else :
							# 交叉回数の制御
							for i3 in range(0, count) :
								# 交叉位置の決定（点の後ろで交叉）
									# 親１の交叉位置
								t11 = int(random() * self.len[p1])
								if t11 > (self.len[p1]-1) :
									t11 = self.len[p1] - 1
								sw = 0
								while sw == 0 :
									t12 = int(random() * self.len[p1])
									if t12 > (self.len[p1]-1) :
										t12 = self.len[p1] - 1
									if t12 != t11 :
										sw = 1
								if t11 > t12 :
									k1  = t11
									t11 = t12
									t12 = k1
									# 親２の交叉位置
								sw  = 0
								t21 = -1
								i4  = 0
								while i4 < self.len[p2] and t21 < 0 :
									i5 = t11
									while i5 <= t12 and t21 < 0 :
										if self.ind[p2][i4] == self.ind[p1][i5] :
											t21 = i4
										i5 += 1
									i4 += 1
								if t21 >= 0 :
									t22 = t21 + t12 - t11
									if t22 < self.len[p2] :
										sw = 1
										i4 = t21 + 1
										while i4 <= t22 and sw > 0 :
											sw = 0
											i5 = t11
											while i5 <= t12 and sw == 0 :
												if self.ind[p2][i4] == self.ind[p1][i5] :
													sw = 1
												i5 += 1
											i4 += 1
									# 交叉の実行
								if sw > 0 :
	
									k1            = self.Position(-1)
									self.pi_w[k1] = 1
									self.len[k1]  = self.len[p1]
									k2            = self.Position(-1)
									self.pi_w[k2] = 1
									self.len[k2]  = self.len[p1]
									k3            = self.Position(-1)
									self.pi_w[k3] = 1
									self.len[k3]  = self.len[p2]
									k4            = self.Position(-1)
									self.pi_w[k4] = 1
									self.len[k4]  = self.len[p2]
	
									for i4 in range(0, t11) :
										self.ind[k1][i4] = self.ind[p1][i4]
										self.ind[k2][i4] = self.ind[p1][i4]
									for i4 in range(t11, t12+1) :
										self.ind[k1][i4] = self.ind[p2][t21+i4-t11]
										self.ind[k2][i4] = self.ind[p2][t22-i4+t11]
									for i4 in range(t12+1, self.len[p1]) :
										self.ind[k1][i4] = self.ind[p1][i4]
										self.ind[k2][i4] = self.ind[p1][i4]
									for i4 in range(0, t21) :
										self.ind[k3][i4] = self.ind[p2][i4]
										self.ind[k4][i4] = self.ind[p2][i4]
									for i4 in range(t21, t22+1) :
										self.ind[k3][i4] = self.ind[p1][t11+i4-t21]
										self.ind[k4][i4] = self.ind[p1][t12-i4+t21]
									for i4 in range(t22+1, self.len[p2]) :
										self.ind[k3][i4] = self.ind[p2][i4]
										self.ind[k4][i4] = self.ind[p2][i4]
	
	#######################################
	# 突然変異（対立遺伝子との置き換え）
	#      pr : 突然変異率
	#######################################

	def M_alle(self, pr) :
	
				# データのチェックと初期設定
		if self.dup_a == 0 :
			print("***error  突然変異方法が不適当 (M_alle)")
				# 実行
		for i1 in range(0, self.size+self.max_ch) :
			if self.pi_w[i1] == 1 :
				for i2 in range(0, self.len[i1]) :
					if random() <= pr :
						lid = int(random() * (self.allele_u - self.allele_l + 1) + self.allele_l)
						if lid > self.allele_u :
							lid = self.allele_u
						if lid != self.ind[i1][i2] :
							self.ind[i1][i2] = lid

	######################################################################
	# 突然変異（移動．２点を選択し，２番目の遺伝子を１番目の遺伝子の前に
	#           移動する）
	#      pr : 突然変異率
	######################################################################

	def M_move(self, pr) :
	
		for i1 in range(0, self.size+self.max_ch) :
	
			if self.pi_w[i1] == 1 and random() <= pr :
				# 位置の決定
						# p1
				p1 = int(random() * self.len[i1])
				if p1 >= self.len[i1] :
					p1 = self.len[i1] - 1
						# p2
				p2 = p1
				sw = 0
				while sw == 0 :
					p2 = int(random() * self.len[i1])
					if p2 >= self.len[i1] :
						p2 = self.len[i1] - 1
					if p2 != p1 :
						sw = 1
				# 実行
				if p2 > p1 :
					ld = self.ind[i1][p2]
					for i2 in range(p2, p1, -1) :
						self.ind[i1][i2] = self.ind[i1][i2-1]
					self.ind[i1][p1] = ld
				else :
					ld = self.ind[i1][p2]
					for i2 in range(p2, p1-1) :
						self.ind[i1][i2] = self.ind[i1][i2+1]
					self.ind[i1][p1-1] = ld
	
	########################################################
	# 突然変異（逆位．２点間の遺伝子順序を逆に並べ替える）
	#      pr : 突然変異率
	#      wd : >0 : 幅を固定
	#           =0 : 幅をランダム(default)
	########################################################

	def M_inv(self, pr, wd = 0) :
	
		for i1 in range(0, self.size+self.max_ch) :
	
			if self.pi_w[i1] == 1 and random() <= pr :
				# 区間の決定
				if wd == 0 :
					p1 = int(random() * self.len[i1])
					if p1 >= self.len[i1] :
						p1 = self.len[i1] - 1
					sw = 0
					p2 = p1
					while sw == 0 :
						p2 = int(random() * self.len[i1])
						if p2 >= self.len[i1] :
							p2 = self.len[i1] - 1
						if p2 != p1 :
							sw = 1
					if p1 > p2 :
						p  = p1
						p1 = p2
						p2 = p
	
				else :
					p1 = self.len[i1]
					while p1 > self.len[i1]-2 :
						p1 = int(random() * self.len[i1])
					p2 = p1 + wd - 1
					if p2 >= self.len[i1] :
						p2 = self.len[i1] - 1
				# 実行
				sw = 0
				while sw == 0 :
					lid              = self.ind[i1][p1]
					self.ind[i1][p1] = self.ind[i1][p2]
					self.ind[i1][p2] = lid
					p1 += 1
					p2 -= 1
					if p1 >= p2 :
						sw = 1
	
	######################################################################
	# 突然変異（スクランブル．２点間の遺伝子順序をランダムに並べ替える）
	#      pr : 突然変異率
	#      wd : >0 : 幅を固定
	#           =0 : 幅をランダム(default)
	######################################################################

	def M_scram(self, pr, wd = 0) :
	
		for i1 in range(0, self.size+self.max_ch) :
	
			if self.pi_w[i1] == 1 and random() <= pr :
				# 区間の決定
				if wd == 0 :
					p1 = int(random() * self.len[i1])
					if p1 >= self.len[i1] :
						p1 = self.len[i1] - 1
					sw = 0
					p2 = p1
					while sw == 0 :
						p2 = int(random() * self.len[i1])
						if p2 >= self.len[i1] :
							p2 = self.len[i1] - 1
						if p2 != p1 :
							sw = 1
					if p1 > p2 :
						p  = p1
						p1 = p2
						p2 = p
	
				else :
					p1 = self.len[i1]
					while p1 > self.len[i1]-2 :
						p1 = int(random() * self.len[i1])
					p2 = p1 + wd - 1
					if p2 >= self.len[i1] :
						p2 = self.len[i1] - 1
				# 実行
				for i2 in range(p1, p2+1) :
					p = int(random() * (p2 - p1 + 1) + p1)
					if p > p2 :
						p = p2
					ld               = self.ind[i1][i2]
					self.ind[i1][i2] = self.ind[i1][p]
					self.ind[i1][p]  = ld
	
	######################################################################
	# 突然変異（転座．２点間の遺伝子を他の位置のものと置き換える．ただし
	#           重複部分はそのままとする）
	#      pr : 突然変異率
	#      wd : >0 : 幅を固定
	#           =0 : 幅をランダム(default)
	######################################################################

	def M_chg(self, pr, wd = 0) :
	
		for i1 in range(0, self.size+self.max_ch) :
	
			if self.pi_w[i1] == 1 and random() <= pr :
				# 区間等の決定（[p1,p2]と[p3,p4]の入れ替え）
						# p1
				p1 = int(random() * self.len[i1])
				if p1 >= self.len[i1] :
					p1 = self.len[i1] - 1
						# p3
				sw = 0
				p3 = p1
				while sw == 0 :
					p3 = int(random() * self.len[i1])
					if p3 >= self.len[i1] :
						p3 = self.len[i1] - 1
					if p3 != p1 :
						sw = 1
						# 小さい方をp1,p2にする
				if p1 > p3 :
					p  = p1
					p1 = p3
					p3 = p
						# p4, p2
				p4 = p1 + wd - 1
				if wd == 0 :
					p4 = int(random() * (self.len[i1] - p3)) + p3
				if p4 >= self.len[i1] :
					p4 = self.len[i1] - 1
				p2 = p1 + (p4 - p3)
						# 重複部分のチェック
				if p2 >= p3 :
					p  = p3 - 1
					p3 = p2 + 1
					p2 = p
					p4 = p3 + (p2 - p1)
				# 実行
				p = p3
				for i2 in range(p1, p2+1) :
					ld                = self.ind[i1][i2]
					self.ind[i1][i2]  = self.ind[i1][p]
					self.ind[i1][p]   = ld
					p                += 1
	
	######################################################################
	# 突然変異（重複．２点間の遺伝子を他の位置にコピーする
	#      pr : 突然変異率
	#      wd : >0 : 幅を固定
	#           =0 : 幅をランダム(deafult)
	######################################################################

	def M_dup(self, pr, wd = 0) :
	
				# データのチェック
		if self.dup_a == 0 :
			print("***error  突然変異方法が不適当 (M_dup)")
				# 実行
		for i1 in range(0, self.size+self.max_ch) :
	
			if self.pi_w[i1] == 1 and random() <= pr :
						# 区間の決定（[p1,p2]を[p3,p4]にコピー）
							# p1
				p1 = int(random() * self.len[i1])
				if p1 >= self.len[i1] :
					p1 = self.len[i1] - 1
							# p3
				sw = 0
				p3 = p1
				while sw == 0 :
					p3 = int(random() * self.len[i1])
					if p3 >= self.len[i1] :
						p3 = self.len[i1] - 1
					if p3 != p1 :
						sw = 1
							# 区間を決める
				p2 = p1
				p4 = p1
				if p3 > p1 :
					p4 = p3 + wd - 1
					if wd == 0 :
						p4 = int(random() * (self.len[i1] - p3)) + p3
					if p4 >= self.len[i1] :
						p4 = self.len[i1] - 1
					p2 = p1 + (p4 - p3)
				else :
					p2 = p1 + wd - 1
					if wd == 0 :
						p2 = int(random() * (self.len[i1] - p1)) + p1
					if p2 >= self.len[i1] :
						p2 = self.len[i1] - 1
					p4 = p3 + (p2 - p1)
						# 実行
				p = p4
				for i2 in range(p2, p1-1, -1) :
					self.ind[i1][p] = self.ind[i1][i2]
					p -= 1
	
	######################################################
	# 突然変異（摂動．値をある量だけ変化させる）
	#      pr : 突然変異率
	#      method : =0 : 正規分布(default)
	#               =1 : 一様分布
	#      m : 平均または一様分布の下限(default=0.0)
	#      s : 標準偏差または一様分布の上限(default=1.0)
	######################################################

	def M_per(self, pr, method = 0, m = 0.0, s = 1.0) :
	
		wd = 0.0
				# データのチェックと初期設定
		if self.dup_a == 0 :
			print("***error  突然変異方法が不適当 (M_per)")
	
		if method > 0 :
			wd = s - m
				# 実行
		for i1 in range(0, self.size+self.max_ch) :
			if self.pi_w[i1] == 1 :
				for i2 in range(0, self.len[i1]) :
					if random() <= pr :
						if method == 0 :
							w = normalvariate(m, s)
						else :
							w = random() * wd
							if random() < 0.5 :
								w = -w
						x1 = float(self.ind[i1][i2]) + w
						if x1 > self.allele_u :
							x1 = self.allele_u
						else :
							if x1 < self.allele_l :
								x1 = self.allele_l
						self.ind[i1][i2] = int(x1)
	
	##############################################
	# 突然変異（挿入．ある長さの遺伝子を挿入する）
	#      pr : 突然変異率
	#      wd : >0 : 幅を固定
	#           =0 : 幅をランダム(default)
	##############################################

	def M_ins(self, pr, wd = 0) :
	
				# データのチェック
		if self.dup_a == 0 or self.min_len < 0 :
			print("***error  突然変異方法が不適当 (M_ins)")
				# 実行
		for i1 in range(0, self.size+self.max_ch) :
	
			if self.pi_w[i1] == 1 and random() <= pr :
						# 挿入位置の決定
				p = int(random() * (self.len[i1] + 1))
				if p > self.len[i1] :
					p = self.len[i1]
						# 挿入する遺伝子長の決定
				l = wd
				if wd == 0 :
					l = int(random() * (self.max_len - self.len[i1] + 1))
				if l > self.max_len-self.len[i1] :
					l = self.max_len - self.len[i1]
				else :
					if l <= 0 :
						l = 1
						# 実行
							# 挿入場所の確保
				if p < self.len[i1] :
					for i2 in range(self.len[i1]+l-1, p-1, -1) :
						self.ind[i1][i2] = self.ind[i1][i2-l]
							# 挿入場所の遺伝子の決定
				for i2 in range(p, p+l) :
					ld = int(random() * (self.allele_u - self.allele_l + 1) + self.allele_l)
					if ld > self.allele_u :
						ld = self.allele_u
					self.ind[i1][i2] = ld
	
				self.len[i1]  += l
	
	##############################################
	# 突然変異（削除．ある長さの遺伝子を削除する）
	#      pr : 突然変異率
	#      wd : >0 : 幅を固定
	#           =0 : 幅をランダム(default)
	##############################################

	def M_del(self, pr, wd = 0) :
	
				# データのチェック
		if self.dup_a == 0 or self.min_len < 0 :
			print("***error  突然変異方法が不適当 (M_del)")
				# 実行
		for i1 in range(0, self.size+self.max_ch) :
	
			if self.pi_w[i1] == 1 and random() <= pr :
						# 削除位置の決定
				p = int(random() * self.len[i1])
				if p >= self.len[i1] :
					p = self.len[i1] - 1
						# 削除する遺伝子長の決定
				max = self.len[i1] - p
				if self.len[i1]-self.min_len < self.len[i1]-p :
					max = self.len[i1] - self.min_len
				l = wd
				if wd == 0 :
					l = int(random() * max + 1)
				if l > max :
					l = max
						# 実行
				for i2 in range(0, self.len[i1]-p-l) :
					self.ind[i1][p+i2] = self.ind[i1][p+i2+l]
	
				self.len[i1]  -= l
	
	######################################################################
	# 淘汰（エリート・ルーレット選択）
	#      elite : エリートで残す個体数(default=0)
	#      s_method : ルーレット板の作成方法(default=1)
	#                   =0 : 適応度をそのまま使用
	#                   =1 : 最小値からの差（ただし，α以下の場合はα）
	#                   =2 : 評価値に順位をつけ，減少率βで線形化
	#      s_bias : α，または，method=2の場合は初期値(default=0)
	#      s_step : β(default=1)
	######################################################################

	def S_roul(self, elite = 0, s_method = 1, s_bias = 0.0, s_step = 1.0) :
	
		count = 0
		k     = 0
		n     = 0
				# 値のチェックと初期設定
		if s_method != 0 and s_method != 2 :
			s_method = 1
	
		if elite > self.size :
			print("***error  エリートで残す数が多すぎる (S_roul)")
	
		if s_method == 2 and s_step <= 0.0 :
			s_step = 1.0
	
		for i1 in range(0, self.size+self.max_ch) :
			self.s_w[i1] = 0
				# 重複個体を削除
		if self.dup_s == 0 :
			for i1 in range(0, self.size+self.max_ch) :
				if self.pi_w[i1] > 0 :
					for i2 in range(i1+1, self.size+self.max_ch) :
						if self.pi_w[i2] > 0 and self.len[i1] == self.len[i2] :
							sw = 0
							for i3 in range(0, self.len[i1]) :
								if self.ind[i1][i3] != self.ind[i2][i3] :
									sw = 1
									break
							if sw == 0 :
								self.pi_w[i2] = 0
	
		for i1 in range(0, self.size+self.max_ch) :
			if self.pi_w[i1] > 1 :
				n += 1
	
		if n < 0 or self.dup_s == 0 and n < self.size :
			print("***error  残す個体がない (S_roul)")
				# 淘汰して残す個体を選ぶ
						# エリートの選択
		sw = 0
	
		while k < elite and k < n and sw == 0 :
			max = -1
			for i1 in range(0, self.size+self.max_ch) :
				if self.pi_w[i1] > 1 and self.s_w[i1] == 0 :
					if max < 0 or self.pi[i1] > self.pi[max] :
						max = i1
			if max < 0 :
				sw = 1
			else :
				self.s_w[max]  = 1
				k             += 1
						# ルーレット選択
		while count < self.size+self.max_ch and k < self.size :
			p = self.Select(s_method, s_bias, s_step)
			if self.dup_s == 0 and self.s_w[p] > 0 :
				count += 1
			else :
				count        = 0
				self.s_w[p] += 1
				k           += 1
							# 選択に失敗した場合の処理
		if self.dup_s == 0 and k < self.size :
			i1 = 0
			while i1 < self.size+self.max_ch and k < self.size :
				if self.pi_w[i1] > 1 and self.s_w[i1] == 0 :
					self.s_w[i1]  = 1
					k            += 1
				i1 += 1
							# 複数回選択されたものの処理
		for i1 in range(0, self.size+self.max_ch) :
			if self.s_w[i1] == 0 :
				self.pi_w[i1] = 0
	
		for i1 in range(0, self.size+self.max_ch) :
			if self.s_w[i1] > 0 :
				if self.s_w[i1] > 1 :
					for i2 in range(2, self.s_w[i1]+1) :
						k            = self.Position(-1)
						self.len[k]  = self.len[i1]
						self.pi_w[k] = 2
						self.pi[k]   = self.pi[i1]
						for i3 in range(0, self.len[i1]) :
							self.ind[k][i3] = self.ind[i1][i3]
	
####################
# クラスTSPの定義
####################

class TSP ( Species ) :

	######################################
	# コンストラクタ
	#      name1 : Species定義ファイル名
	#      name2 : TSP定義ファイル名
	######################################

	def __init__(self, name1, name2) :
	
		Species.__init__(self, name1)   # 親のコンストラクタ
						# 基本データの入力
		inn = open(name2, "r")
	
		s = inn.readline().split()
		self.out_lvl = int(s[1])   # 出力レベル
		                       #   =0 : 最終出力だけ
		                       #   n>0 : ｎ世代毎に出力（負の時はファイル）
		self.out_m   = int(s[3])  # 出力方法
		                      #   =0 : すべてを出力
		                      #   =1 : 最大適応度の個体だけを出力

		s = inn.readline().split()
		self.o_file = s[1]   # 出力ファイル名
		self.out_d  = int(s[3])   # 表示間隔

		s = inn.readline().split()
		self.kosa_m   = int(s[1])   # 交叉方法
		                       #   =-1 : 交叉を使用しない
		                       #   =0 : 親のコピー
		                       #   =1 : 循環交叉
		                       #   =2 : 部分的交叉
		                       #   =3 : 順序交叉
		                       #   =4 : 一様順序交叉
		                       #   =5 : 一様位置交叉
		                       #   =6 : エッジ組み替え交叉
		                       #   =7 : サブツアー交叉
		self.kosa     = float(s[3])   # 交叉確率
		self.k_point  = int(s[5])   # 交差点の数（負の時は，1から-k_point間のランダム）
		self.k_vr     = int(s[7])   # =0 : 両親とも同じ位置で交叉
		                       # =1 : 両親が異なる位置で交叉（遺伝子長は可変）
		self.k_method = int(s[9])   # 交叉の時の親の選択方法
		                       #   =-1 : ランダム
		                       #   =0 : 適応度をそのまま使用
		                       #   =1 : 最小値からの差（ただし，α以下の場合はα）
		                       #   =2 : 評価値に順位をつけ，減少率βで線形化
		self.k_bias   = float(s[11])   # α，または，method=2の場合は初期値
		self.k_step   = float(s[13])   # β

		s = inn.readline().split()
		self.mute_m = int(s[1])   # 突然変異方法
		                       #   =-1 : 突然変異を使用しない
		                       #   =0 : 移動
		                       #   =1 : 逆位
		                       #   =2 : スクランブル
		                       #   =3 : 転座
		self.mute   = float(s[3])   # 突然変異率
		self.wd     = int(s[5])   # 突然変異に使用する部分遺伝子長
		self.m_mean = float(s[7])   # 摂動の平均値
		self.m_std  = float(s[9])   # 摂動の標準偏差

		s = inn.readline().split()
		self.elite    = int(s[1])   # エリート選択で残す数
		self.s_method = int(s[3])   # ルーレット板の作成方法
		                       #   =0 : 適応度をそのまま使用
		                       #   =1 : 最小値からの差（ただし，α以下の場合はα）
		                       #   =2 : 評価値に順位をつけ，減少率βで線形化
		self.s_bias   = float(s[5])   # α，または，s_method=2の場合は初期値
		self.s_step   = float(s[7])   # β

		s = inn.readline().split()
		self.n_city  = int(s[1])   # 都市の数
		self.max_gen = int(s[3])   # 最大世代交代数

		s = inn.readline().split()
		self.kinbo    = int(s[1])   # 近傍探索（０：行わない，１：行う）
		self.neib     = int(s[3])   # 近傍（2 or 3）

		s = inn.readline().split()
		self.sel = int(s[1])   # エッジの選択方法
		                       #   =0 : 最良のものを選択
		                       #   =1 : 最初のものを選択
	
		if self.kinbo > 0 and self.neib != 2 and self.neib != 3 :
			print("***error  近傍の値が不適当")
	
		if self.n_city != self.max_len :
			print("***error  都市数が不適当")
						# 都市の位置データ
		self.city = np.empty((self.n_city, 2), np.int)
		for i1 in range(0, self.n_city) :
			s = inn.readline().split()
			self.city[i1][0] = int(s[0])
			self.city[i1][1] = int(s[1])
						# 距離テーブル
		self.rg = np.empty((self.n_city, self.n_city), np.int)
	
		for i1 in range(0, self.n_city) :
			for i2 in range(i1+1, self.n_city) :
				x               = self.city[i2][0] - self.city[i1][0]
				y               = self.city[i2][1] - self.city[i1][1]
				self.rg[i1][i2] = int(sqrt(x * x + y * y) + 0.5)
	
		for i1 in range(1, self.n_city) :
			for i2 in range(0, i1) :
				self.rg[i1][i2] = self.rg[i2][i1]
	
		inn.close()
	
	###############
	# 全体の制御
	###############

	def Control(self) :
	
		gen = 1
						# 初期集団の発生
		self.Init_std()
						# 評価
		if self.kinbo > 0 :
			self.Kinbo()
		else :
			self.Adap()
						# 出力
		print("***世代 " + str(gen) + " 適応度 max " + str(self.max) + " (" + str(self.max_n) + ") mean " + str(self.mean))
	
		if abs(self.out_lvl) > 0 :
			self.Output(gen)
						# 世代交代
		for gen in range(2, self.max_gen+1) :
							# 交叉
			if self.kosa_m == 0 :
				C_copy()   # 親のコピー
			elif self.kosa_m == 1 :
				self.C_cycle(self.kosa)   # 循環交叉
			elif self.kosa_m == 2 :
				self.C_part(self.kosa)   # 部分的交叉
			elif self.kosa_m == 3 :
				self.C_seq(self.kosa)   # 順序交叉
			elif self.kosa_m == 4 :
				self.C_useq(self.kosa)   # 一様順序交叉
			elif self.kosa_m == 5 :
				self.C_upos(self.kosa)   # 一様位置交叉
			elif self.kosa_m == 6 :
				self.C_edge(self.kosa)   # エッジ組み替え交叉
			elif self.kosa_m == 7 :
				self.C_sub(self.kosa, self.k_point)   # サブツアー交叉
							# 突然変異
			if self.mute_m == 0 :
				self.M_move(self.mute)   # 移動
			elif self.mute_m == 1 :
				self.M_inv(self.mute)   # 逆位
			elif self.mute_m == 2 :
				self.M_scram(self.mute)   # スクランブル
			elif self.mute_m == 3 :
				self.M_chg(self.mute)   # 転座
							# 適応度
			if self.kinbo > 0 :
				self.Kinbo()
			else :
				self.Adap()
							# 淘汰
			self.S_roul(self.elite)
							# 出力
			if gen%self.out_d == 0 :
				print("***世代 " + str(gen) + " 適応度 max " + str(self.max) + " (" + str(self.max_n) + ") mean " + str(self.mean))
	
			if abs(self.out_lvl) > 0 :
				if gen%abs(self.out_lvl) == 0 :
					self.Output(gen)
	
		gen        -= 1
		k1          = self.out_m
		self.out_m  = 0
		print("***世代 " + str(gen) + " 適応度 max " + str(self.max) + " (" + str(self.max_n) + ") mean " + str(self.mean))
		self.Output(gen)
		self.out_m = k1
	
	##########################
	# 距離の計算
	#      n_c : 都市の数
	#      p : 都市番号
	#      return : 距離（負）
	##########################

	def Kyori(self, n_c, p) :
	
		r  = 0
		n1 = p[0]
	
		for i1 in range(1, n_c) :
			n2  = p[i1]
			r  -= self.rg[n1][n2]
			n1  = n2

		n2  = p[0]
		r  -= self.rg[n1][n2]

		return r

	################
	# 適応度の計算
	################

	def Adap(self) :
	
		k          = 0
		self.mean  = 0.0
		self.max   = 0.0
		self.max_n = -1
	
		for i1 in range(0, self.size+self.max_ch) :
			if self.pi_w[i1] == 1 :
				self.pi_w[i1] = 2
				self.pi[i1]   = self.Kyori(self.len[i1], self.ind[i1])
			if self.pi_w[i1] > 0 :
				k         += 1
				self.mean += self.pi[i1]
				if self.max_n < 0 or self.pi[i1] > self.max :
					self.max   = self.pi[i1]
					self.max_n = i1
	
		if k > 0 :
			self.mean /= k

	######################################
	# エッジの入れ替え
	#      n_city : 都市の数
	#      seq : 訪問する順番
	#      r_m : 距離の負値
	#      return : =0 : 改善がなかった
	#               =1 : 改善があった
	######################################

	def Change(self, n_city, seq, r_m) :
	
		ch  = 0
		sw  = 0
		max = r_m[0]
	
		n3  = int(random() * (n_city - 2))
		if n3 > n_city-3 :
			n3 = n_city - 3
	                         # ２近傍
		i1 = 0
		while i1 <= n_city-3 and ch == 0 :
	
			if n3 == 0 :
				n1 = n_city - 2
			else :
				n1 = n_city - 1
	
			i2 = n3 + 2
			while i2 <= n1 and ch == 0 :
	                              # 枝の場所（(n3,n3+1), (k1,k2)）
				k1 = i2
				if i2 == n_city-1 :
					k2 = 0
				else :
					k2 = i2 + 1
	                              # 枝の入れ替え
				self.kou1[0] = seq[n3]
				k            = 1
				for i3 in range(k1, n3, -1) :
					self.kou1[k]  = seq[i3]
					k            += 1
	
				nn = k2
				while nn != n3 :
					self.kou1[k] = seq[nn]
					k  += 1
					nn += 1
					if nn > n_city-1 :
						nn = 0
	                              # 評価
				r = self.Kyori(n_city, self.kou1)
	
				if r > max :
					max = r
					sw  = 1
					for i3 in range(0, n_city) :
						self.kou2[i3] = self.kou1[i3]
					if self.sel > 0 :
						ch = 1
				i2 += 1
	
			n3 += 1
			if n3 > n_city-3 :
				n3 = 0
			i1 += 1
	                         # ３近傍
		if self.neib == 3 and ch == 0 :
	
			i1 = 0
			while i1 <= n_city-3 and ch == 0 :
	
				n1 = n_city - 2
				n2 = n_city - 1
	
				i2 = n3 + 1
				while i2 <= n1 and ch == 0 :
	
					i3 = i2 + 1
					while i3 <= n2 and ch == 0 :
	                              # 枝の場所（(n3,n3+1), (i2,i2+1), (k1,k2)）
						k1 = i3
						k2 = k1
						if i3 == n_city-1 :
							k2 = 0
						else :
							k2 = i3 + 1
	                              # 枝の入れ替えと評価
	                                   # 入れ替え（その１）
						self.kou1[0] = seq[n3]
						k            = 1
						for i4 in range(i2, n3, -1) :
							self.kou1[k]  = seq[i4]
							k            += 1
	
						for i4 in range(k1, i2, -1) :
							self.kou1[k]  = seq[i4]
							k            += 1
	
						nn = k2
						while nn != n3 :
							self.kou1[k] = seq[nn]
							k  += 1
							nn += 1
							if nn > n_city-1 :
								nn = 0
	                                   # 評価（その１）
						r = self.Kyori(n_city, self.kou1)
	
						if r > max :
							max = r
							sw  = 1
							for i3 in range(0, n_city) :
								self.kou2[i3] = self.kou1[i3]
							if self.sel > 0 :
								ch = 1
	                                   # 入れ替え（その２）
						self.kou1[0] = seq[n3]
						k            = 1
						for i4 in range(k1, i2, -1) :
							self.kou1[k]  = seq[i4]
							k            += 1
	
						for i4 in range(n3+1, i2+1) :
							self.kou1[k]  = seq[i4]
							k            += 1
	
						nn = k2
						while nn != n3 :
							self.kou1[k] = seq[nn]
							k  += 1
							nn += 1
							if nn > n_city-1 :
								nn = 0
	                                   # 評価（その２）
						r = self.Kyori(n_city, self.kou1)
	
						if r > max :
							max = r
							sw  = 1
							for i3 in range(0, n_city) :
								self.kou2[i3] = self.kou1[i3]
							if self.sel > 0 :
								ch = 1
	                                   # 入れ替え（その３）
						self.kou1[0] = seq[n3]
						k            = 1
						for i4 in range(i2+1, k1+1) :
							self.kou1[k]  = seq[i4]
							k            += 1
	
						for i4 in range(i2, n3, -1) :
							self.kou1[k]  = seq[i4]
							k            += 1
	
						nn = k2
						while nn != n3 :
							self.kou1[k] = seq[nn]
							k  += 1
							nn += 1
							if nn > n_city-1 :
								nn = 0
	                                   # 評価（その３）
						r = self.Kyori(n_city, self.kou1)
	
						if r > max :
							max = r
							sw  = 1
							for i3 in range(0, n_city) :
								self.kou2[i3] = self.kou1[i3]
							if self.sel > 0 :
								ch = 1
	                                   # 入れ替え（その４）
						self.kou1[0] = seq[n3]
						k            = 1
						for i4 in range(i2+1, k1+1) :
							self.kou1[k]  = seq[i4]
							k            += 1
	
						for i4 in range(n3+1, i2+1) :
							self.kou1[k]  = seq[i4]
							k            += 1
	
						nn = k2
						while nn != n3 :
							self.kou1[k] = seq[nn]
							k  += 1
							nn += 1
							if nn > n_city-1 :
								nn = 0
	                                   # 評価（その４）
						r = self.Kyori(n_city, self.kou1)
	
						if r > max :
							max = r
							sw  = 1
							for i3 in range(0, n_city) :
								self.kou2[i3] = self.kou1[i3]
							if self.sel > 0 :
								ch = 1
						i3 += 1
					i2 += 1
	
				n3 += 1
				if n3 > n_city-3 :
					n3 = 0
				i1 += 1
	                         # 設定
		if sw > 0 :
			r_m[0] = max
			for i1 in range(0, n_city) :
				seq[i1] = self.kou2[i1]
	
		return sw
	
	##############
	# 近傍の探索
	##############

	def Kinbo(self) :
	
		k          = 0
		self.max   = 0.0
		self.max_n = -1
		self.mean  = 0.0
	
		for i1 in range(0, self.size+self.max_ch) :
			if self.pi_w[i1] == 1 :
				self.pi_w[i1] = 2
				sw            = 1
				r             = self.Kyori(self.len[i1], self.ind[i1])
				while sw > 0 :
					r  = np.empty(1, np.int)
					sw = self.Change(self.len[i1], self.ind[i1], r)
				self.pi[i1] = r[0]
			if self.pi_w[i1] > 0 :
				k         += 1
				self.mean += self.pi[i1]
				if self.max_n < 0 or self.pi[i1] > self.max :
					self.max   = self.pi[i1]
					self.max_n = i1
	
		if k > 0 :
			self.mean /= k
	
	#############################
	# 結果の出力
	#      gen : 現在の世代番号
	#############################

	def Output(self, gen) :
	
		k  = 0
		pr = -1

		if self.out_lvl >= 0 :
			print("   出力先は（0:出力なし，n:画面にｎ個づつ，-1:ファイル）？ ", end=" ")
			pr = int(input())

		if pr != 0 :
						# 出力先の決定と評価値の出力
			if pr > 0 :
				out = sys.stdout
				input("")
			else :
				now = datetime.today().time().isoformat()
				out = open(self.o_file, "a")
				out.write("***世代 " + str(gen) + " 適応度 max " + str(self.max) + " (" + str(self.max_n) + ") mean " + str(self.mean) + " 時間 " + now + "\n")
						# 巡回順序の出力
			if self.out_m == 0 :
				for i1 in range(0, self.len[self.max_n]) :
					n = self.ind[self.max_n][i1]
					out.write(str(n) + " " + str(self.city[n][0]) + " " + str(self.city[n][1]) + "\n")
					if pr > 0 :
						k += 1
						if k == pr :
							input("")
							k = 0

			if pr < 0 :
				out.close()

				# 入力ミス
if len(sys.argv) <= 1 :
	print("***error  ファイル名を入力して下さい")
				# 入力ＯＫ
else :
					# データの数と入力データファイル名の入力
	inn     = open(sys.argv[1], "r")

	ss      = inn.readline()
	n       = int(ss)   # データの数
	i_file1 = []
	i_file2 = []

	for i1 in range(0, n) :
		s = inn.readline().split()
		i_file1.append(s[0])
		i_file2.append(s[1])

	inn.close()
					# 実行（乱数の初期値を変える）
	for i1 in range(0, n) :

		print("\n+++++ケース " + str(i1+1) + "+++++")
		seed(1000 * i1 + 1234567);

		tsp = TSP(i_file1[i1], i_file2[i1])

		tsp.Control()

