#!/bin/bash
#--------------------------------------------------
# plot for GMT5
# 1次元時系列グラフ
#--------------------------------------------------

#----- 前処理・設定 ------------------------------

# GMTパラメーター
gmt gmtset FONT_ANNOT_PRIMARY=9p,Helvetica   # 主目盛りのフォント
gmt gmtset FONT_LABEL=10p,Helvetica          # 軸タイトルのフォント
# 日付関係(=前後にスペースを入れない)
gmt gmtset FORMAT_DATE_MAP=mm/dd # 日付目盛り(yy[yy]=年,mm=月,dd=日,o=月名,jjj=日数)
gmt gmtset FORMAT_TIME_PRIMARY_MAP=a # 主目盛りの月・週表示：[f|a|c]

# input
in01=fitness_Tochigi.txt              # 入力ファイル　変更箇所
in01_h=1                        # ヘッダー行数
in01_c=0,1                      # 入力する列(0始まり）

# output
ps=fitness_Tochigi.ps               # 出力psファイル　変更箇所

# グラフのサイズ
mp=X        # 地図の投影法(X:xy, M:メルカトル図法)
gx=10       # グラフの幅[cm]
gy=8        # グラフの高さ[cm]
#log_x=l                          # x軸を対数軸にするならl
#log_y=l                          # y軸を対数軸にするならl
jopt=$mp${gx}c${log_x}/${gy}c${log_y} # -Jオプション

# グラフの範囲(range)（-Rオプション）
rw=2013-02-24T00:00:00             # x軸の最小値(west) 変更箇所
re=2013-03-27T00:00:00             # x軸の最大値(east)　変更箇所
rs=10000                         # y軸の最小値(south)　変更箇所
rn=70000                        # y軸の最大値(north)　変更箇所

# 軸・目盛りの設定
ti_x="Date (2013)"          # x軸タイトル　変更箇所
ti_y="Fitness"       # y軸タイトル　変更箇所
tics_x=a7Df1D               # x軸の目盛り(a=主,f=副)　変更箇所
grid_x=                     # x軸のグリッド(g)(*)　変更箇所
tics_y=a10000f5000                # y軸の目盛り(a=主,f=副)　変更箇所
grid_y=g10000                  # y軸のグリッド(g)(*)　変更箇所
axis=WeSn                   # 軸表示（大文字：軸＋目盛り、小文字：軸）

# 線の設定
line01=0.8p,black               # 線：太さ,色,線種


#----- plot ------------------------------
echo "plot $ps"

# basemap(グリッド線のみ)
gmt psbasemap -J$jopt -R$rw/$re/$rs/$rn -Bx$grid_x -By$grid_y -B+n -P -K > $ps

# psxy(line)
gmt psxy $in01 -h$in01_h -i$in01_c -J -R -W$line01 -K -O >> $ps

# basemap(枠線)
gmt psbasemap -J -R -Bx$tics_x+l"$ti_x" -By$tics_y+l"$ti_y" -B$axis -O >> $ps


#----- 後処理 ------------------------------

# 画像変換(psconvert)の設定
pc_fmt=g            # g=png, b=bmp, j=jpg, t=tiff, e=eps, f=PDF, s=SVG
pc_res=1000          # 画像の解像度(DPI)
pc_mgn=0.1c         # グラフ外側のマージン
pc_aa=2             # アンチエイリアシングの設定(1,2,4)

# png変換
gmt psconvert -A$pc_mgn -E$pc_res -T$pc_fmt -Qt$pc_aa -Qg$pc_aa $ps

# 中間ファイルの削除
rm -f gmt.history gmt.conf
open $ps