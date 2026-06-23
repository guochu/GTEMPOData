set term epslatex standalone size 8.6cm,6.5cm font ",8" header\
'\usepackage{scalerel}'

lw = 3

# beta=5 n: 红色系 -- ED 粉色 + TEMPO 红色
set linestyle 1 lw lw lc rgb "pink"
set linestyle 2 lw lw lc rgb "red"

# beta=Inf n: 青色系 -- ED/实线 浅青 + TEMPO/虚线 深青
set linestyle 3 lw lw lc rgb "#76D7C4"
set linestyle 4 lw lw lc rgb "#009E73"

# beta=5 lt: real 橙色系, imag 蓝色系; ED 实线, TEMPO 虚线
set linestyle 5 lw lw lc rgb "#FF8C00"   # ED real 亮橙
set linestyle 6 lw lw lc rgb "#D2691E"   # TEMPO real 深橙
set linestyle 7 lw lw lc rgb "#87CEFA"   # ED imag 浅蓝
set linestyle 8 lw lw lc rgb "#4169E1"   # TEMPO imag 深蓝

set style arrow 1 noborder lw 2 size 1,20

set out "fig.tex"
set multiplot

# common settings
set format x '\small{%g}'
set format y '\small{%.1f}'
set lmargin 0
set rmargin 0
set tmargin 0
set bmargin 0

w=0.36
h=0.36

y1=0.55
y2=0.10

xL=0.12
xR=0.60

# (a) alpha=0.04, n
set size w,h
set origin xL,y1
set xlabel '$t$' offset 0,1.2
set ylabel '$\langle\hat{n}\rangle$' offset 2,0
set xrange [0:2.5]
set yrange [0.0:1.05]
set xtics 0,0.5,2.5 offset 0,0.5
set ytics 0,0.2,1.0
set mxtics 2
set mytics 2
# 顶部中心 alpha 标签
set label '$\alpha=0.04$' at graph 0.5,1.07 center
# 内部 legend: beta=5 和 beta=Inf 分两行 (各实线+虚线叠加)
set arrow from graph 0.55,0.95 to graph 0.67,0.95 ls 1 nohead
set arrow from graph 0.55,0.91 to graph 0.67,0.91 ls 2 dt 2 nohead
set label '\small{$\beta=5$}' at graph 0.69,0.93
set arrow from graph 0.55,0.84 to graph 0.67,0.84 ls 3 nohead
set arrow from graph 0.55,0.80 to graph 0.67,0.80 ls 4 dt 2 nohead
set label '\small{$\beta=\infty$}' at graph 0.69,0.82
set label '(a)' at graph -0.32,0.95
unset key
plot "../data/freebosons/alpha0.04/panel00_n_ed_beta5.csv"      u 1:2 w l ls 1 notitle,\
     "../data/freebosons/alpha0.04/panel00_n_ed_betaInf.csv"    u 1:2 w l ls 3 notitle,\
     "../data/freebosons/alpha0.04/panel00_n_tempo_beta5.csv"   u 1:2 w l ls 2 dt 2 notitle,\
     "../data/freebosons/alpha0.04/panel00_n_tempo_betaInf.csv" u 1:2 w l ls 4 dt 2 notitle

# (b) alpha=0.08, n
unset label
unset arrow
unset key
set size w,h
set origin xR,y1
set xlabel '$t$' offset 0,1.2
set ylabel '$\langle\hat{n}\rangle$' offset 2,0
set xrange [0:2.5]
set yrange [0.0:1.05]
set xtics 0,0.5,2.5 offset 0,0.5
set ytics 0,0.2,1.0
set mxtics 2
set mytics 2
# 顶部中心 alpha 标签
set label '$\alpha=0.08$' at graph 0.5,1.07 center
set label '(b)' at graph -0.32,0.95
unset key
plot "../data/freebosons/alpha0.08/panel01_n_ed_beta5.csv"      u 1:2 w l ls 1 notitle,\
     "../data/freebosons/alpha0.08/panel01_n_ed_betaInf.csv"    u 1:2 w l ls 3 notitle,\
     "../data/freebosons/alpha0.08/panel01_n_tempo_beta5.csv"   u 1:2 w l ls 2 dt 2 notitle,\
     "../data/freebosons/alpha0.08/panel01_n_tempo_betaInf.csv" u 1:2 w l ls 4 dt 2 notitle

# (c) alpha=0.04, lt
unset label
unset arrow
unset key
set size w,h
set origin xL,y2
set xlabel '$t$' offset 0,1.2
set ylabel '$G^<(t)$' offset 3,0
set xrange [0:2.5]
set yrange [-0.75:1.05]
set xtics 0,0.5,2.5 offset 0,0.5
set ytics -0.5,0.5,1.0
set mxtics 2
set mytics 2
set label '(c)' at graph -0.32,0.95
set arrow from graph 0.35,0.93 to graph 0.47,0.93 ls 5 nohead
set arrow from graph 0.35,0.89 to graph 0.47,0.89 ls 6 dt 2 nohead
set label '\small{real,$\beta=5$}' at graph 0.49,0.91
set arrow from graph 0.35,0.82 to graph 0.47,0.82 ls 7 nohead
set arrow from graph 0.35,0.78 to graph 0.47,0.78 ls 8 dt 2 nohead
set label '\small{imag,$\beta=5$}' at graph 0.49,0.80
plot "../data/freebosons/alpha0.04/panel10_lt_ed_beta5_lt0.csv"      u 1:2 w l ls 7 notitle,\
     "../data/freebosons/alpha0.04/panel10_lt_ed_beta5_lt1.csv"      u 1:2 w l ls 5 notitle,\
     "../data/freebosons/alpha0.04/panel10_lt_tempo_beta5_lt0.csv"   u 1:2 w l ls 8 dt 2 notitle,\
     "../data/freebosons/alpha0.04/panel10_lt_tempo_beta5_lt1.csv"   u 1:2 w l ls 6 dt 2 notitle

# (d) alpha=0.08, lt
unset label
unset arrow
unset key
set size w,h
set origin xR,y2
set xlabel '$t$' offset 0,1.2
set ylabel '$G^<(t)$' offset 3,0
set xrange [0:2.5]
set yrange [-0.75:1.05]
set xtics 0,0.5,2.5 offset 0,0.5
set ytics -0.5,0.5,1.0
set mxtics 2
set mytics 2
set label '(d)' at graph -0.32,0.95
plot "../data/freebosons/alpha0.08/panel11_lt_ed_beta5_lt0.csv"      u 1:2 w l ls 7 notitle,\
     "../data/freebosons/alpha0.08/panel11_lt_ed_beta5_lt1.csv"      u 1:2 w l ls 5 notitle,\
     "../data/freebosons/alpha0.08/panel11_lt_tempo_beta5_lt0.csv"   u 1:2 w l ls 8 dt 2 notitle,\
     "../data/freebosons/alpha0.08/panel11_lt_tempo_beta5_lt1.csv"   u 1:2 w l ls 6 dt 2 notitle

unset multiplot
