set term epslatex standalone size 8.6cm,6.5cm font ",8" header\
'\usepackage{scalerel}'

lw = 3

# beta=5: 红色系 (chi 小 -> 浅, chi 大 -> 深) -- JC 虚线 + Rabi 实线
set linestyle 1 lw lw lc rgb "#FCA5A5"
set linestyle 2 lw lw lc rgb "#EF4444"
set linestyle 3 lw lw lc rgb "#7F1D1D"

# beta=Inf: 蓝色系 -- JC 虚线 + Rabi 实线
set linestyle 4 lw lw lc rgb "#93C5FD"
set linestyle 5 lw lw lc rgb "#2563EB"
set linestyle 6 lw lw lc rgb "#1E3A8A"

# Rabi + beta=5: 同红色系 (实线)
set linestyle 7 lw lw lc rgb "#FCA5A5"
set linestyle 8 lw lw lc rgb "#EF4444"
set linestyle 9 lw lw lc rgb "#7F1D1D"

# Rabi + beta=Inf: 同蓝色系 (实线)
set linestyle 10 lw lw lc rgb "#93C5FD"
set linestyle 11 lw lw lc rgb "#2563EB"
set linestyle 12 lw lw lc rgb "#1E3A8A"

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

y1=0.53
y2=0.08

xL=0.12
xR=0.60

# (a) alpha=0.01
set size w,h
set origin xL,y1
set xlabel '$t$' offset 0,1.2
set ylabel '$\langle\hat{\sigma}_z\rangle$' offset 2.8,0
set xrange [0:5]
set yrange [-0.4:0.5]
set xtics 0,1,5 offset 0,0.5
set ytics -0.4,0.2,0.4
set mxtics 2
set mytics 2
# 顶部图外 legend: 第1行 JC (实线), 第2行 Rabi (虚线), 跨 (a)+(b)
# 第1行: JC, 三个 chi (虚线)
set arrow from graph 0.05,1.22 to graph 0.17,1.22 ls 1 dt 2 nohead
set arrow from graph 0.05,1.18 to graph 0.17,1.18 ls 4 dt 2 nohead
set label '\small{JC,$\chi=30$}'   at graph 0.19,1.20
set arrow from graph 0.75,1.22 to graph 0.87,1.22 ls 2 dt 2 nohead
set arrow from graph 0.75,1.18 to graph 0.87,1.18 ls 5 dt 2 nohead
set label '\small{JC,$\chi=50$}'  at graph 0.89,1.20
set arrow from graph 1.50,1.22 to graph 1.62,1.22 ls 3 dt 2 nohead
set arrow from graph 1.50,1.18 to graph 1.62,1.18 ls 6 dt 2 nohead
set label '\small{JC,$\chi=80$}'  at graph 1.64,1.20
# 第2行: Rabi, 三个 chi (实线)
set arrow from graph 0.05,1.10 to graph 0.17,1.10 ls 7  nohead
set arrow from graph 0.05,1.06 to graph 0.17,1.06 ls 10 nohead
set label '\small{Rabi,$\chi=80$}'   at graph 0.19,1.08
set arrow from graph 0.75,1.10 to graph 0.87,1.10 ls 8  nohead
set arrow from graph 0.75,1.06 to graph 0.87,1.06 ls 11 nohead
set label '\small{Rabi,$\chi=140$}'  at graph 0.89,1.08
set arrow from graph 1.50,1.10 to graph 1.62,1.10 ls 9  nohead
set arrow from graph 1.50,1.06 to graph 1.62,1.06 ls 12 nohead
set label '\small{Rabi,$\chi=200$}' at graph 1.64,1.08
# beta=5 / beta=Inf 合并 legend (放在 (a) 内)
set arrow from graph 0.50,0.93 to graph 0.62,0.93 ls 3 dt 2 nohead
set arrow from graph 0.50,0.89 to graph 0.62,0.89 ls 9 nohead
set label '\small{$\beta=5$}'        at graph 0.64,0.91
set arrow from graph 0.50,0.78 to graph 0.62,0.78 ls 6 dt 2 nohead
set arrow from graph 0.50,0.74 to graph 0.62,0.74 ls 12 nohead
set label '\small{$\beta=\infty$}'   at graph 0.64,0.76
set label '(a)' at graph -0.32,0.95
set label '$\alpha=0.01$' at graph 0.05,0.10
unset key
plot "../data/fig2/panel00_alpha0.01/jc_beta5_chi30.csv"      u 1:2 w l ls 1 dt 2 notitle,\
     "../data/fig2/panel00_alpha0.01/jc_beta5_chi50.csv"      u 1:2 w l ls 2 dt 2 notitle,\
     "../data/fig2/panel00_alpha0.01/jc_beta5_chi80.csv"      u 1:2 w l ls 3 dt 2 notitle,\
     "../data/fig2/panel00_alpha0.01/jc_betaInf_chi30.csv"    u 1:2 w l ls 4 dt 2 notitle,\
     "../data/fig2/panel00_alpha0.01/jc_betaInf_chi50.csv"    u 1:2 w l ls 5 dt 2 notitle,\
     "../data/fig2/panel00_alpha0.01/jc_betaInf_chi80.csv"    u 1:2 w l ls 6 dt 2 notitle,\
     "../data/fig2/panel00_alpha0.01/rabi_beta5_chi80.csv"    u 1:2 w l ls 7  notitle,\
     "../data/fig2/panel00_alpha0.01/rabi_beta5_chi140.csv"   u 1:2 w l ls 8  notitle,\
     "../data/fig2/panel00_alpha0.01/rabi_beta5_chi200.csv"   u 1:2 w l ls 9  notitle,\
     "../data/fig2/panel00_alpha0.01/rabi_betaInf_chi80.csv"  u 1:2 w l ls 10 notitle,\
     "../data/fig2/panel00_alpha0.01/rabi_betaInf_chi140.csv" u 1:2 w l ls 11 notitle,\
     "../data/fig2/panel00_alpha0.01/rabi_betaInf_chi200.csv" u 1:2 w l ls 12 notitle

# (b) alpha=0.04
unset label
unset arrow
unset key
set size w,h
set origin xR,y1
set xlabel '$t$' offset 0,1.2
set ylabel '$\langle\hat{\sigma}_z\rangle$' offset 2.8,0
set xrange [0:5]
set yrange [-0.4:0.5]
set xtics 0,1,5 offset 0,0.5
set ytics -0.4,0.2,0.4
set mxtics 2
set mytics 2
# (b) 顶部 legend 已移到图最上方
set label '(b)' at graph -0.32,0.95
set label '$\alpha=0.04$' at graph 0.05,0.10
plot "../data/fig2/panel01_alpha0.04/jc_beta5_chi30.csv"      u 1:2 w l ls 1 dt 2 notitle,\
     "../data/fig2/panel01_alpha0.04/jc_beta5_chi50.csv"      u 1:2 w l ls 2 dt 2 notitle,\
     "../data/fig2/panel01_alpha0.04/jc_beta5_chi80.csv"      u 1:2 w l ls 3 dt 2 notitle,\
     "../data/fig2/panel01_alpha0.04/jc_betaInf_chi30.csv"    u 1:2 w l ls 4 dt 2 notitle,\
     "../data/fig2/panel01_alpha0.04/jc_betaInf_chi50.csv"    u 1:2 w l ls 5 dt 2 notitle,\
     "../data/fig2/panel01_alpha0.04/jc_betaInf_chi80.csv"    u 1:2 w l ls 6 dt 2 notitle,\
     "../data/fig2/panel01_alpha0.04/rabi_beta5_chi80.csv"    u 1:2 w l ls 7  notitle,\
     "../data/fig2/panel01_alpha0.04/rabi_beta5_chi140.csv"   u 1:2 w l ls 8  notitle,\
     "../data/fig2/panel01_alpha0.04/rabi_beta5_chi200.csv"   u 1:2 w l ls 9  notitle,\
     "../data/fig2/panel01_alpha0.04/rabi_betaInf_chi80.csv"  u 1:2 w l ls 10 notitle,\
     "../data/fig2/panel01_alpha0.04/rabi_betaInf_chi140.csv" u 1:2 w l ls 11 notitle,\
     "../data/fig2/panel01_alpha0.04/rabi_betaInf_chi200.csv" u 1:2 w l ls 12 notitle

# (c) alpha=0.08
unset label
unset arrow
unset key
set size w,h
set origin xL,y2
set xlabel '$t$' offset 0,1.2
set ylabel '$\langle\hat{\sigma}_z\rangle$' offset 1.7,0
set xrange [0:5]
set yrange [0.0:0.5]
set xtics 0,1,5 offset 0,0.5
set ytics 0.0,0.1,0.5
set mxtics 2
set mytics 2
set label '(c)' at graph -0.32,0.95
set label '$\alpha=0.08$' at graph 0.05,0.10
plot "../data/fig2/panel10_alpha0.08/jc_beta5_chi30.csv"      u 1:2 w l ls 1 dt 2 notitle,\
     "../data/fig2/panel10_alpha0.08/jc_beta5_chi50.csv"      u 1:2 w l ls 2 dt 2 notitle,\
     "../data/fig2/panel10_alpha0.08/jc_beta5_chi80.csv"      u 1:2 w l ls 3 dt 2 notitle,\
     "../data/fig2/panel10_alpha0.08/jc_betaInf_chi30.csv"    u 1:2 w l ls 4 dt 2 notitle,\
     "../data/fig2/panel10_alpha0.08/jc_betaInf_chi50.csv"    u 1:2 w l ls 5 dt 2 notitle,\
     "../data/fig2/panel10_alpha0.08/jc_betaInf_chi80.csv"    u 1:2 w l ls 6 dt 2 notitle,\
     "../data/fig2/panel10_alpha0.08/rabi_beta5_chi80.csv"    u 1:2 w l ls 7  notitle,\
     "../data/fig2/panel10_alpha0.08/rabi_beta5_chi140.csv"   u 1:2 w l ls 8  notitle,\
     "../data/fig2/panel10_alpha0.08/rabi_beta5_chi200.csv"   u 1:2 w l ls 9  notitle,\
     "../data/fig2/panel10_alpha0.08/rabi_betaInf_chi80.csv"  u 1:2 w l ls 10 notitle,\
     "../data/fig2/panel10_alpha0.08/rabi_betaInf_chi140.csv" u 1:2 w l ls 11 notitle,\
     "../data/fig2/panel10_alpha0.08/rabi_betaInf_chi200.csv" u 1:2 w l ls 12 notitle

# (d) alpha=0.12
unset label
unset arrow
unset key
set size w,h
set origin xR,y2
set xlabel '$t$' offset 0,1.2
set ylabel '$\langle\hat{\sigma}_z\rangle$' offset 1.7,0
set xrange [0:5]
set yrange [0.0:0.5]
set xtics 0,1,5 offset 0,0.5
set ytics 0.0,0.1,0.5
set mxtics 2
set mytics 2
set label '(d)' at graph -0.32,0.95
set label '$\alpha=0.12$' at graph 0.05,0.10
plot "../data/fig2/panel11_alpha0.12/jc_beta5_chi30.csv"      u 1:2 w l ls 1 dt 2 notitle,\
     "../data/fig2/panel11_alpha0.12/jc_beta5_chi50.csv"      u 1:2 w l ls 2 dt 2 notitle,\
     "../data/fig2/panel11_alpha0.12/jc_beta5_chi80.csv"      u 1:2 w l ls 3 dt 2 notitle,\
     "../data/fig2/panel11_alpha0.12/jc_betaInf_chi30.csv"    u 1:2 w l ls 4 dt 2 notitle,\
     "../data/fig2/panel11_alpha0.12/jc_betaInf_chi50.csv"    u 1:2 w l ls 5 dt 2 notitle,\
     "../data/fig2/panel11_alpha0.12/jc_betaInf_chi80.csv"    u 1:2 w l ls 6 dt 2 notitle,\
     "../data/fig2/panel11_alpha0.12/rabi_beta5_chi80.csv"    u 1:2 w l ls 7  notitle,\
     "../data/fig2/panel11_alpha0.12/rabi_beta5_chi140.csv"   u 1:2 w l ls 8  notitle,\
     "../data/fig2/panel11_alpha0.12/rabi_beta5_chi200.csv"   u 1:2 w l ls 9  notitle,\
     "../data/fig2/panel11_alpha0.12/rabi_betaInf_chi80.csv"  u 1:2 w l ls 10 notitle,\
     "../data/fig2/panel11_alpha0.12/rabi_betaInf_chi140.csv" u 1:2 w l ls 11 notitle,\
     "../data/fig2/panel11_alpha0.12/rabi_betaInf_chi200.csv" u 1:2 w l ls 12 notitle

unset multiplot
