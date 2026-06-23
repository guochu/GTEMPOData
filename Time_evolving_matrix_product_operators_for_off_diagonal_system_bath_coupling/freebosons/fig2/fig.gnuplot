set term epslatex standalone size 8.6cm,6.5cm font ",8" header\
'\usepackage{scalerel}'

lw = 3

# beta=5: 紫色, 空心圆
set linestyle 1 lw lw lc rgb "#8E44AD" pt 6

# beta=Inf: 青色, 空心方
set linestyle 2 lw lw lc rgb "#17A589" pt 4

set style arrow 1 noborder lw 2 size 1,20

set out "fig.tex"
set multiplot

# common settings
set format x '\small{%g}'
set format y '\small{%.2f}'
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

# (a) alpha=0.04, dt error (×100)
set size w,h
set origin xL,y1
set xlabel '$\delta t$' offset 0,1.2
set ylabel '$\mathcal{E}[\hat{n}]$' offset 2.2,0
set format y '\small{%.1f}'
set xrange [0:0.12]
set yrange [0:4.5]
set xtics 0,0.05,0.10 offset 0,0.5
set ytics 0,1,4
set mxtics 2
set mytics 2
set label '$\scriptstyle{\times10^{-2}}$' at graph 0,1.05
# legend
set arrow from graph 0.55,0.95 to graph 0.67,0.95 ls 1 nohead
set label '\small{$\beta=5$}' at graph 0.69,0.93
set arrow from graph 0.55,0.84 to graph 0.67,0.84 ls 2 nohead
set label '\small{$\beta=\infty$}' at graph 0.69,0.82
set label '(a)' at graph -0.32,0.95
unset key
plot "../data/freebosons_error/panel00_dt.csv"      u 1:($3*100) w lp ls 1 notitle,\
     "../data/freebosons_error/panel00_dt.csv"      u 1:($2*100) w lp ls 2 notitle

# (b) alpha=0.08, chi error (×10)
unset label
unset arrow
unset key
set size w,h
set origin xR,y1
set xlabel '$\chi$' offset 0,1.2
set ylabel '$\mathcal{E}[\hat{n}]$' offset 2.2,0
set format y '\small{%.1f}'
set xrange [0:65]
set yrange [0:1.2]
set xtics 0,20,60 offset 0,0.5
set ytics 0,0.4,1.2
set mxtics 2
set mytics 2
set label '$\scriptstyle{\times10^{-1}}$' at graph 0,1.05
set label '(b)' at graph -0.32,0.95
unset key
plot "../data/freebosons_error/panel01_chi.csv"     u 1:($3*10) w lp ls 1 notitle,\
     "../data/freebosons_error/panel01_chi.csv"     u 1:($2*10) w lp ls 2 notitle

# (c) alpha=0.04, k error (×10)
unset label
unset arrow
unset key
set size w,h
set origin xL,y2
set xlabel '$m$' offset 0,1.2
set ylabel '$\mathcal{E}[G^<]$' offset 2.2,0
set format y '\small{%.1f}'
set xrange [1:9]
set yrange [0:1.6]
set xtics 2,2,8 offset 0,0.5
set ytics 0,0.5,1.5
set mxtics 2
set mytics 2
set label '(c)' at graph -0.32,0.95
set label '$\scriptstyle{\times10^{-1}}$' at graph 0,1.05
unset key
plot "../data/freebosons_error/panel10_k.csv"       u 1:($3*10) w lp ls 1 notitle,\
     "../data/freebosons_error/panel10_k.csv"       u 1:($2*10) w lp ls 2 notitle

# (d) alpha=0.08, n error (×100)
unset label
unset arrow
unset key
set size w,h
set origin xR,y2
set xlabel '$n$' offset 0,1.2
set ylabel '$\mathcal{E}[G^<]$' offset 2.2,0
set format y '\small{%.1f}'
set xrange [1:11]
set yrange [0:5.5]
set xtics 2,2,10 offset 0,0.5
set ytics 0,1,5
set mxtics 2
set mytics 2
set label '(d)' at graph -0.32,0.95
set label '$\scriptstyle{\times10^{-2}}$' at graph 0,1.05
unset key
plot "../data/freebosons_error/panel11_n.csv"       u 1:($3*100) w lp ls 1 notitle,\
     "../data/freebosons_error/panel11_n.csv"       u 1:($2*100) w lp ls 2 notitle

unset multiplot
