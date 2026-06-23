set term epslatex standalone size 8.6cm,6.0cm font ",8" header\
'\usepackage{scalerel}'

lw = 4

set linestyle 1 lw lw lc rgb "pink"
set linestyle 2 lw lw lc rgb "red"
set linestyle 3 lw lw lc rgb "#009E73"
set linestyle 4 lw lw lc rgb "#009E73"
set linestyle 5 lw lw lc rgb "#8A2BE2"
set linestyle 6 lw lw lc rgb "#666666"

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

w=0.35
h=0.36

y1=0.55
y2=0.10

xL=0.13
xR=0.62

# (a) sigma_z dynamics
set size w,h
set origin xL,y1
set xlabel '$t$' offset 0,1.2
set ylabel '$\langle\hat{\sigma}_z\rangle$' offset 2.8,0
set xrange [0:5]
set yrange [-0.4:0.05]
set xtics 0,1,5 offset 0,0.5
set ytics -0.4,0.1,0
set mxtics 2
set mytics 2
set label '\small{TEMPO}' at graph 0.50,0.85 right
set arrow from graph 0.53,0.85 to graph 0.73,0.85 ls 2 dt 2 nohead
set label '\small{ED}' at graph 0.50,0.72 right
set arrow from graph 0.53,0.72 to graph 0.73,0.72 ls 1 nohead
set label '(a)' at graph -0.33,0.95
unset key
plot "data/panel00_ed.csv" u 1:2 w l ls 1 notitle,\
     "data/panel00_tempo_chi30.csv" u 1:2 w l ls 2 dt 2 notitle

# (b) error vs chi
unset label
unset arrow
unset key
set size w,h
set origin xR,y1
set xlabel '$\chi$' offset 0,1.2
set ylabel '$\mathcal{E}[\hat{\sigma}_z]$' offset 1.9,0
set xrange [0:42]
set yrange [0:11]
set xtics 0,10,40 offset 0,0.5
set ytics 0,3,11
set format y '\small{%.1f}'
set mxtics 2
set mytics 2
set label '(b)' at graph -0.34,0.95
set label '$\scriptstyle{\times10^{-3}}$' at graph 0,1.05
plot "data/panel01_err_vs_chi.csv" u 1:($2*1000) w lp ls 5 pt 6 lw lw notitle

# (c) error vs k, fit with a*exp(-b*x)+c
unset label
unset arrow
set xrange [2:8]
set yrange [0:4]
set size w,h
set origin xL,y2
set xlabel '$m$' offset 0,1.2
set ylabel '$\mathcal{E}[\hat{\sigma}_z]$' offset 1.9,0
set xtics 2,1,8 offset 0,0.5
set ytics 0,1,4
set format y '\small{%.1f}'
set mxtics 2
set mytics 2
set label '(c)' at graph -0.33,0.95
set label '$\scriptstyle{\times10^{-3}}$' at graph 0,1.05
plot "data/panel10_err_vs_k.csv" u 1:($2*1000) w lp ls 5 pt 6 lw lw notitle

# (d) error vs dt (log-log)
unset label
unset arrow
set xrange [-1.7:-0.5]
set yrange [-4.0:-2.5]
f(x) = a*x + b
set fit logfile "/dev/null"
fit f(x) "data/panel11_err_vs_dt.csv" u (log10($1)):(log10($2)) via a,b
# fit result: log10(E) = 0.96*log10(dt) - 2.24
set size w,h
set origin xR,y2
set xlabel '$\log_{10}(\delta t)$' offset 0,1.2
set ylabel '$\log_{10}(\mathcal{E}[\hat{\sigma}_z])$' offset 2.8,0
set xtics -1.5,0.5,-0.5 offset 0,0.5
set ytics -4.0,0.5,-2.5
set format y '\small{%.1f}'
set format x '\small{%.1f}'
set mxtics 2
set mytics 2
set label '(d)' at graph -0.34,0.95
set arrow from graph 0.05,0.90 to graph 0.20,0.90 nohead ls 6
set label sprintf('{\scriptsize $y = %.2fx %+.2f$}', a, b) at graph 0.55,0.90
plot "data/panel11_err_vs_dt.csv" u (log10($1)):(log10($2)) w lp ls 5 pt 6 lw lw notitle,\
     f(x) w l ls 6 notitle

unset multiplot
