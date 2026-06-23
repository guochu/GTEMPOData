set term epslatex standalone size 8.6cm,13.1cm font ",8" header\
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

w=0.36
h=0.1924

y1=0.7435
y2=0.5134
y3=0.2783
y4=0.0432

xL=0.11
xR=0.60

# (a) alpha=0.1, sigma_z dynamics
set size w,h
set origin xL,y1
set xlabel '$t$' offset 0,1.2
set ylabel '$\langle\hat{\sigma}_z\rangle$' offset 3.0,0
set xrange [0:10]
set yrange [-1:1]
set xtics 0,2,10 offset 0,0.5
set ytics -1,0.5,1
set mxtics 2
set mytics 2
set label '\small{TEMPO}' at graph 0.04,1.07
set arrow from graph 0.43,1.07 to graph 0.63,1.07 ls 2 dt 2 nohead
set label '\small{ED}' at graph 0.8,1.07
set arrow from graph 1.0,1.07 to graph 1.2,1.07 ls 1 nohead
set label '\small{ED, Rabi type}' at graph 1.4,1.07
set arrow from graph 2.0,1.07 to graph 2.2,1.07 ls 4 nohead
set label '(a)' at graph -0.30,0.95
set label '$\lambda^2=0.1$' at graph 0.35,1.22
unset key
plot "data/alpha0.1/panel00_jc_ed.csv" u 1:2 w l ls 1 notitle,\
     "data/alpha0.1/panel00_rabi_ed.csv" u 1:2 w l ls 4 notitle,\
     "data/alpha0.1/panel00_tempo_chi30.csv" u 1:2 w l ls 2 dt 2 notitle

# (b) alpha=0.5, sigma_z dynamics
unset label
unset arrow
unset key
set size w,h
set origin xR,y1
set xrange [0:10]
set yrange [-1:1]
set xtics 0,2,10 offset 0,0.5
set ytics -1,0.5,1
set mxtics 2
set mytics 2
set label '(b)' at graph -0.33,0.95
set label '$\lambda^2=0.5$' at graph 0.35,1.22
plot "data/alpha0.5/panel01_jc_ed.csv" u 1:2 w l ls 1 notitle,\
     "data/alpha0.5/panel01_rabi_ed.csv" u 1:2 w l ls 4 notitle,\
     "data/alpha0.5/panel01_tempo_chi30.csv" u 1:2 w l ls 2 dt 2 notitle

# (c) alpha=0.1, error vs chi
unset label
unset key
set size w,h
set origin xL,y2
set xlabel '$\chi$' offset 0,1.2
set ylabel '$\mathcal{E}[\hat{\sigma}_z]$' offset 2.1,0
set xrange [2:22]
set yrange [0:10]
set xtics 5,5,20 offset 0,0.5
set ytics 0,3,9
set format y '\small{%.1f}'
set mxtics 2
set mytics 2
set label '(c)' at graph -0.30,0.95
set label '$\scriptstyle{\times10^{-4}}$' at graph 0,1.05
plot "data/alpha0.1/panel10_err_vs_chi.csv" u 1:($2*10000) w lp ls 5 pt 6 lw lw notitle

# (d) alpha=0.5, error vs chi
unset label
set size w,h
set origin xR,y2
set xrange [2:22]
set yrange [0:3]
set xtics 5,5,20 offset 0,0.5
set ytics 0,1,3
set format y '\small{%.1f}'
set mxtics 2
set mytics 2
set label '(d)' at graph -0.33,0.95
set label '$\scriptstyle{\times10^{-2}}$' at graph 0,1.05
plot "data/alpha0.5/panel11_err_vs_chi.csv" u 1:($2*100) w lp ls 5 pt 6 lw lw notitle

# (e) alpha=0.1, error vs k
unset label
set size w,h
set origin xL,y3
set xlabel '$m$' offset 0,1.2
set ylabel '$\mathcal{E}[\hat{\sigma}_z]$' offset 2.1,0
set xrange [2:10]
set yrange [0:10]
set xtics 2,2,10 offset 0,0.5
set ytics 0,3,9
set format y '\small{%.1f}'
set mxtics 2
set mytics 2
set label '(e)' at graph -0.30,0.95
set label '$\scriptstyle{\times10^{-4}}$' at graph 0,1.05
plot "data/alpha0.1/panel20_err_vs_k.csv" u 1:($2*10000) w lp ls 5 pt 6 lw lw notitle

# (f) alpha=0.5, error vs k
unset label
set size w,h
set origin xR,y3
set xrange [2:10]
set yrange [0:10]
set xtics 2,2,10 offset 0,0.5
set ytics 0,3,9
set format y '\small{%.1f}'
set mxtics 2
set mytics 2
set label '(f)' at graph -0.33,0.95
set label '$\scriptstyle{\times10^{-3}}$' at graph 0,1.05
plot "data/alpha0.5/panel21_err_vs_k.csv" u 1:($2*1000) w lp ls 5 pt 6 lw lw notitle

# linear fits for (g) and (h) using all points
f1(x) = a1*x + b1
fit [*:*] f1(x) "data/alpha0.1/panel30_err_vs_dt.csv" using (log10($1)):(log10($2)) via a1,b1
f2(x) = a2*x + b2
fit [*:*] f2(x) "data/alpha0.5/panel31_err_vs_dt.csv" using (log10($1)):(log10($2)) via a2,b2
# fit result: log10(E) = 1.58*log10(dt) - 1.30 (alpha=0.1)
# fit result: log10(E) = 1.62*log10(dt) - 0.67 (alpha=0.5)

# (g) alpha=0.1, error vs dt (log-log)
unset label
unset arrow
set size w,h
set origin xL,y4
set xlabel '$\log_{10}(\delta t)$' offset 0,1.2
set ylabel '$\log_{10}(\mathcal{E}[\hat{\sigma}_z])$' offset 1.5,0
set xrange [-2.0:-0.5]
set yrange [-4.0:-2.0]
set xtics -2.0,0.5,-0.5 offset 0,0.5
set ytics -4,1,-2
set format x '\small{%.1f}'
set format y '\small{%g}'
set mxtics 2
set mytics 2
set label '(g)' at graph -0.30,0.95
set arrow from graph 0.05,0.90 to graph 0.18,0.90 nohead ls 6
set label sprintf('{\scriptsize $y = %.2fx %+.2f$}', a1, b1) at graph 0.50,0.90
plot "data/alpha0.1/panel30_err_vs_dt.csv" u (log10($1)):(log10($2)) w lp ls 5 pt 6 lw lw notitle,\
     f1(x) w l ls 6 notitle

# (h) alpha=0.5, error vs dt (log-log)
unset label
unset arrow
set size w,h
set origin xR,y4
set xlabel '$\log_{10}(\delta t)$' offset 0,1.2
set ylabel '$\log_{10}(\mathcal{E}[\hat{\sigma}_z])$' offset 1.2,0
set xrange [-2.0:-0.5]
set yrange [-3.5:-1.5]
set xtics -2.0,0.5,-0.5 offset 0,0.5
set ytics -3,1,-1
set format x '\small{%.1f}'
set format y '\small{%g}'
set mxtics 2
set mytics 2
set label '(h)' at graph -0.33,0.95
set arrow from graph 0.05,0.90 to graph 0.18,0.90 nohead ls 6
set label sprintf('{\scriptsize $y = %.2fx %+.2f$}', a2, b2) at graph 0.50,0.90
plot "data/alpha0.5/panel31_err_vs_dt.csv" u (log10($1)):(log10($2)) w lp ls 5 pt 6 lw lw notitle,\
     f2(x) w l ls 6 notitle

unset multiplot
