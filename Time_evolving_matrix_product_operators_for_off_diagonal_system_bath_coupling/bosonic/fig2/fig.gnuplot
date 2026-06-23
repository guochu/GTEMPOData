set term epslatex standalone size 8.6cm,4.0cm font ",8" header\
'\usepackage{scalerel}'

lw = 3

set linestyle 1 lw lw lc rgb "#8A2BE2" pt 6 dt 2
set linestyle 2 lw lw lc rgb "#666666"

set style arrow 1 noborder lw 2 size 1,20

# linear fit functions (for log10-transformed data)
fa(x) = Aa*x + Ba   # log10(runtime) vs log10(chi)
fb(x) = Ab*x + Bb   # log10(runtime) vs log10(d)
Aa = 2.3; Ba = -0.25
Ab = 2.7; Bb = 1.3

set out "fig.tex"
set multiplot

# common settings
set format x '\small{%g}'
set lmargin 0
set rmargin 0
set tmargin 0
set bmargin 0

w=0.35
h=0.70

xL=0.13
xR=0.62

y1=0.15

# (a) runtime vs chi
set size w,h
set origin xL,y1
set xlabel '$\log_{10}(\chi)$' offset 0,1.2
set ylabel '$\log_{10}(\mathrm{runtime})$' offset 1.5,0
set format y '\small{%.1f}'
set xrange [0:2.0]
set yrange [0:4.5]
set xtics 0,0.5,2.0 offset 0,0.5
set ytics 0,1,4
set mxtics 2
set mytics 2
set label '(a)' at graph -0.35,0.95
fit fa(x) "../data/fig2/panel00_chi_runtime.csv" u (log10($1)):(log10($2)) every ::1 via Aa,Ba
# fit result: log10(runtime) = 3.17*log10(chi) - 1.82
set arrow 2 from graph 0.08,0.90 to graph 0.23,0.90 nohead ls 2
set label 2 '{\scriptsize $y = 3.17x - 1.82$}' at graph 0.55,0.90
plot "../data/fig2/panel00_chi_runtime.csv" u (log10($1)):(log10($2)) every ::1 w p ls 1 notitle, \
     [0:2] fa(x) w l ls 2 notitle

# (b) runtime vs d
unset label
unset arrow

set size w,h
set origin xR,y1
set xlabel '$\log_{10}(d)$' offset 0,1.2
set ylabel '$\log_{10}(\mathrm{runtime})$' offset 1.5,0
set format y '\small{%.1f}'
set xrange [0:1.0]
set yrange [0:5.0]
set xtics 0,0.2,1.0 offset 0,0.5
set ytics 0,1,5
set mxtics 2
set mytics 2
set label '(b)' at graph -0.35,0.95
fit fb(x) "../data/fig2/panel01_d_runtime.csv" u (log10($1)):(log10($2)) via Ab,Bb
# fit result: log10(runtime) = 3.42*log10(d) + 1.22
set arrow 4 from graph 0.10,0.90 to graph 0.25,0.90 nohead ls 2
set label 4 '{\scriptsize $y = 3.42x + 1.22$}' at graph 0.57,0.90
plot "../data/fig2/panel01_d_runtime.csv" u (log10($1)):(log10($2)) w p ls 1 notitle, \
     [0:1] fb(x) w l ls 2 notitle

unset multiplot
set out
