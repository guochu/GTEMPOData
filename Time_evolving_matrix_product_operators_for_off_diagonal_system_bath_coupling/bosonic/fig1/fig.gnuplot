set term epslatex standalone size 18cm,5.8cm font ",8" header\
'\usepackage{scalerel}'

lw = 3

set linestyle 1 lw lw lc rgb "#000080" dt 2
set linestyle 2 lw lw lc rgb "#FF0000" dt 2
set linestyle 3 lw lw lc rgb "#009E73" dt 2
set linestyle 4 lw lw lc rgb "#CC6600" dt 2
set linestyle 5 lw lw lc rgb "#CC79A7" dt 2
set linestyle 6 lw lw lc rgb "#56B4E9" dt 2
set linestyle 7 lw lw lc rgb "#8B4513" dt 2
set linestyle 8 lw lw lc rgb "#006400" dt 2

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

w=0.14
h=0.362

# common key
set key top right maxrows 6 width -6 spacing 1.6 font ",5" samplen 6 opaque

# =====================================================================
# Row 1: data/beta1 - subplots (a)-(e), beta=1
# =====================================================================

# --- (a) <n(t)> for different chi ---
set size w,h
set origin 0.06,0.571
set xlabel '$\tau$' offset 0,1.0
set ylabel '$\mathcal{G}(\tau)$' offset 2.0,0
set xrange [0:1.0]
set yrange [0.5:1.8]
set xtics 0,0.2,1.0 offset 0,0.5
set ytics 0.5,0.3,1.8
set mxtics 2
set mytics 2
set label '(a)' at graph -0.38,0.95
plot "../data/fig1/beta1/panel00_chi.csv" u 1:2 w l ls 1 t '$\chi=10$',\
     "../data/fig1/beta1/panel00_chi.csv" u 1:3 w l ls 2 t '$\chi=20$',\
     "../data/fig1/beta1/panel00_chi.csv" u 1:4 w l ls 3 t '$\chi=30$',\
     "../data/fig1/beta1/panel00_chi.csv" u 1:5 w l ls 4 t '$\chi=40$',\
     "../data/fig1/beta1/panel00_chi.csv" u 1:6 w l ls 5 t '$\chi=50$',\
     "../data/fig1/beta1/panel00_chi.csv" u 1:7 w l ls 6 t '$\chi=60$'

# --- (b) <n(t)> for different d ---
unset label
set size w,h
set origin 0.255,0.571
set xlabel '$\tau$' offset 0,1.0
set ylabel '$\mathcal{G}(\tau)$' offset 2.0,0
set xrange [0:1.0]
set yrange [0.5:1.8]
set xtics 0,0.2,1.0 offset 0,0.5
set ytics 0.5,0.3,1.8
set mxtics 2
set mytics 2
set label '(b)' at graph -0.38,0.95
plot "../data/fig1/beta1/panel01_d.csv" u 1:2 w l ls 1 t '$d=4$',\
     "../data/fig1/beta1/panel01_d.csv" u 1:3 w l ls 2 t '$d=5$',\
     "../data/fig1/beta1/panel01_d.csv" u 1:4 w l ls 3 t '$d=6$',\
     "../data/fig1/beta1/panel01_d.csv" u 1:5 w l ls 4 t '$d=7$',\
     "../data/fig1/beta1/panel01_d.csv" u 1:6 w l ls 5 t '$d=8$'

# --- (c) <n(t)> for different dtau ---
unset label
set size w,h
set origin 0.45,0.571
set xlabel '$\tau$' offset 0,1.0
set ylabel '$\mathcal{G}(\tau)$' offset 2.0,0
set xrange [0:1.0]
set yrange [0.5:1.8]
set xtics 0,0.2,1.0 offset 0,0.5
set ytics 0.5,0.3,1.8
set mxtics 2
set mytics 2
set label '(c)' at graph -0.38,0.95
plot "../data/fig1/beta1/panel02_dtau0.025.csv" u 1:2 w l ls 1 t '$\delta\tau=0.025$',\
     "../data/fig1/beta1/panel02_dtau0.05.csv" u 1:2 w l ls 2 t '$\delta\tau=0.05$',\
     "../data/fig1/beta1/panel02_dtau0.1.csv" u 1:2 w l ls 3 t '$\delta\tau=0.1$',\
     "../data/fig1/beta1/panel02_dtau0.2.csv" u 1:2 w l ls 4 t '$\delta\tau=0.2$'

# --- (d) <n(t)> for different k ---
unset label
set size w,h
set origin 0.645,0.571
set xlabel '$\tau$' offset 0,1.0
set ylabel '$\mathcal{G}(\tau)$' offset 2.0,0
set xrange [0:1.0]
set yrange [0.5:1.8]
set xtics 0,0.2,1.0 offset 0,0.5
set ytics 0.5,0.3,1.8
set mxtics 2
set mytics 2
set label '(d)' at graph -0.38,0.95
plot "../data/fig1/beta1/panel03_k.csv" u 1:2 w l ls 1 t '$m=2$',\
     "../data/fig1/beta1/panel03_k.csv" u 1:4 w l ls 2 t '$m=4$',\
     "../data/fig1/beta1/panel03_k.csv" u 1:6 w l ls 3 t '$m=6$',\
     "../data/fig1/beta1/panel03_k.csv" u 1:7 w l ls 4 t '$m=7$',\
     "../data/fig1/beta1/panel03_k.csv" u 1:8 w l ls 5 t '$m=8$'

# --- (e) <n(t)> for different n ---
unset label
set size w,h
set origin 0.84,0.571
set xlabel '$\tau$' offset 0,1.0
set ylabel '$\mathcal{G}(\tau)$' offset 2.0,0
set xrange [0:1.0]
set yrange [0.5:1.8]
set xtics 0,0.2,1.0 offset 0,0.5
set ytics 0.5,0.3,1.8
set mxtics 2
set mytics 2
set label '(e)' at graph -0.38,0.95
plot "../data/fig1/beta1/panel04_n.csv" u 1:2 w l ls 1 t '$n=1$',\
     "../data/fig1/beta1/panel04_n.csv" u 1:3 w l ls 2 t '$n=2$',\
     "../data/fig1/beta1/panel04_n.csv" u 1:4 w l ls 3 t '$n=4$',\
     "../data/fig1/beta1/panel04_n.csv" u 1:5 w l ls 4 t '$n=6$',\
     "../data/fig1/beta1/panel04_n.csv" u 1:6 w l ls 5 t '$n=8$'

# =====================================================================
# Row 2: data/beta5 - subplots (f)-(j), beta=5
# =====================================================================

unset key

# --- (f) <n(t)> for different chi ---
unset label
set size w,h
set origin 0.06,0.10
set xlabel '$\tau$' offset 0,1.0
set ylabel '$\mathcal{G}(\tau)$' offset 2.0,0
set xrange [0:5]
set yrange [0:1.6]
set xtics 0,1,5 offset 0,0.5
set ytics 0,0.4,1.6
set mxtics 2
set mytics 2
set label '(f)' at graph -0.38,0.95
plot "../data/fig1/beta5/panel10_chi.csv" u 1:2 w l ls 1 t '$\chi=10$',\
     "../data/fig1/beta5/panel10_chi.csv" u 1:3 w l ls 2 t '$\chi=20$',\
     "../data/fig1/beta5/panel10_chi.csv" u 1:4 w l ls 3 t '$\chi=30$',\
     "../data/fig1/beta5/panel10_chi.csv" u 1:5 w l ls 4 t '$\chi=40$',\
     "../data/fig1/beta5/panel10_chi.csv" u 1:6 w l ls 5 t '$\chi=50$',\
     "../data/fig1/beta5/panel10_chi.csv" u 1:7 w l ls 6 t '$\chi=60$'

# --- (g) <n(t)> for different d ---
unset label
set size w,h
set origin 0.255,0.10
set xlabel '$\tau$' offset 0,1.0
set ylabel '$\mathcal{G}(\tau)$' offset 2.0,0
set xrange [0:5]
set yrange [0:1.6]
set xtics 0,1,5 offset 0,0.5
set ytics 0,0.4,1.6
set mxtics 2
set mytics 2
set label '(g)' at graph -0.38,0.95
plot "../data/fig1/beta5/panel11_d.csv" u 1:2 w l ls 1 t '$d=4$',\
     "../data/fig1/beta5/panel11_d.csv" u 1:3 w l ls 2 t '$d=5$',\
     "../data/fig1/beta5/panel11_d.csv" u 1:4 w l ls 3 t '$d=6$',\
     "../data/fig1/beta5/panel11_d.csv" u 1:5 w l ls 4 t '$d=7$',\
     "../data/fig1/beta5/panel11_d.csv" u 1:6 w l ls 5 t '$d=8$'

# --- (h) <n(t)> for different dtau ---
unset label
set size w,h
set origin 0.45,0.10
set xlabel '$\tau$' offset 0,1.0
set ylabel '$\mathcal{G}(\tau)$' offset 2.0,0
set xrange [0:5]
set yrange [0:1.6]
set xtics 0,1,5 offset 0,0.5
set ytics 0,0.4,1.6
set mxtics 2
set mytics 2
set label '(h)' at graph -0.38,0.95
plot "../data/fig1/beta5/panel12_dtau0.025.csv" u 1:2 w l ls 1 t '$\delta\tau=0.025$',\
     "../data/fig1/beta5/panel12_dtau0.05.csv" u 1:2 w l ls 2 t '$\delta\tau=0.05$',\
     "../data/fig1/beta5/panel12_dtau0.1.csv" u 1:2 w l ls 3 t '$\delta\tau=0.1$',\
     "../data/fig1/beta5/panel12_dtau0.2.csv" u 1:2 w l ls 4 t '$\delta\tau=0.2$'

# --- (i) <n(t)> for different k ---
unset label
set size w,h
set origin 0.645,0.10
set xlabel '$\tau$' offset 0,1.0
set ylabel '$\mathcal{G}(\tau)$' offset 2.0,0
set xrange [0:5]
set yrange [0:1.6]
set xtics 0,1,5 offset 0,0.5
set ytics 0,0.4,1.6
set mxtics 2
set mytics 2
set label '(i)' at graph -0.38,0.95
plot "../data/fig1/beta5/panel13_k.csv" u 1:2 w l ls 1 t '$k=2$',\
     "../data/fig1/beta5/panel13_k.csv" u 1:4 w l ls 2 t '$k=4$',\
     "../data/fig1/beta5/panel13_k.csv" u 1:6 w l ls 3 t '$k=6$',\
     "../data/fig1/beta5/panel13_k.csv" u 1:7 w l ls 4 t '$k=7$',\
     "../data/fig1/beta5/panel13_k.csv" u 1:8 w l ls 5 t '$k=8$'

# --- (j) <n(t)> for different n ---
unset label
set size w,h
set origin 0.84,0.10
set xlabel '$\tau$' offset 0,1.0
set ylabel '$\mathcal{G}(\tau)$' offset 2.0,0
set xrange [0:5]
set yrange [0:1.6]
set xtics 0,1,5 offset 0,0.5
set ytics 0,0.4,1.6
set mxtics 2
set mytics 2
set label '(j)' at graph -0.38,0.95
plot "../data/fig1/beta5/panel14_n.csv" u 1:2 w l ls 1 t '$n=1$',\
     "../data/fig1/beta5/panel14_n.csv" u 1:3 w l ls 2 t '$n=2$',\
     "../data/fig1/beta5/panel14_n.csv" u 1:4 w l ls 3 t '$n=4$',\
     "../data/fig1/beta5/panel14_n.csv" u 1:5 w l ls 4 t '$n=6$',\
     "../data/fig1/beta5/panel14_n.csv" u 1:6 w l ls 5 t '$n=8$'

unset multiplot
set out
