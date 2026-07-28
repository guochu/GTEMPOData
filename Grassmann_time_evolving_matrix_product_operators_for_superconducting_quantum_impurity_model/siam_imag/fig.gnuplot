set term epslatex standalone size 8.6cm,7cm font ",8" header\
'\usepackage{scalerel}'

lw = 4

#set linestyle 1 lw lw lc rgb "pink"
#set linestyle 2 lw lw lc rgb "red"
# set linestyle 1 lw lw lc rgb "gray"
# set linestyle 2 lw lw lc rgb "black"
# set linestyle 3 lw lw lc rgb "#90EE90"
# set linestyle 4 lw lw lc rgb "#009E73"
# set linestyle 5 lw lw lc rgb "#1E90FF"
# set linestyle 6 lw lw lc rgb "blue"
# set linestyle 1 lw lw lc rgb "pink"
# set linestyle 2 lw lw lc rgb "red"
# set linestyle 3 lw lw lc rgb "#90EE90"
# set linestyle 4 lw lw lc rgb "#009E73"
# set linestyle 5 lw lw lc rgb "#1E90FF"
# set linestyle 6 lw lw lc rgb "blue"
# set linestyle 7 lw lw lc rgb "violet"
# set linestyle 8 lw lw lc rgb "dark-violet"
# set linestyle 9 lw lw lc rgb "#DEB887"
# set linestyle 10 lw lw lc rgb "#CD853F"
set linestyle 1 lw lw lc rgb "#CCE5FF"
set linestyle 2 lw lw lc rgb "#0066CC"
#set style arrow 1 noborder lw 2 size 1,20

set out "fig.tex"
set multiplot


# main Figure
set xrange [-0.05:1.05]
set yrange [0:0.55]
set ytics 0,0.1,0.5 offset 0.5,0
set format x '\small{%g}'
set format y '\small{%.1f}'
set xlabel '$\Gamma t$' offset 0,1.2
set ylabel '$\mathcal{G}(\tau)$' offset 3,0
#set key samplen 1 maxcolumns 1
set lmargin 0
set rmargin 0
set tmargin 0
set bmargin 0

w=0.38
h=0.4

# (a)
set size w,h
set origin 0.1,0.54
set label '$\scriptstyle\beta=10$' at graph 0.4,1.08
set xlabel '$\tau$' offset 0,1.2
set ylabel '$\mathcal{G}_{\uparrow\uparrow}(\tau)$' offset 3.5,0
set xrange [0:10]
set xtics 0,2,10 offset 0,0.48
set yrange [-0.5:0]
set ytics -0.5,0.1,0
set mxtics 2
set mytics 2
set label '(a)' at graph -0.27,0.95
set key Left maxrows 2 samplen 2.5 width -1 vertical spacing 0.9 offset 0,0.1
plot "Guu10_ctqmc.dat" u 1:2 w l ls 1 t '\tiny{CTQMC}',\
     "Guu10.dat" u 1:2 w l ls 2 dt 2 t '\tiny{GTEMPO}'

# (b)
unset label
unset arrow
set size w,h
set origin 0.59,0.54
set label '$\scriptstyle\beta=100$' at graph 0.4,1.08
set ylabel '$\mathcal{G}_{\uparrow\uparrow}(\tau)$' offset 3.5,0
set xtics 0,20,100 offset 0,0.48
set xrange [0:100]
set yrange [-0.5:0.02]
set label '(b)' at graph -0.27,0.95
plot "Guu100_ctqmc.dat" u 1:2 w l ls 1 notitle,\
     "Guu100.dat" u 1:2 w l ls 2 dt 2 notitle

# (c)
unset label
set size w,h
set origin 0.1,0.07
set ylabel '$\mathcal{G}_{\uparrow\downarrow}(\tau)$' offset 3.5,0
set xrange [0:10]
set xtics 0,2,10 offset 0,0.48
set yrange [-0.4:0.4]
set ytics -0.4,0.2,0.4
set label '(c)' at graph -0.27,0.95
plot "Gud10_ctqmc.dat" u 1:2 w l ls 1 notitle,\
     "Gud10.dat" u 1:2 w l ls 2 dt 2 notitle

# (d)
set size w,h
set origin 0.59,0.07
unset label
set ylabel '$\mathcal{G}_{\uparrow\downarrow}(\tau)$' offset 3.5,0
set xrange [0:100]
set xtics 0,20,100 offset 0,0.48
set yrange [-0.4:0.4]
set ytics -0.4,0.2,0.4
set label '(d)' at graph -0.27,0.95
plot "Gud100_ctqmc.dat" u 1:2 w l ls 1 notitle,\
     "Gud100.dat" u 1:2 w l ls 2 dt 2 notitle
     


# Inset
# (a)
unset label
unset arrow
unset xlabel
unset ylabel
set size 0.1,0.1
set origin 0.18,0.6
set mxtics 1
set mytics 1
set yrange [0:7]
set ytics 0,2,8
set xrange [0.03:0.22]
set xtics 0,0.1,0.2 offset 0,0.6
set ylabel '\scaleto{\mathcal{E}}{4pt}' offset 2.8,0
set xlabel '\scaleto{\delta\tau}{4pt}' offset 0,1.4
set format y '\tiny{%.0f}'
set format x '\tiny{%.1f}'
set label '$\scriptscriptstyle \times10^{-3}$' at graph 0,1.15
plot "Guu10_dtau.dat" u 1:($2*1e3) w lp ls 2 pt 4 lw 2 dt 2 notitle,\

set size 0.1,0.1
set origin 0.33,0.6
set mxtics 1
set mytics 1
set yrange [0:10]
set ytics 0,4,8
set xrange [20:430]
set xtics 0,200,400 offset 0,0.6
set ylabel '\scaleto{\mathcal{E}}{4pt}' offset 2.8,0
set xlabel '\scaleto{\chi}{4pt}' offset 0,1.6
set format y '\tiny{%.0f}'
set format x '\tiny{%.0f}'
plot "Guu10_chi.dat" u 1:($2*1e3) w lp ls 2 pt 4 lw 2 dt 2 notitle,\


# (b)
unset label
unset arrow
unset xlabel
unset ylabel
set size 0.1,0.1
set origin 0.68,0.6
set mxtics 1
set mytics 1
set yrange [0:16]
set ytics 0,4,16
set xrange [0.1:1.1]
set xtics 0,0.5,1 offset 0,0.6
set ylabel '\scaleto{\mathcal{E}}{4pt}' offset 2.8,0
set xlabel '\scaleto{\delta\tau}{4pt}' offset 0,1.4
set format y '\tiny{%.0f}'
set format x '\tiny{%.1f}'
set label '$\scriptscriptstyle \times10^{-3}$' at graph 0,1.15
plot "Guu100_dtau.dat" u 1:($2*1e3) w lp ls 2 pt 4 lw 2 dt 2 notitle,\

set size 0.1,0.1
set origin 0.86,0.6
set mxtics 1
set mytics 1
set yrange [2:5]
set ytics 0,1,5
set xrange [170:430]
set xtics 0,200,400 offset 0,0.6
set ylabel '\scaleto{\mathcal{E}}{4pt}' offset 2.8,0
set xlabel '\scaleto{\chi}{4pt}' offset 0,1.6
set format y '\tiny{%.0f}'
set format x '\tiny{%.0f}'
plot "Guu100_chi.dat" u 1:($2*1e3) w lp ls 2 pt 4 lw 2 dt 2 notitle,\


# # (c)
unset label
unset arrow
unset xlabel
unset ylabel
set size 0.1,0.1
set origin 0.15,0.14
set mxtics 1
set mytics 1
set yrange [2:8]
set ytics 0,2,8
set xrange [0.03:0.22]
set xtics 0,0.1,0.2 offset 0,0.6
set ylabel '\scaleto{\mathcal{E}}{4pt}' offset 2.8,0
set xlabel '\scaleto{\delta\tau}{4pt}' offset 0,1.4
set format y '\tiny{%.0f}'
set format x '\tiny{%.1f}'
set label '$\scriptscriptstyle \times10^{-4}$' at graph 0,1.15
plot "Gud10_dtau.dat" u 1:($2*1e4) w lp ls 2 pt 4 lw 2 dt 2 notitle,\

unset label
set size 0.1,0.1
set origin 0.33,0.32
set mxtics 1
set mytics 1
set yrange [0:10]
set ytics 0,4,8
set xrange [20:430]
set xtics 0,200,400 offset 0,0.6
set ylabel '\scaleto{\mathcal{E}}{4pt}' offset 2.8,0
set xlabel '\scaleto{\chi}{4pt}' offset 0,1.6
set format y '\tiny{%.0f}'
set format x '\tiny{%.0f}'
set label '$\scriptscriptstyle \times10^{-3}$' at graph 0,1.15
plot "Gud10_chi.dat" u 1:($2*1e3) w lp ls 2 pt 4 lw 2 dt 2 notitle,\

# (d)
unset label
unset arrow
unset xlabel
unset ylabel
set size 0.1,0.1
set origin 0.68,0.13
set mxtics 1
set mytics 1
set yrange [0:8]
set ytics 0,4,16
set xrange [0.1:1.1]
set xtics 0,0.5,1 offset 0,0.6
set ylabel '\scaleto{\mathcal{E}}{4pt}' offset 2.8,0
set xlabel '\scaleto{\delta\tau}{4pt}' offset 0,1.4
set format y '\tiny{%.0f}'
set format x '\tiny{%.1f}'
set label '$\scriptscriptstyle \times10^{-3}$' at graph 0,1.15
plot "Gud100_dtau.dat" u 1:($2*1e3) w lp ls 2 pt 4 lw 2 dt 2 notitle,\

set size 0.1,0.1
set origin 0.86,0.33
set mxtics 1
set mytics 1
set yrange [4:16]
set ytics 0,4,16
set xrange [170:430]
set xtics 0,200,400 offset 0,0.6
set ylabel '\scaleto{\mathcal{E}}{4pt}' offset 2.8,0
set xlabel '\scaleto{\chi}{4pt}' offset 0,1.6
set format y '\tiny{%.0f}'
set format x '\tiny{%.0f}'
plot "Gud100_chi.dat" u 1:($2*1e3) w lp ls 2 pt 4 lw 2 dt 2 notitle,\

