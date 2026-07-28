set term epslatex standalone size 8.6cm,7cm font ",8" header\
'\usepackage{scalerel}'

lw = 4



# set linestyle 1 lw lw lc rgb "#000080"
# set linestyle 2 lw lw lc rgb "#DDDDFF"
# set linestyle 3 lw lw lc rgb "#AAAAFF"
# set linestyle 4 lw lw lc rgb "#5555FF"
# set linestyle 5 lw lw lc rgb "#0000FF"

# set linestyle 6 lw lw lc rgb "#800000"
# set linestyle 7 lw lw lc rgb "#FFDDDD"
# set linestyle 8 lw lw lc rgb "#FFAAAA"
# set linestyle 9 lw lw lc rgb "#FF5555"
# set linestyle 10 lw lw lc rgb "#FF0000"

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
set linestyle 3 lw lw lc rgb "#CCFFCC"
set linestyle 4 lw lw lc rgb "#009900"

set out "fig.tex"
set multiplot


# main Figure
set xrange [-0.05:1.05]
set yrange [0:0.55]
set ytics 0,0.1,0.5 offset 0.5,0
set format x '\small{%g}'
set format y '\small{%.1f}'
set xlabel '$\tau$' offset 0,1.2
set ylabel '$\mathcal{G}_{\uparrow\uparrow}(\tau)$' offset 3.5,0
set key Left maxrows 4 samplen 3 width -7 vertical spacing 0.9 offset 0,0.1
#set key samplen 1 maxcolumns 1
set lmargin 0
set rmargin 0
set tmargin 0
set bmargin 0

w=0.38
h=0.4

# (a)
set size w,h
set origin 0.1,0.57
set xrange [0:5]
set xtics 0,1,5 offset 0,0.48
set yrange [-1:0]
set ytics -1,0.2,1
set mxtics 2
set mytics 2
set ylabel '$\mathcal{G}_{\uparrow\uparrow}(\tau)$' offset 3.5,0
set label '(a)' at graph -0.27,0.95
set label '$U=+1$' at graph 0.5,0.3
set arrow from graph 0.3,0.315 to graph 0.45,0.315 ls 1 nohead
set arrow from graph 0.3,0.285 to graph 0.45,0.285 ls 2 dt 2 nohead
set label '$U=-1$' at graph 0.5,0.2
set arrow from graph 0.3,0.215 to graph 0.45,0.215 ls 3 nohead
set arrow from graph 0.3,0.185 to graph 0.45,0.185 ls 4 dt 2 nohead
plot "ed_uup.dat" u 1:2 w l ls 1 notitle,\
     "uup.dat" u 1:2 w l ls 2 dt 2 notitle,\
     "ed_uum.dat" u 1:2 w l ls 3 notitle,\
     "uum.dat" u 1:2 w l ls 4 dt 2 notitle

     

# (b)
unset label
unset arrow
set size w,h
set origin 0.6,0.57
set label '(b)' at graph -0.27,0.95
set yrange [-0.3:0.3]
set ylabel '$\mathcal{G}_{\uparrow\downarrow}(\tau)$' offset 3.5,0
plot "ed_udp.dat" u 1:2 w l ls 1 notitle,\
     "udp.dat" u 1:2 w l ls 2 dt 2 notitle,\
     "ed_udm.dat" u 1:2 w l ls 3 notitle,\
     "udm.dat" u 1:2 w l ls 4 dt 2 notitle
     
# (c)
unset label
set size w,h
set origin 0.1,0.07
set key bottom
set yrange [0:3]
set ytics 0,0.5,3
set xrange [0:0.22]
set xtics 0,0.05,0.2
set format x '\small{%.2f}'
set label '(c)' at graph -0.27,0.95
set ylabel '$\mathcal{E}$' offset 2.5,0
set xlabel '$\delta\tau$'
set label '$\scriptstyle\times10^{-3}$' at graph 0,1.04
plot "uup_dt.dat" u 1:($2*1e3) w lp ls 2 dt 2 pt 4 ps 1.5 notitle,\
     "uum_dt.dat" u 1:($2*1e3) w lp ls 4 dt 2 pt 6 ps 1.5 notitle
     

# (d)
set size w,h
set origin 0.6,0.07
unset label
set label '(d)' at graph -0.27,0.95
set label '$\scriptstyle\times10^{-3}$' at graph 0,1.04
plot "udp_dt.dat" u 1:($2*1e3) w lp ls 2 dt 2 pt 4 ps 1.5 notitle,\
     "udm_dt.dat" u 1:($2*1e3) w lp ls 4 dt 2 pt 6 ps 1.5 notitle
