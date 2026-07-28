set term epslatex standalone size 8.6cm,7cm font ",8" header\
'\usepackage{scalerel}'

lw = 3



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
set linestyle 10 lw lw lc rgb "#CD853F"

set linestyle 1 lw lw lc rgb "#FFA0A0"
set linestyle 2 lw lw lc rgb "#FF0000"
set linestyle 3 lw lw lc rgb "#E5CCFF"
set linestyle 4 lw lw lc rgb "#660099"

set out "fig.tex"
set multiplot


# main Figure
set xrange [-0.05:1.05]
set yrange [0:0.55]
set ytics 0,0.1,0.5 offset 0.5,0
set format x '\small{%g}'
set format y '\small{%.1f}'
set xlabel '$\tau$' offset 0,1.2
set ylabel '$G_{\uparrow\uparrow}(\tau)$' offset 3.5,0
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
set xrange [0:10]
set xtics 0,2,10 offset 0,0.48
set yrange [-1:1]
set ytics -1,0.5,1
set mxtics 2
set mytics 2
set ylabel '$G_{\uparrow\bar{\uparrow}\bullet}(t)$' offset 3.5,0
set xlabel '$t$'
set label '(a)' at graph -0.27,0.95
# set label 'real' at graph 0.3,0.1
# set arrow from graph 0.05,0.115 to graph 0.25,0.115 ls 2 nohead
# set arrow from graph 0.05,0.085 to graph 0.25,0.085 ls 8 dt 2 nohead
# set label 'imag' at graph 0.78,0.1
# set arrow from graph 0.53,0.115 to graph 0.73,0.115 ls 6 nohead
# set arrow from graph 0.53,0.085 to graph 0.73,0.085 ls 4 dt 2 nohead
plot "ed_uu.dat" u 1:2 w l ls 1 notitle,\
     "uu.dat" u 1:2 w l ls 2 dt 2 notitle,\
     "ed_uu.dat" u 1:3 w l ls 3 notitle,\
     "uu.dat" u 1:3 w l ls 4 dt 2 notitle

     

# (b)
unset label
unset arrow
set size w,h
set origin 0.6,0.57
set label '(b)' at graph -0.27,0.95
set yrange [-1:1]
set label '\small{real part}' at graph 0.4,0.25
set arrow from graph 0.15,0.26 to graph 0.35,0.26 ls 1 nohead
set arrow from graph 0.15,0.23 to graph 0.35,0.23 ls 2 dt 2 nohead
set label '\small{imaginary part}' at graph 0.4,0.15
set arrow from graph 0.15,0.16 to graph 0.35,0.16 ls 3 nohead
set arrow from graph 0.15,0.13 to graph 0.35,0.13 ls 4 dt 2 nohead
set ylabel '$G_{\uparrow\bullet\downarrow}(t)$' offset 3.5,0
plot "ed_ud.dat" u 1:2 w l ls 1 notitle,\
     "ud.dat" u 1:2 w l ls 2 dt 2 notitle,\
     "ed_ud.dat" u 1:3 w l ls 3 notitle,\
     "ud.dat" u 1:3 w l ls 4 dt 2 notitle
     
# (c)
unset label
unset arrow
set size w,h
set origin 0.1,0.07
set key bottom
set yrange [1:11]
set ytics 0,2,10
set xrange [0.03:0.42]
set xtics 0.1,0.1,0.4
set format x '\small{%.2f}'
set label '(c)' at graph -0.27,0.95
set ylabel '$\mathcal{E}$' offset 2.5,0
set xlabel '$\delta t$'
set label '$\scriptstyle\times10^{-2}$' at graph 0,1.04
plot "uu_dt.dat" u 1:($2*1e2) w lp ls 10 dt 2 pt 4 ps 1.5 notitle,\
     

# (d)
set size w,h
set origin 0.6,0.07
unset label
set label '(d)' at graph -0.27,0.95
set yrange [1:8]
set label '$\scriptstyle\times10^{-2}$' at graph 0,1.04
plot "ud_dt.dat" u 1:($2*1e2) w lp ls 10 dt 2 pt 4 ps 1.5 notitle,\
