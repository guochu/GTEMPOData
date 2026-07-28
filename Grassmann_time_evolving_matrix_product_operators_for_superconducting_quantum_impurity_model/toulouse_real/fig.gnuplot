set term epslatex standalone size 8.6cm,10cm font ",8" header\
'\usepackage{scalerel}'

lw = 4


set linestyle 1 lw lw lc rgb "pink"
set linestyle 2 lw lw lc rgb "red"
set linestyle 3 lw lw lc rgb "#90EE90"
set linestyle 4 lw lw lc rgb "#009E73"
set linestyle 5 lw lw lc rgb "#1E90FF"
set linestyle 6 lw lw lc rgb "blue"
set linestyle 7 lw lw lc rgb "violet"
set linestyle 8 lw lw lc rgb "dark-violet"
set linestyle 9 lw lw lc rgb "#DEB887"
set linestyle 10 lw lw lc rgb "#CD853F"
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
set linestyle 1 lw lw lc rgb "#FFA0A0"
set linestyle 2 lw lw lc rgb "#FF0000"
set linestyle 3 lw lw lc rgb "#E5CCFF"
set linestyle 4 lw lw lc rgb "#660099"
set style arrow 1 noborder lw 2 size 1,20

set out "fig.tex"
set multiplot


# main Figure
set xrange [0:5]
set yrange [0:1]
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
h=0.27

# (a)
set size w,h
set origin 0.1,0.7
set xlabel '$t$' offset 0,1.2
set ylabel '$\mathrm{Im}\,[G_{\uparrow\bar{\uparrow}\bullet}(t)]$' offset 3.5,0
set xrange [0:10]
set xtics 0,2,10 offset 0,0.48
set yrange [-1:1]
set ytics -1,0.5,1
set mxtics 2
set mytics 2
set label '(a)' at graph -0.27,0.95
set key bottom Left reverse offset 5,0.3
plot "ed_uu.dat" u 1:2 w l ls 1 t '\small{ED}',\
     "uu.dat" u 1:2 w l ls 2 dt 2 t '$\scriptstyle\chi=250,\delta t=0.05$'


# (b)
unset label
unset arrow
set size w,h
set origin 0.6,0.7
set xrange [0:10]
set ylabel '$\mathrm{Re}\,[G_{\uparrow\bullet\downarrow}(t)]$' offset 3.5,0
set label '(b)' at graph -0.27,0.95
plot "ed_ud.dat" u 1:2 w l ls 3 notitle,\
     "ud.dat" u 1:2 w l ls 4 dt 2 notitle



# (c)
unset label
unset arrow
set size w,h
set origin 0.1,0.375
set xrange [0.03:0.22]
set xtics 0.05,0.05,0.2
set yrange [1:9]
set ytics 0,2,10
set format x "\\small{%.2f}"
set key bottom
set label '(c)' at graph -0.27,0.95
set xlabel '$\delta t$'
set ylabel '$\mathcal{E}$' offset 2.5,0
set label '$\scriptstyle{\times10^{-2}}$' at graph 0,1.04
plot "uu_dt.dat" u 1:($2*1e2) w lp ls 2 dt 2 pt 4 ps 1.5 notitle,\

# (d)
set size w,h
set origin 0.6,0.375
unset label
unset arrow
set yrange [1:9]
set ytics 2,2,10
set label '(d)' at graph -0.27,0.95
set xlabel '$\delta t$'
set ylabel '$\mathcal{E}$' offset 2.5,0
set label '$\scriptstyle{\times10^{-2}}$' at graph 0,1.04
plot "ud_dt.dat" u 1:($2*1e2) w lp ls 4 dt 2 pt 4 ps 1.5 notitle,\

# (e)
unset label
unset arrow
set size w,h
set origin 0.1,0.05
set yrange [3.5:6.5]
set xrange [40:260]
unset format x
set xtics 50,50,250
set ytics 3,0.5,7
#set label '$\scriptstyle\epsilon_d=5,d=1,\alpha=1$' at graph 0.25,1.08
set key Left reverse left top maxrows 5 samplen 3 width 0 vertical spacing 0.9 offset 0,-0.2
set ylabel '$\mathcal{E}$' offset 2.5,0
set xlabel '$\chi$'
set label '(e)' at graph -0.27,0.95
set label '$\scriptstyle{\times10^{-2}}$' at graph 0,1.04
plot "uu_chi.dat" u 1:($2*1e2) w lp ls 2 dt 2 pt 4 ps 1.5 notitle


# (f)
set size w,h
set origin 0.6,0.05
unset label
unset arrow
set yrange [3.5:6.5]
#set label '$\scriptstyle\epsilon_d=5,d=3,\alpha=1$' at graph 0.25,1.08
set label '(f)' at graph -0.27,0.95
set label '$\scriptstyle{\times10^{-2}}$' at graph 0,1.04
plot "ud_chi.dat" u 1:($2*1e2) w lp ls 4 dt 2 pt 4 ps 1.5 notitle,\


# Inset
# (a)
# unset label
# unset arrow
# unset xlabel
# unset ylabel
# set size 0.18,0.1
# set origin 0.21,0.755
# set mxtics 1
# set mytics 1
# set yrange [0:7]
# set xrange [0:0.2]
# set xtics 0,0.1,0.2
# set ytics 0,2,7
# set ylabel '\scaleto{\mathcal{E}}{6pt}' offset 6,0
# set xlabel '\scaleto{\delta t}{5pt}' offset 0,1.2
# set format x '\scaleto{%g}{4pt}'
# set format y '\scaleto{%.0f}{4pt}'
# set label '$\scriptstyle\times10^{-2}$' at graph 0,1.1
# plot "keldysh-error.dat" u 1:($2*10**2) w lp ls 4 pt 6 lw 2 notitle,\

# # (b)
# unset label
# unset arrow
# set size 0.18,0.1
# set origin 0.7,0.755
# set mxtics 1
# set mytics 1
# set format x '\scaleto{%g}{4pt}'
# set format y '\scaleto{%.0f}{4pt}'
# set label '$\scriptstyle\times10^{-2}$' at graph 0,1.1
# plot "kadanoff-error.dat" u 1:($2*10**2) w lp ls 4 pt 6 lw 2 notitle,\

# # (c)
# unset label
# unset arrow
# unset ylabel
# set size 0.18,0.1
# set origin 0.2,0.51
# set yrange [0:6]
# set ytics 0,2,6
# set ylabel '\scaleto{\mathcal{E}}{6pt}' offset 6,0
# set format x '\scaleto{%g}{4pt}'
# set format y '\scaleto{%.0f}{4pt}'
# set label '$\scriptstyle\times10^{-2}$' at graph 0,1.1
# plot "keldysh-error.dat" u 1:($3*10**2) w lp ls 4 pt 6 lw 2 notitle,\

# # (d)
# unset label
# unset arrow
# set size 0.18,0.1
# set origin 0.7,0.52
# set yrange [0:9]
# set ytics 0,4,10
# set format x '\scaleto{%g}{4pt}'
# set format y '\scaleto{%.0f}{4pt}'
# set label '$\scriptstyle\times10^{-2}$' at graph 0,1.1
# plot "kadanoff-error.dat" u 1:($3*10**2) w lp ls 4 pt 6 lw 2 notitle,\

# # (e)
# unset label
# unset arrow
# unset ylabel
# set size 0.18,0.1
# set origin 0.2,0.11
# set yrange [0:6]
# set ytics 0,2,6
# set ylabel '\scaleto{\mathcal{E}}{6pt}' offset 6,0
# set format x '\scaleto{%g}{4pt}'
# set format y '\scaleto{%.0f}{4pt}'
# set label '$\scriptstyle\times10^{-2}$' at graph 0,1.1
# plot "keldysh-error.dat" u 1:($4*10**2) w lp ls 4 pt 6 lw 2 notitle,\

# # (f)
# unset label
# unset arrow
# set size 0.18,0.1
# set origin 0.7,0.11
# set yrange [0:2]
# set ytics 0,1,10
# set format x '\scaleto{%g}{4pt}'
# set format y '\scaleto{%.0f}{4pt}'
# set label '$\scriptstyle\times10^{-2}$' at graph 0,1.1
# plot "kadanoff-error.dat" u 1:($4*10**2) w lp ls 4 pt 6 lw 2 notitle,\
