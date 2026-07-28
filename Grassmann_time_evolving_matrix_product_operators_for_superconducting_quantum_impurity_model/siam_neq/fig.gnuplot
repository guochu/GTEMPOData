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


#set linestyle 1 lw lw lc rgb "#FFC8C8"
#set linestyle 2 lw lw lc rgb "#FFA0A0"
#set linestyle 3 lw lw lc rgb "#FF7F7F"
#set linestyle 4 lw lw lc rgb "#FF5050"
#set linestyle 5 lw lw lc rgb "#FF0000"
set linestyle 1 lw lw lc rgb "#FFF2E6"
set linestyle 2 lw lw lc rgb "#FFD8B1"
set linestyle 3 lw lw lc rgb "#FFA500"
set linestyle 4 lw lw lc rgb "#E67300"
set linestyle 5 lw lw lc rgb "#CC5500"
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
#set label '$\scriptstyle\beta=10$' at graph 0.4,1.08
set xlabel '$t$' offset 0,1.2
set ylabel '$P_0(t)$' offset 2.5,0
set xrange [0:10]
set xtics 0,2,10 offset 0,0.48
set yrange [-0.03:1]
set ytics -1,0.2,1
set mxtics 2
set mytics 2
# set label '$\chi=200$' at graph 0.4,0.9
# set arrow from graph 0.75,0.915 to graph 0.9,0.915 ls 1 nohead
# set arrow from graph 0.75,0.885 to graph 0.9,0.885 ls 5 dt 2 nohead
# set label '$\chi=250$' at graph 0.4,0.8
# set arrow from graph 0.75,0.815 to graph 0.9,0.815 ls 2 nohead
# set arrow from graph 0.75,0.785 to graph 0.9,0.785 ls 6 dt 2 nohead
# set label '$\chi=300$' at graph 0.4,0.7
# set arrow from graph 0.75,0.715 to graph 0.9,0.715 ls 3 nohead
# set arrow from graph 0.75,0.685 to graph 0.9,0.685 ls 7 dt 2 nohead
# set label '$\chi=400$' at graph 0.4,0.6
# set arrow from graph 0.75,0.615 to graph 0.9,0.615 ls 4 nohead
# set arrow from graph 0.75,0.585 to graph 0.9,0.585 ls 8 dt 2 nohead
set label '(a)' at graph -0.27,0.95
set key Left maxrows 1 samplen 2.5 width -2 vertical spacing 0.9 offset 27.4,1.6
plot "vacuum_Nt200_chi150.dat" u 1:2 w l ls 1 t '$\scriptstyle\chi=150$',\
     "vacuum_Nt200_chi200.dat" u 1:2 w l ls 2 t '$\scriptstyle\chi=200$',\
     "vacuum_Nt200_chi250.dat" u 1:2 w l ls 3 t '$\scriptstyle\chi=250$',\
     "vacuum_Nt200_chi300.dat" u 1:2 w l ls 4 t '$\scriptstyle\chi=300$',\
     "vacuum_Nt200_chi400.dat" u 1:2 w l ls 5 t '$\scriptstyle\chi=400$',\
#     "vacuum_Nt100_chi200.dat" u 1:2 w l ls 5 dt 2 notitle,\
     "vacuum_Nt100_chi250.dat" u 1:2 w l ls 6 dt 2 notitle,\
     "vacuum_Nt100_chi300.dat" u 1:2 w l ls 7 dt 2 notitle,\
     "vacuum_Nt100_chi400.dat" u 1:2 w l ls 8 dt 2 notitle

# (b)
unset label
unset arrow
set size w,h
set origin 0.6,0.54
#set label '$\scriptstyle\beta=100$' at graph 0.4,1.08
set ylabel '$P_{\uparrow}(t)$' offset 2.5,0
set xtics 0,2,100 offset 0,0.48
set xrange [0:10]
set label '(b)' at graph -0.27,0.95
plot "nup_Nt200_chi150.dat" u 1:2 w l ls 1 notitle,\
     "nup_Nt200_chi200.dat" u 1:2 w l ls 2 notitle,\
     "nup_Nt200_chi250.dat" u 1:2 w l ls 3 notitle,\
     "nup_Nt200_chi300.dat" u 1:2 w l ls 4 notitle,\
     "nup_Nt200_chi400.dat" u 1:2 w l ls 5 notitle,\
     #"nup_Nt100_chi200.dat" u 1:2 w l ls 5 dt 2 notitle,\
     "nup_Nt100_chi250.dat" u 1:2 w l ls 6 dt 2 notitle,\
     "nup_Nt100_chi300.dat" u 1:2 w l ls 7 dt 2 notitle,\
     "nup_Nt100_chi400.dat" u 1:2 w l ls 8 dt 2 notitle,\
     


# (c)
unset label
set size w,h
set origin 0.1,0.07
set ylabel '$P_{\downarrow}(t)$' offset 2.5,0
set xrange [0:10]
set xtics 0,2,10 offset 0,0.48
set yrange [0:1]
set ytics 0,0.2,1
set label '(c)' at graph -0.27,0.95
plot "ndown_Nt200_chi150.dat" u 1:2 w l ls 1 notitle,\
     "ndown_Nt200_chi200.dat" u 1:2 w l ls 2 notitle,\
     "ndown_Nt200_chi250.dat" u 1:2 w l ls 3 notitle,\
     "ndown_Nt200_chi300.dat" u 1:2 w l ls 4 notitle,\
     "ndown_Nt200_chi400.dat" u 1:2 w l ls 5 notitle,\
     #"ndown_Nt100_chi200.dat" u 1:2 w l ls 5 dt 2 notitle,\
     "ndown_Nt100_chi250.dat" u 1:2 w l ls 6 dt 2 notitle,\
     "ndown_Nt100_chi300.dat" u 1:2 w l ls 7 dt 2 notitle,\
     "ndown_Nt100_chi400.dat" u 1:2 w l ls 8 dt 2 notitle,\

# (d)
set size w,h
set origin 0.6,0.07
unset label
set ylabel '$P_{\uparrow\downarrow}(t)$' offset 2.5,0
set xrange [0:10]
set xtics 0,2,10 offset 0,0.48
set yrange [0:1]
set ytics 0,0.2,1
set label '(d)' at graph -0.27,0.95
plot "nn_Nt200_chi150.dat" u 1:2 w l ls 1 notitle,\
     "nn_Nt200_chi200.dat" u 1:2 w l ls 2 notitle,\
     "nn_Nt200_chi250.dat" u 1:2 w l ls 3 notitle,\
     "nn_Nt200_chi300.dat" u 1:2 w l ls 4 notitle,\
     "nn_Nt200_chi400.dat" u 1:2 w l ls 5 notitle,\
     #"nn_Nt100_chi200.dat" u 1:2 w l ls 5 dt 2 notitle,\
     "nn_Nt100_chi250.dat" u 1:2 w l ls 6 dt 2 notitle,\
     "nn_Nt100_chi300.dat" u 1:2 w l ls 7 dt 2 notitle,\
     "nn_Nt100_chi400.dat" u 1:2 w l ls 8 dt 2 notitle



# Inset
# (a)
unset label
unset arrow
unset xlabel
unset ylabel
set size 0.18,0.15
set origin 0.19,0.75
set mxtics 1
set mytics 1
set yrange [-3:2]
set ytics -4,2,8
set xrange [0:10]
set xtics 0,2,10 offset 0,0.6
set ylabel '\scaleto{\mathrm{error}}{2.5pt}' offset 3.2,0
set xlabel '\scaleto{t}{4pt}' offset 0,1.4
set format y '\tiny{%.0f}'
set format x '\tiny{%.0f}'
set label '$\scriptstyle\chi=400$' at graph 0.4,0.3
set label '\scaleto{\times10^{-2}}{4pt}' at graph 0,1.1
plot "vacuum_diff_chi400.dat" u 1:($2*1e2) w l ls 5 lw 1 notitle,\


# (b)
unset label
unset arrow
unset xlabel
unset ylabel
set size 0.18,0.15
set origin 0.69,0.75
set mxtics 1
set mytics 1
set yrange [-3:2]
#set ylabel '\scaleto{\mathcal{E}}{4pt}' offset 2.8,0
#set xlabel '\scaleto{\delta\tau}{4pt}' offset 0,1.4
#set format y '\tiny{%.0f}'
#set format x '\tiny{%.1f}'
set label '\scaleto{\times10^{-2}}{4pt}' at graph 0,1.1
set ylabel '\scaleto{\mathrm{error}}{2.5pt}' offset 3.2,0
set xlabel '\scaleto{t}{4pt}' offset 0,1.4
plot "nup_diff_chi400.dat" u 1:($2*1e2) w l ls 5 lw 1 notitle,\

# (c)
unset label
unset arrow
unset xlabel
unset ylabel
set size 0.18,0.15
set origin 0.19,0.28
set mxtics 1
set mytics 1
set yrange [-2:2]
#set ylabel '\scaleto{\mathcal{E}}{4pt}' offset 2.8,0
#set xlabel '\scaleto{\delta\tau}{4pt}' offset 0,1.4
#set format y '\tiny{%.0f}'
#set format x '\tiny{%.1f}'
set label '\scaleto{\times10^{-2}}{4pt}' at graph 0,1.1
set ylabel '\scaleto{\mathrm{error}}{2.5pt}' offset 3.2,0
set xlabel '\scaleto{t}{4pt}' offset 0,1.4
plot "ndown_diff_chi400.dat" u 1:($2*1e2) w l ls 5 lw 1 notitle,\

# (d)
unset label
unset arrow
unset xlabel
unset ylabel
set size 0.18,0.15
set origin 0.69,0.28
set mxtics 1
set mytics 1
set yrange [-2:4]
set ytics -10,2,8
#set ylabel '\scaleto{\mathcal{E}}{4pt}' offset 2.8,0
#set xlabel '\scaleto{\delta\tau}{4pt}' offset 0,1.4
#set format y '\tiny{%.0f}'
#set format x '\tiny{%.1f}'
set label '\scaleto{\times10^{-2}}{4pt}' at graph 0,1.1
set ylabel '\scaleto{\mathrm{error}}{2.5pt}' offset 3.2,0
set xlabel '\scaleto{t}{4pt}' offset 0,1.4
plot "nn_diff_chi400.dat" u 1:($2*1e2) w l ls 5 lw 1 notitle,\
