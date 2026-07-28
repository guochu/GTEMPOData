set term epslatex standalone size 8.6cm,7.4cm font ",8" header\
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

set linestyle 1 lw lw lc rgb "#FFCCCC"
set linestyle 2 lw lw lc rgb "#CCE5FF"
set linestyle 3 lw lw lc rgb "#CCFFCC"
set linestyle 4 lw lw lc rgb "#FFFFCC"
set linestyle 5 lw lw lc rgb "#E5CCFF"
set linestyle 6 lw lw lc rgb "#B30000"
set linestyle 7 lw lw lc rgb "#0066CC"
set linestyle 8 lw lw lc rgb "#009900"
set linestyle 9 lw lw lc rgb "#FFCC00"
set linestyle 10 lw lw lc rgb "#660099"

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
h=0.39

# (a)
set size w,h
set origin 0.1,0.53
set label '$\scriptstyle\beta=10$' at graph 0.4,1.07
set xlabel '$\tau$' offset 0,1.2
set ylabel '$\mathcal{G}_{\uparrow\uparrow}(\tau)$' offset 3.5,0
set xrange [0:10]
set xtics 0,2,10 offset 0,0.48
set yrange [-0.5:-0.1]
set ytics -0.5,0.1,0
set mxtics 2
set mytics 2
set label '\small{Iter[1]}' at graph 0.4,1.16
set arrow from graph 0.22,1.175 to graph 0.37,1.175 ls 1 nohead
set arrow from graph 0.22,1.145 to graph 0.37,1.145 ls 6 dt 2 nohead
set label '\small{Iter[3]}' at graph 0.9,1.16
set arrow from graph 0.72,1.175 to graph 0.87,1.175 ls 2 nohead
set arrow from graph 0.72,1.145 to graph 0.87,1.145 ls 7 dt 2 nohead
#set label '\small{Iter[8]}' at graph 0.5,0.4
#set arrow from graph 0.3,0.415 to graph 0.45,0.415 ls 3 nohead
#set arrow from graph 0.3,0.385 to graph 0.45,0.385 ls 8 dt 2 nohead
set label '\small{Iter[8]}' at graph 1.35,1.16
set arrow from graph 1.17,1.175 to graph 1.32,1.175 ls 3 nohead
set arrow from graph 1.17,1.145 to graph 1.32,1.145 ls 8 dt 2 nohead
set label '\small{Iter[10]}' at graph 1.8,1.16
set arrow from graph 1.62,1.175 to graph 1.77,1.175 ls 5 nohead
set arrow from graph 1.62,1.145 to graph 1.77,1.145 ls 10 dt 2 nohead
set label '(a)' at graph -0.27,0.95
set key Left maxrows 2 samplen 2.5 width -1 vertical spacing 0.9 offset 0,0.1
plot "Guu10_iter1_ctqmc.dat" u 1:2 every 2 w l ls 1 notitle,\
     "Guu10_iter3_ctqmc.dat" u 1:2 every 2 w l ls 2 notitle,\
     "Guu10_iter8_ctqmc.dat" u 1:2 every 2 w l ls 3 notitle,\
     "Guu10_iter10_ctqmc.dat" u 1:2 every 2 w l ls 5 notitle,\
     "Guu10_iter1.dat" u 1:2 w l ls 6 dt 2 notitle,\
     "Guu10_iter3.dat" u 1:2 w l ls 7 dt 2 notitle,\
     "Guu10_iter8.dat" u 1:2 w l ls 8 dt 2 notitle,\
     "Guu10_iter10.dat" u 1:2 w l ls 10 dt 2 notitle


# (b)
unset label
unset arrow
set size w,h
set origin 0.59,0.53
set label '$\scriptstyle\beta=100$' at graph 0.4,1.07
set ylabel '$\mathcal{G}_{\uparrow\uparrow}(\tau)$' offset 3.5,0
set yrange [-0.5:0.02]
set xtics 0,20,100 offset 0,0.48
set xrange [0:100]
set label '(b)' at graph -0.27,0.95
plot "Guu100_iter1_ctqmc.dat" u 1:2 every 2 w l ls 1 notitle,\
     "Guu100_iter10_ctqmc.dat" u 1:2 every 2 w l ls 5 notitle,\
     "Guu100_iter1.dat" u 1:2 every 2 w l ls 6 dt 2 notitle,\
     "Guu100_iter10.dat" u 1:2 every 2 w l ls 10 dt 2 notitle,\
     #"Guu100_iter3_ctqmc.dat" u 1:2 every 2 w l ls 2 notitle,\
     "Guu100_iter8_ctqmc.dat" u 1:2 every 2 w l ls 3 notitle,\
     "Guu100_iter9_ctqmc.dat" u 1:2 every 2 w l ls 4 notitle,\
     "Guu100_iter10_ctqmc.dat" u 1:2 every 10 w l ls 5 notitle,\
     "Guu100_iter1.dat" u 1:2 w l ls 6 dt 2 notitle,\
     "Guu100_iter3.dat" u 1:2 w l ls 7 dt 2 notitle,\
     "Guu100_iter8.dat" u 1:2 w l ls 8 dt 2 notitle,\
     "Guu100_iter9.dat" u 1:2 w l ls 9 dt 2 notitle,\
     "Guu100_iter10.dat" u 1:2 w l ls 10 dt 2 notitle

# (c)
unset label
unset arrow
set size w,h
set origin 0.1,0.07
set ylabel '$\mathcal{G}_{\uparrow\downarrow}(\tau)$' offset 3.5,0
set xrange [0:10]
set xtics 0,2,10 offset 0,0.48
set yrange [-0.2:0.2]
set ytics -0.4,0.1,0.4
set label '(c)' at graph -0.27,0.95
plot "Gud10_iter1_ctqmc.dat" u 1:2 every 2 w l ls 1 notitle,\
     "Gud10_iter3_ctqmc.dat" u 1:2 every 2 w l ls 2 notitle,\
     "Gud10_iter8_ctqmc.dat" u 1:2 every 2 w l ls 3 notitle,\
     "Gud10_iter10_ctqmc.dat" u 1:2 every 2 w l ls 5 notitle,\
     "Gud10_iter1.dat" u 1:2 w l ls 6 dt 2 notitle,\
     "Gud10_iter3.dat" u 1:2 w l ls 7 dt 2 notitle,\
     "Gud10_iter8.dat" u 1:2 w l ls 8 dt 2 notitle,\
     "Gud10_iter10.dat" u 1:2 w l ls 10 dt 2 notitle
     

# (d)
set size w,h
set origin 0.59,0.07
unset label
unset arrow
set ylabel '$\mathcal{G}_{\uparrow\downarrow}(\tau)$' offset 3.5,0
set xrange [0:100]
set xtics 0,20,100 offset 0,0.48
set yrange [-0.3:0.3]
set ytics -0.4,0.1,0.4
set label '(d)' at graph -0.27,0.95
plot "Gud100_iter1_ctqmc.dat" u 1:2 every 2 w l ls 1 notitle,\
     "Gud100_iter10_ctqmc.dat" u 1:2 every 2 w l ls 5 notitle,\
     "Gud100_iter1.dat" u 1:2 w l ls 6 dt 2 notitle,\
     "Gud100_iter10.dat" u 1:2 w l ls 10 dt 2 notitle,\
     #"Gud100_iter3_ctqmc.dat" u 1:2 every 10 w l ls 2 notitle,\
     "Gud100_iter8_ctqmc.dat" u 1:2 every 10 w l ls 3 notitle,\
     "Gud100_iter9_ctqmc.dat" u 1:2 every 10 w l ls 4 notitle,\
     "Gud100_iter10_ctqmc.dat" u 1:2 every 10 w l ls 5 notitle,\
     "Gud100_iter1.dat" u 1:2 w l ls 1 dt 2 notitle,\
     "Gud100_iter3.dat" u 1:2 w l ls 2 dt 2 notitle,\
     "Gud100_iter8.dat" u 1:2 w l ls 3 dt 2 notitle,\
     "Gud100_iter9.dat" u 1:2 w l ls 4 dt 2 notitle,\
     "Gud100_iter10.dat" u 1:2 w l ls 5 dt 2 notitle



# Inset
# (a)
unset label
unset arrow
unset xlabel
unset ylabel
set size 0.15,0.13
set origin 0.23,0.62
set mxtics 1
set mytics 1
set yrange [-8:8]
set ytics -8,4,8
set xrange [0:10]
set xtics 0,2,10 offset 0,0.6
set ylabel '\scaleto{\mathrm{error}}{2.5pt}' offset 3,0
set xlabel '\scaleto{\tau}{2.5pt}' offset 0,1.7
set format y '\tiny{%.0f}'
set format x '\tiny{%.0f}'
set label '\scaleto{\times10^{-3}}{4pt}' at graph 0,1.1
plot "Guu10_iter1_diff.dat" u 1:($2*1e3) w l ls 6 lw 2 notitle,\
     "Guu10_iter3_diff.dat" u 1:($2*1e3) w l ls 7 lw 2 notitle,\
     "Guu10_iter8_diff.dat" u 1:($2*1e3) w l ls 8 lw 2 notitle,\
     "Guu10_iter10_diff.dat" u 1:($2*1e3) w l ls 10 lw 2 notitle


# (b)
unset label
unset arrow
unset xlabel
unset ylabel
set size 0.15,0.13
set origin 0.73,0.62
set mxtics 1
set mytics 1
set yrange [-5:5]
set ytics -4,2,4
set xrange [0:100]
set xtics 0,50,100 offset 0,0.6
set ylabel '\scaleto{\mathrm{error}}{2.5pt}' offset 3,0
set xlabel '\scaleto{\tau}{2.5pt}' offset 0,1.7
set format y '\tiny{%.0f}'
set format x '\tiny{%.0f}'
set label '\scaleto{\times10^{-2}}{4pt}' at graph 0,1.1
plot "Guu100_iter1_diff.dat" u 1:($2*1e2) w l ls 6 lw 2 notitle,\
     "Guu100_iter10_diff.dat" u 1:($2*1e2) w l ls 10 lw 2 notitle,\

# # (c)
unset label
unset arrow
unset xlabel
unset ylabel
set size 0.15,0.13
set origin 0.17,0.112
set mxtics 1
set mytics 1
set yrange [-8:8]
set ytics -8,4,8
set xrange [0:10]
set xtics 0,2,10 offset 0,0.6
set ylabel '\scaleto{\mathrm{error}}{2.5pt}' offset 3,0
set xlabel '\scaleto{\tau}{2.5pt}' offset 0,1.7
set format y '\tiny{%.0f}'
set format x '\tiny{%.0f}'
set label '\scaleto{\times10^{-4}}{4pt}' at graph 0,1.1
plot "Gud10_iter1_diff.dat" u 1:($2*1e3) w l ls 6 lw 2 notitle,\
     "Gud10_iter3_diff.dat" u 1:($2*1e3) w l ls 7 lw 2 notitle,\
     "Gud10_iter8_diff.dat" u 1:($2*1e3) w l ls 8 lw 2 notitle,\
     "Gud10_iter10_diff.dat" u 1:($2*1e3) w l ls 10 lw 2 notitle,\


# (d)
unset label
unset arrow
unset xlabel
unset ylabel
set size 0.15,0.13
set origin 0.67,0.112
set mxtics 1
set mytics 1
set yrange [-4:4]
set ytics -8,2,16
set xrange [0:100]
set xtics 0,50,100 offset 0,0.6
set ylabel '\scaleto{\mathrm{error}}{2.5pt}' offset 3,0
set xlabel '\scaleto{\tau}{2.5pt}' offset 0,1.7
set format y '\tiny{%.0f}'
set format x '\tiny{%.0f}'
set label '\scaleto{\times10^{-2}}{4pt}' at graph 0,1.1
plot "Gud100_iter1_diff.dat" u 1:($2*1e2) w l ls 6 lw 2 notitle,\
     "Gud100_iter10_diff.dat" u 1:($2*1e2) w l ls 10 lw 2 notitle,\

# set size 0.1,0.1
# set origin 0.86,0.33
# set mxtics 1
# set mytics 1
# set yrange [4:16]
# set ytics 0,4,16
# set xrange [170:430]
# set xtics 0,200,400 offset 0,0.6
# set ylabel '\scaleto{\mathcal{E}}{4pt}' offset 2.8,0
# set xlabel '\scaleto{\chi}{4pt}' offset 0,1.6
# set format y '\tiny{%.0f}'
# set format x '\tiny{%.0f}'
# plot "Gud100_chi.dat" u 1:($2*1e3) w lp ls 8 pt 4 lw 2 dt 2 notitle,\

