#!/usr/bin/gnuplot

w = 17
h = 8
mu = 0.0
kRange = 32
w_min = -8.0
w_max =  8.0

Nw_data = 1500
w_min_data = -13.0
w_max_data =  13.0
w_range_data = w_max_data - w_min_data


#set terminal x11 size 600,800


set yrange [w_min:w_max]
set ytics ('-16' -16,'-12' -12,'-8' -8,'-4' -4, '0' 0.0, '4' +4.0, '8' +8.0, '12' +12.0,'-16' +16)

plot 'output/dos_uhb.dat' u ($2):($1-mu) with filledcurves notitle lc '#bb6666',\
     'output/dos_p.dat'   u ($2):($1-mu) with filledcurves notitle lc '#ffaaaa',\
     'output/dos_lhb.dat' u ($2):($1-mu) with filledcurves notitle lc '#6666bb',\

#MAX=GPVAL_X_MAX
#MIN=GPVAL_X_MIN
#set arrow from MIN,0.0 to MAX,0.0 nohead lw 1 lc 'black' front

#plot 'output/dos_uhb.dat' u ($2):($1-mu) with filledcurves notitle lc '#bb6666',\
#     'output/dos_p.dat'   u ($2):($1-mu) with filledcurves notitle lc '#ffaaaa',\
#     'output/dos_lhb.dat' u ($2):($1-mu) with filledcurves notitle lc '#6666bb',\

pause -1
