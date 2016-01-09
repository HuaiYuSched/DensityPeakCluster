set autoscale
set xtic auto
set ytic auto 
set title "CLUSTERING BY FAST SEARCH AND FIND OF DENSITY PEAKS"
set xlabel "X"
set ylabel "Y"
set term jpeg
set output strftime('%F_%H-%M-%S.jpeg', time(0))
set palette model RGB defined (0 "red",1 "blue", 2 "green", 3 "black")
plot 'Clusterfile' using 1:2:3 notitle with points pt 2 palette
