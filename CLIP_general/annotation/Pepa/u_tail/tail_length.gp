#!/usr/bin/env gnuplot
set term postscript eps size 7,5 color 'Arial,14'
#set term pdf color enhanced size 7,5
set datafile separator '\t'
set pointintervalbox 3

set style line 1 lw 2 lt 1 lc rgb "#5b4d41" pt 7 pi -1 ps 1.5
set style line 2 lw 2 lt 1 lc rgb "#a7bf42" pt 7 pi -1 ps 1.5
set style line 3 lw 2 lt 1 lc rgb "#524e86" pt 7 pi -1 ps 1.5
set style line 4 lw 2 lt 1 lc rgb "#c22e13" pt 7 pi -1 ps 1.5
set style line 5 lw 2 lt 1 lc rgb "#515c73" pt 7 pi -1 ps 1.5
set style line 6 lw 2 lt 1 lc rgb "#899096" pt 7 pi -1 ps 1.5
set style line 7 lw 2 lt 1 lc rgb "#000000" pt 7 pi -1 ps 1.5
set style line 8 lw 2 lt 1 lc rgb "#333333" pt 7 pi -1 ps 1.5
set style line 9 lw 2 lt 1 lc rgb "#cccccc" pt 7 pi -1 ps 1.5

set xrange[1:25]
#set yrange[0:20]
set output "tail_length_histogram.eps"
set xtics 1; set xlabel 'Tail length [nt]'; set format x '%.0f'
set ytics; set format y '%g%%'; set ylabel 'Percentage of all reads'
set title "Histogram of uridylated tail length in various RNA categories"

plot "results_histogram_u_tail.txt" index 0 using 1:2 with lines ls 1 ti columnhead,\
     "" index 0 using 1:3 with lines ls 2 ti columnhead,\
     "" index 0 using 1:4 with lines ls 3 ti columnhead,\
     "" index 0 using 1:5 with lines ls 4 ti columnhead,\
     "" index 0 using 1:6 with lines ls 5 ti columnhead,\
     "" index 0 using 1:7 with lines ls 6 ti columnhead,\
     "" index 0 using 1:8 with lines ls 7 ti columnhead,\
     "" index 0 using 1:9 with lines ls 8 ti columnhead,\
     "" index 0 using 1:10 with lines ls 9 ti columnhead

# with linespoints
