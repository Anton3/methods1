set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 0 ps 1.5
set style line 2 lc rgb '#dd181f' lt 1 lw 2 pt 0 ps 1.5
plot './runge-kutta4-system-4' index 0 with linespoints ls 1, \
     '' index 1 with linespoints ls 2
