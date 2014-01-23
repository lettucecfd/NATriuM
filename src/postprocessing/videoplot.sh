for (( i=1; i<=$1; i++ )); do gnuplot -e "set terminal jpeg; set title 't = $i'; set view 90,0,1; set zlabel 'z'; set ylabel 'x'; unset ytics; splot 't_$i.dat' with pm3d" > t_$i.jpg; done; 
avconv -f image2 -i t_%d.jpg video.mpg
