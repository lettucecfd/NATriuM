for i in {0..499}; do gnuplot -e "set terminal jpeg; set title 't = $i'; set view 0,0,1; set xlabel 'x'; set ylabel 'x'; unset ztics; splot 't_$i.dat' with pm3d" > t_$i.jpg; done; 
avconv -f image2 -i t_%d.jpg video.mpg
