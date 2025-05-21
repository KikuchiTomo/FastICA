all:
	g++ -std=c++11 src/main.cpp -o main -I/usr/local/include/eigen3

run:
	./main

plot:
	gnuplot plot.gpl

clean:
	-rm *.dat
	-rm *.png
	-rm main
