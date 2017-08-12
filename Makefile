all:
	#g++ -m64 main.cpp libGL.so.1 libglut.so libGLU.so libarmadillo.so
	#g++ -g3 -O0 -m64 main.cpp -lGL -lglut -lGLU -larmadillo
	g++ -g0 -O3 -m64 main.cpp -lGL -lglut -lGLU -larmadillo

clean:
	rm -rf a.out

