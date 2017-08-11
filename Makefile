all:
	g++ -m64 main.cpp libGL.so.1 libglut.so libGLU.so libarmadillo.so

clean:
	rm -rf a.out

