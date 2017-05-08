INCLUDE = ./

COMPILATEUR = g++

OPTIONS = -g -lm

SOURCES = main.cpp

OBJECTS = main.o matrix.o fichiers.o pred.o dct.o

EXECUTABLE = main

$(EXECUTABLE): $(OBJECTS)
	$(COMPILATEUR) $(OPTIONS) $(OBJECTS)  -o $(EXECUTABLE)

main.o: main.cpp
	$(COMPILATEUR) $(OPTIONS) -c main.cpp

