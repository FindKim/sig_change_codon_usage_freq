all: StopLight

StopLight: main.o Parse_pvalue_ortholog_map.o StopLight.o TTest.o
	g++ main.o Parse_pvalue_ortholog_map.o StopLight.o TTest.o -o StopLight

main.o: main.cpp
	g++ -c main.cpp

Parse_pvalue_ortholog_map.o: Parse_pvalue_ortholog_map.cpp Parse_pvalue_ortholog_map.h
	g++ -c Parse_pvalue_ortholog_map.cpp

StopLight.o: StopLight.cpp StopLight.h
	g++ -c StopLight.cpp

TTest.o: TTest.cpp TTest.h
	g++ -c TTest.cpp

clean:
	rm -f *.o *.pyc *~ StopLight
