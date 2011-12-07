test.o: simulator.o data.o
	g++ simulator.o data.o main.cpp -o transitionMatrix

data.o: communityData.cpp communityData.h simulator.o
	g++ -c communityData.cpp -o data.o

simulator.o : communitySimulator.cpp communitySimulator.h 
	g++ -c communitySimulator.cpp -o simulator.o

clean:
	rm data.o simulator.o