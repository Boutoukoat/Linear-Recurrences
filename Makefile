


all: ./lr2 ./lr3 ./lr4

./lr2: lr2.cpp
	g++ -O3 -o ./lr2 lr2.cpp -lgmp

./lr3: lr3.cpp
	g++ -O3 -o ./lr3 lr3.cpp -lgmp

./lr4: lr4.cpp
	g++ -O3 -o ./lr4 lr4.cpp -lgmp

clear:
	rm -f ./lr2 ./lr3 ./lr4

tests:
	./lr2 1 -2 1000000
	./lr3 1 -2 1 1000000
	./lr4 1 -2 1 2 1000000


