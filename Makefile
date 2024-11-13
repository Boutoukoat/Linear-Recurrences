


all: ./lr2_psp_display ./lr2_psp_count \
	./lr3_psp_display ./lr3_psp_count ./lr3_primality ./lr3_sqr_disc


./lr2_psp_display: lr2.cpp lnrcr.h lr2_psp_display.cpp 
	g++ -O3 -o ./lr2_psp_display lr2_psp_display.cpp lr2.cpp -lgmp

./lr2_psp_count: lr2.cpp lnrcr.h lr2_psp_count.cpp 
	g++ -O3 -o ./lr2_psp_count lr2_psp_count.cpp lr2.cpp -lgmp

./lr3_psp_display: lr3.cpp lnrcr.h lr3_psp_display.cpp 
	g++ -O3 -o ./lr3_psp_display lr3_psp_display.cpp lr3.cpp -lgmp

./lr3_psp_count: lr3.cpp lnrcr.h lr3_psp_count.cpp 
	g++ -O3 -o ./lr3_psp_count lr3_psp_count.cpp lr3.cpp -lgmp

./lr3_primality: lr3.cpp lnrcr.h lr3_primality.cpp 
	g++ -O3 -o ./lr3_primality lr3_primality.cpp lr3.cpp -lgmp

./lr3_sqr_disc: lr3.cpp lnrcr.h lr3_sqr_disc.cpp 
	g++ -O3 -o ./lr3_sqr_disc lr3_sqr_disc.cpp lr3.cpp -lgmp

./lr4_psp_display: lr4.cpp lnrcr.h lr4_psp_display.cpp 
	g++ -O3 -o ./lr4_psp_display l43_psp_display.cpp lr4.cpp -lgmp

./lr4_psp_count: lr4.cpp lnrcr.h lr4_psp_count.cpp 
	g++ -O3 -o ./lr4_psp_count lr4_psp_count.cpp lr4.cpp -lgmp

clear:
	rm -f ./lr2_psp_display ./lr2_psp_count \
	./lr3_psp_display ./lr3_psp_count \
	./lr4_psp_display ./lr4_psp_count

tests:
	./lr2_psp_display 1 -2 1000000
	./lr3_psp_display 1 -2 1 1000000
	./lr3_primality 0 1 1 271441 904631 16532714 24658561 27422714 27664033 46672291 102690901 130944133


