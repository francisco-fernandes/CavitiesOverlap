all: CavitiesOverlap.c
	gcc -Wall -Wextra -O3 -o CavitiesOverlap CavitiesOverlap.c -lm
clean:
	rm CavitiesOverlap
