objects = main.o matrice.o grad.o fonctions.o

programme: $(objects)
cc -Wall -std=c99 -o programme $(objects) -lm

main.o : main.c
	gcc -Wall -std=c99 -c main.c

matrice.o : matrice.c
	gcc -Wall -std=c99 -c matrice.c

grad.o : grad.c
	gcc -Wall -std=c99 -c grad.c

fonctions.o : fonctions.c
	gcc -Wall -std=c99 -c fonctions.c


clean :
	rm $(objects)