#Robert Grumbine
NAME ?= test

$(NAME) : $(NAME).C makefile
	g++ -Wall -O2 -DLINUX -std=c++17 $< -o $@ ../*.a -I ../include
	#g++ -Wall -O2 -DLINUX -std=c++20 $< -o $@ ../*.a -I ../include # issues with == on mvector<mvector<T>>[i]
	./$@
