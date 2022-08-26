
CC 			:= gcc
CCFLAGS		:= -Wall -std=c11
CPATH		:= ./src/
HPATH		:= -I ./include/
NAME		:= fp

doclear: clear do

clear:
	clear

do:
	$(CC) $(CCFLAGS) $(HPATH) -o $(NAME) $(CPATH)fp.c $(CPATH)fp_helper.c \
	$(CPATH)montgomory_redc.c $(CPATH)fp2.c fp_test.c




run :
	./$(NAME)