
CC 			:= gcc
CCFLAGS		:= -Wall -std=c11
CPATH		:= ./src/
HPATH		:= -I ./include/
NAME		:= fp

doclear: clear clear do

clear:
	clear

do:
	$(CC) $(CCFLAGS) $(HPATH) -o $(NAME) $(CPATH)fp.c $(CPATH)fp_helper.c \
	$(CPATH)montgomory_redc.c $(CPATH)fp2.c $(CPATH)sike_encodings_helper.c \
	$(CPATH)sike_encodings.c $(CPATH)montgomery_curve.c $(CPATH)point_arithmetic.c \
	$(CPATH)point_arithmetic_helper.c $(CPATH)params.c $(CPATH)isogeny_strat.c \
	$(CPATH)sidh.c fp_test.c -lpthread

test:
	python3 Extern_Implementation/tester.py



run :
	./$(NAME)