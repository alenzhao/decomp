.SUFFIXES: .o .c
CC = gcc 
#DEBUG = -g -Wall
DEBUG = -Wall
FLAGS = $(DEBUG) -O3
LIBS = -lm

CORE = ovc.o parse_option.o

CORE4 = grn_cmp_hcc.o parse_option.o

EXEC = ovc

EXEC4 = grn_cmp_hcc

all: $(EXEC) $(EXEC4)

.c.o: 
	$(CC) $(FLAGS) -c $< -o $@

core : $(CORE)

$(EXEC) : $(CORE)
	$(CC) $(FLAGS) $(CORE) -o $(EXEC) $(LIBS)
	@echo 'Made '

$(EXEC4) : $(CORE4)
	$(CC) $(FLAGS) $(CORE4) -o $(EXEC4) -lm
	@echo 'Made '

clean:
	rm -f *.o *.bak *.*~
	rm -f $(EXEC)
	@echo 'Made'
