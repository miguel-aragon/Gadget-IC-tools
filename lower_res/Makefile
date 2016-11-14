#
#  Change the compiler if neccesary


CC       =  gcc            # sets the C-compiler

OPTIMIZE =   -O2  #-Wall 

MATH_LIB =   -lm  #-Wall 

OPTIONS =  $(OPTIMIZE) $(MATH_LIB) $(OPT1) $(OPT2) $(OPT3) $(OPT4) $(OPT5) $(OPT6)\
           $(OPT7) $(OPT8) $(OPT9) 

EXEC   = lower_resolution

OBJS   = lower_resolution.o io.o

INCL   = lower_resolution.h io.h

CFLAGS =   $(OPTIONS)

$(EXEC): $(OBJS) 
	$(CC) $(OBJS) $(LIBS)   -o  $(EXEC)  

$(OBJS): $(INCL) 

.PHONY : clean
clean:
	rm -f $(OBJS) $(EXEC)

