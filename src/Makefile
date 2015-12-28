NAME 	:= poisson
OBDIR	:= ../OBJ/
RUNDIR	:= ../RUN/

OBJECTS =  $(OBDIR)output.o   $(OBDIR)poisson.o $(OBDIR)parallel.o

HEADER 	:= poisson.h

.DEFAULT:

INCPATH	=	-I ~/ 

LIBPATH	= 	-L ~/

LIBS	=	-lm 

CXX	:=	mpic++
CXXFLAGS:= 	-g  $(INCPATH) -Wall -O3
#FLAGS	:=	-g -pg $(INCPATH) -Wall -fopenmp -O3
#LFLAGS	:=	$(LIBPATH) $(LIBS) -lm -lgsl -lblas
LFLAGS	:=	$(LIBPATH) $(LIBS) 

####################################################################

$(NAME): $(OBJECTS) $(RUNDIR)
	$(CXX) -o $(NAME) $(OBJECTS) $(FLAGS) $(LFLAGS)
	@mv $(NAME) $(RUNDIR) 
#	@cp iwrite $(RUNDIR)

$(RUNDIR): 
	@test -d $(RUNDIR) || mkdir $(RUNDIR)

$(OBDIR)poisson.o: poisson.cpp $(HEADER)
	@test -d $(OBDIR) || mkdir $(OBDIR)        
	$(CXX) -c $(FLAGS) -o $(OBDIR)poisson.o poisson.cpp

$(OBDIR)output.o: output.cpp $(HEADER)
	@test -d $(OBDIR) || mkdir $(OBDIR)        
	$(CXX) -c $(FLAGS) -o $(OBDIR)output.o output.cpp

$(OBDIR)parallel.o: parallel.cpp $(HEADER)
	@test -d $(OBDIR) || mkdir $(OBDIR)        
	$(CXX) -c $(FLAGS) -o $(OBDIR)parallel.o parallel.cpp


copy:
	cp -f clean $(RUNDIR)/.	
clean: 		
	$(RM) *~
	$(RM) $(OBDIR)*.o
