SOURCE  := $(wildcard lib/*.c) $(wildcard *.c) $(wildcard *.cpp)
OBJS    := $(patsubst %.c,%.o,$(patsubst %.cpp,%.o,$(SOURCE)))
  
TARGET  := deSAMBA
  
#compile and lib parameter
CC      := gcc -std=c99
LIBS    := -lm -lz -lpthread
INCLUDE := -I ./lib
CFLAGS  := -g -Wall -O3 $(INCLUDE)

.PHONY : everything objs clean veryclean rebuild
  
everything : $(TARGET)
  
all : $(TARGET)
  
objs : $(OBJS)
  
rebuild: veryclean everything
                
clean :
	rm -fr *.so
	rm -fr *.o
	rm -fr lib/*.so
	rm -fr lib/*.o
	rm -fr $(TARGET) 
    
veryclean : clean
	rm -fr $(TARGET)
  
$(TARGET) : $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LIBS)
