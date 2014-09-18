NVCC		= nvcc
FLAGS		= -arch=sm_35 -I./include
DFLAGS	= -G -g -lineinfo
RFLAGS 	= -O3 -DNDEBUG 

#SILO_INCLUDES := /usr/local/silo-4.8/include
#SILO_LIBS := /usr/local/silo-4.8/lib

LINKFLAGS = 
#LINKFLAGS += -L$(SILO_LIBS) -lsilo

#INC_SILO:= -I$(SILO_INCLUDES)

all: release 

debug: LINKFLAGS += -G -g

release: 	FLAGS += $(RFLAGS)
debug: 		FLAGS += $(DFLAGS)

release: lulesh
debug: lulesh

lulesh: allocator.o lulesh.o
	$(NVCC) $(LINKFLAGS) allocator.o lulesh.o -o lulesh

allocator.o: src/allocator.cu include/vector.h
	$(NVCC) $(FLAGS) src/allocator.cu -I ./ -c -o allocator.o

lulesh.o: src/lulesh.cu include/util.h include/vector.h include/texture_objAPI.h include/allocator.h
	$(NVCC) $(FLAGS) src/lulesh.cu -I ./  $(INC_SILO) -c -o lulesh.o

clean: 
	rm -rf allocator.o  lulesh.o lulesh




