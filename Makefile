OBJS := main.cpp.o timer.cpp.o matrix_op.cpp.o linear_solvers.cpp.o
TARGET := solve

CXX := mpicxx
DEBUG := 0

ifeq ($(DEBUG),1)
CFLAGS := -c -Wall -Weffc++ -g -pg -fno-inline
LDFLAGS := -g -pg -fno-inline
else
#CFLAGS := -c -O3 -Wall -Winline -march=native -ftree-vectorize -ffast-math
#LDFLAGS := -O3 -Wall -Winline -march=native -ftree-vectorize -ffast-math
CFLAGS := -c -O2 -Wall -Winline
LDFLAGS := -O2 -Wall -Winline
endif

all: $(OBJS) $(TARGET)

.SUFFIXES: .cpp .o .h

$(TARGET): $(OBJS)
	$(CXX) $(LDFLAGS) $(OBJS) -o $@

%.cpp.o: %.cpp
	$(CXX) $(CFLAGS) $< -o $@

clean:
	rm -f $(OBJS) $(TARGET) gmon.out leak.out*
	
