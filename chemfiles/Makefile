CPPFLAGS=-std=c++11 -O2 -lrt
INCLUDE=-I${CONDA_PREFIX}/include
LIBS=-L${CONDA_PREFIX}/lib -lchemfiles -Wl,-rpath,${CONDA_PREFIX}/lib
TARGETS=count distance parse_mmcif parse_mmtf parse_pdb ramachandran

all: ${TARGETS}

%: %.cpp
	${CXX} ${CPPFLAGS} ${INCLUDE} $< -o $@ ${LIBS}

.PHONY: clean
clean:
	rm -rf ${TARGETS}
