TARGETS=parse_pdb parse_mmcif count distance
all: ${TARGETS}

%: %.cc
	${CXX} -std=c++11 -lrt -O2 -Igemmi/include $< -o $@

clean:
	rm -rf ${TARGETS}
