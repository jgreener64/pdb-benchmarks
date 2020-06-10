TARGETS=parse_pdb parse_mmcif count distance #ramachandran
all: ${TARGETS}

%: %.cc
	${CXX} -O2 -Igemmi/include $< -o $@

clean:
	rm -rf ${TARGETS}
