CC=g++
CFLAGS=-O3 -ffast-math
LDFLAGS=-static# -lm
PROGRAM=RTMalign se pdb2xyz xyz_sfetch pdb2fasta pdb2ss

all: ${PROGRAM}

RTMalign: RTMalign.cpp param_set.h basic_fun.h Kabsch.h NW.h RTMalign.h pstream.h NWalign.h BLOSUM.h
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

# se: se.cpp se.h param_set.h basic_fun.h Kabsch.h NW.h RNAalign.h pstream.h
# 	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

# pdb2xyz: pdb2xyz.cpp basic_fun.h pstream.h
# 	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

# xyz_sfetch: xyz_sfetch.cpp
# 	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

# pdb2fasta: pdb2fasta.cpp basic_fun.h pstream.h
# 	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

# pdb2ss: pdb2ss.cpp RNAalign.h basic_fun.h pstream.h
# 	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

# clean:
# 	rm -f ${PROGRAM}
