#Robert Grumbine
# 1998

all : libombf_4.a libombc_4.a

libombc_4.a :
	sorc/makelibombC.sh

libombf_4.a :
	sorc/makelibombF.sh

clean : 
	rm libombf_4.a libombc_4.a sorc/*.o

