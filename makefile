#----------------------------------------------------
# GCC LINUX make file for cnd_frac3 directory
#----------------------------------------------------


#----- make NDEBUG=1 for nodebugging ------#

ifeq ($(NDEBUG), 1)
CL = gcc -O2 -c -DNDEBUG
LINK = gcc

#----- DEBUG case is below -----#
else   
CL = gcc -g -c -Wall -DDEBUG
LINK = gcc -g -Wall -DDEBUG
endif


#----- compiling and linking of program -----#

findpts : findpts.o 
	$(LINK) findpts.o -lm
#	$(LINK) -o findpts findpts.o -lm

findpts.o : findpts.c 
	$(CL) findpts.c 


#----- organization of files -----#
zoo :
	zoo ahP: cntd_frac3.zoo *

gzip :
	tar cvf cntd_frac3.tar *
	gzip -9 cntd_frac3.tar
	mv cntd_frac3.tar.gz  cntd_frac3.tgz 

bzip :
	tar cvf cntd_frac3.tar *
	bzip2 -9 cntd_frac3.tar


clean :
	for f in *.o;   do rm -f $$f; done
	for f in *.bak; do rm -f $$f; done	
	for f in *.zoo; do rm -f $$f; done
	for f in a.out; do rm -f $$f; done
	for f in core;  do rm -f $$f; done
	for f in *.tgz; do rm -f $$f; done
	for f in *.tar; do rm -f $$f; done
	for f in *.bz2; do rm -f $$f; done
	for f in *.pyc; do rm -f $$f; done
	for f in *~;    do rm -f $$f; done
	for f in *.dvi; do rm -f $$f; done
	for f in *.log; do rm -f $$f; done

     
 
