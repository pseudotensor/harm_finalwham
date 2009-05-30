
# for long doubles
CFLAGSGCCOTHER= -O2 -Wall -Wuninitialized -fno-builtin
CFLAGSOTHER= -O3 -tpp7 -axKW -unroll -long_double -Wall -Wuninitialized -fno-builtin


ifeq ($(USEMPI),1)
CC=$(MCC) -cc=$(COMP) $(CFLAGS3) $(CFLAGS2) $(CDEBUG)
CC2=$(MCC) -cc=gcc $(GCCCFLAGS)
endif # endif usempich==1

ifeq ($(USEMPI),0)
CC=$(COMP)  $(CFLAGS3) $(CFLAGS2) $(CDEBUG)
CC2=gcc  $(GCCCFLAGS)
endif # endif usempich==0

include maketail.harm.inc
include maketail.ldouble.inc

# "make prep" must be called first so these are set correctly
#MYCOMMANDS1 := $(shell cat objsonepoint.cmd.list ) 
#MYCOMMANDS2 := $(shell cat objsmanypoints.cmd.list ) 


#all:	$(PREP) $(CMD) $(FINISH)
all:	$(CMD) $(FINISH)



$(PREP):
	( sh ./makedecs.h.sh )
	( sh ./makempidecs.h.sh )
	#
	( echo $(OBJSONEPOINT) > objsonepoint )
	# remove and touch (so exists) files to create
	( rm -rf objsonepoint.cmd.list ; rm -rf onepoint ; touch objsonepoint.cmd.list ; touch onepoint)
	# create list of prefixes
	( for fil in `cat objsonepoint` ; do prefa=`echo $$fil | sed 's/\.o/ /'` ; pref=`echo $$prefa | awk '{print $$1}'` ; echo -e $$pref >> onepoint ; done )
	# create list of commands
	( for fil in `cat onepoint` ; do echo -e $$fil.o :\\n\\t$(CC) $(CFLAGSNONPRECISE) -c $$fil.c $(LDFLAGS) >> objsonepoint.cmd.list ; done )
	# remove old depenedencies file and create
	( rm -rf deps.onepoint.list ; touch deps.onepoint.list )
	# create new dependencies file
	( for fil in `cat onepoint` ; do  echo $$fil.c ; $(CC) -MM $$fil.c >> deps.onepoint.list ; done )
	#
	# below is only different in CFLAGS part
	( echo $(OBJSMANYPOINTS) > objsmanypoints )
	# remove and touch (so exists) files to create
	( rm -rf objsmanypoints.cmd.list ; rm -rf manypoints ; touch objsmanypoints.cmd.list ; touch manypoints)
	# create list of prefixes
	( for fil in `cat objsmanypoints` ; do prefa=`echo $$fil | sed 's/\.o/ /'` ; pref=`echo $$prefa | awk '{print $$1}'` ; echo -e $$pref >> manypoints ; done )
	# create list of commands
	( for fil in `cat manypoints` ; do echo -e $$fil.o :\\n\\t$(CC) $(CFLAGS) -c $$fil.c $(LDFLAGS) >> objsmanypoints.cmd.list ; done )
	# remove old depenedencies file and create
	( rm -rf deps.manypoints.list ; touch deps.manypoints.list )
	# create new dependencies file
	( for fil in `cat manypoints` ; do  echo $$fil.c ; $(CC) -MM $$fil.c >> deps.manypoints.list ; done )
	#
	#
#	( echo $(OBJSMANYPOINTS) > objsmanypoints )
#	( rm -rf objsmanypoints.cmd.list ; touch objsmanypoints.cmd.list )
#	( for fil in `cat objsmanypoints` ; do echo $(CC) $(CFLAGS) -c$$fil $(LDFLAGS) >> objsmanypoints.cmd.list ; done )
#	( for fil in `cat objsmanypoints` ; do echo $$fil.c ; $(CC) -MM $$fil >> deps.list ; done )
	
$(FINISH):
#	rm -rf objsonepoint objsmanypoints

#$(CMD):   	GROUP1 GROUP2 makefile freespace.o
#	echo objsmanypoints: $(OBJSMANYPOINTS)
#	echo objsonepoint: $(OBJSONEPOINT)
#	$(MYCOMMANDS2)
#	$(MYCOMMANDS1)


$(CMD):	$(OBJSMANYPOINTS) $(OBJSONEPOINT) makefile freespace.o
	$(CC)  $(CFLAGSNONPRECISE) -o $(CMD) $(OBJSMANYPOINTS) $(OBJSONEPOINT) freespace.o $(LDFLAGS)


$(CMD2):       $(OBJSMANYPOINTS) $(OBJSONEPOINT) makefile libml.a $(OBJSOTHER)
	       $(CC2) $(CFLAGSGCCOTHER) -c polevll.c $(LDFLAGSOTHER)
	       $(CC)  $(CFLAGS) -o $(CMD2) $(OBJSMANYPOINTS) $(OBJSONEPOINT) freespace.o polevll.o libml.a $(LDFLAGSOTHER)



-include maketail.ldouble.deps.inc

-include objsonepoint.cmd.list

-include objsmanypoints.cmd.list


-include deps.manypoints.list

-include deps.onepoint.list

freespace.o : 
	    $(CC2) $(GCCCFLAGS) -c freespace.c $(LDFLAGS)







#$(OBJSMANYPOINTS) :

#$(MYCOMMANDS2)

#$(OBJSONEPOINT) :

#$(MYCOMMANDS1)





#GROUP2:	$(OBJSMANYPOINTS) makefile
#	$(CC) $(CFLAGS) -c $(SRCSONEPOINT) $(LDFLAGS)
#	CFLAGS=$(CFLAGS)
#	echo $(MYCOMMANDS2)
#	$(MYCOMMANDS2)

#GROUP1:	$(OBJSONEPOINT) makefile
#	$(CC) $(CFLAGSNONPRECISE) -c $(SRCSONEPOINT) $(LDFLAGS)
#	CFLAGS=$(CFLAGSNONPRECISE)
#	echo $(MYCOMMANDS1)
#	$(MYCOMMANDS1)

clean:
	rm *.o

cleani:
	rm *.o *.il

cleanall:
	rm *.o *.il *~

cleanbackup:
	rm *~


# below file not used anymore since generated automagically
#include maketail.harm.deps.inc

