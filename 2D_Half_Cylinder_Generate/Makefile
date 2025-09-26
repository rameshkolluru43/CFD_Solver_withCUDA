EXEC := Genenerate_Half_Cylinder	
EXEC_DIRECTORY := ./
EXECUTABLES := Half_Cylinder

BASE_PATH := ../
Euler_Cube_Path   := 2D_Half_Cylinder_Generate/
BASIC_PATH:= Basic_Function_Files/
IPATH  := MYLIBRARY/include/
SPATH  := MYLIBRARY/src/
LPATH  := ../
IFLAGS := -I $(IPATH)  -I$(LPATH)
LFLAGS :=  #-lgsl -lgslcblas -lm -l fnccheck -L$(CFDRC)/lib -lDTF
DFLAGS := -ggdb3 # -g
OFLAGS :=  -O3  -O2 -O
WFLAGS := -Wall -Wno-long-long
EFLAGS := # -pg
AFLAGS := # -m64
CFLAGS := $(EFLAGS) -c $(AFLAGS) $(IFLAGS) $(DFLAGS) $(OFLAGS) $(WFLAGS) 
CC     := g++ # i586-mingw32msvc-g++
STRIP  := strip # i586-mingw32msvc-strip

OBJS          = $(BASE_PATH)$(BASIC_PATH)Point.o \
			$(BASE_PATH)$(BASIC_PATH)Line.o \
			$(BASE_PATH)$(BASIC_PATH)Circle.o \
			$(BASE_PATH)$(BASIC_PATH)Grid_Gen_Functions.o \
			$(BASE_PATH)$(BASIC_PATH)gridgen_Elliptic_TFI.o \
			$(BASE_PATH)$(BASIC_PATH)Utility_Files.o \
			$(BASE_PATH)$(Euler_Cube_Path)Half_Cylinder.o\

default:@echo "Type one of the following "
		@echo "  make  pc"
		@echo "  make  pc_debug"
		@echo "  make  pc_static"
		@echo "  make  pc_64"
		@echo "  make  pc_64_debug"
		@echo "  make  pc_64_static"

$(EXEC): $(OBJS)
		$(CC)  $(EFLAGS) $(AFLAGS) $(OBJS)  $(LFLAGS) -o $(EXEC_DIRECTORY)$@

%.o: %.cc
		$(CC) $(CFLAGS) `pwd`/$< -o $@
		@-makedepend -I$(IPATH) -I$(LPATH) -a -fMake_dependencies $< > /dev/null 2>&1

%.o: %.cpp
		$(CC)$(CFLAGS) `pwd`/$< -o $@
		@-makedepend -I$(IPATH) -I$(LPATH) -a -fMake_dependencies $< > /dev/null 2>&1

clean:
		rm -f $(OBJS)
		rm -f $(EXEC_DIRECTORY)$(EXEC)
		@-touch Make_dependencies
		@-makedepend -fMake_dependencies __arbitrary_file > /dev/null 2>&1
archive:
		ar rcs $(EXEC_DIRECTORY)$(EXEC).a $(OBJS)

#linux 32 bit

pc:
	@for i in ${EXECUTABLES}; do make  EXEC="$$i" "$$i" \
	"AFLAGS=-m32" "OFLAGS=-O3"; done
pc_debug:
		@for i in ${EXECUTABLES}; do make  EXEC="$$i" "$$i" \
		"AFLAGS=-m32" "DFLAGS=-ggdb3"; done
pc_static:
		@for i in ${EXECUTABLES}; do make  EXEC="$$i" "$$i" \
		"AFLAGS=-m32" "OFLAGS=-O2" "LFLAGS=-static"; done
pc64:
		@for i in ${EXECUTABLES}; do make  EXEC="$$i" "$$i" \
		"AFLAGS=-m64" "OFLAGS=-O3"; done
pc64_debug:
		@for i in ${EXECUTABLES}; do make  EXEC="$$i" "$$i" \
		"AFLAGS=-m64" "DFLAGS=-ggdb3"; done
pc64_static:
		@for i in ${EXECUTABLES}; do make  EXEC="$$i" "$$i" \
		"AFLAGS=-m64" "OFLAGS=-O2" "LFLAGS=-static"; done
-include Make_dependencies