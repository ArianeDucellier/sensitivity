# Declaration de l'extension
.SUFFIXES:.f90

# Compilateur et editeur de liens
FC      = ifort -traceback
FLD     = $(FC)
RM      = rm -f

# Option de profilage
PROF    = -O3 -openmp

# Option de compilation
FFLAGS  = $(PROF)

# Option edition des liens
LDFLAGS = $(PROF)

# Objets Fortran
FOBJS   = $(FSRCS:.f90=.o)

# Sources Fortran
FSRCS   = $(wildcard *.f90)

# Compilation Fortran
.f90.o:
		$(FC) -c $(FFLAGS) $<

EXEC = ./MC.exe
DEST =

all:		$(FOBJS)
		$(FLD) $(LDFLAGS) $(FOBJS) $(LIBS) -o $(EXEC)

# Nettoyage
clean:
		$(RM) $(FOBJS) $(EXEC) *.mod
