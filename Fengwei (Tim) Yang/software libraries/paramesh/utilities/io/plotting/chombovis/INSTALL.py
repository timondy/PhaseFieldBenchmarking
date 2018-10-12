import os
import sys

#determine if using double or single precision
file = open("../../../../headers/paramesh_preprocessor.fh","r")
real8 = 0
while 1:
        line = file.readline()
        if line == '':
                break
        if '\n' in line:
                line = line[:-1]
        if '\r' in line:
                line = line[:-1]
        if line[0:13] == "#define REAL8":
                real8 = 1
                break
file.close()

file = open("underscore.h","w")

# determine the system we will be compiling and set preprocessor flag
# for adding an underscore or not

uname = os.uname()
if (uname[0] == 'Darwin' or uname[0] == 'IBM'):
	file.write("#undef UNDERSCORE\n")
elif (uname[1] == 'bohr'):
	file.write("#define DOUBLE_UNDERSCORE\n")
else:
	file.write("#define UNDERSCORE\n")
if (real8 == 1):
        file.write("#define REAL8\n")
else:
        file.write("#define REAL4\n")
file.close()

# copy the files to paramesh/mpi_source

os.system('cp *.c ../../../../mpi_source/.')
os.system('cp *.F90 ../../../../mpi_source/.')
os.system('cp *.h ../../../../headers/.')

