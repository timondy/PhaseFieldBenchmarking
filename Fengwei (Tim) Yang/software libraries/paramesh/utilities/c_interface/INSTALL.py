import os
import sys

#determine if using double or single precision
file = open("../../headers/paramesh_preprocessor.fh","r")
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

#determine if running REAL8 or REAL4

# copy the files to paramesh/mpi_source

os.system('cp *.c ../../source/.')
os.system('cp *.F90 ../../source/.')
os.system('cp *.h ../../headers/.')
os.system('cp test_c_interface.c ../../Tests/.')

# edit the Makefile.gnu in source and add the file names which need to be
# compiled in

#first check to see if Makefile has already been modified for c_interfaces
file = open("../../source/Makefile.gnu","r")
modified = 0
while 1:
	line = file.readline()
	if line == 'c_amr_get_pointer.c \\\n':
		modified = 1
	if line == '':
		break
file.close()

if (modified == 0):
	file = open("../../source/Makefile.gnu","r")
	file2 = open("Makefile.gnu","w")
	while 1:
		line = file.readline()
		file2.write(line)
		if line == 'sources := \\\n':
			file2.write("c_amr_get_pointer.c \\\n")
			file2.write("fortran_get_pointer.F90 \\\n")
			file2.write("fortran_wrappers.F90 \\\n")
		if line == '':
			break
	file.close()
	file2.close()

	# copy the new Makefile.gnu to source direstory
	os.system('cp Makefile.gnu ../../source/Makefile.gnu')

# edit the Makefile.gnu in Tests and add the file names which need to be
# compiled in

#first check to see if Makefile has already been modified for c_interfaces
file = open("../../Tests/Makefile.gnu","r")
modified = 0
while 1:
	line = file.readline()
	if line == 'c_interface \\\n':
		modified = 1
	if line == '':
		break
file.close()

if (modified == 0):
	file = open("../../Tests/Makefile.gnu","r")
	file2 = open("Makefile.gnu","w")
	while 1:
		line = file.readline()
		file2.write(line)
		if line == 'tests := \\\n':
			file2.write("c_interface \\\n")
		if line == '':
			break
	file.close()
	file2.close()

	# copy the new Makefile.gnu to source direstory
	os.system('cp Makefile.gnu ../../Tests/Makefile.gnu')




