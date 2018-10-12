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

#determine if running REAL8 or REAL4

# copy the files to paramesh/mpi_source

os.system('cp *.c ../../../../mpi_source/.')
os.system('cp mpi_amr_checkpoint_re_mpiio.F90 ../../../../mpi_source/.')
os.system('cp mpi_amr_checkpoint_wr_mpiio.F90 ../../../../mpi_source/.')
os.system('cp *.h ../../../../headers/.')
os.system('cp test_checkpoint_mpiio.F90 ../../../../Tests/.')

# edit the Makefile.gnu in mpi_source and add the file names which need to be
# compiled in

#first check to see if Makefile has already been modified for c_interfaces
file = open("../../../../mpi_source/Makefile.gnu","r")
modified = 0
while 1:
	line = file.readline()
	if line == 'read_blocks_mpiio_r4.c \\\n':
		modified = 1
	if line == '':
		break
file.close()

if (modified == 0):
	file = open("../../../../mpi_source/Makefile.gnu","r")
	file2 = open("Makefile.gnu","w")
	while 1:
		line = file.readline()
		file2.write(line)
		if line == 'sources := \\\n':
			file2.write("read_blocks_mpiio_r4.c \\\n")
			file2.write("read_blocks_mpiio_r8.c \\\n")
			file2.write("write_blocks_mpiio_r4.c \\\n")
			file2.write("write_blocks_mpiio_r8.c \\\n")
		if line == '':
			break
	file.close()
	file2.close()

	# copy the new Makefile.gnu to source direstory
	os.system('cp Makefile.gnu ../../../../mpi_source/Makefile.gnu')

# edit the Makefile.gnu in Tests and add the file names which need to be
# compiled in

#first check to see if Makefile has already been modified for c_interfaces
file = open("../../../../Tests/Makefile.gnu","r")
modified = 0
while 1:
	line = file.readline()
	if line == 'checkpoint_mpiio \\\n':
		modified = 1
	if line == '':
		break
file.close()

if (modified == 0):
	file = open("../../../../Tests/Makefile.gnu","r")
	file2 = open("Makefile.gnu","w")
	while 1:
		line = file.readline()
		file2.write(line)
		if line == 'tests := \\\n':
			file2.write("checkpoint_mpiio \\\n")
		if line == '':
			break
	file.close()
	file2.close()

	# copy the new Makefile.gnu to source direstory
	os.system('cp Makefile.gnu ../../../../Tests/Makefile.gnu')




