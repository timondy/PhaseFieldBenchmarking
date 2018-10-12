import os
import sys

# copy the files to paramesh/source

os.system('cp amr_get_new_nodetypes.F90 ../../source/.')
os.system('cp amr_mg_common.F90 ../../headers/.')
os.system('cp amr_mg_init.F90 ../../source/.')
os.system('cp amr_mg_morton_process.F90 ../../source/.')
os.system('cp amr_mg_prolong.F90 ../../source/.')
os.system('cp amr_mg_restrict.F90 ../../source/.')
os.system('cp mpi_amr_store_comm_info_mg.F90 ../../source/.')

os.system('cp test_multigrid.F90 ../../Tests/.')

# edit the Makefile.gnu in source and add the file names which need to be
# compiled in

#first check to see if Makefile has already been modified for multigrid in source

file = open("../../source/Makefile.gnu","r")
modified = 0
while 1:
	line = file.readline()
	if line == 'amr_get_new_nodetypes.F90 \\\n':
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
			file2.write("amr_get_new_nodetypes.F90 \\\n")
			file2.write("amr_mg_init.F90 \\\n")
			file2.write("amr_mg_morton_process.F90 \\\n")
			file2.write("amr_mg_prolong.F90 \\\n")
			file2.write("amr_mg_restrict.F90 \\\n")
			file2.write("mpi_amr_store_comm_info_mg.F90 \\\n")
		if line == '':
			break
	file.close()
	file2.close()

	# copy the new Makefile.gnu to source direstory
	os.system('cp Makefile.gnu ../../source/Makefile.gnu')

# edit the Makefile.gnu in headers and add the file names which need to be
# compiled in

#first check to see if Makefile has already been modified for multigrid in headers

file = open("../../headers/Makefile.gnu","r")
modified = 0
while 1:
	line = file.readline()
	if line == 'amr_mg_common.F90 \\\n':
		modified = 1
	if line == '':
		break
file.close()

if (modified == 0):
	file = open("../../headers/Makefile.gnu","r")
	file2 = open("Makefile.gnu","w")
	while 1:
		line = file.readline()
		file2.write(line)
		if line == 'tree.F90 \\\n':
			file2.write("amr_mg_common.F90 \\\n")
		if line == '':
			break
	file.close()
	file2.close()

	# copy the new Makefile.gnu to source direstory
	os.system('cp Makefile.gnu ../../headers/Makefile.gnu')

#first check to see if Makefile has already been modified for multigrid
file = open("../../Tests/Makefile.gnu","r")
modified = 0
while 1:
        line = file.readline()
        if line == 'multigrid \\\n':
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
                        file2.write("multigrid \\\n")
                if line == '':
                        break
        file.close()
        file2.close()

        # copy the new Makefile.gnu to source direstory
        os.system('cp Makefile.gnu ../../Tests/Makefile.gnu')
