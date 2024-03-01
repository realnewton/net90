path = src/

#FC = ifort -c -traceback        #compile
#LK = ifort -o runnet90  #link

FC = gfortran -c -O3          #compile
LK = gfortran -o runnet90  #link

OBJS =\
	eose.o\
	eosi.o\
	eosr.o\
	net90.o\
	const_eos_mod.o\
	helm_table_storage_mod.o\
	read_helm_table.o\
	helmeos.o\
	interp.o\
	read_rates_table.o\
	nuclear90_module.o\
	testnuclear.o

nuclearnet:\
  $(OBJS); $(LK) $(OBJS); rm *.o *.mod

eose.o:\
  $(path)$(@:.o=.f); $(FC) $(path)$(@:.o=.f)

eosi.o:\
  $(path)$(@:.o=.f); $(FC) $(path)$(@:.o=.f)

eosr.o:\
  $(path)$(@:.o=.f); $(FC) $(path)$(@:.o=.f)

const_eos_mod.o:\
  $(path)$(@:.o=.f90); $(FC) $(path)$(@:.o=.f90)

helm_table_storage_mod.o:\
  $(path)$(@:.o=.f90); $(FC) $(path)$(@:.o=.f90)

read_helm_table.o:\
	helm_table_storage_mod.o\
  $(path)$(@:.o=.f90); $(FC) $(path)$(@:.o=.f90)

helmeos.o:\
	const_eos_mod.o\
	helm_table_storage_mod.o\
  $(path)$(@:.o=.f90); $(FC) $(path)$(@:.o=.f90)

read_rates_table.o:\
	nuclear90_module.o\
  $(path)$(@:.o=.f90); $(FC) $(path)$(@:.o=.f90)

interp.o:\
	nuclear90_module.o\
  $(path)$(@:.o=.f90); $(FC) $(path)$(@:.o=.f90)

net90.o:\
    nuclear90_module.o\
	interp.o\
  $(path)$(@:.o=.f90); $(FC) $(path)$(@:.o=.f90)

nuclear90_module.o:\
  $(path)$(@:.o=.f90); $(FC) $(path)$(@:.o=.f90)

testnuclear.o:\
  nuclear90_module.o\
  eose.o\
  eosi.o\
  eosr.o\
  read_helm_table.o\
  helmeos.o\
  read_rates_table.o\
  $(@:.o=.f90); $(FC) $(@:.o=.f90)

clean:
	rm *.o *.mod
