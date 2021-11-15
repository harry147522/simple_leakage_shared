mod_usr.o: mod_map_profile.o mod_yeofrho.o
gmunu: mod_map_profile.o mod_yeofrho.o

%.o: %.f
	$(F90) $(F90FLAGS) -c $< -o $@ $(addprefix -I,$(INC_DIRS))
