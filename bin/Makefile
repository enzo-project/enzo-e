clean:
	-rm -f test_* *.hdf5 *.png *.unit cello *.out *~
run-tests:
	./test_array
	./test_data
	./test_disk_hdf5
	./test_disk_ifrit
	./test_enzo_method_ppm
	./test_error
	./test_memory
	./test_monitor
	./test_parameters
	./test_performance
	./test_schedule
	./test_tree
	./test_tree_k
	./test_parse < implosion.in > implosion.out
	./test_parse < shock-pool.in > shock-pool.out
	./test_parse < test_groups.in > test_groups.out
	./test_parse < test.in > test.out
	./test_parse < test_list.in > test_list.out
	./test_parse < test_logical.in > test_logical.out
	./test_parse < test_scalar.in > test_scalar.out
	-grep -i -l error *.out

#	./test_ppm
#	./test_ppml

