#!/bin/bash
# Created by Artur Nowicki on 03.02.2018.
ok_status=0
err_missing_program_input=100
err_f_open=101
err_f_read=102
err_f_write=103
err_f_close=104
err_memory_alloc=105

in_dir='../../data/boundary_conditions/tmp_bin_data/'
out_dir='../../data/boundary_conditions/spread_data/'

source ./../common_code/assertions.sh
total_tests=0
failed_tests=0

echo "Compile program."
gfortran ../common_code/messages.f90 ../common_code/error_codes.f90 poisson_solver.f90 -o poisson_solver
if [[ $? -ne 0 ]]; then
	exit
fi

test_in_file=${in_dir}'2018-01-01-46800_SU_0600_0640_0001_0001.ieeer8'
test_out_file=${out_dir}'2018-01-01-46800_SU_0600_0640_0001_0001.ieeer8'
x_in=600
y_in=640
z_in=1

expected_error_code=${ok_status}
echo "-------------------------"
echo "Test missing parameters"
./poisson_solver ${test_in_file}
assertNotEquals ${expected_error_code} $?
failed_tests=$((failed_tests+$?))
total_tests=$((total_tests+1))

echo "-------------------------"
echo "Test wrong dim parameter"
./poisson_solver 'ban_in_file.bin' ${test_out_file} blabla ${y_in} ${z_in}
assertNotEquals ${expected_error_code} $?
failed_tests=$((failed_tests+$?))
total_tests=$((total_tests+1))

echo "-------------------------"
echo "Test wrong input file"
./poisson_solver 'ban_in_file.bin' ${test_out_file} ${x_in} ${y_in} ${z_in}
assertNotEquals ${expected_error_code} $?
failed_tests=$((failed_tests+$?))
total_tests=$((total_tests+1))

echo "-------------------------"
echo "Test wrong output file"
./poisson_solver ${test_in_file} 'badpath/bad_out_file.bin' ${x_in} ${y_in} ${z_in}
assertNotEquals ${expected_error_code} $?
failed_tests=$((failed_tests+$?))
total_tests=$((total_tests+1))

echo "-------------------------"
echo "Test all ok"
expected_error_code=${ok_status}
./poisson_solver ${test_in_file} ${test_out_file} ${x_in} ${y_in} ${z_in}
assertEquals ${expected_error_code} $?
failed_tests=$((failed_tests+$?))
total_tests=$((total_tests+1))

echo
echo "-------------------------"
echo "TESTING RESULTS:"
echo "Tests failed: ${failed_tests} out of ${total_tests}"

if [[ ${failed_tests} -ne 0 ]]; then
	exit
fi

echo "-------------------------"
echo "Start actual script:"

for in_file in ${in_dir}; do
	# echo ${in_f1:0:23}
	out_file=${in_file/${in_dir}}
	echo ${out_file}
	./poisson_solver ${in_file} ${out_file}
	# ./poisson_solver ${in_path} ${in_file2} ${out_file2}
done