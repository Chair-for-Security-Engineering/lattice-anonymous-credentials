cmake_minimum_required(VERSION 3.11)

project(lac)

# Use C11 and crash if that is not available/possible
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED True)

add_executable(test)
add_executable(bench)

set(SRC_DIR "${CMAKE_SOURCE_DIR}/src")
set(TEST_DIR "${CMAKE_SOURCE_DIR}/test")
set(BENCH_DIR "${CMAKE_SOURCE_DIR}/bench")
set(INCLUDE_DIR
  "${CMAKE_SOURCE_DIR}/include"
  "${CMAKE_SOURCE_DIR}/include/arith"
  "${CMAKE_SOURCE_DIR}/include/arith/arith_q"
  "${CMAKE_SOURCE_DIR}/include/arith/arith_qiss"
  "${CMAKE_SOURCE_DIR}/include/arith/arith_qshow"
  "${CMAKE_SOURCE_DIR}/include/arith/arith_real"
  "${CMAKE_SOURCE_DIR}/include/arith/arith_z")

set(SOURCES_LAC
  "${SRC_DIR}/arith/arith.c"
  "${SRC_DIR}/arith/arith_q/poly_q.c"
  "${SRC_DIR}/arith/arith_q/poly_q_vec_d.c"
  "${SRC_DIR}/arith/arith_q/poly_q_vec_k.c"
  "${SRC_DIR}/arith/arith_q/poly_q_vec_m.c"
  "${SRC_DIR}/arith/arith_q/poly_q_mat_d_d.c"
  "${SRC_DIR}/arith/arith_q/poly_q_mat_d_k.c"
  "${SRC_DIR}/arith/arith_q/poly_q_mat_d_m.c"
  "${SRC_DIR}/arith/arith_q/poly_q_sampling.c"
  "${SRC_DIR}/arith/arith_qiss/poly_qiss.c"
  "${SRC_DIR}/arith/arith_qiss/poly_qiss_vec_d.c"
  "${SRC_DIR}/arith/arith_qiss/poly_qiss_vec_k.c"
  "${SRC_DIR}/arith/arith_qiss/poly_qiss_vec_l.c"
  "${SRC_DIR}/arith/arith_qiss/poly_qiss_vec_256.c"
  "${SRC_DIR}/arith/arith_qiss/poly_qiss_vec_256_l.c"
  "${SRC_DIR}/arith/arith_qiss/poly_qiss_vec_m2.c"
  "${SRC_DIR}/arith/arith_qiss/poly_qiss_mat_k_k.c"
  "${SRC_DIR}/arith/arith_qiss/poly_qiss_mat_d_k.c"
  "${SRC_DIR}/arith/arith_qiss/poly_qiss_mat_d_m2.c"
  "${SRC_DIR}/arith/arith_qiss/poly_qiss_mat_256l_m2.c"
  "${SRC_DIR}/arith/arith_qiss/poly_qiss_sampling.c"
  "${SRC_DIR}/arith/arith_qshow/poly_qshow.c"
  "${SRC_DIR}/arith/arith_qshow/poly_qshow_vec_d.c"
  "${SRC_DIR}/arith/arith_qshow/poly_qshow_vec_k.c"
  "${SRC_DIR}/arith/arith_qshow/poly_qshow_vec_l.c"
  "${SRC_DIR}/arith/arith_qshow/poly_qshow_vec_256.c"
  "${SRC_DIR}/arith/arith_qshow/poly_qshow_vec_256_l.c"
  "${SRC_DIR}/arith/arith_qshow/poly_qshow_vec_m1.c"
  "${SRC_DIR}/arith/arith_qshow/poly_qshow_vec_m2.c"
  "${SRC_DIR}/arith/arith_qshow/poly_qshow_mat_k_k.c"
  "${SRC_DIR}/arith/arith_qshow/poly_qshow_mat_d_k.c"
  "${SRC_DIR}/arith/arith_qshow/poly_qshow_mat_d_m1.c"
  "${SRC_DIR}/arith/arith_qshow/poly_qshow_mat_d_m2.c"
  "${SRC_DIR}/arith/arith_qshow/poly_qshow_mat_256l_m2.c"
  "${SRC_DIR}/arith/arith_qshow/poly_qshow_sampling.c"
  "${SRC_DIR}/arith/arith_real/poly_real.c"
  "${SRC_DIR}/arith/arith_real/poly_real_vec_2d.c"
  "${SRC_DIR}/arith/arith_real/poly_real_mat_2d_2d.c"
  "${SRC_DIR}/arith/arith_z/poly_z.c"
  "${SRC_DIR}/arith/arith_z/poly_z_vec_d.c"
  "${SRC_DIR}/arith/arith_z/poly_z_mat_d_d.c"
  "${SRC_DIR}/covariance.c"
  "${SRC_DIR}/four_squares.c"
  "${SRC_DIR}/fips202.c"
  "${SRC_DIR}/random.c"
  "${SRC_DIR}/randombytes.c"
  "${SRC_DIR}/sep.c"
  "${SRC_DIR}/osig.c"
  "${SRC_DIR}/osig_prove.c"
  "${SRC_DIR}/osig_proof_verify.c"
  "${SRC_DIR}/show.c"
  "${SRC_DIR}/show_prove.c"
  "${SRC_DIR}/show_proof_verify.c")

target_sources(test PUBLIC
  "${TEST_DIR}/main.c"
  "${SOURCES_LAC}")

target_sources(bench PUBLIC
  "${BENCH_DIR}/main.c"
  "${BENCH_DIR}/benchmark.c"
  "${BENCH_DIR}/bench_sep.c"
  "${BENCH_DIR}/bench_osig.c"
  "${BENCH_DIR}/bench_show.c"
  "${SOURCES_LAC}")

target_include_directories(test PUBLIC ${INCLUDE_DIR})
target_include_directories(bench PUBLIC ${INCLUDE_DIR} ${BENCH_DIR})

set(WARN_ERROR_FLAGS
  -Wall -Wextra -Wundef -Wshadow -Wpointer-arith
  -Wcast-align -Wstrict-prototypes -Wwrite-strings -Waggregate-return
  -Wcast-qual -Wswitch-default -Wswitch-enum -Wunreachable-code)

set(SANITIZER_FLAGS
  -fsanitize=address -fsanitize=leak -fno-omit-frame-pointer
  -fsanitize=undefined -fsanitize=float-divide-by-zero
  -fsanitize=float-cast-overflow)

set(TEST_OPTS
  -march=native
  -g
  -O3
  )

set(BENCH_OPTS
  -march=native
  -O3
  )

target_compile_options(test PUBLIC ${WARN_ERROR_FLAGS} ${TEST_OPTS} ${SANITIZER_FLAGS})
target_link_options(test PUBLIC ${TEST_OPTS} ${SANITIZER_FLAGS})

target_compile_definitions(bench PUBLIC REMOVE_ASSERTS=1)
target_compile_options(bench PUBLIC ${WARN_ERROR_FLAGS} ${BENCH_OPTS})
target_link_options(bench PUBLIC ${BENCH_OPTS})

target_link_libraries(test PUBLIC "${CMAKE_SOURCE_DIR}/libs/flint/libflint.so")
target_link_libraries(test PUBLIC m)

target_link_libraries(bench PUBLIC "${CMAKE_SOURCE_DIR}/libs/flint/libflint.so")
target_link_libraries(bench PUBLIC m)
