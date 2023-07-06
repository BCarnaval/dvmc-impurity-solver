# Include guard
cmake_minimum_required(VERSION 3.12 FATAL_ERROR)
if(${CMAKE_PROJECT_NAME} STREQUAL "Project")
  message(
    FATAL_ERROR
      "CMake shouldn't be executed for 'src' subdirectory, but for the top directory of mVMC."
  )
endif(${CMAKE_PROJECT_NAME} STREQUAL "Project")

add_executable(dvmc_gnuplot src/dvmc_gnuplot.c)
add_library(dvmc_speedup SHARED src/dvmc_speedup.c)
set_property(TARGET dvmc_speedup PROPERTY POSITION_INDEPENDENT_CODE ON)

install(TARGETS dvmc_gnuplot RUNTIME DESTINATION $ENV{HOME}/.local/share/dvmc)
install(TARGETS dvmc_speedup DESTINATION $ENV{HOME}/.local/share/dvmc)

set(SCRIPTS
    filter_bin.sh
    dvmc.py
    fij_sym.py
    averageOhgoe.py
    cluster.py
    dos_analysis.py
    dvmc_spectrum.py
    dvmc_spectrum_eigh_w_sqrtS.py
    dvmc_spectrum_w_SVD.py
    dvmc_spectrum_w_cond_number.py
    generate_template_gnuplot.py
    init_groundstate.py
    init_params.py
    makeDH4.py
    makeExcitation_from_hopping.py
    makeExcitation_from_hopping_only_t.py
    makeExcitation.py
    convertOutputBin.py
    mergeOutputBin.py
    postprocess_dVMC.sh
    preprocess_dVMC.sh
    selectExcitation.py
    select_kPath.py
    template_plot_Akw.gp
    template_plot_allAkw.gp
    template_plot_dos.gp
    hopping.py)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/${SCRIPTS} DESTINATION .)

install(
  FILES ${CMAKE_CURRENT_BINARY_DIR}/${SCRIPTS}
  DESTINATION $ENV{HOME}/.local/share/dvmc
  PERMISSIONS
    OWNER_READ
    OWNER_WRITE
    OWNER_EXECUTE
    GROUP_READ
    GROUP_EXECUTE
    WORLD_READ
    WORLD_EXECUTE)
