cmake_minimum_required(VERSION 3.20)

project(libomb)

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -ansi ")

message("mmablib sorc " CMAKE_INSTALL_LIBDIR="${CMAKE_INSTALL_LIBDIR}")

#set(CMAKE_INSTALL_DIR "lib")

list(APPEND c_src
mapll.c
mapxy.c
)
add_library(ombc_4 ${c_src})

list(APPEND fortran_src
GRIBIT.f90
VCC.f90
arcdis.f90
bsslz1.f90
cfsread.f90
ciagree.f90
correl.f90
darcdis.f90
fit.f90
gaulat.f90
iagree.f90
mapxy.f90
new_wdir.f90
ssanaly2.f90
sumx.f90
sumx2.f90
sumxy.f90
tfreez.f90
tvap.f90
vectorize.f90
w3ft01.f90
waldrop_mse.f90
wdir.f90
wmoout.f90
)

add_library(ombf_4 ${fortran_src})

#install(TARGETS ${exe_name} RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(TARGETS ombf_4  RUNTIME DESTINATION ${CMAKE_INSTALL_LIBDIR} ) 
install(TARGETS ombc_4  RUNTIME DESTINATION ${CMAKE_INSTALL_LIBDIR} ) 

