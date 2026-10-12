create_standard_test(NAME aerosol_property SOURCES aerosol_property.cpp)
create_standard_test(NAME condensation_rate SOURCES condensation_rate.cpp)
create_standard_test(NAME model SOURCES model.cpp)

add_subdirectory(processes)
add_subdirectory(constraints)
add_subdirectory(representations)


