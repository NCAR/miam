create_standard_test(NAME aerosol_property SOURCES aerosol_property.cpp)
create_standard_test(NAME condensation_rate SOURCES condensation_rate.cpp)
create_standard_test(NAME model SOURCES model.cpp)
create_standard_test(NAME process_set SOURCES process_set.cpp)

add_subdirectory(processes)
add_subdirectory(constraints)
add_subdirectory(representations)


