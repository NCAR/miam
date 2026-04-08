################################################################################
# Tests

create_standard_test(NAME miam_api SOURCES test_miam_api.cpp)
create_standard_test(NAME dissolved_reaction_integration SOURCES test_dissolved_reaction.cpp)
create_standard_test(NAME dissolved_reversible_reaction_integration SOURCES test_dissolved_reversible_reaction.cpp)
create_standard_test(NAME henry_law_phase_transfer_integration SOURCES test_henry_law_phase_transfer.cpp)
create_standard_test(NAME readme_example SOURCES test_readme_example.cpp)
create_standard_test(NAME equilibrium_constraints_integration SOURCES test_equilibrium_constraints.cpp)
create_standard_test(NAME henry_law_equilibrium_constraint_integration SOURCES test_henry_law_equilibrium_constraint.cpp)
create_standard_test(NAME kinetic_vs_constrained SOURCES test_kinetic_vs_constrained.cpp)
create_standard_test(NAME jacobian_verification SOURCES test_jacobian_verification.cpp)
