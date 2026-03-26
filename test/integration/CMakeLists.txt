################################################################################
# Tests

create_standard_test(NAME miam_api SOURCES test_miam_api.cpp)
create_standard_test(NAME dissolved_reversible_reaction_integration SOURCES test_dissolved_reversible_reaction.cpp)
create_standard_test(NAME henry_law_phase_transfer_integration SOURCES test_henry_law_phase_transfer.cpp)
