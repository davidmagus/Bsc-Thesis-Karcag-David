INCLUDE_DIRECTORIES(
  ${PROJECT_SOURCE_DIR}
  ${PROJECT_BINARY_DIR}
)

LINK_DIRECTORIES(
  ${PROJECT_BINARY_DIR}/lemon
)

ADD_EXECUTABLE(test_run 
    test_run.cc
    TSPrandom.cc
    Heuristic.cc
)
TARGET_LINK_LIBRARIES(test_run lemon)

ADD_EXECUTABLE(test_runs 
    test_runs.cc
    TSPrandom.cc
    Heuristic.cc
)
TARGET_LINK_LIBRARIES(test_runs lemon)