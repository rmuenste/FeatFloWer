## This file should be placed in the root directory of your project.
## Then modify the CMakeLists.txt file in the root directory of your
## project to incorporate the testing dashboard.
## # The following are required to uses Dart and the Cdash dashboard
##   ENABLE_TESTING()
##   INCLUDE(CTest)
set(CTEST_PROJECT_NAME "FlowAroundACylinder")
#set(CTEST_NIGHTLY_START_TIME "01:00:00 UTC")

set(CTEST_DROP_METHOD "http")
set(CTEST_DROP_SITE "129.217.165.10/CDash/public")
set(CTEST_DROP_LOCATION "/submit.php?project=FEATFLOW_Q2P1")
set(CTEST_DROP_SITE_CDASH TRUE)
