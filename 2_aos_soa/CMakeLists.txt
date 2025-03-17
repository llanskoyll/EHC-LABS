
file(GLOB SOURCE_V0 *.cpp)

project(aos_soa LANGUAGES CXX)

add_executable(aos_soa ${SOURCE_V0})

target_compile_options(aos_soa PRIVATE /ZI /DEBUG /Qipo- /Qopt-report:max /QxCORE-AVX2 /Qopt-zmm-usage:high /O3 /Qopenmp /I"D:\\Program Files\\oneAPI\\advisor\\latest\\include")

message(STATUS ${CMAKE_CURRENT_SOURCE_DIR})
