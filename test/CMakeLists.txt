file ( GLOB code "*.cpp" "*.h" )

add_executable ( geometry3_tests ${code}  )
target_link_libraries ( geometry3_tests geometry3 )
