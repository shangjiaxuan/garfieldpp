
# Enables all the relevant warnings for the different compilers
function(garfield_enable_default_warnings target_name)
target_compile_options(${target_name} PRIVATE
$<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:GNU>>:
    -Wall -Wextra -pedantic -Wfatal-errors -Wno-missing-braces -Wshadow>
$<$<CXX_COMPILER_ID:MSVC>:/W3>)
endfunction()


##########################################################################################################
# Forces the color output for GNU and Clang in all cases
# Useful to enable that for Ninja build
function(force_color_output)
	option(force_colored_output "Always produce ANSI-colored output (GNU/Clang only)." TRUE )
	mark_as_advanced(force_colored_output)
	if (force_colored_output)
	add_compile_options( $<$<CXX_COMPILER_ID:GNU>:-fdiagnostics-color=always>
		$<$<AND:$<CXX_COMPILER_ID:Clang>,$<VERSION_GREATER:$<CXX_COMPILER_VERSION>,4.0.0>>:-fcolor-diagnostics> )
	endif ()
endfunction()