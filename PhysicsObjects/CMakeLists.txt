
include_directories(include)
include_directories(.)
include_directories(${EXTRA_INCLUDE_DIRS})

file(GLOB EXHH "include/*.hh")

EXCLUDEEMACS(EXHH)

SET (USERPOLIBS "")
FOREACH (lib ${EXHH})
	GET_FILENAME_COMPONENT(libName ${lib} NAME_WE)
	IF (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/include/LinkDef_${libName}.h)
		ROOT_GENERATE_DICTIONARY(${libName}DICT include/${libName}.hh LINKDEF include/LinkDef_${libName}.h MODULE ${libName})
		
		#Create fake library object to avoid compiling sources twice
        add_library(${libName}ObjLib OBJECT include/${libName}.hh src/${libName}.cc ${libName}DICT.cxx)
        set_property(TARGET ${${libName}ObjLib} PROPERTY POSITION_INDEPENDENT_CODE 1)
        
		add_library(${libName}-static STATIC $<TARGET_OBJECTS:${libName}ObjLib>)
		add_library(${libName} SHARED $<TARGET_OBJECTS:${libName}ObjLib>)
	ELSE()
		add_library(${libName}${LIBTYPEPOSTFIX} ${LIBTYPE} include/${libName}.hh src/${libName}.cc)
	ENDIF()
	SET (USERPOLIBS ${USERPOLIBS} ${libName}${LIBTYPEPOSTFIX})
ENDFOREACH(lib)

SET (USERPOLIBS ${USERPOLIBS} PARENT_SCOPE)
