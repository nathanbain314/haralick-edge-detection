# Copyright (C) 2011-2012, Haldo Spontón <haldos@fing.edu.uy>
# Copyright (C) 2011-2012, Juan Cardelino <juanc@fing.edu.uy>

PROJECT(edge_detectors)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

SET(CMAKE_C_FLAGS "-std=c99 -lrt -msse3 -mavx -lpthread -O3 -ftree-vectorize -ftree-vectorizer-verbose=1")
SET(PROJECT_SOURCE_DIR src)
SET(PROJECT_BINARY_DIR bin)
 
FIND_PACKAGE(TIFF)
IF(TIFF_FOUND)
   include_directories(${TIFF_INCLUDE_DIR})
   #  link_directories(${TIFF_LIBRARY})
   SET( LIBS ${LIBS} ${TIFF_LIBRARIES})
ENDIF()

FIND_PACKAGE(JPEG)
if(JPEG_FOUND)
   include_directories(${JPEG_INCLUDE_DIR})
   # link_directories(${JPEG_LIBRARY})
   SET( LIBS ${LIBS} ${JPEG_LIBRARIES})
ENDIF()

FIND_PACKAGE(PNG)
IF(PNG_FOUND)
   include_directories(${PNG_INCLUDE_DIR})
   #  link_directories(${PNG_LIBRARY})
   SET( LIBS ${LIBS} ${PNG_LIBRARIES})
ENDIF()

MESSAGE("libraries: ${LIBS}")

ADD_SUBDIRECTORY(src)
MESSAGE("standard libraries: ${LIBS}")
MESSAGE("custom libraries: ${E_LIBS}")

ADD_EXECUTABLE(${PROJECT_BINARY_DIR}/test_mh ${PROJECT_SOURCE_DIR}/test_mh.c)
ADD_EXECUTABLE(${PROJECT_BINARY_DIR}/test_mh_log ${PROJECT_SOURCE_DIR}/test_mh_log.c)
ADD_EXECUTABLE(${PROJECT_BINARY_DIR}/test_fded ${PROJECT_SOURCE_DIR}/test_fded.c)
ADD_EXECUTABLE(${PROJECT_BINARY_DIR}/test_haralick ${PROJECT_SOURCE_DIR}/test_haralick.c)
TARGET_LINK_LIBRARIES(${PROJECT_BINARY_DIR}/test_mh ${E_LIBS} ${LIBS})
TARGET_LINK_LIBRARIES(${PROJECT_BINARY_DIR}/test_mh_log ${E_LIBS} ${LIBS})
TARGET_LINK_LIBRARIES(${PROJECT_BINARY_DIR}/test_fded ${E_LIBS} ${LIBS})
TARGET_LINK_LIBRARIES(${PROJECT_BINARY_DIR}/test_haralick ${E_LIBS} ${LIBS})

INSTALL(TARGETS ${PROJECT_BINARY_DIR}/test_mh RUNTIME DESTINATION bin)
INSTALL(TARGETS ${PROJECT_BINARY_DIR}/test_mh_log RUNTIME DESTINATION bin)
INSTALL(TARGETS ${PROJECT_BINARY_DIR}/test_fded RUNTIME DESTINATION bin)
INSTALL(TARGETS ${PROJECT_BINARY_DIR}/test_haralick RUNTIME DESTINATION bin)
