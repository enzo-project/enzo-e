
add_compile_definitions(NO_FREETYPE)

add_library(external_library
    ${CMAKE_CURRENT_SOURCE_DIR}/ExprEval/defs.h
    ${CMAKE_CURRENT_SOURCE_DIR}/ExprEval/except.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ExprEval/except.h
    ${CMAKE_CURRENT_SOURCE_DIR}/ExprEval/expr.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ExprEval/expreval.h
    ${CMAKE_CURRENT_SOURCE_DIR}/ExprEval/expr.h
    ${CMAKE_CURRENT_SOURCE_DIR}/ExprEval/func.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ExprEval/funclist.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ExprEval/funclist.h
    ${CMAKE_CURRENT_SOURCE_DIR}/ExprEval/node.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ExprEval/node.h
    ${CMAKE_CURRENT_SOURCE_DIR}/ExprEval/parser.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ExprEval/parser.h
    ${CMAKE_CURRENT_SOURCE_DIR}/ExprEval/vallist.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ExprEval/vallist.h
    ${CMAKE_CURRENT_SOURCE_DIR}/pngconf.h
    ${CMAKE_CURRENT_SOURCE_DIR}/pngwriter.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/pngwriter.h
    )
