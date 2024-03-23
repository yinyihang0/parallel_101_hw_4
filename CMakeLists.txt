cmake_minimum_required(VERSION 3.12)
project(intrin LANGUAGES CXX)

add_executable(ori_sse_compare main.cpp)

set_target_properties(ori_sse_compare PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED ON
    CMAKE_BUILD_TYPE "Release"
)
# 检查是否是GCC编译器
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    # 对于GCC，设置对应的编译选项
    target_compile_options(ori_sse_compare PRIVATE
        -O3        # 最大优化
        -ffast-math # 优先考虑速度而不是精确度的数学运算
        -msse3     # 使用SSE3指令集
        -fopenmp   # 启用OpenMP支持
        -Wall        # 在GCC中没有直接对应的-W4选项，可能需要使用-Wall或者具体调整
    )
    target_link_options(ori_sse_compare PRIVATE
        -fopenmp   # 启用OpenMP支持
    )
else()
    # 如果不是GCC，打印一条警告消息
    message(WARNING "这个配置是为GCC编译器准备的，当前编译器不是GCC。请检查并修改编译选项。")
endif()