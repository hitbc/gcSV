################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/cpp_lib/cpp_utils.cpp \
../src/cpp_lib/graph.cpp 

CPP_DEPS += \
./src/cpp_lib/cpp_utils.d \
./src/cpp_lib/graph.d 

OBJS += \
./src/cpp_lib/cpp_utils.o \
./src/cpp_lib/graph.o 


# Each subdirectory must supply rules for building sources it contributes
src/cpp_lib/%.o: ../src/cpp_lib/%.cpp src/cpp_lib/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -I../src/htslib -I../src -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


clean: clean-src-2f-cpp_lib

clean-src-2f-cpp_lib:
	-$(RM) ./src/cpp_lib/cpp_utils.d ./src/cpp_lib/cpp_utils.o ./src/cpp_lib/graph.d ./src/cpp_lib/graph.o

.PHONY: clean-src-2f-cpp_lib

