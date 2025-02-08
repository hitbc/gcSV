################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/cpp_lib/statistics/StatsManager.cpp \
../src/cpp_lib/statistics/StatsTracker.cpp 

CPP_DEPS += \
./src/cpp_lib/statistics/StatsManager.d \
./src/cpp_lib/statistics/StatsTracker.d 

OBJS += \
./src/cpp_lib/statistics/StatsManager.o \
./src/cpp_lib/statistics/StatsTracker.o 


# Each subdirectory must supply rules for building sources it contributes
src/cpp_lib/statistics/%.o: ../src/cpp_lib/statistics/%.cpp src/cpp_lib/statistics/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -I../src/htslib -I../src/ -O3 -Wall -c -fmessage-length=0 -static -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


clean: clean-src-2f-cpp_lib-2f-statistics

clean-src-2f-cpp_lib-2f-statistics:
	-$(RM) ./src/cpp_lib/statistics/StatsManager.d ./src/cpp_lib/statistics/StatsManager.o ./src/cpp_lib/statistics/StatsTracker.d ./src/cpp_lib/statistics/StatsTracker.o

.PHONY: clean-src-2f-cpp_lib-2f-statistics

