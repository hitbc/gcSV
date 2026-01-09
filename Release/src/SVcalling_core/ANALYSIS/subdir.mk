################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/SVcalling_core/ANALYSIS/analysis.cpp 

CPP_DEPS += \
./src/SVcalling_core/ANALYSIS/analysis.d 

OBJS += \
./src/SVcalling_core/ANALYSIS/analysis.o 


# Each subdirectory must supply rules for building sources it contributes
src/SVcalling_core/ANALYSIS/%.o: ../src/SVcalling_core/ANALYSIS/%.cpp src/SVcalling_core/ANALYSIS/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -std=c++17 -I../src/htslib -I../src -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


clean: clean-src-2f-SVcalling_core-2f-ANALYSIS

clean-src-2f-SVcalling_core-2f-ANALYSIS:
	-$(RM) ./src/SVcalling_core/ANALYSIS/analysis.d ./src/SVcalling_core/ANALYSIS/analysis.o

.PHONY: clean-src-2f-SVcalling_core-2f-ANALYSIS

