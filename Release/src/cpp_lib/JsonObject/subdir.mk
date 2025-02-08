################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/cpp_lib/JsonObject/CJsonObject.cpp \
../src/cpp_lib/JsonObject/demo.cpp 

C_SRCS += \
../src/cpp_lib/JsonObject/cJSON.c 

CPP_DEPS += \
./src/cpp_lib/JsonObject/CJsonObject.d \
./src/cpp_lib/JsonObject/demo.d 

C_DEPS += \
./src/cpp_lib/JsonObject/cJSON.d 

OBJS += \
./src/cpp_lib/JsonObject/CJsonObject.o \
./src/cpp_lib/JsonObject/cJSON.o \
./src/cpp_lib/JsonObject/demo.o 


# Each subdirectory must supply rules for building sources it contributes
src/cpp_lib/JsonObject/%.o: ../src/cpp_lib/JsonObject/%.cpp src/cpp_lib/JsonObject/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -I../src/htslib -I../src/ -O3 -Wall -c -fmessage-length=0 -static -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

src/cpp_lib/JsonObject/%.o: ../src/cpp_lib/JsonObject/%.c src/cpp_lib/JsonObject/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: Cross GCC Compiler'
	gcc -I../src/htslib -I../src/ -O3 -Wall -c -fmessage-length=0 -static -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


clean: clean-src-2f-cpp_lib-2f-JsonObject

clean-src-2f-cpp_lib-2f-JsonObject:
	-$(RM) ./src/cpp_lib/JsonObject/CJsonObject.d ./src/cpp_lib/JsonObject/CJsonObject.o ./src/cpp_lib/JsonObject/cJSON.d ./src/cpp_lib/JsonObject/cJSON.o ./src/cpp_lib/JsonObject/demo.d ./src/cpp_lib/JsonObject/demo.o

.PHONY: clean-src-2f-cpp_lib-2f-JsonObject

