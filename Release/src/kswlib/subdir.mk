################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../src/kswlib/kalloc.c \
../src/kswlib/ksw2_extd2_sse.c 

C_DEPS += \
./src/kswlib/kalloc.d \
./src/kswlib/ksw2_extd2_sse.d 

OBJS += \
./src/kswlib/kalloc.o \
./src/kswlib/ksw2_extd2_sse.o 


# Each subdirectory must supply rules for building sources it contributes
src/kswlib/%.o: ../src/kswlib/%.c src/kswlib/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: Cross GCC Compiler'
	gcc -I../src/htslib -I../src -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


clean: clean-src-2f-kswlib

clean-src-2f-kswlib:
	-$(RM) ./src/kswlib/kalloc.d ./src/kswlib/kalloc.o ./src/kswlib/ksw2_extd2_sse.d ./src/kswlib/ksw2_extd2_sse.o

.PHONY: clean-src-2f-kswlib

