################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/SVcalling_core/GT_cal.cpp \
../src/SVcalling_core/ReadHandler.cpp \
../src/SVcalling_core/SveHandler_combine_contigs.cpp \
../src/SVcalling_core/SveHandler_ins_del.cpp \
../src/SVcalling_core/SveHandler_inv_tra.cpp \
../src/SVcalling_core/analysis.cpp \
../src/SVcalling_core/forceCalling.cpp \
../src/SVcalling_core/sve.cpp 

CPP_DEPS += \
./src/SVcalling_core/GT_cal.d \
./src/SVcalling_core/ReadHandler.d \
./src/SVcalling_core/SveHandler_combine_contigs.d \
./src/SVcalling_core/SveHandler_ins_del.d \
./src/SVcalling_core/SveHandler_inv_tra.d \
./src/SVcalling_core/analysis.d \
./src/SVcalling_core/forceCalling.d \
./src/SVcalling_core/sve.d 

OBJS += \
./src/SVcalling_core/GT_cal.o \
./src/SVcalling_core/ReadHandler.o \
./src/SVcalling_core/SveHandler_combine_contigs.o \
./src/SVcalling_core/SveHandler_ins_del.o \
./src/SVcalling_core/SveHandler_inv_tra.o \
./src/SVcalling_core/analysis.o \
./src/SVcalling_core/forceCalling.o \
./src/SVcalling_core/sve.o 


# Each subdirectory must supply rules for building sources it contributes
src/SVcalling_core/%.o: ../src/SVcalling_core/%.cpp src/SVcalling_core/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -I../src/htslib -I../src/ -O3 -Wall -c -fmessage-length=0 -static -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


clean: clean-src-2f-SVcalling_core

clean-src-2f-SVcalling_core:
	-$(RM) ./src/SVcalling_core/GT_cal.d ./src/SVcalling_core/GT_cal.o ./src/SVcalling_core/ReadHandler.d ./src/SVcalling_core/ReadHandler.o ./src/SVcalling_core/SveHandler_combine_contigs.d ./src/SVcalling_core/SveHandler_combine_contigs.o ./src/SVcalling_core/SveHandler_ins_del.d ./src/SVcalling_core/SveHandler_ins_del.o ./src/SVcalling_core/SveHandler_inv_tra.d ./src/SVcalling_core/SveHandler_inv_tra.o ./src/SVcalling_core/analysis.d ./src/SVcalling_core/analysis.o ./src/SVcalling_core/forceCalling.d ./src/SVcalling_core/forceCalling.o ./src/SVcalling_core/sve.d ./src/SVcalling_core/sve.o

.PHONY: clean-src-2f-SVcalling_core

