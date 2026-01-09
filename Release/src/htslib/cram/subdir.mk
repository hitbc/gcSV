################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../src/htslib/cram/cram_codecs.c \
../src/htslib/cram/cram_decode.c \
../src/htslib/cram/cram_encode.c \
../src/htslib/cram/cram_external.c \
../src/htslib/cram/cram_index.c \
../src/htslib/cram/cram_io.c \
../src/htslib/cram/cram_samtools.c \
../src/htslib/cram/cram_stats.c \
../src/htslib/cram/files.c \
../src/htslib/cram/mFILE.c \
../src/htslib/cram/open_trace_file.c \
../src/htslib/cram/pooled_alloc.c \
../src/htslib/cram/rANS_static.c \
../src/htslib/cram/sam_header.c \
../src/htslib/cram/string_alloc.c 

C_DEPS += \
./src/htslib/cram/cram_codecs.d \
./src/htslib/cram/cram_decode.d \
./src/htslib/cram/cram_encode.d \
./src/htslib/cram/cram_external.d \
./src/htslib/cram/cram_index.d \
./src/htslib/cram/cram_io.d \
./src/htslib/cram/cram_samtools.d \
./src/htslib/cram/cram_stats.d \
./src/htslib/cram/files.d \
./src/htslib/cram/mFILE.d \
./src/htslib/cram/open_trace_file.d \
./src/htslib/cram/pooled_alloc.d \
./src/htslib/cram/rANS_static.d \
./src/htslib/cram/sam_header.d \
./src/htslib/cram/string_alloc.d 

OBJS += \
./src/htslib/cram/cram_codecs.o \
./src/htslib/cram/cram_decode.o \
./src/htslib/cram/cram_encode.o \
./src/htslib/cram/cram_external.o \
./src/htslib/cram/cram_index.o \
./src/htslib/cram/cram_io.o \
./src/htslib/cram/cram_samtools.o \
./src/htslib/cram/cram_stats.o \
./src/htslib/cram/files.o \
./src/htslib/cram/mFILE.o \
./src/htslib/cram/open_trace_file.o \
./src/htslib/cram/pooled_alloc.o \
./src/htslib/cram/rANS_static.o \
./src/htslib/cram/sam_header.o \
./src/htslib/cram/string_alloc.o 


# Each subdirectory must supply rules for building sources it contributes
src/htslib/cram/%.o: ../src/htslib/cram/%.c src/htslib/cram/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: Cross GCC Compiler'
	gcc -I../src/htslib -I../src -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


clean: clean-src-2f-htslib-2f-cram

clean-src-2f-htslib-2f-cram:
	-$(RM) ./src/htslib/cram/cram_codecs.d ./src/htslib/cram/cram_codecs.o ./src/htslib/cram/cram_decode.d ./src/htslib/cram/cram_decode.o ./src/htslib/cram/cram_encode.d ./src/htslib/cram/cram_encode.o ./src/htslib/cram/cram_external.d ./src/htslib/cram/cram_external.o ./src/htslib/cram/cram_index.d ./src/htslib/cram/cram_index.o ./src/htslib/cram/cram_io.d ./src/htslib/cram/cram_io.o ./src/htslib/cram/cram_samtools.d ./src/htslib/cram/cram_samtools.o ./src/htslib/cram/cram_stats.d ./src/htslib/cram/cram_stats.o ./src/htslib/cram/files.d ./src/htslib/cram/files.o ./src/htslib/cram/mFILE.d ./src/htslib/cram/mFILE.o ./src/htslib/cram/open_trace_file.d ./src/htslib/cram/open_trace_file.o ./src/htslib/cram/pooled_alloc.d ./src/htslib/cram/pooled_alloc.o ./src/htslib/cram/rANS_static.d ./src/htslib/cram/rANS_static.o ./src/htslib/cram/sam_header.d ./src/htslib/cram/sam_header.o ./src/htslib/cram/string_alloc.d ./src/htslib/cram/string_alloc.o

.PHONY: clean-src-2f-htslib-2f-cram

