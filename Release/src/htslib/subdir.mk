################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../src/htslib/bcf_sr_sort.c \
../src/htslib/bgzf.c \
../src/htslib/bgzip.c \
../src/htslib/errmod.c \
../src/htslib/faidx.c \
../src/htslib/hfile.c \
../src/htslib/hfile_gcs.c \
../src/htslib/hfile_net.c \
../src/htslib/hts.c \
../src/htslib/hts_os.c \
../src/htslib/htsfile.c \
../src/htslib/kfunc.c \
../src/htslib/knetfile.c \
../src/htslib/kstring.c \
../src/htslib/md5.c \
../src/htslib/multipart.c \
../src/htslib/probaln.c \
../src/htslib/realn.c \
../src/htslib/regidx.c \
../src/htslib/sam.c \
../src/htslib/synced_bcf_reader.c \
../src/htslib/tabix.c \
../src/htslib/tbx.c \
../src/htslib/textutils.c \
../src/htslib/thread_pool.c \
../src/htslib/vcf.c \
../src/htslib/vcf_sweep.c \
../src/htslib/vcfutils.c 

C_DEPS += \
./src/htslib/bcf_sr_sort.d \
./src/htslib/bgzf.d \
./src/htslib/bgzip.d \
./src/htslib/errmod.d \
./src/htslib/faidx.d \
./src/htslib/hfile.d \
./src/htslib/hfile_gcs.d \
./src/htslib/hfile_net.d \
./src/htslib/hts.d \
./src/htslib/hts_os.d \
./src/htslib/htsfile.d \
./src/htslib/kfunc.d \
./src/htslib/knetfile.d \
./src/htslib/kstring.d \
./src/htslib/md5.d \
./src/htslib/multipart.d \
./src/htslib/probaln.d \
./src/htslib/realn.d \
./src/htslib/regidx.d \
./src/htslib/sam.d \
./src/htslib/synced_bcf_reader.d \
./src/htslib/tabix.d \
./src/htslib/tbx.d \
./src/htslib/textutils.d \
./src/htslib/thread_pool.d \
./src/htslib/vcf.d \
./src/htslib/vcf_sweep.d \
./src/htslib/vcfutils.d 

OBJS += \
./src/htslib/bcf_sr_sort.o \
./src/htslib/bgzf.o \
./src/htslib/bgzip.o \
./src/htslib/errmod.o \
./src/htslib/faidx.o \
./src/htslib/hfile.o \
./src/htslib/hfile_gcs.o \
./src/htslib/hfile_net.o \
./src/htslib/hts.o \
./src/htslib/hts_os.o \
./src/htslib/htsfile.o \
./src/htslib/kfunc.o \
./src/htslib/knetfile.o \
./src/htslib/kstring.o \
./src/htslib/md5.o \
./src/htslib/multipart.o \
./src/htslib/probaln.o \
./src/htslib/realn.o \
./src/htslib/regidx.o \
./src/htslib/sam.o \
./src/htslib/synced_bcf_reader.o \
./src/htslib/tabix.o \
./src/htslib/tbx.o \
./src/htslib/textutils.o \
./src/htslib/thread_pool.o \
./src/htslib/vcf.o \
./src/htslib/vcf_sweep.o \
./src/htslib/vcfutils.o 


# Each subdirectory must supply rules for building sources it contributes
src/htslib/%.o: ../src/htslib/%.c src/htslib/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: Cross GCC Compiler'
	gcc -I../src/htslib -I../src -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


clean: clean-src-2f-htslib

clean-src-2f-htslib:
	-$(RM) ./src/htslib/bcf_sr_sort.d ./src/htslib/bcf_sr_sort.o ./src/htslib/bgzf.d ./src/htslib/bgzf.o ./src/htslib/bgzip.d ./src/htslib/bgzip.o ./src/htslib/errmod.d ./src/htslib/errmod.o ./src/htslib/faidx.d ./src/htslib/faidx.o ./src/htslib/hfile.d ./src/htslib/hfile.o ./src/htslib/hfile_gcs.d ./src/htslib/hfile_gcs.o ./src/htslib/hfile_net.d ./src/htslib/hfile_net.o ./src/htslib/hts.d ./src/htslib/hts.o ./src/htslib/hts_os.d ./src/htslib/hts_os.o ./src/htslib/htsfile.d ./src/htslib/htsfile.o ./src/htslib/kfunc.d ./src/htslib/kfunc.o ./src/htslib/knetfile.d ./src/htslib/knetfile.o ./src/htslib/kstring.d ./src/htslib/kstring.o ./src/htslib/md5.d ./src/htslib/md5.o ./src/htslib/multipart.d ./src/htslib/multipart.o ./src/htslib/probaln.d ./src/htslib/probaln.o ./src/htslib/realn.d ./src/htslib/realn.o ./src/htslib/regidx.d ./src/htslib/regidx.o ./src/htslib/sam.d ./src/htslib/sam.o ./src/htslib/synced_bcf_reader.d ./src/htslib/synced_bcf_reader.o ./src/htslib/tabix.d ./src/htslib/tabix.o ./src/htslib/tbx.d ./src/htslib/tbx.o ./src/htslib/textutils.d ./src/htslib/textutils.o ./src/htslib/thread_pool.d ./src/htslib/thread_pool.o ./src/htslib/vcf.d ./src/htslib/vcf.o ./src/htslib/vcf_sweep.d ./src/htslib/vcf_sweep.o ./src/htslib/vcfutils.d ./src/htslib/vcfutils.o

.PHONY: clean-src-2f-htslib

