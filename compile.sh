
#!/usr/bin/env bash


if [ $1 == "all" ]; then
	for i in *.c
	do
	    echo "arm-linux-gnueabihf-gcc -static -o ${i%.c}.out -mfpu=neon -mcpu=cortex-a9 $i"
	    arm-linux-gnueabihf-gcc -static -o ${i%.c}.out -mfpu=neon -mcpu=cortex-a9 $i -lm
	done
    exit 1

fi


