
#!/usr/bin/env bash

if [ $# -eq 0 ]; then
    echo "Podaj jako argumenty plik wejsciowy i wyjsciowy lub all by skompilowac wszystko"
    exit 1
fi
if [ $# -eq 2 ]; then
    echo "arm-linux-gnueabihf-gcc -static -o $2 -mfpu=neon -mcpu=cortex-a9 $1"
    arm-linux-gnueabihf-gcc -static -o $2 -mfpu=neon -mcpu=cortex-a9 $1
    exit 1
fi

if [ $1 == "all" ]; then
    echo "Kompilowanie wszystkich plikow C i C++"
	for i in *.c
	do
	    echo "arm-linux-gnueabihf-gcc -static -o ${i%.c}.out -mfpu=neon -mcpu=cortex-a9 $i"
	    arm-linux-gnueabihf-gcc -static -o ${i%.c}.out -mfpu=neon -mcpu=cortex-a9 $i -lm
	done


	for i in *.cpp
	do
	    echo "arm-linux-gnueabihf-gcc -static -o ${i%.cpp}.out -mfpu=neon -mcpu=cortex-a9 $i"
	    arm-linux-gnueabihf-g++	 -static -o ${i%.cpp}+.out -mfpu=neon -mcpu=cortex-a9 $i
	done
    exit 1

fi


