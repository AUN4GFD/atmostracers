gcc -c -fPIC addcomp.c -l /lib/surface/libsurface.so -Wl,-rpath=/lib -Wall -o addcomp.o
gcc addcomp.o -shared -o libaddcomp.so
rm addcomp.o
