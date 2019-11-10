# matrix-inversion

### Usage
1. Clone repository
2. cd matrix-inversion
3. mkdir build && cd build
4. cmake .. && make -j4
5. ./matrix_inversion [DIMENSION]\
   Use /usr/bin/time to measure execution time\
   **Example:**\
   /usr/bin/time ./matrix-inversion 1000

### Debugging
1. Build program with symbols - Debug version in Cmake file
2. Install Valgrind and kCachegrind

For Debian distros:
```console
sudo apt install valgrind kcachegrind dot
```

For Redhat distros:
```console
sudo dnf install valgrind kcachegrind dot
```

3. Run valgrind:
```console
valgrind --tool=cachegrind ./matrix-inversion 100
```
4. Open cachegrind.out.XXX file.
```console
kcachegrind cachegrind.out.16901
```
