Compile with:
```
git clone https://github.com/project-gemmi/gemmi.git
make
```
or

```
c++ -std=c++11 -Igemmi/include -O2 parse_pdb.cc   -o parse_pdb
c++ -std=c++11 -Igemmi/include -O2 parse_mmcif.cc -o parse_mmcif
# etc.
```
