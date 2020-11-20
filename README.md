<table>
<tr><td>SuperMoon.rmd</td><td> </td><td>R code, intended for use in RStudio.</td></tr>
<tr><td>SuperMoon.md</td><td> </td><td>Github-browsable code with example results, including 2001-2100 supermoon listing.</td></tr>
</table>

The R program is supported by files implementing ELP-2000/82 (the Fortran program elp82b_2.f and associated data files ELP1 to ELP36, see Chapront-Touzé, M., & Chapront, J. 1983, A&A, 124, 50). It's assumed elp82b_2.f has been compiled into "elp82b_2" (either a shared object on Unix/Linux or dll on Windows), in the same directory as the R program. An elp82b_2.dll, compiled with MinGW-w64, has been included in this repository and may save some Windows users from compiling elp82b_2.f themselves.

Reference: Stenborg, T.N., "21st Century Supermoon Estimation in R", in JE Ruiz and F Pierfederici (eds), Astronomical Data Analysis Software and Systems XXX, Granada, Spain.
