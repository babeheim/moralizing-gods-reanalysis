```markdown
## System requirements

This project uses `renv` to reproduce the R package environment. Some R packages also require system-level compilers that are not managed by `renv`.

### macOS

Install the Apple Command Line Tools:

```bash
xcode-select --install
```

For Apple Silicon Macs using R 4.6.x, install GNU Fortran 14.2 from the official R for macOS tools page:

https://mac.r-project.org/tools/

The relevant installer is:

```text
gfortran-14.2-universal.pkg
```

After installation, verify that R can find the Fortran compiler:

```r
system("R CMD config FC")
system("R CMD config F77")
```

These should point to a `gfortran` installation, typically under:

```text
/opt/gfortran/bin/gfortran
```

Then restore the project R environment:

```r
renv::restore()
```

Note that `renv.lock` records R and R-package dependencies, but does not install external system tools such as `clang`, `make`, or `gfortran`.
```